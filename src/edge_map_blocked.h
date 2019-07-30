#pragma once

#include "bridge.h"
#include "edge_map_utils.h"
#include "graph.h"
#include "vertex_subset.h"

template <
    class data /* data associated with vertices in the output vertex_subset */,
    class G /* graph type */, class VS /* vertex_subset type */,
    class F /* edgeMap struct */>
inline vertexSubsetData<data> edgeMapSparseNoOutput(G& GA, VS& indices, F& f,
                                                    const flags fl) {
  size_t m = indices.numNonzeros();
#ifdef NVM
  bool inner_parallel = false;
#else
  bool inner_parallel = true;
#endif
  auto n = GA.n;
  auto g = get_emsparse_nooutput_gen<data>();
  auto h = get_emsparse_nooutput_gen_empty<data>();
  par_for(0, m, 1, [&](size_t i) {
    uintT v = indices.vtx(i);
    auto vert = GA.get_vertex(v);
    (fl & in_edges) ? vert.decodeInNghSparse(v, 0, f, g, h, inner_parallel)
                    : vert.decodeOutNghSparse(v, 0, f, g, h, inner_parallel);
  });
  return vertexSubsetData<data>(n);
}

struct block {
  uintE id;
  uintE block_num;
  block(uintE _id, uintE _b) : id(_id), block_num(_b) {}
  block() {}
  void print() { std::cout << id << " " << block_num << "\n"; }
};

template <
    class data /* data associated with vertices in the output vertex_subset */,
    class G /* graph type */, class VS /* vertex_subset type */,
    class F /* edgeMap struct */>
inline vertexSubsetData<data> edgeMapBlocked(G& GA, VS& indices, F& f,
                                             const flags fl) {
  if (fl & no_output) {
    return edgeMapSparseNoOutput<data, G, VS, F>(GA, indices, f, fl);
  }
  using S = std::tuple<uintE, data>;
  size_t n = indices.n;

  auto block_f = [&](size_t i) -> size_t {
    return (fl & in_edges) ? GA.get_vertex(indices.vtx(i)).getNumInBlocks()
                           : GA.get_vertex(indices.vtx(i)).getNumOutBlocks();
  };
  auto block_imap = pbbslib::make_sequence<uintE>(indices.size(), block_f);

  // 1. Compute the number of blocks each vertex is subdivided into.
  auto vertex_offs = sequence<uintE>(indices.size() + 1);
  par_for(0, indices.size(), pbbslib::kSequentialForThreshold,
          [&](size_t i) { vertex_offs[i] = block_imap[i]; });
  vertex_offs[indices.size()] = 0;
  size_t num_blocks = pbbslib::scan_add_inplace(vertex_offs.slice());
  cout << "num_blocks = " << num_blocks << endl;

  auto blocks = sequence<block>(num_blocks);
  auto degrees = sequence<uintT>(num_blocks);

  // 2. Write each block to blocks and scan degree array.
  par_for(0, indices.size(), pbbslib::kSequentialForThreshold, [&](size_t i) {
    size_t vtx_off = vertex_offs[i];
    size_t num_blocks = vertex_offs[i + 1] - vtx_off;
    uintE vtx_id = indices.vtx(i);
    assert(vtx_id < n);
    if (vtx_id == 978407842) {
      cout << "i = " << i << " indices.size = " << indices.size() << endl;
    }
    auto vtx = GA.get_vertex(vtx_id);
    par_for(0, num_blocks, pbbslib::kSequentialForThreshold, [&](size_t j) {
      size_t block_deg = (fl & in_edges)
                             ? vtx.in_block_degree(j)
                             : vtx.out_block_degree(j);
      // assert(block_deg <= PARALLEL_DEGREE); // only for compressed
      blocks[vtx_off + j] = block(i, j);  // j-th block of the i-th vertex.
      degrees[vtx_off + j] = block_deg;
    });
  });
  pbbslib::scan_add_inplace(degrees.slice(), pbbslib::fl_scan_inclusive);
  size_t outEdgeCount = degrees[num_blocks - 1];

  // 3. Compute the number of threads, binary search for offsets.
  size_t n_threads = pbbs::num_blocks(outEdgeCount, kEMBlockSize);
  size_t* thread_offs = pbbslib::new_array_no_init<size_t>(n_threads + 1);
  auto lt = [](const uintT& l, const uintT& r) { return l < r; };
  par_for(0, n_threads, 1, [&](size_t i) {  // TODO: granularity of 1?
    size_t start_off = i * kEMBlockSize;
    thread_offs[i] = pbbslib::binary_search(degrees, start_off, lt);
  });
  thread_offs[n_threads] = num_blocks;

  // 4. Run each thread in parallel
  auto cts = sequence<uintE>(n_threads + 1);
  S* outEdges = pbbslib::new_array_no_init<S>(outEdgeCount);
  auto g = get_emsparse_blocked_gen<data>(outEdges);
  par_for(0, n_threads, 1, [&](size_t i) {
    size_t start = thread_offs[i];
    size_t end = thread_offs[i + 1];
    // <= kEMBlockSize edges in this range, sequentially process
    if (start != end && start != num_blocks) {
      size_t start_offset = (start == 0) ? 0 : degrees[start - 1];
      size_t k = start_offset;
      for (size_t j = start; j < end; j++) {
        auto& block = blocks[j];
        uintE id = block.id;  // id in vset
        uintE block_num = block.block_num;

        uintE vtx_id = indices.vtx(id);  // actual vtx_id corresponding to id

        auto our_vtx = GA.get_vertex(vtx_id);
        size_t num_in =
            (fl & in_edges)
                ? our_vtx.decodeInBlock(vtx_id, k, block_num, f, g)
                : our_vtx.decodeOutBlock(vtx_id, k, block_num, f, g);
        k += num_in;
      }
      cts[i] = k - start_offset;
    } else {
      cts[i] = 0;
    }
  });
  cts[n_threads] = 0;
  size_t out_size = pbbslib::scan_add_inplace(cts.slice());

  // 5. Use cts to get
  S* out = pbbslib::new_array_no_init<S>(out_size);
  cout << "outEdgeCount (blocked) = " << (indices.numNonzeros() + outEdgeCount) << " bytes allocated = " << (sizeof(S) * outEdgeCount) << " only needed: " << (sizeof(S)*out_size) << endl;
  par_for(0, n_threads, 1, [&](size_t i) {
    size_t start = thread_offs[i];
    size_t end = thread_offs[i + 1];
    if (start != end) {
      size_t start_offset = (start == 0) ? 0 : degrees[start - 1];
      size_t out_offset = cts[i];
      size_t num_live = cts[i + 1] - out_offset;
      for (size_t j = 0; j < num_live; j++) {
        out[out_offset + j] = outEdges[start_offset + j];
      }
    }
  });
  pbbslib::free_array(outEdges);
  pbbslib::free_array(thread_offs);
  cts.clear();
  vertex_offs.clear();
  blocks.clear();
  degrees.clear();

  return vertexSubsetData<data>(n, out_size, out);
}
