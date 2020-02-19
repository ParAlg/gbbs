#pragma once

#include "bridge.h"
#include "edge_map_utils.h"
#include "graph.h"
#include "vertex_subset.h"
#include "pbbslib/list_allocator.h"
#include "pbbslib/binary_search.h"

#include <vector>

template <class data  /* data associated with vertices in the output vertex_subset */,
          class Graph /* graph type */,
          class VS    /* vertex_subset type */,
          class F     /* edgeMap struct */>
inline vertexSubsetData<data> edgeMapSparse(Graph& G,
                                            VS& indices, F& f,
                                            const flags fl) {
  using S = std::tuple<uintE, data>;
  size_t n = indices.n;
  size_t m = indices.size();
  S* outEdges;

  if (should_output(fl)) {
    auto offsets = sequence<uintT>(indices.size(), [&](size_t i) {
      return (fl & in_edges) ? G.get_vertex(indices.vtx(i)).getInDegree()
                             : G.get_vertex(indices.vtx(i)).getOutDegree();
    });
    size_t outEdgeCount = pbbslib::scan_add_inplace(offsets.slice());

    outEdges = pbbslib::new_array_no_init<S>(outEdgeCount);
    auto g = get_emsparse_gen_full<data>(outEdges);
    auto h = get_emsparse_gen_empty<data>(outEdges);
    par_for(0, m, 1, [&] (size_t i) {
      uintT v = indices.vtx(i);
      uintT o = offsets[i];
      auto vert = G.get_vertex(v);
      (fl & in_edges) ? vert.decodeInNghSparse(v, o, f, g, h)
                      : vert.decodeOutNghSparse(v, o, f, g, h);
    });

    S* nextIndices = pbbslib::new_array_no_init<S>(outEdgeCount);
    auto p = [](std::tuple<uintE, data>& v) {
      return std::get<0>(v) != UINT_E_MAX;
    };
    size_t nextM = pbbslib::filterf(outEdges, nextIndices, outEdgeCount, p);
    pbbslib::free_array(outEdges);
    return vertexSubsetData<data>(n, nextM, nextIndices);
  }

  auto g = get_emsparse_nooutput_gen<data>();
  auto h = get_emsparse_nooutput_gen_empty<data>();
  par_for(0, m, 1, [&] (size_t i) {
    uintT v = indices.vtx(i);
    auto vert = G.get_vertex(v);
    (fl & in_edges) ? vert.decodeInNghSparse(v, 0, f, g, h)
                    : vert.decodeOutNghSparse(v, 0, f, g, h);
  });
  return vertexSubsetData<data>(n);
}

template <class data  /* data associated with vertices in the output vertex_subset */,
          class Graph /* graph type */,
          class VS    /* vertex_subset type */,
          class F     /* edgeMap struct */>
inline vertexSubsetData<data> edgeMapSparseNoOutput(Graph& G, VS& indices, F& f,
                                                    const flags fl) {
  size_t m = indices.numNonzeros();
#ifdef NVM
  bool inner_parallel = false;
#else
  bool inner_parallel = true;
#endif
  auto n = G.n;
  auto g = get_emsparse_nooutput_gen<data>();
  auto h = get_emsparse_nooutput_gen_empty<data>();
  par_for(0, m, 1, [&](size_t i) {
    uintT v = indices.vtx(i);
    auto vert = G.get_vertex(v);
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


template <class data  /* data associated with vertices in the output vertex_subset */,
          class Graph /* graph type */,
          class VS    /* vertex_subset type */,
          class F     /* edgeMap struct */>
inline vertexSubsetData<data> edgeMapBlocked(Graph& G, VS& indices, F& f,
                                             const flags fl) {
  if (fl & no_output) {
    return edgeMapSparseNoOutput<data, Graph, VS, F>(G, indices, f, fl);
  }
  using S = std::tuple<uintE, data>;
  size_t n = indices.n;

  auto block_f = [&](size_t i) -> size_t {
    return (fl & in_edges) ? G.get_vertex(indices.vtx(i)).getNumInBlocks()
                           : G.get_vertex(indices.vtx(i)).getNumOutBlocks();
  };
  auto block_imap = pbbslib::make_sequence<uintE>(indices.size(), block_f);

  // 1. Compute the number of blocks each vertex is subdivided into.
  auto vertex_offs = sequence<uintE>(indices.size() + 1);
  par_for(0, indices.size(), pbbslib::kSequentialForThreshold,
          [&](size_t i) { vertex_offs[i] = block_imap[i]; });
  vertex_offs[indices.size()] = 0;
  size_t num_blocks = pbbslib::scan_add_inplace(vertex_offs.slice());
  cout << "# num_blocks = " << num_blocks << endl;

  auto blocks = sequence<block>(num_blocks);
  auto degrees = sequence<uintT>(num_blocks);

  // 2. Write each block to blocks and scan degree array.
  par_for(0, indices.size(), pbbslib::kSequentialForThreshold, [&](size_t i) {
    size_t vtx_off = vertex_offs[i];
    size_t num_vertex_blocks = vertex_offs[i + 1] - vtx_off;
    uintE vtx_id = indices.vtx(i);
    assert(vtx_id < n);
    auto vtx = G.get_vertex(vtx_id);
    par_for(0, num_vertex_blocks, pbbslib::kSequentialForThreshold, [&](size_t j) {
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

        auto our_vtx = G.get_vertex(vtx_id);
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

constexpr size_t kDataBlockSizeBytes = 16384;
struct em_data_block {
  size_t block_size;
  uint8_t data[kDataBlockSizeBytes];
};
using data_block_allocator = list_allocator<em_data_block>;


// block format:
// size_t block_size (8 bytes for alignment)
// remainder is used as std::tuple<uintE, data>

template <class data, class Graph>
struct emhelper {
  using thread_blocks = std::vector<em_data_block*>;
  using ngh_data = std::tuple<uintE, data>;
  static constexpr size_t max_block_size = kDataBlockSizeBytes/sizeof(ngh_data);

  using vertex = typename Graph::vertex;
  static constexpr size_t work_block_size = vertex::getInternalBlockSize();

  bool alloc;
  uintE n_groups;
  thread_blocks* perthread_blocks;
  size_t* perthread_counts;

  static constexpr size_t kPerThreadStride = 128/sizeof(size_t);
  static constexpr size_t kThreadBlockStride = 1; //128/sizeof(thread_blocks);

  emhelper(size_t n_groups) : n_groups(n_groups) {
    perthread_blocks = pbbs::new_array<thread_blocks>(n_groups*kThreadBlockStride);
    perthread_counts = pbbs::new_array<size_t>(n_groups);
    alloc = true;
  }

  void del() {
    if (alloc) {
      pbbs::free_array(perthread_blocks);
      pbbs::free_array(perthread_counts);
    }
  }

  inline size_t scan_perthread_blocks() {
    size_t ct = 0;
    for (size_t i=0; i<n_groups; i++) {
      size_t i_sz = perthread_blocks[i*kThreadBlockStride].size();
      perthread_counts[i] = ct;
      ct += i_sz;
    }
    return ct;
  }

  auto get_all_blocks() {
    size_t total_blocks = scan_perthread_blocks();
    auto all_blocks = pbbs::sequence<em_data_block*>(total_blocks);
    if (total_blocks < 1000) { // handle sequentially
      size_t k=0;
      for (size_t i=0; i<n_groups; i++) {
        auto& vec = perthread_blocks[i*kThreadBlockStride];
        size_t this_thread_size = vec.size();
        for (size_t j=0; j<this_thread_size; j++) {
          all_blocks[k++] = vec[j];
        }
        vec.clear();
      }
    } else {
      parallel_for(0, n_groups, [&] (size_t thread_id) {
        size_t this_thread_offset = perthread_counts[thread_id];
        auto& vec = perthread_blocks[thread_id*kThreadBlockStride];
        size_t this_thread_size = vec.size();
        for (size_t j=0; j<this_thread_size; j++) {
          all_blocks[this_thread_offset + j] = vec[j];
        }
        vec.clear();
      }, 1);
    }
    return std::move(all_blocks);
  }


  // returns the next block for this group, reallocates if nec.
  em_data_block* get_block_and_offset_for_group(size_t group_id) {
     // fetch current block for thread
    em_data_block* block_ptr;
    size_t offset = 0;
    auto& vec = perthread_blocks[group_id*kThreadBlockStride];
    if (vec.size() == 0) { // alloc new
      block_ptr = data_block_allocator::alloc();
      vec.emplace_back(block_ptr);
      block_ptr->block_size = 0;
    } else {
      block_ptr = vec.back();
      offset = block_ptr->block_size;
      if (offset + work_block_size > max_block_size) { // realloc
        block_ptr = data_block_allocator::alloc();
        vec.emplace_back(block_ptr);
        block_ptr->block_size = 0;
        offset = 0;
      }
    }
    return block_ptr;
  }
};

template <class Graph>
void alloc_init(Graph& G) {
  size_t uintes_per_block = kDataBlockSizeBytes/sizeof(uintE);
  size_t list_alloc_init_blocks = std::max(static_cast<size_t>(0.5 * (G.n/uintes_per_block)),
      static_cast<size_t>(1000));
  cout << "# list_alloc init_blocks: " << list_alloc_init_blocks << endl;
  data_block_allocator::reserve(list_alloc_init_blocks);
  cout << "# after init: " << endl;
  data_block_allocator::print_stats();
}

template <class data  /* data associated with vertices in the output vertex_subset */,
          class Graph /* graph type */,
          class VS    /* vertex_subset type */,
          class F     /* edgeMap struct */>
inline vertexSubsetData<data> edgeMapChunked(Graph& G, VS& indices, F& f,
                                             const flags fl) {
  if (fl & no_output) {
    return edgeMapSparseNoOutput<data, Graph, VS, F>(G, indices, f, fl);
  }
  using S = std::tuple<uintE, data>;
  size_t n = indices.n;

  auto block_f = [&](size_t i) -> size_t {
    return (fl & in_edges) ? G.get_vertex(indices.vtx(i)).getNumInBlocks()
                           : G.get_vertex(indices.vtx(i)).getNumOutBlocks();
  };
  auto block_imap = pbbslib::make_sequence<uintE>(indices.size(), block_f);

  // 1. Compute the number of blocks each vertex is subdivided into.
  auto vertex_offs = sequence<uintE>(indices.size() + 1);
  par_for(0, indices.size(), pbbslib::kSequentialForThreshold,
          [&](size_t i) { vertex_offs[i] = block_imap[i]; });
  vertex_offs[indices.size()] = 0;
  size_t num_blocks = pbbslib::scan_add_inplace(vertex_offs.slice());

  auto blocks = sequence<block>(num_blocks);
  auto degrees = sequence<uintT>(num_blocks);

  // 2. Write each block to blocks and scan degree array.
  par_for(0, indices.size(), pbbslib::kSequentialForThreshold, [&](size_t i) {
    size_t vtx_off = vertex_offs[i];
    size_t num_vertex_blocks = vertex_offs[i + 1] - vtx_off;
    uintE vtx_id = indices.vtx(i);
    assert(vtx_id < n);
    auto vtx = G.get_vertex(vtx_id);
    par_for(0, num_vertex_blocks, pbbslib::kSequentialForThreshold, [&](size_t j) {
      size_t block_deg = (fl & in_edges)
                             ? vtx.in_block_degree(j)
                             : vtx.out_block_degree(j);
      // assert(block_deg <= PARALLEL_DEGREE); // only for compressed
      blocks[vtx_off + j] = block(i, j);  // j-th block of the i-th vertex.
      degrees[vtx_off + j] = block_deg;
    });
  });
  pbbslib::scan_add_inplace(degrees.slice(), pbbslib::fl_scan_inclusive);
  vertex_offs.clear();
  size_t outEdgeCount = degrees[num_blocks - 1];

  // 3. Compute the number of threads, binary search for offsets.
  // try to use 8*p threads, less only if guess'd blocksize is smaller than kEMBlockSize
  size_t edge_block_size_guess = pbbs::num_blocks(outEdgeCount, num_workers() << 3);
  size_t edge_block_size = std::max(kEMBlockSize, edge_block_size_guess);
  size_t n_groups = pbbs::num_blocks(outEdgeCount, edge_block_size);

//  cout << "outEdgeCount = " << outEdgeCount << endl;
//  cout << "n_blocks = " << num_blocks << endl;
//  cout << "n_groups = " << n_groups << endl;

  auto our_emhelper = emhelper<data, Graph>(n_groups);

  // Run each thread in parallel
  auto lt = [](const uintT& l, const uintT& r) { return l < r; };
  parallel_for(0, n_groups, [&](size_t group_id) {
    size_t start_off = group_id * edge_block_size;
    size_t our_start = pbbslib::binary_search(degrees, start_off, lt);
    size_t our_end;
    if (group_id < (n_groups - 1)) {
      size_t next_start_off = (group_id+1) * edge_block_size;
      our_end = pbbslib::binary_search(degrees, next_start_off, lt);
    } else {
      our_end = num_blocks;
    }

    // <= block_size edges in this range, sequentially process
    if (our_start != our_end && our_start != num_blocks) {
      for (size_t work_id = our_start; work_id < our_end; work_id++) {
        // 1. before starting next work block check whether we need to reallocate
        // the output block. This guarantees that there is enough space in the
        // output block even if all items in the work block are written out
        em_data_block* out_block = our_emhelper.get_block_and_offset_for_group(group_id);
        size_t offset = out_block->block_size;
        auto out_block_data = (std::tuple<uintE, data>*)out_block->data;

        auto g = get_emblock_gen<data>(out_block_data);

        // 2. process the work block
        auto& block = blocks[work_id];
        uintE id = block.id;  // id in vset
        uintE block_num = block.block_num;
        uintE vtx_id = indices.vtx(id);  // actual vtx_id corresponding to id
        auto vtx = G.get_vertex(vtx_id);
        size_t num_in = (fl & in_edges)
          ? vtx.decodeInBlock(vtx_id, offset, block_num, f, g)
          : vtx.decodeOutBlock(vtx_id, offset, block_num, f, g);
        out_block->block_size += num_in;
      }
    }
  }, 1);

  // scan the #output blocks/thread
  sequence<em_data_block*> all_blocks = our_emhelper.get_all_blocks();
  auto block_offsets = pbbs::sequence<size_t>(all_blocks.size(), [&] (size_t i) {
    return all_blocks[i]->block_size;
  });
  size_t output_size = pbbslib::scan_add_inplace(block_offsets.slice());
  vertexSubsetData<data> ret(n);
  if (output_size > 0) {
    S* out = pbbslib::new_array_no_init<S>(output_size);

    parallel_for(0, all_blocks.size(), [&] (size_t block_id) {
      em_data_block* block = all_blocks[block_id];
      size_t block_size = block->block_size;
      std::tuple<uintE, data>* block_data = (std::tuple<uintE, data>*)block->data;
      size_t block_offset = block_offsets[block_id];
      for (size_t i=0; i<block_size; i++) {
        out[block_offset + i] = block_data[i];
      }
      // deallocate block to list_alloc
      data_block_allocator::free(block);
    }, 1);
    ret = vertexSubsetData<data>(n, output_size, out);
  } else {
    parallel_for(0, all_blocks.size(), [&] (size_t block_id) {
      em_data_block* block = all_blocks[block_id];
      data_block_allocator::free(block);
    });
  }

  all_blocks.clear();
  block_offsets.clear();
  our_emhelper.del();

//  our_em_block.reset(); (handled by get_all_blocks)
  return std::move(ret);
}
