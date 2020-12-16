// This code is part of the project "Theoretically Efficient Parallel Graph
// Algorithms Can Be Fast and Scalable", presented at Symposium on Parallelism
// in Algorithms and Architectures, 2018.
// Copyright (c) 2018 Laxman Dhulipala, Guy Blelloch, and Julian Shun
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all  copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#pragma once

#include <math.h>
#include <limits>

// Library dependencies
#include "gbbs/bucket.h"
#include "gbbs/edge_map_reduce.h"
#include "gbbs/gbbs.h"
#include "gbbs/pbbslib/dyn_arr.h"
#include "gbbs/pbbslib/sparse_table.h"
#include "gbbs/pbbslib/sparse_additive_map.h"
#include "pbbslib/assert.h"
#include "pbbslib/list_allocator.h"
#include "pbbslib/integer_sort.h"

// Ordering files
#include "benchmarks/DegeneracyOrder/BarenboimElkin08/DegeneracyOrder.h"
#include "benchmarks/DegeneracyOrder/GoodrichPszona11/DegeneracyOrder.h"
#include "benchmarks/KCore/JulienneDBS17/KCore.h"
#include "benchmarks/TriangleCounting/ShunTangwongsan15/Triangle.h"
#include "benchmarks/CliqueCounting/Clique.h"

// Clique files
#include "benchmarks/CliqueCounting/intersect.h"
#include "benchmarks/CliqueCounting/induced_intersection.h"
#include "benchmarks/CliqueCounting/induced_neighborhood.h"
#include "benchmarks/CliqueCounting/induced_hybrid.h"
#include "benchmarks/CliqueCounting/induced_split.h"
#include "benchmarks/CliqueCounting/relabel.h"
namespace gbbs {
  struct hash128 {
    inline size_t operator () (unsigned __int128 t) const {
      size_t l = t >> 64;
      unsigned __int128 mask = l << 64;
      size_t r = t & (~mask);
      return pbbs::hash_combine(pbbslib::hash64_2(l), pbbslib::hash64_2(r));
    }
  };

  struct EndTable {
    pbbslib::sparse_table<unsigned __int128, long, hash128> table;
    uintE vtx;
    //MidTable* up_table;
  };

  struct MidTable {
    pbbslib::sparse_table<uintE, EndTable*, std::hash<uintE>> table;
    //MidTable* up_table;
  };
  
  class TwolevelHash {
    public:
      using T = pbbslib::sparse_table<unsigned __int128, long, hash128>;
      using X = std::tuple<unsigned __int128, long>;
      MidTable top_table;
      sequence<long> top_table_sizes;
      int rr;
      std::size_t total = 0;
      X* space = nullptr;
  
      template <class Graph>
      TwolevelHash(int r, Graph& DG, size_t max_deg) {
        rr = r;
        //top_table.up_table = nullptr;
        // How many vert in top level?
        // For each vert in top level, how many pairs of vert finish it?
        auto tmp_table = pbbslib::sparse_additive_map<uintE, long>(
          DG.n, std::make_tuple(UINT_E_MAX, long{0}));
        auto base_f = [&](sequence<uintE>& base){
          auto min_vert = pbbslib::reduce_min(base);
          auto tmp = std::make_tuple<uintE, long>(static_cast<uintE>(min_vert), long{1});
          tmp_table.insert(tmp);
        };
        auto init_induced = [&](HybridSpace_lw* induced) { induced->alloc(max_deg, r-1, DG.n, true, true); };
        auto finish_induced = [&](HybridSpace_lw* induced) { if (induced != nullptr) { delete induced; } };
        parallel_for_alloc<HybridSpace_lw>(init_induced, finish_induced, 0, DG.n, [&](size_t i, HybridSpace_lw* induced) {
          if (DG.get_vertex(i).getOutDegree() != 0) {
            induced->setup(DG, r-1, i);
            auto base = sequence<uintE>(r);
            base[0] = i;
            NKCliqueDir_fast_hybrid_rec(DG, 1, r-1, induced, base_f, base);
          }
        }, 1, false);
        auto top_table_sizes2 = tmp_table.entries();
        // Modify top_table_sizes2 to be appropriately oversized
        parallel_for(0, top_table_sizes2.size(), [&](std::size_t i) {
          auto m = (size_t)1 << pbbslib::log2_up((size_t)(1.1 * std::get<1>(top_table_sizes2[i])) + 1);
          top_table_sizes2[i] = std::make_tuple(std::get<0>(top_table_sizes2[i]), m);
        });
        // Do a scan inplace
        auto add_tuple_monoid = pbbs::make_monoid([](std::tuple<uintE, long> a, std::tuple<uintE, long> b){
          return std::make_tuple(std::get<0>(b), std::get<1>(a) + std::get<1>(b));
        }, std::make_tuple(UINT_E_MAX, 0));
        std::tuple<uintE, long> total_top_table_sizes2 = scan_inplace(top_table_sizes2.slice(), add_tuple_monoid);
        // Allocate space for the second level tables
        space = pbbslib::new_array_no_init<X>(std::get<1>(total_top_table_sizes2));
        tmp_table.del();
        top_table.table = pbbslib::sparse_table<uintE, EndTable*, std::hash<uintE>>(
          top_table_sizes2.size(),
          std::make_tuple<uintE, EndTable*>(UINT_E_MAX, static_cast<EndTable*>(nullptr)),
          std::hash<uintE>());
        top_table_sizes = sequence<long>(top_table.table.m + 1, long{0});
        parallel_for(0, top_table_sizes2.size(), [&](std::size_t i){
          auto vtx = i == top_table_sizes2.size() - 1 ? std::get<0>(total_top_table_sizes2) : 
            std::get<0>(top_table_sizes2[i + 1]);
          auto upper_size = i == top_table_sizes2.size() - 1 ? std::get<1>(total_top_table_sizes2) : 
            std::get<1>(top_table_sizes2[i + 1]);
          auto size = upper_size - std::get<1>(top_table_sizes2[i]);
          EndTable* end_table = new EndTable();
          end_table->vtx = vtx;
          end_table->table = pbbslib::sparse_table<unsigned __int128, long, hash128>(
            size, 
            std::make_tuple<unsigned __int128, long>(static_cast<unsigned __int128>(0), static_cast<long>(0)),
            hash128{},
            space + std::get<1>(top_table_sizes2[i])
            );
            //pbbslib::sparse_table<unsigned __int128, long, hash128>(
            //size,
            //std::make_tuple<unsigned __int128, long>(static_cast<unsigned __int128>(0), static_cast<long>(0)),
            //hash128{});
          top_table.table.insert(std::make_tuple(vtx, end_table));
          std::size_t l = top_table.table.find_index(vtx);
          top_table_sizes[l] = end_table->table.m;
        });
        total = scan_inplace(top_table_sizes.slice(), pbbs::addm<long>());
        /*for (std::size_t i = 1; i < top_table_sizes.size(); i++) {
          assert(top_table_sizes[i - 1] <= top_table_sizes[i]);
          EndTable* end_table = std::get<1>(top_table.table.table[i - 1]);
          if (end_table != nullptr) {
            if (top_table_sizes[i] - top_table_sizes[i - 1] != (end_table->table).m) {
              std::cout << "M: " << (end_table->table).m << std::endl;
              std::cout << "Diff: " << top_table_sizes[i] - top_table_sizes[i - 1] << std::endl;
              fflush(stdout);
            }
            assert(top_table_sizes[i] - top_table_sizes[i - 1] == (end_table->table).m);
          }
        }
        assert(top_table_sizes[0] == 0);*/
      }
  
      void insert(sequence<uintE>& base, int r, int k) {
        auto add_f = [&] (long* ct, const std::tuple<unsigned __int128, long>& tup) {
          pbbs::fetch_and_add(ct, (long)1);
        };
        // Sort base
        sequence<uintE> base2(base);
        pbbs::sample_sort_inplace(base2.slice(), std::less<uintE>());

        std::string bitmask(r+1, 1); // K leading 1's
        bitmask.resize(k+1, 0); // N-K trailing 0's

        do {
          bool use_vtx = false;
          uintE vtx = 0;
          unsigned __int128 key = 0;
          for (int i = 0; i < static_cast<int>(k)+1; ++i) {
            if (bitmask[i]) {
              if (!use_vtx) {
                use_vtx = true;
                vtx = base2[i];
              } else {
                key = key << 32;
                key |= static_cast<int>(base2[i]);
              }
            }
          }
          EndTable* end_table = top_table.table.find(vtx, nullptr);
          //assert(end_table != nullptr);
          (end_table->table).insert_f(std::make_tuple(key, (long) 1), add_f);
        } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
      }

      std::size_t return_total() {return total;}

      size_t get_top_index(std::size_t index) {
        //assert(top_table_sizes[0] == 0);
        for (std::size_t i = 0; i < top_table_sizes.size(); i++) {
          //if (i > 0) {
          //  assert(top_table_sizes[i - 1] <= top_table_sizes[i]);
          //  EndTable* end_table = std::get<1>(top_table.table.table[i - 1]);
          //  if (end_table != nullptr) assert(top_table_sizes[i] - top_table_sizes[i - 1] == (end_table->table).m);
          //}
          if (top_table_sizes[i] == index) {
            return i;
          } else if (top_table_sizes[i] > index) {
            if (i == 0) return i;
            return i - 1;
          }
        }
        return top_table_sizes.size() - 1;
        //return pbbslib::binary_search(top_table_sizes, long{index}, std::greater<long>());
      }

      long get_count(std::size_t index) {
        //assert(index < total);
        size_t top_index = get_top_index(index);
        //assert(top_index != top_table_sizes.size());
        /*if (top_table_sizes[top_index] > index) {
          std::cout << "Size: " << top_table_sizes[top_index] << std::endl;
          std::cout << "Idx: " << index << std::endl;
          fflush(stdout);
        }
        assert(top_table_sizes[top_index] <= index);
        assert(top_index < top_table.table.m);
        if (top_index < top_table_sizes.size() - 1) {
          if ((top_table_sizes[top_index + 1] < index)) {
            std::cout << "Top: " << top_table_sizes[top_index + 1] << ", Idx: " << index << std::endl;
            fflush(stdout);
          }
          assert(top_table_sizes[top_index + 1] >= index);
        }*/
        EndTable* end_table = std::get<1>(top_table.table.table[top_index]);
        if (end_table == nullptr) return 0;
        size_t bottom_index = index - top_table_sizes[top_index];
        /*if (bottom_index >= (end_table->table).m) {
          std::cout << "Bottom: " << bottom_index << std::endl;
          std::cout << "M: " << (end_table->table).m << std::endl;
          fflush(stdout);
        }
        assert(bottom_index < (end_table->table).m);*/
        return std::get<1>((end_table->table).table[bottom_index]);
      }

      size_t update_count(std::size_t index, size_t update){
        size_t top_index = get_top_index(index);
        EndTable* end_table = std::get<1>(top_table.table.table[top_index]);
        size_t bottom_index = index - top_table_sizes[top_index];
        auto val = std::get<1>((end_table->table).table[bottom_index]) - update;
        (end_table->table).table[bottom_index] = std::make_tuple(
          std::get<0>((end_table->table).table[bottom_index]), val
        );
        return val;
      }

      void clear_count(std::size_t index) {
        size_t top_index = get_top_index(index);
        EndTable* end_table = std::get<1>(top_table.table.table[top_index]);
        size_t bottom_index = index - top_table_sizes[top_index];
        (end_table->table).table[bottom_index] = std::make_tuple(
          std::get<0>((end_table->table).table[bottom_index]), 0
        );
      }

      template<class I>
      void extract_indices(sequence<uintE>& base, I func, int r, int k) {
        sequence<uintE> base2(base);
        pbbs::sample_sort_inplace(base2.slice(), std::less<uintE>());
        std::string bitmask(r+1, 1); // K leading 1's
        bitmask.resize(k+1, 0); // N-K trailing 0's
        do {
          bool use_vtx = false;
          uintE vtx = 0;
          unsigned __int128 key = 0;
          for (int i = 0; i < static_cast<int>(k)+1; ++i) {
            if (bitmask[i]) {
              if (!use_vtx) {
                use_vtx = true;
                vtx = base2[i];
              } else {
                key = key << 32;
                key |= static_cast<int>(base2[i]);
              }
            }
          }
          // First, find index in top_table
          // This should populate into a prefix sum of sizes
          auto top_index = top_table.table.find_index(vtx);
          auto prefix = top_table_sizes[top_index];
          EndTable* end_table = std::get<1>(top_table.table.table[top_index]);
          auto index = (end_table->table).find_index(key);
          func(prefix + index);
        } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
      }

      // Given an index, get the clique
      template<class S, class Graph>
      void extract_clique(S index, sequence<uintE>& base, Graph& G, int k) {
        // First, do a binary search for index in prefix
        size_t top_index = get_top_index(index);
        base[0] = std::get<0>(top_table.table.table[top_index]);
        EndTable* end_table = std::get<1>(top_table.table.table[top_index]);
        size_t bottom_index = index - top_table_sizes[top_index];
        auto vert = std::get<0>((end_table->table).table[bottom_index]);
        for (int j = 0; j < rr - 1; ++j) {
          int extract = (int) vert; // vert & mask
          /*if (static_cast<uintE>(extract) >= G.n) {
            std::cout << "Vert: " << static_cast<uintE>(extract) << ", n: " << G.n << std::endl;
          }*/
          assert(static_cast<uintE>(extract) < G.n);
          base[k - j] = static_cast<uintE>(extract);
          vert = vert >> 32;
        }
      }
  };

  class OnelevelHash {
    public:
      using T = pbbslib::sparse_table<unsigned __int128, long, hash128>;
      T table;
      int rr;
      template <class Graph>
      OnelevelHash(int r, Graph& DG, size_t max_deg) {
        rr = r;
        // count cliques
        timer t_pre; t_pre.start();
        size_t pre_count = 0;
        // Clique counting
        if (r == 3) pre_count = TriClique_count(DG, false, nullptr);
        else pre_count = Clique_count(DG, r, 5, true, false, 0, nullptr);
        std::cout << "Pre count " << r << ": " << pre_count << std::endl;
        double tt_pre = t_pre.stop();
        std::cout << "### Pre count: " << tt_pre << std::endl;

        std::cout << "Start table" << std::endl;
        table = pbbslib::sparse_table<unsigned __int128, long, hash128>(
          pre_count,
          std::make_tuple(static_cast<unsigned __int128>(0), long{0}), hash128{});
      }

      void insert(sequence<uintE>& base, int r, int k) {
        auto add_f = [&] (long* ct, const std::tuple<unsigned __int128, long>& tup) {
          pbbs::fetch_and_add(ct, (long)1);
        };
        // Sort base
        sequence<uintE> base2(base);
        pbbs::sample_sort_inplace(base2.slice(), std::less<uintE>());

        std::string bitmask(r+1, 1); // K leading 1's
        bitmask.resize(k+1, 0); // N-K trailing 0's

        do {
          unsigned __int128 key = 0;
          for (int i = 0; i < static_cast<int>(k)+1; ++i) {
            if (bitmask[i]) {
              key = key << 32;
              key |= static_cast<int>(base2[i]);
            }
          }
          table.insert_f(std::make_tuple(key, (long) 1), add_f);
          // TESTING EXTRACTION OF KEy
          /*for (int j = 0; j < static_cast<int>(r)+1; ++j) {
            int extract = (int) key;
            if (static_cast<uintE>(extract) >= DG.n) {
              std::cout << "Extract: " << static_cast<uintE>(extract) << ", n: " << DG.n << std::endl;
            }
            assert(static_cast<uintE>(extract) < DG.n);
          key = key >> 32;
        }*/
        } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
      }

      std::size_t return_total() {return table.m;}
      long get_count(std::size_t i) {
        if (std::get<0>((table.table)[i]) == 0) return 0;
        return std::get<1>((table.table)[i]);
      }
      size_t update_count(std::size_t i, size_t update){
        auto val = std::get<1>(table.table[i]) - update;
        table.table[i] = std::make_tuple(std::get<0>(table.table[i]), val);
        return val;
      }
      void clear_count(std::size_t index) {
        table.table[index] = std::make_tuple(std::get<0>(table.table[index]),0);
      }

      template<class I>
      void extract_indices(sequence<uintE>& base, I func, int r, int k) {
        // Sort base
        sequence<uintE> base2(base);
        pbbs::sample_sort_inplace(base2.slice(), std::less<uintE>());
        std::string bitmask(r+1, 1); // K leading 1's
        bitmask.resize(k+1, 0); // N-K trailing 0's

        do {
          unsigned __int128 key = 0;
          for (int i = 0; i < static_cast<int>(k)+1; ++i) {
            if (bitmask[i]) {
              key = key << 32;
              key |= static_cast<int>(base2[i]);
            }
          }
          auto index = table.find_index(key);
          func(index);
        } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
      }
    
    template<class S, class Graph>
    void extract_clique(S index, sequence<uintE>& base, Graph& G, int k) {
      auto vert = std::get<0>(table.table[index]);
      for (int j = 0; j < rr; ++j) {
        int extract = (int) vert; // vert & mask
        //if (static_cast<uintE>(extract) >= G.n) {
        //  std::cout << "Vert: " << static_cast<uintE>(extract) << ", n: " << G.n << std::endl;
        //}
        assert(static_cast<uintE>(extract) < G.n);
        if (j == rr - 1) base[0] = static_cast<uintE>(extract);
        else base[k - j] = static_cast<uintE>(extract);
        vert = vert >> 32;
      }
    }
  };

  template<class Graph>
  size_t get_max_deg3(Graph& DG) {
    size_t max_deg = 0;
    parallel_for(0, DG.n, [&] (size_t i) {
      size_t deg = DG.get_vertex(i).getOutDegree();
      pbbs::write_min(&max_deg, deg, std::greater<size_t>());
    });
    return max_deg;
  }

  template <class Graph, class F>
  inline size_t NKCliqueDir_fast_hybrid_rec(Graph& DG, size_t k_idx, size_t k, HybridSpace_lw* induced, F base_f,
  sequence<uintE>& base) {
    size_t num_induced = induced->num_induced[k_idx-1];
    if (num_induced == 0) return 0;
    uintE* prev_induced = induced->induced + induced->nn * (k_idx - 1);

    for (size_t i=0; i < num_induced; i++) { induced->labels[prev_induced[i]] = k_idx; }

    if (k_idx + 1 == k) {
      size_t counts = 0;
      for (size_t i=0; i < num_induced; i++) {
        uintE vtx = prev_induced[i];
        //  get neighbors of vtx
        uintE* intersect = induced->induced_edges + vtx * induced->nn;
        size_t tmp_counts = 0;
        base[k_idx] = induced->relabel[vtx];
        for (size_t j=0; j < induced->induced_degs[vtx]; j++) {
          if (static_cast<size_t>(induced->labels[intersect[j]]) == k_idx) {
            base[k] = induced->relabel[intersect[j]];
            tmp_counts++;
            base_f(base);
          }
        }
        counts += tmp_counts;
      }
      for (size_t i=0; i < num_induced; i++) { induced->labels[prev_induced[i]] = k_idx - 1; }
      return counts;
    }

size_t total_ct = 0;
    for (size_t i=0; i < num_induced; ++i) {
      uintE vtx = prev_induced[i];
      uintE* intersect = induced->induced_edges + vtx * induced->nn;
      uintE* out = induced->induced + induced->nn * k_idx;
      uintE count = 0;
      for (size_t j=0; j < induced->induced_degs[vtx]; j++) {
        if (static_cast<size_t>(induced->labels[intersect[j]]) == k_idx) {
          out[count] = intersect[j];
          count++;
        }
      }
      induced->num_induced[k_idx] = count;
      if (induced->num_induced[k_idx] > k - k_idx - 1) {
        base[k_idx] = induced->relabel[vtx];
        auto curr_counts = NKCliqueDir_fast_hybrid_rec(DG, k_idx + 1, k, induced, base_f, base);
        total_ct += curr_counts;
      }
    }
    for (size_t i=0; i < num_induced; i++) { induced->labels[prev_induced[i]] = k_idx - 1; }
    return total_ct;
  }


  template <class Graph, class T>
  inline size_t CountCliquesNuc(Graph& DG, size_t k, size_t r, size_t max_deg, T* table) {
    k--; r--;
    timer t2; t2.start();

    auto add_f = [&] (long* ct, const std::tuple<unsigned __int128, long>& tup) {
      pbbs::fetch_and_add(ct, (long)1);
    };
// Nested hash tables where the first level you hash the first vertex, second level you hash the second, etc.
// Saves space if cliques share vertices
// Renumber vert by which is in the most r-cliques
    auto base_f = [&](sequence<uintE>& base){
      table->insert(base, r, k);
    };
    auto tots = sequence<size_t>(DG.n, size_t{0});

    auto init_induced = [&](HybridSpace_lw* induced) { induced->alloc(max_deg, k, DG.n, true, true); };
    auto finish_induced = [&](HybridSpace_lw* induced) { if (induced != nullptr) { delete induced; } }; //induced->del();
    parallel_for_alloc<HybridSpace_lw>(init_induced, finish_induced, 0, DG.n, [&](size_t i, HybridSpace_lw* induced) {
        if (DG.get_vertex(i).getOutDegree() != 0) {
          induced->setup(DG, k, i);
          auto base = sequence<uintE>(k + 1);
          base[0] = i;
          //auto base_f2 = [&](uintE vtx, size_t _count) {};
          //tots[i] = induced_hybrid::KCliqueDir_fast_hybrid_rec(DG, 1, k, induced, base_f2, 0);
          tots[i] = NKCliqueDir_fast_hybrid_rec(DG, 1, k, induced, base_f, base);
        } else tots[i] = 0;
    }, 1, false);
    double tt2 = t2.stop();
    std::cout << "##### Actual counting: " << tt2 << std::endl;

    return pbbslib::reduce_add(tots);
  }

class list_buffer {
  public:
    int buffer;
    sequence<size_t> list;
    sequence<size_t> starts;
    sequence<bool> to_pack;
    size_t next;
    size_t num_workers2;
    list_buffer(size_t s){
      num_workers2 = num_workers();
      buffer = 64;
      int buffer2 = 64;
      list = sequence<size_t>(s + buffer2 * num_workers2, static_cast<size_t>(0));
      starts = sequence<size_t>(num_workers2, [&](size_t i){return i * buffer2;});
      next = num_workers2 * buffer2;
      to_pack = sequence<bool>(s + buffer2 * num_workers2, true);
    }
    void add(size_t index) {
      size_t worker = worker_id();
      list[starts[worker]] = index;
      starts[worker]++;
      if (starts[worker] % buffer == 0) {
        size_t use_next = pbbs::fetch_and_add(&next, buffer);
        starts[worker] = use_next;
      }
    }

    template <class I>
    size_t filter(I update_changed, sequence<size_t>& per_processor_counts) {
      parallel_for(0, num_workers2, [&](size_t worker) {
        size_t divide = starts[worker] / buffer;
        for (size_t j = starts[worker]; j < (divide + 1) * buffer; j++) {
          to_pack[j] = false;
        }
      });
      // Pack out 0 to next of list into pack
      parallel_for(0, next, [&] (size_t i) {
        if (to_pack[i])
          update_changed(per_processor_counts, i, list[i]);
        else
          update_changed(per_processor_counts, i, UINT_E_MAX);
      });
      parallel_for(0, num_workers2, [&](size_t worker) {
        size_t divide = starts[worker] / buffer;
        for (size_t j = starts[worker]; j < (divide + 1) * buffer; j++) {
          to_pack[j] = true;
        }
      });
      return next;
    }

    void reset() {
      parallel_for (0, num_workers2, [&] (size_t j) {
        starts[j] = j * buffer;
      });
      next = num_workers2 * buffer;
    }
};

template <class Graph, class Graph2, class F, class H, class I, class T>
inline size_t cliqueUpdate(Graph& G, Graph2& DG, size_t r, 
size_t k, size_t max_deg, bool label, F get_active, size_t active_size,
  size_t granularity, char* still_active, sequence<uintE> &rank, 
  sequence<size_t>& per_processor_counts, H update,
  bool do_update_changed, I update_changed,
  T* cliques, size_t n, list_buffer& count_idxs, timer& t1) {

  // Set up space for clique counting
  auto init_induced = [&](HybridSpace_lw* induced) { induced->alloc(max_deg, k-r+1, G.n, true, true); };
  auto finish_induced = [&](HybridSpace_lw* induced) { if (induced != nullptr) { delete induced; } };

  // Mark every vertex in the active set
  parallel_for (0, active_size, [&] (size_t j) {
    auto index = get_active(j); //cliques->find_index(get_active(j));
    still_active[index] = 1;
    }, 2048);

  // Hash table to contain clique count updates
  //size_t edge_table_size = (size_t) (n);
  //auto edge_table = pbbslib::sparse_table<uintE, bool, hashtup>(edge_table_size, std::make_tuple(0, false), hashtup());

  // Function that dictates which edges to consider in first level of recursion
  //auto ignore_f = [&](const uintE& u, const uintE& v) {
  //  return true;
    //auto index_u = cliques->find_index(u);
    //auto index_v = cliques->find_index(v);
    //auto status_u = still_active[index_u]; auto status_v = still_active[index_v];
    //if (status_u == 2 || status_v == 2) return false; // deleted edge
    //if (status_v == 0) return true; // non-deleted, non-active edge
    //return rank[u] < rank[v]; // orient edges if in active set
  //};

  // Collate clique counts by processor
  //count_idxs[0] = 0;
  auto update_d = [&](sequence<uintE>& base){
    cliques->extract_indices(base, [&](std::size_t index){
      size_t ct = pbbs::fetch_and_add(&(per_processor_counts[index]), 1);
      if (ct == 0) count_idxs.add(index);
    }, r, k);
  };

t1.start();
  // Clique count updates
  parallel_for_alloc<HybridSpace_lw>(init_induced, finish_induced, 0, active_size,
                                     [&](size_t i, HybridSpace_lw* induced) {
    // TODO: THIS PART IS WRONG
    // you wanna start from the clique given by vert
    auto x = get_active(i);
    auto base = sequence<uintE>(k + 1);
    cliques->extract_clique(x, base, G, k);
      // Fill base[k] ... base[k-r+1] and base[0]
      induced->setup_nucleus(G, DG, k, base, r);
      // Need to fix so that k_idx is 1, but ends as if it was r
      NKCliqueDir_fast_hybrid_rec(G, 1, k-r, induced, update_d, base);
  }, granularity, false);
t1.stop();

  // Extract all vertices with changed clique counts
  //auto changed_vtxs = edge_table.entries();
  //edge_table.del();

  // Aggregate the updated counts across all worker's local arrays, as specified by update
  /*parallel_for(0, changed_vtxs.size(), [&] (size_t i) {
    size_t nthreads = num_workers();
    uintE v = std::get<0>(changed_vtxs[i]);
    auto index = cliques->find_index(v);
    for (size_t j=0; j<nthreads; j++) {
      update(per_processor_counts, j, index);
    }
  }, 128);*/

  // Perform update_changed on each vertex with changed clique counts
  std::size_t num_count_idxs = 0;
  if (do_update_changed) {
    /*parallel_for(0, changed_vtxs.size(), [&] (size_t i) {
      auto index = cliques->find_index(std::get<0>(changed_vtxs[i]));
      update_changed(per_processor_counts, i, index);
    });*/
      num_count_idxs = count_idxs.filter(update_changed, per_processor_counts);
      count_idxs.reset();
    /*
    parallel_for(0, num_count_idxs, [&] (size_t i) {//count_idxs[0]
      //assert(count_idxs[i+1] < n);
      //assert(per_processor_counts[count_idxs[i+1]] > 0);
      update_changed(per_processor_counts, i, count_idxs.pack[i]);//count_idxs[i + 1]
    });*/
    
  }

  // Mark every vertex in the active set as deleted
  parallel_for (0, active_size, [&] (size_t j) {
    auto index = get_active(j); //cliques->find_index(get_active(j));
    still_active[index] = 2;}, 2048);

  return num_count_idxs; //count_idxs[0];
}

template <typename bucket_t, class Graph, class Graph2, class T>
sequence<bucket_t> Peel(Graph& G, Graph2& DG, size_t r, size_t k, 
  T* cliques, sequence<uintE> &rank,
  size_t num_buckets=16) {
    k--; r--;
  timer t2; t2.start();
// Here's the mistake: You're thinking that get_active should return the key,
// which is the concatenation of vertices. That's wrong. get_active should
// return an index into cliques->table, which we can then convert to the
// appropriate key. When we bucket, we only maintain indices into
// cliques->table, and we transfer to the key as necessary.

  //sequence<std::tuple<__int128 unsigned, long int> > entries = cliques->entries;
  size_t num_entries = cliques->return_total();
  auto D = sequence<bucket_t>(num_entries, [&](size_t i) -> bucket_t { 
    return cliques->get_count(i);
  });

  auto D_filter = sequence<std::tuple<uintE, bucket_t>>(num_entries);

  auto b = make_vertex_custom_buckets<bucket_t>(num_entries, D, increasing, num_buckets);

  auto per_processor_counts = sequence<size_t>(num_entries , static_cast<size_t>(0));
  
  list_buffer count_idxs(num_entries);
  //auto count_idxs = sequence<size_t>(num_entries + 1 + 
  //  count_buffer * num_workers(), static_cast<size_t>(0));

  char* still_active = (char*) calloc(num_entries, sizeof(char));
  size_t max_deg = induced_hybrid::get_max_deg(G); // could instead do max_deg of active?
  //auto update_idxs = sequence<uintE>(max_deg);
  size_t n = num_entries;

  timer t_extract;
  timer t_count;
  timer t_update;
  timer t_x;

  size_t rounds = 0;
  size_t finished = 0;
  bucket_t cur_bkt = 0;
  bucket_t max_bkt = 0;
  double max_density = 0;
  bool use_max_density = false;
  // Peel each bucket
  auto update_clique = [&](sequence<size_t>& ppc, size_t j, uintE v) {
    if (j == 0) return;
    ppc[v] += ppc[j*n + v];
    ppc[j*n + v] = 0;
  };
  while (finished != n) {
    t_extract.start();
    // Retrieve next bucket
    auto bkt = b.next_bucket();
    //auto active = vertexSubset(n, bkt.identifiers);
    auto active_size = (bkt.identifiers).size();
    cur_bkt = bkt.id;
    t_extract.stop();

    finished += active_size;
    max_bkt = std::max(cur_bkt, max_bkt);

    auto get_active = [&](size_t j) -> unsigned __int128 { return (bkt.identifiers)[j]; };
      //return active.vtx(j); };

    if (active_size == 0 || D[get_active(0)] == 0) continue;
    //std::cout << "PEEL" << std::endl;
    //fflush(stdout);

    // TESTING EXTRACTION OF KEy
    /*
    for (int i = 0; i < static_cast<int>(active_size); i++){
    auto x = get_active(i);
    auto key = std::get<0>((cliques->table)[x]);
       for (int j = 0; j < static_cast<int>(r)+1; ++j) {
          int extract = (int) key;
          if (static_cast<uintE>(extract) >= DG.n) {
            std::cout << "Extract: " << static_cast<uintE>(extract) << ", n: " << DG.n << std::endl;
          }
          assert(static_cast<uintE>(extract) < DG.n);
        key = key >> 32;
      }
    }*/

    size_t granularity = (cur_bkt * active_size < 10000) ? 1024 : 1;

    size_t filter_size = 0;

      auto update_changed = [&](sequence<size_t>& ppc, size_t i, uintE v){
        /* Update the clique count for v, and zero out first worker's count */
        //auto index = cliques->find_index(v);
        //auto val = std::get<1>((cliques->table)[index]) - ppc[v];
    //(cliques->table)[index] = std::make_tuple(std::get<0>((cliques->table)[index]), val);
        //cliques[v] -= ppc[v];
        if (v == UINT_E_MAX) {
          D_filter[i] = std::make_tuple(num_entries + 1, 0);
          return;
        }
        if (ppc[v] == 0) D_filter[i] = std::make_tuple(num_entries + 1, 0);
        else {
          auto val = cliques->update_count(v, ppc[v]);
        ppc[v] = 0;
        bucket_t deg = D[v];
        if (deg > cur_bkt) {
          bucket_t new_deg = std::max((bucket_t) val, (bucket_t) cur_bkt);
          D[v] = new_deg;
          // store (v, bkt) in an array now, pass it to apply_f below instead of what's there right now -- maybe just store it in D_filter?
          D_filter[i] = std::make_tuple(v, b.get_bucket(deg, new_deg));
        } else D_filter[i] = std::make_tuple(num_entries + 1, 0);
        }
      };
    t_count.start();
     filter_size = cliqueUpdate(G, DG, r, k, max_deg, true, get_active, active_size, 
     granularity, still_active, rank, per_processor_counts,
      update_clique, true, update_changed, cliques, num_entries, count_idxs, t_x);
      t_count.stop();

    auto apply_f = [&](size_t i) -> std::optional<std::tuple<unsigned __int128, bucket_t>> {
      auto v = std::get<0>(D_filter[i]);
      bucket_t bucket = std::get<1>(D_filter[i]);
      if (v != num_entries + 1) {
        if (still_active[v] != 2) return wrap(v, bucket);
      }
      return std::nullopt;
    };

t_update.start();
    b.update_buckets(apply_f, filter_size);

    parallel_for (0, active_size, [&] (size_t j) {
      auto index = get_active(j);
      //auto index = cliques->find_index(v);
      cliques->clear_count(index);
      //cliques[active.vtx(j)] = 0;
      }, 2048);
      t_update.stop();

    //active.del();

    rounds++;
  }

  t_extract.reportTotal("### Peel Extract time: ");
  t_count.reportTotal("### Peel Count time: ");
  t_update.reportTotal("### Peel Update time: ");
  t_x.reportTotal("Inner counting: ");

  double tt2 = t2.stop();
  std::cout << "### Peel Running Time: " << tt2 << std::endl;

  std::cout.precision(17);
  std::cout << "rho: " << rounds << std::endl;
  std::cout << "clique core: " << static_cast<uintE>(max_bkt) << std::endl;
  if (use_max_density) std::cout << "max density: " << max_density << std::endl;

  b.del();
  free(still_active);

  return D;
}

template <class Graph>
inline size_t NucleusDecomposition(Graph& GA, size_t r, size_t s) {
  // TODO: if r = 2
  using W = typename Graph::weight_type;

  // Obtain vertex ordering
  timer t_rank; t_rank.start();
  sequence<uintE> rank = get_ordering(GA, 3, 0.1); // in clique counting
  double tt_rank = t_rank.stop();
  std::cout << "### Rank Running Time: " << tt_rank << std::endl;

  // Direct the graph based on ordering
  timer t_filter; t_filter.start();
  auto pack_predicate = [&](const uintE& u, const uintE& v, const W& wgh) {
    return (rank[u] < rank[v]) && GA.get_vertex(u).getOutDegree() >= r-1 && GA.get_vertex(v).getOutDegree() >= r-1;
  };
  auto DG = filterGraph(GA, pack_predicate);
  double tt_filter = t_filter.stop();
  std::cout << "### Filter Graph Running Time: " << tt_filter << std::endl;

  timer t; t.start();

  auto max_deg = get_max_deg3(DG);
  TwolevelHash table(r, DG, max_deg);

  std::cout << "End table + start count" << std::endl;
  size_t count = CountCliquesNuc(DG, s, r, max_deg, &table);
  std::cout << "End count" << std::endl;
  fflush(stdout);

  double tt = t.stop();
  std::cout << "### Count Running Time: " << tt << std::endl;
  std::cout << "### Num " << s << " cliques = " << count << "\n";

timer t2; t2.start();
  auto peel = Peel<std::size_t>(GA, DG, r, s, &table, rank);
  double tt2 = t2.stop();
  std::cout << "### Peel Running Time: " << tt2 << std::endl;

//table.del();
  DG.del();

  return count;
}

}