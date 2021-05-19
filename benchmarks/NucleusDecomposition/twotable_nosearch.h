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

// Clique files
#include "benchmarks/CliqueCounting/intersect.h"
#include "benchmarks/CliqueCounting/induced_intersection.h"
#include "benchmarks/CliqueCounting/induced_neighborhood.h"
#include "benchmarks/CliqueCounting/induced_hybrid.h"
#include "benchmarks/CliqueCounting/induced_split.h"
#include "benchmarks/CliqueCounting/relabel.h"

#include "commontable.h"
  
namespace gbbs {

namespace twotable_nosearch {

  template <class Y, class H>
  struct EndTable {
    pbbslib::sparse_table<Y, long, H> table;
    uintE vtx;
    //MidTable* up_table;
  };

  template <class Y, class H>
  struct MidTable {
    using EndTableY = EndTable<Y, H>;
    pbbslib::sparse_table<uintE, EndTableY*, std::hash<uintE>> table;
    sequence<EndTableY*> arr;
  };

  template<class Y, class S, class EndSpace>
  uintE get_mtable(S index, EndSpace* end_space) {
    using X = std::tuple<Y, long>;
    while (true) {
      auto max_val = std::get<0>(static_cast<X>(end_space[index]));
      std::size_t max_bit = sizeof(Y) * 8;
      Y one = 1;
      Y check_bit = (max_val >> (max_bit - 1)) & 1U;
      if (check_bit != 0) {
        max_val &= ~(one << (max_bit - 1));
        return static_cast<uintE>(max_val);
      }
      index++;
      //index = index % mtable.m;
    }
  }
  
  template <class Y, class H>
  class TwolevelHash {
    public:
      using T = pbbslib::sparse_table<Y, long, H>;
      using X = std::tuple<Y, long>;
      using EndTableY = EndTable<Y, H>;
      using MidTableY = MidTable<Y, H>;
      MidTableY top_table;
      sequence<long> top_table_sizes;
      int rr;
      std::size_t total = 0;
      X* space = nullptr;
      int shift_factor;
  
      template <class Graph>
      TwolevelHash(int r, Graph& DG, size_t max_deg, bool relabel, int _shift_factor) {
        shift_factor = _shift_factor;
        rr = r;
        //top_table.up_table = nullptr;
        // How many vert in top level?
        // For each vert in top level, how many pairs of vert finish it?
        auto tmp_table = pbbslib::sparse_additive_map<uintE, long>(
          DG.n, std::make_tuple(UINT_E_MAX, long{0}));
        auto base_f = [&](sequence<uintE>& base){
          auto min_vert = relabel ? base[0] : pbbslib::reduce_min(base);
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
        sequence<long> actual_sizes(top_table_sizes2.size() + 1);
        // Modify top_table_sizes2 to be appropriately oversized
        parallel_for(0, top_table_sizes2.size(), [&](std::size_t i) {
          auto m = 1 + ((size_t)1 << pbbslib::log2_up((size_t)(1.1 * std::get<1>(top_table_sizes2[i])) + 1));
          actual_sizes[i] = m;
        });
        actual_sizes[top_table_sizes2.size()] = 0;
        // Do a scan inplace
        //auto add_tuple_monoid = pbbs::make_monoid([](std::tuple<uintE, long> a, std::tuple<uintE, long> b){
        //  return std::make_tuple(std::get<0>(b), std::get<1>(a) + std::get<1>(b));
        //}, std::make_tuple(UINT_E_MAX, 0));
        long total_top_table_sizes2 = scan_inplace(actual_sizes.slice(), pbbs::addm<long>());
        // Allocate space for the second level tables
        space = pbbslib::new_array_no_init<X>(total_top_table_sizes2);
        tmp_table.del();
  
        //*** for arr
        top_table.arr = sequence<EndTableY*>(DG.n, [](std::size_t i){return nullptr;});
        top_table_sizes = sequence<long>(DG.n + 1, long{0});
        /*top_table.table = pbbslib::sparse_table<uintE, EndTable*, std::hash<uintE>>(
          top_table_sizes2.size(),
          std::make_tuple<uintE, EndTable*>(UINT_E_MAX, static_cast<EndTable*>(nullptr)),
          std::hash<uintE>());
        top_table_sizes = sequence<long>(top_table.table.m + 1, long{0});*/
  
        parallel_for(0, top_table_sizes2.size(), [&](std::size_t i){
          auto vtx = std::get<0>(top_table_sizes2[i]);
          auto upper_size = actual_sizes[i + 1];
          auto size = upper_size - actual_sizes[i];
          EndTableY* end_table = new EndTableY();

          Y max_val = static_cast<Y>(vtx); 
          std::size_t max_bit = sizeof(Y) * 8;
          Y one = 1;
          max_val |= one << (max_bit - 1);

          end_table->vtx = vtx;
          end_table->table = pbbslib::sparse_table<Y, long, H>(
            size - 1, 
            std::make_tuple<Y, long>(static_cast<Y>(max_val), static_cast<long>(0)),
            H{},
            space + actual_sizes[i]
            );
          space[size - 1] = std::make_tuple<Y, long>(static_cast<Y>(max_val), static_cast<long>(0));
          /*top_table.table.insert(std::make_tuple(vtx, end_table));
          std::size_t l = top_table.table.find_index(vtx);
          top_table_sizes[l] = end_table->table.m;*/
          //***for arr
          top_table.arr[vtx] = end_table;
          top_table_sizes[vtx] = end_table->table.m;
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
  
      void insert(sequence<uintE>& base2, int r, int k) {
        auto add_f = [&] (long* ct, const std::tuple<Y, long>& tup) {
          pbbs::fetch_and_add(ct, (long)1);
        };
        // Sort base
        uintE base[10];
        assert(10 > k);
        for(std::size_t i = 0; i < k + 1; i++) {
          base[i] = base2[i];
        }
        std::sort(base, base + k + 1,std::less<uintE>());

        std::string bitmask(r+1, 1); // K leading 1's
        bitmask.resize(k+1, 0); // N-K trailing 0's

        do {
          bool use_vtx = false;
          uintE vtx = 0;
          Y key = 0;
          for (int i = 0; i < static_cast<int>(k)+1; ++i) {
            if (bitmask[i]) {
              if (!use_vtx) {
                use_vtx = true;
                vtx = base[i];
              } else {
                key = key << shift_factor;
                key |= static_cast<int>(base[i]);
              }
            }
          }
          //EndTable* end_table = top_table.table.find(vtx, nullptr);
          //***for arr
          EndTableY* end_table = top_table.arr[vtx];
          //assert(end_table != nullptr);
          (end_table->table).insert_f(std::make_tuple(key, (long) 1), add_f);
        } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
      }

      std::size_t return_total() {return total;}

      size_t get_top_index(std::size_t index) {
        // This gives the first i such that top_table_sizes[i] >= index
        auto idx = pbbslib::binary_search(top_table_sizes, long{index}, std::less<long>());
        if (idx >= top_table_sizes.size()) return top_table_sizes.size() - 1;
        if (top_table_sizes[idx] == index) {
          while(idx < top_table_sizes.size() && top_table_sizes[idx] == index) {
            idx++;
          }
          return idx-1;
        }
        if (idx == 0) return idx;
        return idx - 1;
      }

      long get_count(std::size_t index) {
        return std::get<1>(space[index]);
      }

      size_t update_count(std::size_t index, size_t update){
        auto val = std::get<1>(space[index]) - update;
        space[index] =
          std::make_tuple(std::get<0>(space[index]), val);
        return val;
      }

      void clear_count(std::size_t index) {
        space[index] = std::make_tuple(std::get<0>(space[index]), 0);
      }

      template<class I>
      void extract_indices(sequence<uintE>& base2, I func, int r, int k) {
        // Sort base
        uintE base[10];
        assert(10 > k);
        for(std::size_t i = 0; i < k + 1; i++) {
          base[i] = base2[i];
        }
        std::sort(base, base + k + 1,std::less<uintE>());

        std::string bitmask(r+1, 1); // K leading 1's
        bitmask.resize(k+1, 0); // N-K trailing 0's
        do {
          bool use_vtx = false;
          uintE vtx = 0;
          Y key = 0;
          for (int i = 0; i < static_cast<int>(k)+1; ++i) {
            if (bitmask[i]) {
              if (!use_vtx) {
                use_vtx = true;
                vtx = base[i];
              } else {
                key = key << shift_factor;
                key |= static_cast<int>(base[i]);
              }
            }
          }
          // First, find index in top_table
          // This should populate into a prefix sum of sizes
          /*auto top_index = top_table.table.find_index(vtx);
          auto prefix = top_table_sizes[top_index];
          EndTable* end_table = std::get<1>(top_table.table.table[top_index]);*/
          //***for arr
          EndTableY* end_table = top_table.arr[vtx];
          auto prefix = top_table_sizes[vtx];
          auto index = (end_table->table).find_index(key);
          func(prefix + index);
        } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
      }

      // Given an index, get the clique
      template<class S, class Graph>
      void extract_clique(S index, sequence<uintE>& base, Graph& G, int k) {
        Y vert;
        uintE v = get_mtable<Y>(index, space);
        base[0] = v;
        vert = std::get<0>(space[index]);
        for (int j = 0; j < rr - 1; ++j) {
          int mask = (1UL << shift_factor) - 1;
          int extract = (int) vert & mask; // vert & mask
          /*if (static_cast<uintE>(extract) >= G.n) {
            std::cout << "Vert: " << static_cast<uintE>(extract) << ", n: " << G.n << std::endl;
          }*/
          assert(static_cast<uintE>(extract) < G.n);
          base[k - j] = static_cast<uintE>(extract);
          vert = vert >> shift_factor;
        }
      }
  };

} // end namespace twotable

} //end namespace gbbs