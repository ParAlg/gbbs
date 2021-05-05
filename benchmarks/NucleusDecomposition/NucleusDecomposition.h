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

  struct MultiSizeTable {
    using NextSizeTable = sequence<MultiSizeTable*>;
    // TODO: This could be optimized, because only the level == max_level - 2
    // needs sizetable to be an additive map (otherwise, we just want the # of 
    // distinct entries on that level)
    using SizeTable = pbbslib::sparse_additive_map<uintE, long>;
    using NextSizeTableSE = pbbslib::sparse_table<uintE, MultiSizeTable*, std::hash<uintE>>;
    using EndTable = pbbslib::sparse_table<unsigned __int128, long, hash128>;
    SizeTable size_table;
    NextSizeTable next_size_table;
    NextSizeTableSE next_size_table_se;
    uintE level;
    uintE max_level;
    bool is_hash;
    bool use_se = false;

    sequence<long> table_sizes;
    uintE vtx;
    sequence<EndTable*> end_tables;

    MultiSizeTable* prev_size_table = nullptr;
    std::tuple<unsigned __int128, long>* end_space;
    long total_size;

    void constructFinalTables(std::tuple<unsigned __int128, long>* space) {
      if (level == max_level - 2) {
        end_tables = sequence<EndTable*>(table_sizes.size() - 1, [](std::size_t i){return nullptr;});
        parallel_for(0, table_sizes.size() - 1, [&](std::size_t i){
          end_tables[i] = new EndTable(table_sizes[i+1] - table_sizes[i],
            std::make_tuple<unsigned __int128, long>(static_cast<unsigned __int128>(0), static_cast<long>(0)),
            hash128{},
            space + table_sizes[i]);
        });
        return;
      }

      if (!use_se) {
        table_sizes = sequence<long>(next_size_table.size(), [&](std::size_t i){ return 0; });
        for (std::size_t i = 0; i < next_size_table.size(); i++) {
          // Here pass only the space for this section
          next_size_table[i]->constructFinalTables(space + table_sizes[i]);
        }
      } else {
        table_sizes = sequence<long>(next_size_table_se.m, [&](std::size_t i){ return 0; });
        for (std::size_t i = 0; i < next_size_table_se.m; i++) {
          if (std::get<0>(next_size_table_se.table[i]) != next_size_table_se.empty_key) {
            std::get<1>(next_size_table_se.table[i])->constructFinalTables(space + table_sizes[i]);
          }
        }
      }

    }

    void constructFullTable() {
      total_size = computeSizes();

      // Allocate space for the last level tables
      using X = std::tuple<unsigned __int128, long>;
      end_space = pbbslib::new_array_no_init<X>(total_size);

      constructFinalTables(end_space);
    }

    long computeSizes() {
      if (level == max_level - 2) {
        auto size_table_entries = size_table.entries();
        table_sizes = sequence<long>(size_table_entries.size() + 1, [](std::size_t i){ return 0; });
        parallel_for(0, size_table_entries.size(), [&](std::size_t i) {
          auto m = (size_t)1 << pbbslib::log2_up((size_t)(1.1 * std::get<1>(size_table_entries[i])) + 1);
          vtx = std::get<0>(size_table_entries[i]);
          table_sizes[i] = m;
        });
        // Do a scan inplace
        auto add_tuple_monoid = pbbs::make_monoid([](long a, long b){
          return a + b;
        }, long{0});
        long total_size_table_entries = scan_inplace(table_sizes.slice(), add_tuple_monoid);
        return total_size_table_entries;
      }
      if (!use_se) {
        table_sizes = sequence<long>(next_size_table.size(), [&](std::size_t i){ return 0; });
        for (std::size_t i = 0; i < next_size_table.size(); i++) {
          table_sizes[i] = next_size_table[i]->computeSizes();
        }
      } else {
        table_sizes = sequence<long>(next_size_table_se.m, [&](std::size_t i){ return 0; });
        for (std::size_t i = 0; i < next_size_table_se.m; i++) {
          if (std::get<0>(next_size_table_se.table[i]) != next_size_table_se.empty_key) {
            table_sizes[i] = std::get<1>(next_size_table_se.table[i])->computeSizes();
          }
        }
      }
      auto add_tuple_monoid = pbbs::make_monoid([](long a, long b){
        return a + b;
      }, long{0});
      long total_size_table_entries = scan_inplace(table_sizes.slice(), add_tuple_monoid);
      return total_size_table_entries;
    }

    // Doing non space efficient is easy because we just call this, and then run
    // r clique counting once, calling insertSize
    void createMultiSizeTable(uintE _level, uintE _max_level, uintE n, sequence<bool>& is_array_level,
      MultiSizeTable* _prev_size_table = nullptr) {
      prev_size_table = _prev_size_table;
      level = _level;
      max_level = _max_level;
      is_hash = !is_array_level[level] || level == max_level - 2;
      if (is_hash) size_table = SizeTable(n, std::make_tuple(UINT_E_MAX, long{0}));
      if (level == max_level - 2) return;
      next_size_table = NextSizeTable(n, [](std::size_t i){
        return new MultiSizeTable();
      });
      // TODO: This is space inefficient because we may not need size n, but I think
      // it's computationally hard to do better (practically hard?)
      // We need to do something like, for each hash level, re-do the level-clique counting,
      // then iteratively start constructing the full multi level hash table, so that
      // we don't use too much space
      for (std::size_t i = 0; i < n; i++) {
        next_size_table[i]->createMultiSizeTable(level + 1, max_level, n, is_array_level, this);
      }
    }

    void insertSize(uintE* base, int r) {
      if (level >= r) return;
      if (is_hash) {
        auto entry = std::make_tuple<uintE, long>(static_cast<uintE>(base[level]), long{1});
        size_table.insert(entry);
      }
      if (level == max_level - 2) return;
      next_size_table[base[level]]->insertSize(base, r);
    }

    // Zeroth level should always be an array level
    // The idea here is for every level, if it's a hash level, run
    // insertCurrSizeSE with r = level + 1
    // Then run createNextMultiSizeTableSE if it's a hash level, and else
    // run createNextMultiSizeTable
    // Terminate after creating and inserting for level = max_level - 2
    // TODO: modify clique counting to work for k = 2
    void createMultiSizeTableSE(uintE _level, uintE _max_level, uintE n, sequence<bool>& is_array_level,
      MultiSizeTable* _prev_size_table = nullptr) {
      prev_size_table = _prev_size_table;
      level = _level;
      max_level = _max_level;
      is_hash = !is_array_level[level] || level == max_level - 2;
      use_se = is_hash;
      if (is_hash) size_table = SizeTable(n, std::make_tuple(UINT_E_MAX, long{0}));
    }
    void insertCurrSizeSE(uintE* base, int r) {
      if (level >= r) return;
      if (is_hash && level + 1 == r) {
        //std::cout << "Flag 1" << std::endl; fflush(stdout);
        auto entry = std::make_tuple<uintE, long>(static_cast<uintE>(base[level]), long{1});
        size_table.insert(entry);
        return;
      }
      //std::cout << "Flag 2" << std::endl; fflush(stdout);
      // TODO: figure out why this is not working
      //if (!use_se) {
      //  next_size_table[base[level]]->insertCurrSizeSE(base, r);
      //  return;
      //}
      MultiSizeTable* next = next_size_table_se.find(base[level], static_cast<MultiSizeTable*>(nullptr));
      assert(next != nullptr);
      next->insertCurrSizeSE(base, r);
    }
    void createNextMultiSizeTable(uintE n, sequence<bool>& is_array_level) {
      next_size_table = NextSizeTable(n, [](std::size_t i){
        return new MultiSizeTable();
      });
      for (std::size_t i = 0; i < n; i++) {
        next_size_table[i]->createMultiSizeTableSE(level + 1, max_level, n, is_array_level, this);
      }
    }
    void createNextMultiSizeTableSE(uintE nn, sequence<bool>& is_array_level) {
      auto entries = size_table.entries();
      uintE n = entries.size();
      use_se = true;
      next_size_table_se = NextSizeTableSE(
          n,
          std::make_tuple<uintE, MultiSizeTable*>(UINT_E_MAX, static_cast<MultiSizeTable*>(nullptr)),
          std::hash<uintE>()
      );
      for (std::size_t i = 0; i < n; i++) {
        auto entry = entries[i];
        MultiSizeTable* table = new MultiSizeTable();
        table->createMultiSizeTableSE(level + 1, max_level, nn, is_array_level, this);
        next_size_table_se.insert(std::make_tuple(std::get<0>(entry), table));
      }
    }

    void insert(uintE* base, int curr_idx, int r, int k, std::string& bitmask) {
      if (level == max_level - 2) {
        auto add_f = [&] (long* ct, const std::tuple<unsigned __int128, long>& tup) {
          pbbs::fetch_and_add(ct, (long)1);
        };
        unsigned __int128 key = 0;
        for (int i = curr_idx + 1; i < static_cast<int>(k)+1; ++i) {
          if (bitmask[i]) {
            key = key << 32;
            key |= static_cast<int>(base[i]);
          }
        }
        auto end_table = end_tables[base[curr_idx]];
        end_table->insert_f(std::make_tuple(key, (long) 1), add_f);
        return;
      }
      int next_idx = curr_idx;
      for (int i = curr_idx + 1; i < static_cast<int>(k)+1; ++i) {
        if (bitmask[i]) {
          next_idx = i;
          break;
        }
      }
      if (!use_se) next_size_table[base[curr_idx]]->insert(base, next_idx, r, k, bitmask);
      else {
        auto next = next_size_table_se.find(base[curr_idx], nullptr);
        next->insert(base, next_idx, r, k, bitmask);
      }
    }

    void insert(sequence<uintE>& base, int r, int k) {
        // Sort base
        uintE base2[10];
        assert(10 > k);
        for(std::size_t i = 0; i < k + 1; i++) {
          base2[i] = base[i];
        }
        std::sort(base2, base2 + k + 1,std::less<uintE>());

        std::string bitmask(r+1, 1); // K leading 1's
        bitmask.resize(k+1, 0); // N-K trailing 0's

        do {
          unsigned __int128 key = 0;
          for (int i = 0; i < static_cast<int>(k)+1; ++i) {
            if (bitmask[i]) {
              insert(base2, i, r, k, bitmask);
            }
          }
        } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
    }

    std::tuple<unsigned __int128, long>* extract_indices(uintE* base, int curr_idx,
      int r, int k, std::string& bitmask) {
      if (level == max_level - 2) {
        unsigned __int128 key = 0;
        for (int i = curr_idx + 1; i < static_cast<int>(k)+1; ++i) {
          if (bitmask[i]) {
            key = key << 32;
            key |= static_cast<int>(base[i]);
          }
        }
        auto end_table = end_tables[base[curr_idx]];
        return end_table->table + end_table->find_index(key);
      }
      int next_idx = curr_idx;
      for (int i = curr_idx + 1; i < static_cast<int>(k)+1; ++i) {
        if (bitmask[i]) {
          next_idx = i;
          break;
        }
      }
      if (!use_se) return next_size_table[base[curr_idx]]->extract_indices(base, next_idx, r, k, bitmask);
      else {
        auto next = next_size_table_se.find(base[curr_idx], nullptr);
        return next->extract_indices(base, next_idx, r, k, bitmask);
      }
    }

    template<class I>
    void extract_indices(sequence<uintE>& base, I func, int r, int k) {
        uintE base2[10];
        assert(10 > k);
        for(std::size_t i = 0; i < k + 1; i++) {
          base2[i] = base[i];
        }
        std::sort(base2, base2 + k + 1,std::less<uintE>());

        std::string bitmask(r+1, 1); // K leading 1's
        bitmask.resize(k+1, 0); // N-K trailing 0's
        do {
          unsigned __int128 key = 0;
          for (int i = 0; i < static_cast<int>(k)+1; ++i) {
            if (bitmask[i]) {
              // TODO: make sure this pointer arithmetic is fine
              auto idx = extract_indices(base2, i, r, k, bitmask) - end_space;
              func(idx);
            }
          }
        } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
      }
  
    void extract_clique(sequence<uintE>& base, int base_idx) {
      base[base_idx] = vtx;
      if (prev_size_table != nullptr) prev_size_table->extract_clique(base, base_idx + 1);
    }

    template<class S, class Graph>
    void extract_clique(S index, sequence<uintE>& base, Graph& G, int k, int rr) {
      auto vert = std::get<0>(end_space[index]);
      // TOOD: make sure this calc is correct
      for (int j = 0; j < rr - max_level + 1; ++j) {
        int extract = (int) vert;
        assert(static_cast<uintE>(extract) < G.n);
        base[j] = static_cast<uintE>(extract);
        vert = vert >> 32;
      }
      base[rr - max_level + 1] = vtx;
      if (prev_size_table != nullptr) prev_size_table->extract_clique(base, rr - max_level + 2);
    }
  };

  template <class EndKey, class HashKey>
  class MultilevelHash {
    public:
      //FullMultiTable<EndKey, HashKey> table;
      int rr = 0;
      uintE max_level;
      MultiSizeTable multisize_table;

      template<class Graph>
      MultilevelHash(int r, Graph& DG, size_t max_deg, uintE _max_level, sequence<bool> is_array_level){
        rr = r;
        max_level = _max_level;
        using W = typename Graph::weight_type;
        // First, compute sizes for every level 0 to max_level - 1 (assume max_level is a size, e.g., twolevel sets max_level = 2)
        // We only need sizes on non-array levels, indexed by the key given to the previous level
        // First, create a MultiSizeTable
        //MultiSizeTable multisize_table;
        // TODO: If the last level is not a hash, max_level here can be the
        // last index of a false in is_array_level
        std::cout << "FLAG: start create" << std::endl; fflush(stdout);
        multisize_table.createMultiSizeTableSE(0, max_level, DG.n, is_array_level);
        std::cout << "FLAG: end create" << std::endl; fflush(stdout);
        for (std::size_t i = 0; i < max_level - 1; i++) {
          auto tmp_r = i + 1;
          if (tmp_r == 1) {
            std::cout << "FLAG: start tmp_r = 2" << std::endl; fflush(stdout);
            parallel_for(0, DG.n, [&](std::size_t j) {
              uintE base[2];
              auto map_label_f = [&] (const uintE& u, const uintE& v, const W& wgh) {
                base[0] = std::min(static_cast<uintE>(j), v);
                base[1] = std::max(static_cast<uintE>(j), v);
                multisize_table.insertCurrSizeSE(base, tmp_r);
              };
              DG.get_vertex(j).mapOutNgh(j, map_label_f, false);
            });
            std::cout << "FLAG: end tmp_r = 2" << std::endl; fflush(stdout);
            std::cout << "FLAG: start next table" << std::endl; fflush(stdout);
            if (is_array_level[i]) multisize_table.createNextMultiSizeTable(DG.n, is_array_level);
            else multisize_table.createNextMultiSizeTableSE(DG.n, is_array_level);
            std::cout << "FLAG: end next table" << std::endl; fflush(stdout);
            continue;
          }
          std::cout << "FLAG: start tmp_r = " << tmp_r << std::endl; fflush(stdout);
          auto base_f = [&](sequence<uintE>& base){
            uintE base2[10];
            for(std::size_t j = 0; j < rr; j++) { base2[j] = base[j]; }
            std::sort(base2, base2 + rr,std::less<uintE>());
            multisize_table.insertCurrSizeSE(base2, tmp_r);
          };
          auto init_induced = [&](HybridSpace_lw* induced) { induced->alloc(max_deg, tmp_r, DG.n, true, true); };
          auto finish_induced = [&](HybridSpace_lw* induced) { if (induced != nullptr) { delete induced; } };
          parallel_for_alloc<HybridSpace_lw>(init_induced, finish_induced, 0, DG.n, [&](size_t i, HybridSpace_lw* induced) {
            if (DG.get_vertex(i).getOutDegree() != 0) {
              induced->setup(DG, tmp_r, i);
              auto base = sequence<uintE>(r);
              base[0] = i;
              NKCliqueDir_fast_hybrid_rec(DG, 1, tmp_r, induced, base_f, base);
            }
          }, 1, false);
          std::cout << "FLAG: end tmp_r" << std::endl; fflush(stdout);
          std::cout << "FLAG: start next table" << std::endl; fflush(stdout);
          if (is_array_level[i]) multisize_table.createNextMultiSizeTable(DG.n, is_array_level);
          else multisize_table.createNextMultiSizeTableSE(DG.n, is_array_level);
          std::cout << "FLAG: end next table" << std::endl; fflush(stdout);
        }
        std::cout << "FLAG: end create loop" << std::endl; fflush(stdout);
        // Now that the last layer has a hash table with the correct vertex counts,
        // we must do prefix sums to figure out indices
        // And allocate space
        std::cout << "FLAG: start create full" << std::endl; fflush(stdout);
        multisize_table.constructFullTable();
        std::cout << "FLAG: end create full" << std::endl; fflush(stdout);

      }

      void insert(sequence<uintE>& base, int r, int k) {
        multisize_table.insert(base, r, k);
      }

      std::size_t return_total() { return multisize_table.total_size; }

      long get_count(std::size_t index) {
        return std::get<1>(multisize_table.end_space[index]);
      }

      size_t update_count(std::size_t index, size_t update){
        auto val = std::get<1>(multisize_table.end_space[index]) - update;
        multisize_table.end_space[index] =
          std::make_tuple(std::get<0>(multisize_table.end_space[index]), val);
        return val;
      }

      void clear_count(std::size_t index) {
        multisize_table.end_space[index] =
          std::make_tuple(std::get<0>(multisize_table.end_space[index]), 0);
      }

      template<class I>
      void extract_indices(sequence<uintE>& base, I func, int r, int k) {
        multisize_table.extract_indices(base, func, r, k);
      }

      template<class S, class Graph>
      void extract_clique(S index, sequence<uintE>& base, Graph& G, int k) {
        multisize_table.extract_clique(index, base, G, k, rr);
      }
  };

  // max_lvl should be set to # levels - 2
  // two level hash is equiv to setting max_level to 0
  struct MTable {
    using NextMTable = pbbslib::sparse_table<uintE, MTable*, std::hash<uintE>>;
    using EndTable = pbbslib::sparse_table<unsigned __int128, long, hash128>;
    NextMTable mtable;
    EndTable end_table;
    uintE vtx;
    uintE lvl;
    uintE max_lvl;
    MTable* prev_mtable = nullptr;
    long total_size = 0;
    sequence<long> table_sizes;
    std::tuple<unsigned __int128, long>* end_space = nullptr;
    bool set_table_size_flag = false;

    void set_end_table_rec(std::tuple<unsigned __int128, long>* _end_space) {
      end_space = _end_space;

      if (lvl == max_lvl) {
        // Allocate end_table here using end_space
        end_table = EndTable(
          total_size,
          std::make_tuple<unsigned __int128, long>(std::numeric_limits<unsigned __int128>::max(), static_cast<long>(0)),
          hash128{},
          end_space
        );

      } else {
        for (std::size_t i = 0; i < mtable.m; i++) {
          if (std::get<0>(mtable.table[i]) != UINT_E_MAX) {
            std::get<1>(mtable.table[i])->set_end_table_rec(end_space + table_sizes[i]);
          }
        }
      }
    }

    void increment_size() {
      assert(lvl == max_lvl);
      total_size++;
    }

    long set_table_sizes() {
      if (lvl != max_lvl) {
        if (lvl + 1 == max_lvl) {
          table_sizes = sequence<long>(mtable.m, [](std::size_t i){ return 0; });
          parallel_for(0, mtable.m, [&](std::size_t i){
            if (std::get<0>(mtable.table[i]) != UINT_E_MAX) {
              auto tbl = std::get<1>(mtable.table[i]);
              tbl->total_size = (size_t)1 << pbbslib::log2_up((size_t)(1.5 * tbl->total_size) + 1);
              table_sizes[i] = tbl->total_size;
            }
          });
          total_size = scan_inplace(table_sizes.slice(), pbbs::addm<long>());
          return total_size;
        }
        table_sizes = sequence<long>(mtable.m, [&](std::size_t i){
          if (std::get<0>(mtable.table[i]) == UINT_E_MAX) return long{0};
          return std::get<1>(mtable.table[i])->total_size;
        });
        total_size = scan_inplace(table_sizes.slice(), pbbs::addm<long>());
        return total_size;
      }
      //else if (set_table_size_flag == false) {
      //  set_table_size_flag = true;
      //  total_size = (size_t)1 << pbbslib::log2_up((size_t)(1.5 * total_size) + 1);
      //  return total_size;
      //}
      return total_size;
    }

    void initialize(uintE _v, uintE _lvl, uintE _max_lvl, MTable* _prev_mtable = nullptr) {
      vtx = _v;
      lvl = _lvl;
      max_lvl = _max_lvl;
      prev_mtable = _prev_mtable;
      total_size = 0;
    }

    // k_idx goes from 1 to k - 1 (inclusive)
    void allocate(uintE count, size_t k_idx, size_t k) {
      //lvl = k_idx - 1;
      if (lvl != max_lvl) {
        mtable = NextMTable(
          count * 1.5,
          std::make_tuple<uintE, MTable*>(UINT_E_MAX, static_cast<MTable*>(nullptr)),
          std::hash<uintE>()
        );
      }
    }
    
    MTable* next(uintE v, size_t k_idx, size_t k) {
      if (lvl == max_lvl) return nullptr;
      MTable* next_table = new MTable();
      next_table->initialize(v, lvl + 1, max_lvl, this);
      mtable.insert(std::make_tuple(v, next_table));
      return next_table;
    }

    void insert(sequence<uintE>& base, int curr_idx, int r, int k, std::string& bitmask) {
      if (lvl == max_lvl) {
        assert(end_table.m > 0);
        auto add_f = [&] (long* ct, const std::tuple<unsigned __int128, long>& tup) {
          pbbs::fetch_and_add(ct, (long)1);
        };
        unsigned __int128 key = 0;
        for (int i = curr_idx; i < static_cast<int>(k)+1; ++i) {
          if (bitmask[i]) {
            key = key << 32;
            key |= static_cast<int>(base[i]);
          }
        }
        end_table.insert_f(std::make_tuple(key, (long) 1), add_f);
        return;
      }
      int next_idx = curr_idx;
      for (int i = curr_idx + 1; i < static_cast<int>(k)+1; ++i) {
        if (bitmask[i]) {
          next_idx = i;
          break;
        }
      }
      auto next = mtable.find(base[curr_idx], nullptr);
      assert(next != nullptr);
      next->insert(base, next_idx, r, k, bitmask);
    }

    void insert(sequence<uintE>& base, int r, int k) {
        std::string bitmask(r+1, 1); // K leading 1's
        bitmask.resize(k+1, 0); // N-K trailing 0's

        do {
          unsigned __int128 key = 0;
          for (int i = 0; i < static_cast<int>(k)+1; ++i) {
            if (bitmask[i]) {
              insert(base, i, r, k, bitmask);
              break;
            }
          }
        } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
    }

//std::tuple<unsigned __int128, long>*
    size_t extract_indices(sequence<uintE>& base, int curr_idx,
      int r, int k, std::string& bitmask) {
      if (lvl == max_lvl) {
        unsigned __int128 key = 0;
        for (int i = curr_idx; i < static_cast<int>(k)+1; ++i) {
          if (bitmask[i]) {
            key = key << 32;
            key |= static_cast<int>(base[i]);
          }
        }
        return end_table.find_index(key); // end_table.table + 
      }
      int next_idx = curr_idx;
      for (int i = curr_idx + 1; i < static_cast<int>(k)+1; ++i) {
        if (bitmask[i]) {
          next_idx = i;
          break;
        }
      }
      auto next_mtable_index = mtable.find_index(base[curr_idx]);
      auto next = std::get<1>(mtable.table[next_mtable_index]);
      assert(next != nullptr);
      return table_sizes[next_mtable_index] + next->extract_indices(base, next_idx, r, k, bitmask);
    }

    template<class I>
    void extract_indices(sequence<uintE>& base, I func, int r, int k) {
      std::string bitmask(r+1, 1); // K leading 1's
      bitmask.resize(k+1, 0); // N-K trailing 0's
      do {
        //unsigned __int128 key = 0;
        for (int i = 0; i < static_cast<int>(k)+1; ++i) {
          if (bitmask[i]) {
            // TODO: make sure this pointer arithmetic is fine
            size_t idx = extract_indices(base, i, r, k, bitmask); // - end_space;
            func(idx);
            break;
          }
        }
      } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
    }

    template <class S>
    long get_top_index(S index) {
      // This gives the first i such that top_table_sizes[i] >= index
      auto idx = pbbslib::binary_search(table_sizes.slice(), long{index}, std::less<long>());
      if (idx >= table_sizes.size()) return table_sizes.size() - 1;
      if (idx == 0) return idx;
      if (table_sizes[idx] == index) {
        while(table_sizes[idx] == index) {
          idx--;
        }
        return idx++;
      }
      return idx - 1;
    }
  
    //Fill base[k] ... base[k-r+1] and base[0]
    template<class S>
    void extract_clique(S index, sequence<uintE>& base, int base_idx, int rr, int k) {
      if (lvl == max_lvl) {
        if (lvl != 0) {
          base[base_idx] = vtx;
          if (base_idx == 0) base_idx = k - rr + 1;
          else base_idx++;
        }
        assert(end_space != nullptr);
        auto vert = std::get<0>(end_space[index]);
        // TOOD: make sure this calc is correct
        for (int j = k; j >= base_idx; --j) { //rr - 1, base_idx
          int extract = (int) vert;
          //assert(static_cast<uintE>(extract) < G.n);
          base[j] = static_cast<uintE>(extract);
          vert = vert >> 32;
        }
        return;
      }
      auto next_mtable_idx = get_top_index(index);
      if (next_mtable_idx >= mtable.m) {
        std::cout << "Idx: " << next_mtable_idx << std::endl;
        std::cout << "m: " << mtable.m << std::endl;
        fflush(stdout);
      }
      assert(next_mtable_idx < mtable.m);
      assert(index >= table_sizes[next_mtable_idx]);
      S next_index = index - table_sizes[next_mtable_idx];
      if (lvl != 0) {
        base[base_idx] = vtx;
        if (base_idx == 0) base_idx = k - rr + 1;
        else base_idx++;
      }
      if (std::get<1>(mtable.table[next_mtable_idx]) == nullptr) {
        std::cout << "rr: " << rr << std::endl; fflush(stdout);
        std::cout << "base_idx: " << base_idx << std::endl; fflush(stdout);
        std::cout << "index: " << long{index} << std::endl; fflush(stdout);
        std::cout << "top index: " << long{next_mtable_idx} << std::endl; fflush(stdout);
        std::cout << "size: " << table_sizes[next_mtable_idx] << std::endl; fflush(stdout);
      }
      assert(std::get<1>(mtable.table[next_mtable_idx]) != nullptr);
      std::get<1>(mtable.table[next_mtable_idx])->extract_clique(next_index, base, base_idx, rr, k);
    }

  };

  template <class Graph, class Space>
  inline size_t NKCliqueDir_fast_hybrid_rec_multi(Graph& DG, size_t k_idx, size_t k,
  HybridSpace_lw* induced, Space* space) {
    size_t num_induced = induced->num_induced[k_idx-1];
    if (num_induced == 0) return 0;
    uintE* prev_induced = induced->induced + induced->nn * (k_idx - 1);

    for (size_t i=0; i < num_induced; i++) { induced->labels[prev_induced[i]] = k_idx; }

    if (k_idx + 1 == k) {
      size_t counts = 0;
      space->allocate(num_induced, k_idx, k);
      for (size_t i=0; i < num_induced; i++) {
        uintE vtx = prev_induced[i];
        //  get neighbors of vtx
        uintE* intersect = induced->induced_edges + vtx * induced->nn;
        size_t tmp_counts = 0;
        //base[k_idx] = induced->relabel[vtx];
        Space* next_space = space->next(induced->relabel[vtx], k_idx, k);
        if (next_space == nullptr) next_space = space;
        for (size_t j=0; j < induced->induced_degs[vtx]; j++) {
          if (static_cast<size_t>(induced->labels[intersect[j]]) == k_idx) {
            //base[k] = induced->relabel[intersect[j]];
            tmp_counts++;
            //base_f(base, space);
            // TODO: This could be improved by collating on count before
            // returning to recursive level of max_lvl
            next_space->increment_size();
          }
        }
        next_space->set_table_sizes();
        counts += tmp_counts;
      }
      for (size_t i=0; i < num_induced; i++) { induced->labels[prev_induced[i]] = k_idx - 1; }
      space->set_table_sizes();
      return counts;
    }

    size_t total_ct = 0;
    space->allocate(num_induced, k_idx, k);
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
        //base[k_idx] = induced->relabel[vtx];
        Space* next_space = space->next(induced->relabel[vtx], k_idx, k);
        if (next_space == nullptr) next_space = space;
        auto curr_counts = NKCliqueDir_fast_hybrid_rec_multi(DG, k_idx + 1, k, induced, next_space);
        total_ct += curr_counts;
      }
    }
    for (size_t i=0; i < num_induced; i++) { induced->labels[prev_induced[i]] = k_idx - 1; }
    space->set_table_sizes();
    return total_ct;
  }

  class MHash {
    public:
      using X = std::tuple<unsigned __int128, long>;
      int rr = 0;
      uintE max_lvl;
      MTable mtable;
      X* space;

      template<class Graph>
      MHash(int r, Graph& DG, size_t max_deg, uintE _max_level){
        std::cout << "Init MHash" << std::endl; fflush(stdout);
        rr = r;
        max_lvl = _max_level;
        using W = typename Graph::weight_type;
        
        mtable.initialize(0, 0, max_lvl);
        mtable.allocate(DG.n, 0, r-1);
  
        auto init_induced = [&](HybridSpace_lw* induced) { induced->alloc(max_deg, r-1, DG.n, true, true); };
        auto finish_induced = [&](HybridSpace_lw* induced) { if (induced != nullptr) { delete induced; } };
        parallel_for_alloc<HybridSpace_lw>(init_induced, finish_induced, 0, DG.n, [&](size_t i, HybridSpace_lw* induced) {
          if (DG.get_vertex(i).getOutDegree() != 0) {
            auto next_space = mtable.next(i, 0, r-1);
            if (next_space == nullptr) next_space = &mtable;
            induced->setup(DG, r-1, i);
            //auto base = sequence<uintE>(r);
            //base[0] = i;
            NKCliqueDir_fast_hybrid_rec_multi(DG, 1, r-1, induced, next_space);
          }
        }, 1, false);

        std::cout << "End MHash Count" << std::endl; fflush(stdout);

        long total = mtable.set_table_sizes();
        space = pbbslib::new_array_no_init<X>(total);
        mtable.set_end_table_rec(space);

        std::cout << "End MHash" << std::endl; fflush(stdout);
      }

      void insert(sequence<uintE>& base, int r, int k) {
        mtable.insert(base, r, k);
      }

      std::size_t return_total() { return mtable.total_size; }

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
      void extract_indices(sequence<uintE>& base, I func, int r, int k) {
        mtable.extract_indices(base, func, r, k);
      }

      //Fill base[k] ... base[k-r+1] and base[0]
      template<class S, class Graph>
      void extract_clique(S index, sequence<uintE>& base, Graph& G, int k) {
        mtable.extract_clique(index, base, 0, rr, k);
      }
  };

  struct EndTable {
    pbbslib::sparse_table<unsigned __int128, long, hash128> table;
    uintE vtx;
    //MidTable* up_table;
  };

  struct MidTable {
    pbbslib::sparse_table<uintE, EndTable*, std::hash<uintE>> table;
    sequence<EndTable*> arr;
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
  
        //*** for arr
        top_table.arr = sequence<EndTable*>(DG.n, [](std::size_t i){return nullptr;});
        top_table_sizes = sequence<long>(DG.n + 1, long{0});
        /*top_table.table = pbbslib::sparse_table<uintE, EndTable*, std::hash<uintE>>(
          top_table_sizes2.size(),
          std::make_tuple<uintE, EndTable*>(UINT_E_MAX, static_cast<EndTable*>(nullptr)),
          std::hash<uintE>());
        top_table_sizes = sequence<long>(top_table.table.m + 1, long{0});*/
  
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
  
      void insert(sequence<uintE>& base, int r, int k) {
        auto add_f = [&] (long* ct, const std::tuple<unsigned __int128, long>& tup) {
          pbbs::fetch_and_add(ct, (long)1);
        };
        // Sort base
        uintE base2[10];
        assert(10 > k);
        for(std::size_t i = 0; i < k + 1; i++) {
          base2[i] = base[i];
        }
        std::sort(base2, base2 + k + 1,std::less<uintE>());
        //sequence<uintE> base2(base);
        //pbbs::sample_sort_inplace(base2.slice(), std::less<uintE>());

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
          //EndTable* end_table = top_table.table.find(vtx, nullptr);
          //***for arr
          EndTable* end_table = top_table.arr[vtx];
          //assert(end_table != nullptr);
          (end_table->table).insert_f(std::make_tuple(key, (long) 1), add_f);
        } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
      }

      std::size_t return_total() {return total;}

      size_t get_top_index(std::size_t index) {
        // This gives the first i such that top_table_sizes[i] >= index
        auto idx = pbbslib::binary_search(top_table_sizes, long{index}, std::less<long>());
        if (idx >= top_table_sizes.size()) return top_table_sizes.size() - 1;
        if (top_table_sizes[idx] == index) return idx;
        //if (top_table_sizes[idx] > index) {
          if (idx == 0) return idx;
          return idx - 1;
        //}
        //assert(top_table_sizes[0] == 0);
        /*
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
        return top_table_sizes.size() - 1;*/
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
        //***for arr
        //EndTable* end_table = std::get<1>(top_table.table.table[top_index]);
        EndTable* end_table = top_table.arr[top_index];
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
        //EndTable* end_table = std::get<1>(top_table.table.table[top_index]);
        //***for arr
        EndTable* end_table = top_table.arr[top_index];
        size_t bottom_index = index - top_table_sizes[top_index];
        auto val = std::get<1>((end_table->table).table[bottom_index]) - update;
        (end_table->table).table[bottom_index] = std::make_tuple(
          std::get<0>((end_table->table).table[bottom_index]), val
        );
        return val;
      }

      void clear_count(std::size_t index) {
        size_t top_index = get_top_index(index);
        //***for arr
        //EndTable* end_table = std::get<1>(top_table.table.table[top_index]);
        EndTable* end_table = top_table.arr[top_index];
        size_t bottom_index = index - top_table_sizes[top_index];
        (end_table->table).table[bottom_index] = std::make_tuple(
          std::get<0>((end_table->table).table[bottom_index]), 0
        );
      }

      template<class I>
      void extract_indices(sequence<uintE>& base, I func, int r, int k) {
        uintE base2[10];
        assert(10 > k);
        for(std::size_t i = 0; i < k + 1; i++) {
          base2[i] = base[i];
        }
        std::sort(base2, base2 + k + 1,std::less<uintE>());
        //sequence<uintE> base2(base);
        //pbbs::sample_sort_inplace(base2.slice(), std::less<uintE>());
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
          /*auto top_index = top_table.table.find_index(vtx);
          auto prefix = top_table_sizes[top_index];
          EndTable* end_table = std::get<1>(top_table.table.table[top_index]);*/
          //***for arr
          EndTable* end_table = top_table.arr[vtx];
          auto prefix = top_table_sizes[vtx];
          auto index = (end_table->table).find_index(key);
          func(prefix + index);
        } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
      }

      // Given an index, get the clique
      template<class S, class Graph>
      void extract_clique(S index, sequence<uintE>& base, Graph& G, int k) {
        // First, do a binary search for index in prefix
        size_t top_index = get_top_index(index);
        /*base[0] = std::get<0>(top_table.table.table[top_index]);
        EndTable* end_table = std::get<1>(top_table.table.table[top_index]);*/
        //***for arr
        base[0] = top_index;
        EndTable* end_table = top_table.arr[top_index];
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
        uintE base2[10];
        assert(10 > k);
        for(std::size_t i = 0; i < k + 1; i++) {
          base2[i] = base[i];
        }
        std::sort(base2, base2 + k + 1,std::less<uintE>());
        //sequence<uintE> base2(base);
        //pbbs::sample_sort_inplace(base2.slice(), std::less<uintE>());

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
        uintE base2[10];
        assert(10 > k);
        for(std::size_t i = 0; i < k + 1; i++) {
          base2[i] = base[i];
        }
        std::sort(base2, base2 + k + 1,std::less<uintE>());
        //sequence<uintE> base2(base);
        //pbbs::sample_sort_inplace(base2.slice(), std::less<uintE>());
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
  //TwolevelHash table(r, DG, max_deg);
  MHash table(r, DG, max_deg, 2);

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