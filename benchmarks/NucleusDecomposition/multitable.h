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

namespace multitable {

  inline bool is_uint_e_max(uintE max_val) {
    return max_val == UINT_E_MAX;
  }

  // max_lvl should be set to # levels - 2
  // two level hash is equiv to setting max_level to 0
  template <class Y, class H>
  struct MTable {
    using MTableY = MTable<Y, H>;
    using NextMTable = pbbslib::sparse_table<uintE, MTableY*, std::hash<uintE>>;
    using EndTable = pbbslib::sparse_table<Y, long, H>;
    NextMTable mtable;
    EndTable end_table;
    uintE lvl;
    uintE max_lvl;
    long total_size = 0;
    sequence<long> table_sizes;
    std::tuple<Y, long>* end_space = nullptr;

//#ifdef NUCLEUS_USE_VERTEX
    uintE vtx;
    uintE get_vtx(long idx) {
      return vtx;
    }
/*#else
    uintE vtx;
    uintE get_vtx(long idx) {
      if (lvl != max_lvl) {
        while (true) {
          uintE max_val = std::get<0>(mtable.table[idx]);
          std::size_t max_bit = sizeof(uintE) * 8;
          auto check_bit = (max_val >> (max_bit - 1)) & 1U;
          if (check_bit) {
            max_val &= ~(1UL << (max_bit - 1));
            return max_val;
          }
          idx++;
          idx = idx % mtable.m;
        }
      }
      return vtx;
    }
#endif*/

    void set_end_table_rec(std::tuple<Y, long>* _end_space) {
      end_space = _end_space;

      if (lvl == max_lvl) {
        // Allocate end_table here using end_space
        end_table = EndTable(
          total_size,
          std::make_tuple<Y, long>(std::numeric_limits<Y>::max(), static_cast<long>(0)),
          H{},
          end_space
        );

      } else {
        for (std::size_t i = 0; i < mtable.m; i++) {
          if (!is_uint_e_max(std::get<0>(mtable.table[i]))) {
            std::get<1>(mtable.table[i])->set_end_table_rec(end_space + table_sizes[i]);
          }
        }
      }
    }

    void set_end_table_rec() {
      if (lvl == max_lvl) {
        // Allocate end_table here
        end_table = EndTable(
          total_size,
          std::make_tuple<Y, long>(std::numeric_limits<Y>::max(), static_cast<long>(0)),
          H{}, 1, true
        );

      } else {
        for (std::size_t i = 0; i < mtable.m; i++) {
          if (!is_uint_e_max(std::get<0>(mtable.table[i]))) {
            std::get<1>(mtable.table[i])->set_end_table_rec();
          }
        }
      }
    }

    void increment_size(long s) {
      if (lvl == max_lvl) total_size += s;
    }

    long set_table_sizes() {
      if (lvl != max_lvl) {
        if (lvl + 1 == max_lvl) {
          table_sizes = sequence<long>(mtable.m, [](std::size_t i){ return 0; });
          parallel_for(0, mtable.m, [&](std::size_t i){
            if (!is_uint_e_max(std::get<0>(mtable.table[i]))) {
              auto tbl = std::get<1>(mtable.table[i]);
              tbl->total_size = (size_t)1 << pbbslib::log2_up((size_t)(1.1 * tbl->total_size) + 1);
              table_sizes[i] = tbl->total_size;
            }
          });
          total_size = scan_inplace(table_sizes.slice(), pbbs::addm<long>());
          return total_size;
        }
        table_sizes = sequence<long>(mtable.m, [&](std::size_t i){
          if (is_uint_e_max(std::get<0>(mtable.table[i]))) return long{0};
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

    void initialize(uintE _v, uintE _lvl, uintE _max_lvl, MTableY* _prev_mtable = nullptr) {
//#ifdef NUCLEUS_USE_VERTEX
      vtx = _v;
//#endif
      lvl = _lvl;
      max_lvl = _max_lvl;
      total_size = 0;
    }

    // k_idx goes from 1 to k - 1 (inclusive)
    void allocate(uintE count, size_t k_idx, size_t k, uintE prev_vtx) {
      //lvl = k_idx - 1;
      if (lvl != max_lvl) {
//#ifdef NUCLEUS_USE_VERTEX
        uintE max_val = UINT_E_MAX;
/*#else
        uintE max_val = prev_vtx;
        std::size_t max_bit = sizeof(uintE) * 8;
        max_val |= 1ULL << (max_bit - 1);
        uintE x = max_val & ~(1UL << (max_bit - 1));
#endif*/
        mtable = NextMTable(
          count * 1.1,
          std::make_tuple<uintE, MTableY*>(uintE{max_val}, static_cast<MTableY*>(nullptr)),
          std::hash<uintE>()
        );
      }
    }
    
    MTableY* next(uintE v, size_t k_idx, size_t k) {
      if (lvl == max_lvl) return nullptr;
      MTableY* next_table = new MTableY();
      next_table->initialize(v, lvl + 1, max_lvl, this);
      mtable.insert(std::make_tuple(v, next_table));
      return next_table;
    }

    void insert(uintE* base, int curr_idx, int r, int k, std::string& bitmask) {
      if (lvl == max_lvl) {
        assert(end_table.m > 0);
        auto add_f = [&] (long* ct, const std::tuple<Y, long>& tup) {
          pbbs::fetch_and_add(ct, (long)1);
        };
        Y key = 0;
        unsigned __int128 mask = (1ULL << (nd_global_shift_factor)) - 1;
        for (int i = curr_idx; i < static_cast<int>(k)+1; ++i) {
          if (bitmask[i]) {
            key = key << nd_global_shift_factor;
            key |= (base[i] & mask);
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

    void insert(uintE* base, int r, int k) {
        std::string bitmask(r+1, 1); // K leading 1's
        bitmask.resize(k+1, 0); // N-K trailing 0's

        do {
          Y key = 0;
          for (int i = 0; i < static_cast<int>(k)+1; ++i) {
            if (bitmask[i]) {
              insert(base, i, r, k, bitmask);
              break;
            }
          }
        } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
    }

    size_t extract_indices(uintE* base, int curr_idx,
      int r, int k, std::string& bitmask) {
      if (lvl == max_lvl) {
        Y key = 0;
        unsigned __int128 mask = (1ULL << (nd_global_shift_factor)) - 1;
        for (int i = curr_idx; i < static_cast<int>(k)+1; ++i) {
          if (bitmask[i]) {
            key = key << nd_global_shift_factor;
            key |= (base[i] & mask);
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

    Y extract_indices_check(uintE* base, int curr_idx, int r) {
      if (lvl == max_lvl) {
        Y key = 0;
        unsigned __int128 mask = (1ULL << (nd_global_shift_factor)) - 1;
        for (int i = curr_idx; i < static_cast<int>(r)+1; ++i) {
          key = key << nd_global_shift_factor;
          key |= (base[i] & mask);
        }
        return end_table.find_index(key);
      }
      int next_idx = curr_idx + 1;
      auto next_mtable_index = mtable.find_index(base[curr_idx]);
      auto next = std::get<1>(mtable.table[next_mtable_index]);
      assert(next != nullptr);
      return table_sizes[next_mtable_index] + next->extract_indices_check(base, next_idx, r);
    }

    template<class I, class HH>
    void extract_indices(uintE* base, I func, int r, int k, HH& h_func) {
      std::string bitmask(r+1, 1); // K leading 1's
      bitmask.resize(k+1, 0); // N-K trailing 0's
      do {
        for (int i = 0; i < static_cast<int>(k)+1; ++i) {
          if (bitmask[i]) {
            // TODO: make sure this pointer arithmetic is fine
            size_t idx = extract_indices(base, i, r, k, bitmask); // - end_space;
            h_func(idx);
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
      if (table_sizes[idx] == index) {
        while(idx < table_sizes.size() && table_sizes[idx] == index) {
          idx++;
        }
        return idx - 1;
      }
      if (idx == 0) return idx;
      return idx - 1;
    }
  
    //Fill base[k] ... base[k-r+1] and base[0]
    template<class S>
    void extract_clique(S index, uintE* base, int base_idx, int rr, int k) {
      if (lvl == max_lvl) {
        if (lvl != 0) {
          // TODO: not sure if we should be doing 0...
          base[base_idx] = get_vtx(0);
          if (base_idx == 0) base_idx = k - rr + 2;
          else base_idx++;
        }
        //assert(end_space != nullptr);
        auto vert = std::get<0>(end_table.table[index]); //end_space[index]
        // TOOD: make sure this calc is correct
        unsigned __int128 mask = (1ULL << (nd_global_shift_factor)) - 1;
        for (int j = k; j >= base_idx; --j) { //rr - 1, base_idx
          uintE extract = (uintE) (vert & mask); // vert & mask
          //assert(static_cast<uintE>(extract) < G.n);
          base[j] = static_cast<uintE>(extract);
          vert = vert >> nd_global_shift_factor;
        }
        return;
      }
      auto next_mtable_idx = get_top_index(index);
      /*if (next_mtable_idx >= mtable.m) {
        std::cout << "Idx: " << next_mtable_idx << std::endl;
        std::cout << "m: " << mtable.m << std::endl;
        fflush(stdout);
      }
      assert(next_mtable_idx < mtable.m);
      assert(index >= table_sizes[next_mtable_idx]);*/
      S next_index = index - table_sizes[next_mtable_idx];
      if (lvl != 0) {
        base[base_idx] = get_vtx(next_mtable_idx);
        if (base_idx == 0) base_idx = k - rr + 2;
        else base_idx++;
      }
      /*if (std::get<1>(mtable.table[next_mtable_idx]) == nullptr) {
        std::cout << "rr: " << rr << std::endl; fflush(stdout);
        std::cout << "base_idx: " << base_idx << std::endl; fflush(stdout);
        std::cout << "index: " << long{index} << std::endl; fflush(stdout);
        std::cout << "top index: " << long{next_mtable_idx} << std::endl; fflush(stdout);
        std::cout << "size: " << table_sizes[next_mtable_idx] << std::endl; fflush(stdout);
      }
      assert(std::get<1>(mtable.table[next_mtable_idx]) != nullptr);*/
      std::get<1>(mtable.table[next_mtable_idx])->extract_clique(next_index, base, base_idx, rr, k);
    }

    template<class S, class F>
    void find_table_loc(S index, F func) {
      if (lvl == max_lvl) {
        func(&(end_table.table[index]));
        return;
      }
      auto next_mtable_idx = get_top_index(index);
      /*if (next_mtable_idx >= mtable.m) {
        std::cout << "Idx: " << next_mtable_idx << std::endl;
        std::cout << "m: " << mtable.m << std::endl;
        fflush(stdout);
      }
      assert(next_mtable_idx < mtable.m);
      assert(index >= table_sizes[next_mtable_idx]);*/
      S next_index = index - table_sizes[next_mtable_idx];
      if (std::get<1>(mtable.table[next_mtable_idx]) == nullptr) return;
      /*if (std::get<1>(mtable.table[next_mtable_idx]) == nullptr) {
        std::cout << "rr: " << rr << std::endl; fflush(stdout);
        std::cout << "base_idx: " << base_idx << std::endl; fflush(stdout);
        std::cout << "index: " << long{index} << std::endl; fflush(stdout);
        std::cout << "top index: " << long{next_mtable_idx} << std::endl; fflush(stdout);
        std::cout << "size: " << table_sizes[next_mtable_idx] << std::endl; fflush(stdout);
      }
      assert(std::get<1>(mtable.table[next_mtable_idx]) != nullptr);*/
      std::get<1>(mtable.table[next_mtable_idx])->find_table_loc(next_index, func);
    }

  };

  template <class Graph, class Space>
  inline size_t NKCliqueDir_fast_hybrid_rec_multi(Graph& DG, size_t k_idx, size_t k,
  HybridSpace_lw* induced, Space* space, uintE prev_vtx, bool valid_space = true) {
    size_t num_induced = induced->num_induced[k_idx-1];
    if (num_induced == 0) return 0;
    uintE* prev_induced = induced->induced + induced->nn * (k_idx - 1);

    for (size_t i=0; i < num_induced; i++) { induced->labels[prev_induced[i]] = k_idx; }

    if (k_idx + 1 == k) {
      size_t counts = 0;
      space->allocate(num_induced, k_idx, k, prev_vtx);
      for (size_t i=0; i < num_induced; i++) {
        uintE vtx = prev_induced[i];
        //  get neighbors of vtx
        uintE* intersect = induced->induced_edges + vtx * induced->nn;
        size_t tmp_counts = 0;
        //base[k_idx] = induced->relabel[vtx];
        Space* next_space = space->next(induced->relabel[vtx], k_idx, k);
        bool next_valid_space = true;
        if (next_space == nullptr) {
          next_space = space;
          next_valid_space = false;
        }
        for (size_t j=0; j < induced->induced_degs[vtx]; j++) {
          if (static_cast<size_t>(induced->labels[intersect[j]]) == k_idx) {
            //base[k] = induced->relabel[intersect[j]];
            tmp_counts++;
            //base_f(base, space);
            // TODO: This could be improved by collating on count before
            // returning to recursive level of max_lvl
          }
        }
        next_space->increment_size(tmp_counts);
        if (next_valid_space) {
          next_space->set_table_sizes();
        }
        counts += tmp_counts;
      }
      for (size_t i=0; i < num_induced; i++) { induced->labels[prev_induced[i]] = k_idx - 1; }
      if (valid_space) {
        //space->increment_size(counts);
        space->set_table_sizes();
      }
      return counts;
    }

    size_t total_ct = 0;
    space->allocate(num_induced, k_idx, k, prev_vtx);
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
        bool next_valid_space = true;
        if (next_space == nullptr) {
          next_space = space;
          next_valid_space = false;
        }
        auto curr_counts = NKCliqueDir_fast_hybrid_rec_multi(DG, k_idx + 1, k, induced, next_space, induced->relabel[vtx], next_valid_space);
        total_ct += curr_counts;
      }
    }
    for (size_t i=0; i < num_induced; i++) { induced->labels[prev_induced[i]] = k_idx - 1; }
    if (valid_space) {
      //space->increment_size(total_ct);
      space->set_table_sizes();
    }
    return total_ct;
  }

  template <class Y, class H, class F>
  class MHash {
    public:
      using X = std::tuple<Y, long>;
      using MTableY = MTable<Y, H>;
      int rr = 0;
      uintE max_lvl;
      MTableY mtable;
      X* space = nullptr;
      bool contiguous_space;
      bool relabel;
      F sort_func;

      template<class Graph>
      MHash(int r, Graph& DG, size_t max_deg, uintE _max_level, bool _contiguous_space,
        F _rank_func) : sort_func(_rank_func) {
        using W = typename Graph::weight_type;
        contiguous_space = _contiguous_space;
        //std::cout << "Init MHash" << std::endl; fflush(stdout);
        rr = r;
        max_lvl = _max_level;
        
        mtable.initialize(0, 0, max_lvl);
        mtable.allocate(DG.n, 0, r-1, 0);

        if (r == 2) {
          parallel_for(0, DG.n, [&](std::size_t i){
            auto next_space = mtable.next(i, 0, r-1);
            bool valid_space = true;
            if (next_space == nullptr) {
              valid_space = false;
              next_space = &mtable;
            }
            auto map_f = [&](const uintE& src, const uintE& ngh, const W& wgh) {
              next_space->allocate(DG.get_vertex(i).getOutDegree(), r-2, r-1, i);
              auto next_next_space = next_space->next(ngh, r-2, r-1);
              bool next_next_valid_space = true;
              if (next_next_space == nullptr) {
                next_next_space = next_space;
                next_next_valid_space = false;
              }
              next_next_space->increment_size(1);
              if (next_next_valid_space) {
                next_next_space->set_table_sizes();
              }
            };
            DG.get_vertex(i).mapOutNgh(i, map_f, true);
            if (valid_space) {
              next_space->set_table_sizes();
            }
          });
        } else {
  
        auto init_induced = [&](HybridSpace_lw* induced) { induced->alloc(max_deg, r-1, DG.n, true, true); };
        auto finish_induced = [&](HybridSpace_lw* induced) { if (induced != nullptr) { delete induced; } };
        parallel_for_alloc<HybridSpace_lw>(init_induced, finish_induced, 0, DG.n, [&](size_t i, HybridSpace_lw* induced) {
          if (DG.get_vertex(i).getOutDegree() != 0) {
            auto next_space = mtable.next(i, 0, r-1);
            if (next_space == nullptr) next_space = &mtable;
            induced->setup(DG, r-1, i);
            //auto base = sequence<uintE>(r);
            //base[0] = i;
            NKCliqueDir_fast_hybrid_rec_multi(DG, 1, r-1, induced, next_space, i);
          }
        }, 1, false);

        }

        //std::cout << "End MHash Count" << std::endl; fflush(stdout);

        long total = mtable.set_table_sizes();
        if (contiguous_space) {
          space = pbbslib::new_array_no_init<X>(total);
          mtable.set_end_table_rec(space);
        } else {
          mtable.set_end_table_rec();
        }

        //std::cout << "End MHash" << std::endl; fflush(stdout);
      }

      void insert(sequence<uintE>& base2, int r, int k) {
        uintE base[10];
        assert(10 > k);
        for(std::size_t i = 0; i < k + 1; i++) {
          base[i] = base2[i];
        }
        std::sort(base, base + k + 1,sort_func);
        mtable.insert(base, r, k);
      }

      void insert_twothree(uintE v1, uintE v2, uintE v3, int r, int k) {
        auto add_f = [&] (long* ct, const std::tuple<Y, long>& tup) {
          pbbs::fetch_and_add(ct, (long)1);
        };
        unsigned __int128 mask = (1ULL << (nd_global_shift_factor)) - 1;

        uintE min13 = sort_func(v1, v3) ? v1 : v3;
        uintE max13 = sort_func(v1, v3) ? v3 : v1;
        // Level 0
        auto next13 = mtable.mtable.find(min13, nullptr);
        // Level 1
        Y key13 = max13 & mask;
        next13->end_table.insert_f(std::make_tuple(key13, (long) 1), add_f);

        uintE min23 = sort_func(v2, v3) ? v2 : v3;
        uintE max23 = sort_func(v2, v3) ? v3 : v2;
        // Level 0
        auto next23 = mtable.mtable.find(min23, nullptr);
        // Level 1
        Y key23 = max23 & mask;
        next23->end_table.insert_f(std::make_tuple(key23, (long) 1), add_f);

        uintE min12 = sort_func(v1, v2) ? v1 : v2;
        uintE max12 = sort_func(v1, v2) ? v2 : v1;
        // Level 0
        auto next12 = mtable.mtable.find(min12, nullptr);
        // Level 1
        Y key12 = max12 & mask;
        next12->end_table.insert_f(std::make_tuple(key12, (long) 1), add_f);
      }

      std::size_t return_total() { return mtable.total_size; }

      long get_count(std::size_t index) {
        if (contiguous_space) {
          if (std::get<0>(space[index]) == std::numeric_limits<Y>::max()) return 0;
          return std::get<1>(space[index]);
        }

        long count = 0;
        auto func = [&](std::tuple<Y, long>* loc){
          if (std::get<0>(*loc) == std::numeric_limits<Y>::max()) count = 0;
          else count = std::get<1>(*loc);
        };
        mtable.find_table_loc(index, func);
        return count;
      }

      size_t update_count(std::size_t index, size_t update){
        if (contiguous_space) {
          auto val = std::get<1>(space[index]) - update;
          space[index] =
            std::make_tuple(std::get<0>(space[index]), val);
          return val;
        }
        size_t val = 0;
        auto func = [&](std::tuple<Y, long>* loc){
          val = std::get<1>(*loc) - update;
          *loc = std::make_tuple(std::get<0>(*loc), val);
        };
        mtable.find_table_loc(index, func);
        return val;
      }

      void clear_count(std::size_t index) {
        if (contiguous_space) {
          space[index] = std::make_tuple(std::get<0>(space[index]), 0);
          return;
        }
        auto func = [&](std::tuple<Y, long>* loc){
          *loc = std::make_tuple(std::get<0>(*loc), 0);
        };
        mtable.find_table_loc(index, func);
      }

      void set_count(std::size_t index, size_t update) {
        if (contiguous_space) {
          space[index] = std::make_tuple(std::get<0>(space[index]), update);
          return;
        }
        auto func = [&](std::tuple<Y, long>* loc){
          *loc = std::make_tuple(std::get<0>(*loc), update);
        };
        mtable.find_table_loc(index, func);
      }

      Y extract_indices_check(uintE* base2, int r) {
        // Size of base2 should be r + 1
        uintE base[10];
        assert(10 > r + 1);
        for(std::size_t i = 0; i < r + 1; i++) {
          base[i] = base2[i];
        }
        std::sort(base, base + r + 1,std::less<uintE>());

        return mtable.extract_indices_check(base, 0, r);
      }

      template<class HH, class HG, class I>
      void extract_indices(uintE* base2, HH is_active, HG is_inactive, I func, int r, int k) {
        uintE base[10];
        assert(10 > k);
        for(std::size_t i = 0; i < k + 1; i++) {
          base[i] = base2[i];
        }
        std::sort(base, base + k + 1,sort_func);

        std::vector<size_t> indices;
        size_t num_active = 0;
        bool use_func = true;

        auto h_func = [&](std::size_t idx){
          indices.push_back(idx);
          if (is_active(idx)) num_active++;
          if (is_inactive(idx)) use_func = false;
        };

        mtable.extract_indices(base, func, r, k, h_func);

        assert(num_active != 0);

        if (use_func) {
          for (std::size_t i = 0; i < indices.size(); i++) {
            if (!is_active(indices[i]) && !is_inactive(indices[i]))
              func(indices[i], 1.0 / (double) num_active);
          }
        }
      }

      Y extract_indices_two(uintE v1, uintE v3) {
        unsigned __int128 mask = (1ULL << (nd_global_shift_factor)) - 1;
        size_t index13 = 0;
        // Level 0
        uintE min13 = sort_func(v1, v3) ? v1 : v3;
        uintE max13 = sort_func(v1, v3) ? v3 : v1;
        auto next_mtable_index13 = mtable.mtable.find_index(min13);
        auto next13 = std::get<1>(mtable.mtable.table[next_mtable_index13]);
        index13 += mtable.table_sizes[next_mtable_index13];

        // Level 1
        Y key13 = max13 & mask;
        index13 += next13->end_table.find_index(key13);
        return index13;
      }

      template<class HH, class HG, class I>
      void extract_indices_twothree(uintE v1, uintE v2, uintE v3, HH is_active, 
        HG is_inactive, I func, int r, int k) {
        unsigned __int128 mask = (1ULL << (nd_global_shift_factor)) - 1;
        size_t index13 = 0;
        size_t index23 = 0;
        size_t num_active = 1;
        // Level 0
        uintE min13 = sort_func(v1, v3) ? v1 : v3;
        uintE max13 = sort_func(v1, v3) ? v3 : v1;
        auto next_mtable_index13 = mtable.mtable.find_index(min13);
        auto next13 = std::get<1>(mtable.mtable.table[next_mtable_index13]);
        index13 += mtable.table_sizes[next_mtable_index13];

        // Level 1
        Y key13 = max13 & mask;
        index13 += next13->end_table.find_index(key13);
        if (is_inactive(index13)) return;

        // Level 0
        uintE min23 = sort_func(v2, v3) ? v2 : v3;
        uintE max23 = sort_func(v2, v3) ? v3 : v2;
        auto next_mtable_index23 = mtable.mtable.find_index(min23);
        auto next23 = std::get<1>(mtable.mtable.table[next_mtable_index23]);
        index23 += mtable.table_sizes[next_mtable_index23];

        // Level 1
        Y key23 = max23 & mask;
        index23 += next23->end_table.find_index(key23);
        if (is_inactive(index23)) return;

        bool is_active_13 = is_active(index13);
        bool is_active_23 = is_active(index23);
        if (is_active_13 && is_active_23) return;
        if (is_active_13 || is_active_23) {
          // Level 0
        size_t index12 = 0;
        uintE min12 = sort_func(v1, v2) ? v1 : v2;
        uintE max12 = sort_func(v1, v2) ? v2 : v1;
        auto next_mtable_index12 = mtable.mtable.find_index(min12);
        auto next12 = std::get<1>(mtable.mtable.table[next_mtable_index12]);
        index12 += mtable.table_sizes[next_mtable_index12];

        // Level 1
        Y key12 = max12 & mask;
        index12 += next12->end_table.find_index(key12);
          if (is_active_13 && index12 < index13) func(index23, 1.0);
          else if (is_active_23 && index12 < index23) func(index13, 1.0);
          return;
        }
        func(index23, 1.0);
        func(index13, 1.0);
      }


      template<class HH, class HG, class I>
      void extract_indices_threefour(uintE v1, uintE v2, uintE v3, uintE v4,
        HH is_active, HG is_inactive, I func, int r, int k) {
        unsigned __int128 mask = (1ULL << (nd_global_shift_factor)) - 1;
        size_t index134 = 0;
        size_t index234 = 0;
        size_t index124 = 0;
        size_t num_active = 1;
        uintE min134 = minOfThree(v1, v3, v4, sort_func);
        uintE mid134 = middleOfThree(v1, v3, v4, sort_func);
        uintE max134 = maxOfThree(v1, v3, v4, sort_func);

        uintE min234 = minOfThree(v2, v3, v4, sort_func);
        uintE mid234 = middleOfThree(v2, v3, v4, sort_func);
        uintE max234 = maxOfThree(v2, v3, v4, sort_func);

        uintE min124 = minOfThree(v1, v2, v4, sort_func);
        uintE mid124 = middleOfThree(v1, v2, v4, sort_func);
        uintE max124 = maxOfThree(v1, v2, v4, sort_func);

        if (max_lvl == 2) { // three levels
        // Level 0
        auto next_mtable_index134 = mtable.mtable.find_index(min134);
        auto next134 = std::get<1>(mtable.mtable.table[next_mtable_index134]);
        index134 += mtable.table_sizes[next_mtable_index134];
        // Level 1
        auto next_next_mtable_index134 = next134->mtable.find_index(mid134);
        auto next_next134 = std::get<1>(next134->mtable.table[next_next_mtable_index134]);
        index134 += next134->table_sizes[next_next_mtable_index134];
        // Level 2
        Y key134 = max134 & mask;
        index134 += next_next134->end_table.find_index(key134);
        if (is_inactive(index134)) return;

        // Level 0
        auto next_mtable_index234 = mtable.mtable.find_index(min234);
        auto next234 = std::get<1>(mtable.mtable.table[next_mtable_index234]);
        index234 += mtable.table_sizes[next_mtable_index234];
        // Level 1
        auto next_next_mtable_index234 = next234->mtable.find_index(mid234);
        auto next_next234 = std::get<1>(next234->mtable.table[next_next_mtable_index234]);
        index234 += next234->table_sizes[next_next_mtable_index234];
        // Level 2
        Y key234 = max234 & mask;
        index234 += next_next234->end_table.find_index(key234);
        if (is_inactive(index234)) return;

        // Level 0
        auto next_mtable_index124 = mtable.mtable.find_index(min124);
        auto next124 = std::get<1>(mtable.mtable.table[next_mtable_index124]);
        index124 += mtable.table_sizes[next_mtable_index124];
        // Level 1
        auto next_next_mtable_index124 = next124->mtable.find_index(mid124);
        auto next_next124 = std::get<1>(next124->mtable.table[next_next_mtable_index124]);
        index124 += next124->table_sizes[next_next_mtable_index124];
        // Level 2
        Y key124 = max124 & mask;
        index124 += next_next124->end_table.find_index(key124);
        if (is_inactive(index124)) return;



        bool is_active_124 = is_active(index124);
        bool is_active_134 = is_active(index134);
        bool is_active_234 = is_active(index234);
        if (is_active_124 && is_active_134 && is_active_234) return;
        if (is_active_124 || is_active_134 || is_active_234) {
          uintE min123 = minOfThree(v1, v3, v2, sort_func);
          uintE mid123 = middleOfThree(v1, v3, v2, sort_func);
          uintE max123 = maxOfThree(v1, v3, v2, sort_func);
          size_t index123 = 0;
          // Level 0
          auto next_mtable_index123 = mtable.mtable.find_index(min123);
          auto next123 = std::get<1>(mtable.mtable.table[next_mtable_index123]);
          index123 += mtable.table_sizes[next_mtable_index123];
          // Level 1
          auto next_next_mtable_index123 = next123->mtable.find_index(mid123);
          auto next_next123 = std::get<1>(next123->mtable.table[next_next_mtable_index123]);
          index123 += next123->table_sizes[next_next_mtable_index123];
          // Level 2
          Y key123 = max123 & mask;
          index123 += next_next123->end_table.find_index(key123);
          if (is_active_124 && index123 > index124) return;
          if (is_active_134 && index123 > index134) return;
          if (is_active_234 && index123 > index234) return;
          if (!is_active(index124)) func(index124, 1.0);
          if (!is_active(index134)) func(index134, 1.0);
          if (!is_active(index234)) func(index234, 1.0);
          return;
        }
        func(index124, 1.0);
        func(index134, 1.0);
        func(index234, 1.0);
        return;
        }
        // max_lvl = 1
        // Level 0
        auto next_mtable_index134 = mtable.mtable.find_index(min134);
        auto next134 = std::get<1>(mtable.mtable.table[next_mtable_index134]);
        index134 += mtable.table_sizes[next_mtable_index134];
        // Level 1
        Y key134 = mid134 & mask;
        key134 = key134 << nd_global_shift_factor;
        key134 |= (max134 & mask);
        index134 += next134->end_table.find_index(key134);
        if (is_inactive(index134)) return;

        // Level 0
        auto next_mtable_index234 = mtable.mtable.find_index(min234);
        auto next234 = std::get<1>(mtable.mtable.table[next_mtable_index234]);
        index234 += mtable.table_sizes[next_mtable_index234];
        // Level 1
        Y key234 = mid234 & mask;
        key234 = key234 << nd_global_shift_factor;
        key234 |= (max234 & mask);
        index234 += next234->end_table.find_index(key234);
        if (is_inactive(index234)) return;

        // Level 0
        auto next_mtable_index124 = mtable.mtable.find_index(min124);
        auto next124 = std::get<1>(mtable.mtable.table[next_mtable_index124]);
        index124 += mtable.table_sizes[next_mtable_index124];
        // Level 1
        Y key124 = mid124 & mask;
        key124 = key124 << nd_global_shift_factor;
        key124 |= (max124 & mask);
        index124 += next124->end_table.find_index(key124);
        if (is_inactive(index124)) return;

        bool is_active_124 = is_active(index124);
        bool is_active_134 = is_active(index134);
        bool is_active_234 = is_active(index234);
        if (is_active_124 && is_active_134 && is_active_234) return;
        if (is_active_124 || is_active_134 || is_active_234) {
          uintE min123 = minOfThree(v1, v3, v2, sort_func);
          uintE mid123 = middleOfThree(v1, v3, v2, sort_func);
          uintE max123 = maxOfThree(v1, v3, v2, sort_func);
          size_t index123 = 0;
          // Level 0
          auto next_mtable_index123 = mtable.mtable.find_index(min123);
          auto next123 = std::get<1>(mtable.mtable.table[next_mtable_index123]);
          index123 += mtable.table_sizes[next_mtable_index123];
          // Level 1
          Y key123 = mid123 & mask;
          key123 = key123 << nd_global_shift_factor;
          key123 |= (max123 & mask);
          index123 += next123->end_table.find_index(key123);
          if (is_active_124 && index123 > index124) return;
          if (is_active_134 && index123 > index134) return;
          if (is_active_234 && index123 > index234) return;
          if (!is_active(index124)) func(index124, 1.0);
          if (!is_active(index134)) func(index134, 1.0);
          if (!is_active(index234)) func(index234, 1.0);
          return;
        }
        func(index124, 1.0);
        func(index134, 1.0);
        func(index234, 1.0);
        return;
      }

      //Fill base[k] ... base[k-r+2] and base[0]
      template<class S, class Graph>
      void extract_clique(S index, uintE* base, Graph& G, int k) {
        mtable.extract_clique(index, base, 0, rr, k);
      }
    
    template<class S>
    std::tuple<uintE, uintE> extract_clique_two(S index, int k) {
      // If max level is 1 (2 levels)
      S next_index = index;
      // Level 0
      auto next_mtable_idx = mtable.get_top_index(next_index);
      next_index -= mtable.table_sizes[next_mtable_idx];
      auto mtable_ptr_1 = std::get<1>(mtable.mtable.table[next_mtable_idx]);

      // Level 1
      uintE v1 = mtable_ptr_1->get_vtx(0);
      auto vert = std::get<0>(mtable_ptr_1->end_table.table[next_index]);
      unsigned __int128 mask = (1ULL << (nd_global_shift_factor)) - 1;
      uintE extract = (uintE) (vert & mask); // vert & mask
      uintE v2 = static_cast<uintE>(extract);
      return std::make_tuple(v1, v2);
    }

    template<class S>
    std::tuple<uintE, uintE, uintE> extract_clique_three(S index, int k) {
      if (max_lvl == 2) { // 3 levels
        S next_index = index;
        // Level 0
        auto next_mtable_idx = mtable.get_top_index(next_index);
        next_index -= mtable.table_sizes[next_mtable_idx];
        auto mtable_ptr_1 = std::get<1>(mtable.mtable.table[next_mtable_idx]);

        // Level 1
        auto next_mtable_idx_1 = mtable_ptr_1->get_top_index(next_index);
        next_index -= mtable_ptr_1->table_sizes[next_mtable_idx_1];
        uintE v1 = mtable_ptr_1->get_vtx(0);
        auto mtable_ptr_2 = std::get<1>(mtable_ptr_1->mtable.table[next_mtable_idx_1]);

        // Level 2
        uintE v2 = mtable_ptr_2->get_vtx(0);
        auto vert = std::get<0>(mtable_ptr_2->end_table.table[next_index]);
        unsigned __int128 mask = (1ULL << (nd_global_shift_factor)) - 1;
        uintE extract = (uintE) (vert & mask);
        uintE v3 = static_cast<uintE>(extract);
        return std::make_tuple(v1, v2, v3);
      }
        // If max level is 1 (2 levels)
        S next_index = index;
        // Level 0
        auto next_mtable_idx = mtable.get_top_index(next_index);
        next_index -= mtable.table_sizes[next_mtable_idx];
        auto mtable_ptr_1 = std::get<1>(mtable.mtable.table[next_mtable_idx]);

        // Level 1
        uintE v1 = mtable_ptr_1->get_vtx(0);
        auto vert = std::get<0>(mtable_ptr_1->end_table.table[next_index]);
        unsigned __int128 mask = (1ULL << (nd_global_shift_factor)) - 1;
        uintE v2 = UINT_E_MAX; uintE v3 = UINT_E_MAX;
        for (int j = 0; j < 2; j++) {
          uintE extract = (uintE) (vert & mask); // vert & mask
          if (j == 0) v2 = static_cast<uintE>(extract);
          else v3 = static_cast<uintE>(extract);
          vert = vert >> nd_global_shift_factor;
        }
        return std::make_tuple(v1, v2, v3);
    }
  };

} // end namespace multitable

} //end namespace gbbs