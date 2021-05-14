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
//#ifdef NUCLEUS_USE_VERTEX
    return max_val == UINT_E_MAX;
/*#else
    std::size_t max_bit = sizeof(uintE) * 8;
    auto check_bit = (max_val >> (max_bit - 1)) & 1U;
    return check_bit;
#endif*/
  }

  // max_lvl should be set to # levels - 2
  // two level hash is equiv to setting max_level to 0
  struct MTable {
    using NextMTable = pbbslib::sparse_table<uintE, MTable*, std::hash<uintE>>;
    using EndTable = pbbslib::sparse_table<unsigned __int128, long, hash128>;
    NextMTable mtable;
    EndTable end_table;
    uintE lvl;
    uintE max_lvl;
    long total_size = 0;
    sequence<long> table_sizes;
    std::tuple<unsigned __int128, long>* end_space = nullptr;

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
          std::make_tuple<unsigned __int128, long>(std::numeric_limits<unsigned __int128>::max(), static_cast<long>(0)),
          hash128{}, 1
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
          table_sizes = sequence<long>(mtable.m + 1, [](std::size_t i){ return 0; });
          parallel_for(0, mtable.m, [&](std::size_t i){
            if (!is_uint_e_max(std::get<0>(mtable.table[i]))) {
              auto tbl = std::get<1>(mtable.table[i]);
              tbl->total_size = (size_t)1 << pbbslib::log2_up((size_t)(1.5 * tbl->total_size) + 1);
              table_sizes[i] = tbl->total_size;
            }
          });
          table_sizes[mtable.m] = 0;
          total_size = scan_inplace(table_sizes.slice(), pbbs::addm<long>());
          return total_size;
        }
        table_sizes = sequence<long>(mtable.m + 1, [&](std::size_t i){
          if (is_uint_e_max(std::get<0>(mtable.table[i]))) return long{0};
          return std::get<1>(mtable.table[i])->total_size;
        });
        table_sizes[mtable.m] = 0;
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
          count * 1.5,
          std::make_tuple<uintE, MTable*>(uintE{max_val}, static_cast<MTable*>(nullptr)),
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

    void insert(uintE* base, int curr_idx, int r, int k, std::string& bitmask) {
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

    void insert(uintE* base, int r, int k) {
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
    size_t extract_indices(uintE* base, int curr_idx,
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
    void extract_indices(uintE* base, I func, int r, int k) {
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
      if (idx >= table_sizes.size()) return table_sizes.size() - 2;
      if (idx == 0) return idx;
      if (table_sizes[idx] == index) {
        while(table_sizes[idx] == index) {
          idx--;
        }
        return idx + 1;
      }
      return idx - 1;
    }
  
    //Fill base[k] ... base[k-r+1] and base[0]
    template<class S>
    void extract_clique(S index, sequence<uintE>& base, int base_idx, int rr, int k) {
      if (lvl == max_lvl) {
        if (lvl != 0) {
          // TODO: not sure if we should be doing 0...
          base[base_idx] = get_vtx(0);
          if (base_idx == 0) base_idx = k - rr + 1;
          else base_idx++;
        }
        //assert(end_space != nullptr);
        auto vert = std::get<0>(end_table.table[index]); //end_space[index]
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
        base[base_idx] = get_vtx(next_mtable_idx);
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

  template <class F>
  class MHash {
    public:
      using X = std::tuple<unsigned __int128, long>;
      int rr = 0;
      uintE max_lvl;
      MTable mtable;
      X* space = nullptr;
      bool contiguous_space;
      bool relabel;
      F sort_func;

      template<class Graph>
      MHash(int r, Graph& DG, size_t max_deg, uintE _max_level, bool _contiguous_space,
        F _rank_func) : sort_func(_rank_func) {
        contiguous_space = _contiguous_space;
        std::cout << "Init MHash" << std::endl; fflush(stdout);
        rr = r;
        max_lvl = _max_level;
        using W = typename Graph::weight_type;
        
        mtable.initialize(0, 0, max_lvl);
        mtable.allocate(DG.n, 0, r-1, 0);
  
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

        std::cout << "End MHash Count" << std::endl; fflush(stdout);

        long total = mtable.set_table_sizes();
        if (contiguous_space) {
          space = pbbslib::new_array_no_init<X>(total);
          mtable.set_end_table_rec(space);
        } else {
          mtable.set_end_table_rec();
        }

        std::cout << "End MHash" << std::endl; fflush(stdout);
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

      std::size_t return_total() { return mtable.total_size; }

      long get_count(std::size_t index) {
        if (contiguous_space) return std::get<1>(space[index]);

        long count = 0;
        auto func = [&](std::tuple<unsigned __int128, long>* loc){
          count = std::get<1>(*loc);
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
        auto func = [&](std::tuple<unsigned __int128, long>* loc){
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
        auto func = [&](std::tuple<unsigned __int128, long>* loc){
          *loc = std::make_tuple(std::get<0>(*loc), 0);
        };
        mtable.find_table_loc(index, func);
      }

      template<class I>
      void extract_indices(sequence<uintE>& base2, I func, int r, int k) {
        uintE base[10];
        assert(10 > k);
        for(std::size_t i = 0; i < k + 1; i++) {
          base[i] = base2[i];
        }
        std::sort(base, base + k + 1,sort_func);
        mtable.extract_indices(base, func, r, k);
      }

      //Fill base[k] ... base[k-r+1] and base[0]
      template<class S, class Graph>
      void extract_clique(S index, sequence<uintE>& base, Graph& G, int k) {
        mtable.extract_clique(index, base, 0, rr, k);
      }
  };

} // end namespace multitable

} //end namespace gbbs