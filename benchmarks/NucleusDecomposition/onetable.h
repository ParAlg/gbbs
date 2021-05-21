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

namespace onetable {

  template <class Y, class H>
  class OnelevelHash {
    public:
      using T = pbbslib::sparse_table<Y, long, H>;
      T table;
      int rr;
      int shift_factor;

      template <class Graph>
      OnelevelHash(int r, Graph& DG, size_t max_deg, int _shift_factor) {
        shift_factor = _shift_factor;
        rr = r;
        // count cliques
        timer t_pre; t_pre.start();
        size_t pre_count = 0;
        // Clique counting
        if (r == 2) pre_count = DG.m;
        else if (r == 3) pre_count = TriClique_count(DG, false, nullptr);
        else pre_count = Clique_count(DG, r, 5, true, false, 0, nullptr);
        //std::cout << "Pre count " << r << ": " << pre_count << std::endl;
        double tt_pre = t_pre.stop();
        //std::cout << "### Pre count: " << tt_pre << std::endl;

        //std::cout << "Start table" << std::endl;
        table = pbbslib::sparse_table<Y, long, H>(
          pre_count,
          std::make_tuple(std::numeric_limits<Y>::max(), long{0}), H{});
        
        size_t data_structure_size = sizeof(*this) + table.m * sizeof(std::tuple<Y, long>);
        std::cout << "Data Structure Size: " << data_structure_size << std::endl;
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
          Y key = 0;
          for (int i = 0; i < static_cast<int>(k)+1; ++i) {
            if (bitmask[i]) {
              key = key << shift_factor;
              key |= static_cast<uintE>(base[i]);
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

/*auto index = table.find_index(key);
assert(index < table.m);
sequence<uintE> base3 = sequence<uintE>(k + 1, [](std::size_t i){return UINT_E_MAX;});
extract_clique(index, base3, k, k);
int base3_idx = 0;
for (int i = 0; i < static_cast<int>(k)+1; ++i) {
  if (bitmask[i]) {
    if (base[i] != base3[base3_idx]) {
      std::cout << "base3idx: " << base3_idx << ", base3: " << base3[base3_idx] << std::endl;
      std::cout << "base: " << base[i] << std::endl;
      fflush(stdout);
    }
    assert(base[i] == base3[base3_idx]);
    base3_idx++;
    if (base3_idx == 1) base3_idx = k - rr + 2;
  }
}*/

        } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
      }

      std::size_t return_total() {return table.m;}
      long get_count(std::size_t i) {
        if (std::get<0>((table.table)[i]) == std::numeric_limits<Y>::max()) return UINT_E_MAX;
        return std::get<1>((table.table)[i]);
      }
      size_t update_count(std::size_t i, size_t update){
        /*if (std::get<0>((table.table)[i]) == std::numeric_limits<Y>::max()) return 0
        }*/
        if (get_count(i) < update) {
          std::cout << "i: " << i << ", count: " << get_count(i) << ", update: " << update << std::endl;
          fflush(stdout);
          exit(0);
        }
        auto val = std::get<1>(table.table[i]) - update;
        table.table[i] = std::make_tuple(std::get<0>(table.table[i]), val);
        return val;
      }
      void clear_count(std::size_t index) {
        table.table[index] = std::make_tuple(std::get<0>(table.table[index]),0);
      }

      void set_count(std::size_t index, size_t update) {
        table.table[index] = std::make_tuple(std::get<0>(table.table[index]),update);
      }

      Y extract_indices_check(sequence<uintE>& base2, int r) {
        // Size of base2 should be r + 1
        uintE base[10];
        assert(10 > r + 1);
        for(std::size_t i = 0; i < r + 1; i++) {
          base[i] = base2[i];
        }
        std::sort(base, base + r + 1,std::less<uintE>());

        Y key = 0;
        for (int i = 0; i < static_cast<int>(r)+1; ++i) {
          key = key << shift_factor;
          key |= static_cast<uintE>(base[i]);
        }
        auto index = table.find_index(key);
        assert(index < table.m);
        return index;
      }

      template<class HH, class HG, class I>
      void extract_indices(sequence<uintE>& base2, HH is_active, HG is_inactive, I func, int r, int k) {
        // Sort base
        // Sort base
        uintE base[10];
        assert(10 > k);
        for(std::size_t i = 0; i < k + 1; i++) {
          base[i] = base2[i];
        }
        std::sort(base, base + k + 1,std::less<uintE>());

        std::string bitmask(r+1, 1); // K leading 1's
        bitmask.resize(k+1, 0); // N-K trailing 0's

        std::vector<size_t> indices;
        size_t num_active = 0;
        bool use_func = true;

        do {
          Y key = 0;
          for (int i = 0; i < static_cast<int>(k)+1; ++i) {
            if (bitmask[i]) {
              key = key << shift_factor;
              key |= static_cast<uintE(base[i]);
            }
          }
          auto index = table.find_index(key);
          assert(index < table.m);
 
 
sequence<uintE> base3 = sequence<uintE>(k + 1, [](std::size_t i){return UINT_E_MAX;});
extract_clique(index, base3, k, k);
int base3_idx = 0;
for (int i = 0; i < static_cast<int>(k)+1; ++i) {
  if (bitmask[i]) {
    if (base[i] != base3[base3_idx]) {
      std::cout << "base3idx: " << base3_idx << ", base3: " << base3[base3_idx] << std::endl;
      std::cout << "base: " << base[i] << std::endl;
      fflush(stdout);
    }
    assert(base[i] == base3[base3_idx]);
    base3_idx++;
    if (base3_idx == 1) base3_idx = k - rr + 2;
  }
}


          indices.push_back(index);
          if (is_active(index)) num_active++;
          if (is_inactive(index)) use_func = false;
          //func(index);
        } while (std::prev_permutation(bitmask.begin(), bitmask.end()));

        assert(num_active != 0);
        if (use_func) {
          for (std::size_t i = 0; i < indices.size(); i++) {
            if (!is_active(indices[i]) && !is_inactive(indices[i]))
              func(indices[i], 1.0 / (double) num_active);
          }
        }
      }
    
    template<class S, class Graph>
    void extract_clique(S index, sequence<uintE>& base, Graph& G, int k) {
      auto vert = std::get<0>(table.table[index]);
      for (int j = 0; j < rr; ++j) {
        unsigned __int128 mask = (1ULL << (shift_factor + 1)) - 1;
        uintE extract = (uintE) vert & mask; // vert & mask
        /*if (static_cast<uintE>(extract) >= G.n) {
          std::cout << "Vert: " << static_cast<uintE>(extract) << ", n: " << G.n << std::endl;
        }*/
        //assert(static_cast<uintE>(extract) < G.n);
        if (j == rr - 1) base[0] = static_cast<uintE>(extract);
        else base[k - j] = static_cast<uintE>(extract);
        vert = vert >> shift_factor;
      }
    }
  };

} // end namespace onetable

} //end namespace gbbs