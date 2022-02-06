#pragma once

#include <math.h>
#include <limits>

// Library dependencies
#include "gbbs/bucket.h"
#include "gbbs/edge_map_reduce.h"
#include "gbbs/gbbs.h"
#include "gbbs/helpers/dyn_arr.h"
#include "gbbs/helpers/sparse_table.h"
#include "gbbs/helpers/sparse_additive_map.h"

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

  template <class Y, class H, class C>
  class OnelevelHash {
    public:
      using T = gbbs::sparse_table<Y, C, H>;
      T table;
      int rr;
      int shift_factor;

      template <class Graph>
      OnelevelHash(int r, Graph& DG, size_t max_deg, int _shift_factor) {
        shift_factor = _shift_factor;
        rr = r;
        // count cliques
        //timer t_pre; t_pre.start();
        size_t pre_count = 0;
        // Clique counting
        if (r == 2) pre_count = DG.m;
        else if (r == 3) pre_count = TriClique_count(DG, false, nullptr);
        else pre_count = Clique_count(DG, r, 5, true, false, 0, nullptr);
        //std::cout << "Pre count " << r << ": " << pre_count << std::endl;
        //double tt_pre = t_pre.stop();
        //std::cout << "### Pre count: " << tt_pre << std::endl;

        //std::cout << "Start table" << std::endl;
        table = gbbs::sparse_table<Y, C, H>(
          pre_count,
          std::make_tuple(std::numeric_limits<Y>::max(), C{0}), H{});
        
        size_t data_structure_size = sizeof(*this) + table.m * sizeof(std::tuple<Y, C>);
        std::cout << "Data Structure Size: " << data_structure_size << std::endl;
      }

      void insert(sequence<uintE>& base2, int r, int k) {
        auto add_f = [&] (C* ct, const std::tuple<Y, C>& tup) {
          gbbs::fetch_and_add(ct, (C)1);
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
          unsigned __int128 mask = (1ULL << (shift_factor)) - 1;
          for (int i = 0; i < static_cast<int>(k)+1; ++i) {
            if (bitmask[i]) {
              key = key << shift_factor;
              assert((base[i] & mask) == base[i]);
              key |= (base[i] & mask);
            }
          }
          table.insert_f(std::make_tuple(key, (C) 1), add_f);
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

      void insert_twothree(uintE v1, uintE v2, uintE v3, int r, int k) {
        auto add_f = [&] (C* ct, const std::tuple<Y, C>& tup) {
          gbbs::fetch_and_add(ct, (C)1);
        };
        unsigned __int128 mask = (1ULL << (shift_factor)) - 1;

        Y key12 = std::min(v1, v2) & mask;
        key12  = key12 << shift_factor;
        key12 |= (std::max(v1, v2) & mask);
        table.insert_f(std::make_tuple(key12, (C) 1), add_f);

        Y key13 = std::min(v1, v3) & mask;
        key13 = key13 << shift_factor;
        key13 |= (std::max(v1, v3) & mask);
        table.insert_f(std::make_tuple(key13, (C) 1), add_f);

        Y key23 = std::min(v2, v3) & mask;
        key23 = key23 << shift_factor;
        key23 |= (std::max(v2, v3) & mask);
        table.insert_f(std::make_tuple(key23, (C) 1), add_f);
      }

      std::size_t return_total() {return table.m;}
      C get_count(std::size_t i) {
        if (std::get<0>((table.table)[i]) == std::numeric_limits<Y>::max()) return 0;
        return std::get<1>((table.table)[i]);
      }
      C update_count(std::size_t i, C update){
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

      void set_count(std::size_t index, C update) {
        table.table[index] = std::make_tuple(std::get<0>(table.table[index]),update);
      }

      C update_count_atomic(std::size_t index, C update) {
        if (get_count(index) < update) {
          std::cout << "i: " << index << ", count: " << get_count(index) << ", update: " << update << std::endl;
          fflush(stdout);
          exit(0);
        }
        gbbs::write_add(&std::get<1>(table.table[index]), -1 * update);
        return std::get<1>(table.table[index]);
      }

      Y extract_indices_check(uintE* base2, int r) {
        // Size of base2 should be r + 1
        uintE base[10];
        assert(10 > r + 1);
        for(std::size_t i = 0; i < r + 1; i++) {
          base[i] = base2[i];
        }
        std::sort(base, base + r + 1,std::less<uintE>());

        Y key = 0;
        unsigned __int128 mask = (1ULL << (shift_factor)) - 1;
        for (int i = 0; i < static_cast<int>(r)+1; ++i) {
          key = key << shift_factor;
          assert((base[i] & mask) == base[i]);
          key |= (base[i] & mask);
        }
        auto index = table.find_index(key);
        assert(index < table.m);
        return index;
      }

      template<class HH, class HG, class I>
      void extract_indices(uintE* base2, HH is_active, HG is_inactive, I func, int r, int k) {
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
          unsigned __int128 mask = (1ULL << (shift_factor)) - 1;
          for (int i = 0; i < static_cast<int>(k)+1; ++i) {
            if (bitmask[i]) {
              key = key << shift_factor;
              key |= (base[i] & mask);
            }
          }
          auto index = table.find_index(key);
          assert(index < table.m);
 
 /*
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


          indices.push_back(index);
          if (is_active(index)) num_active++;
          if (is_inactive(index)) {
            use_func = false;
            return; //TODO: check that this is ok
          }
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

      Y extract_indices_two(uintE v1, uintE v3) {
        unsigned __int128 mask = (1ULL << (shift_factor)) - 1;
        //size_t num_active = 1;
        // Assume v1, v2 is the active edge
        // Need to get indices for v1, v3 and v2, v3
        Y key13 = (std::min(v1, v3) & mask);
        key13 = key13 << shift_factor;
        key13 |= (std::max(v1, v3) & mask);
        return table.find_index(key13);
      }

      template<class HH, class HG, class I>
      void extract_indices_twothree(uintE v1, uintE v2, uintE v3, HH is_active, 
        HG is_inactive, I func, int r, int k) {
        unsigned __int128 mask = (1ULL << (shift_factor)) - 1;
        //size_t num_active = 1;
        // Assume v1, v2 is the active edge
        // Need to get indices for v1, v3 and v2, v3
        Y key13 = (std::min(v1, v3) & mask);
        Y key23 = (std::min(v2, v3) & mask);
        key13 = key13 << shift_factor;
        key23 = key23 << shift_factor;
        key13 |= (std::max(v1, v3) & mask);
        key23 |= (std::max(v2, v3) & mask);
        auto index13 = table.find_index(key13);
        if (is_inactive(index13)) return;
        auto index23 = table.find_index(key23);
        if (is_inactive(index23)) return;

        //if (is_active(index13)) num_active++;
        //if (is_active(index23)) num_active++;
        //if (!is_active(index13)) func(index13, 1.0 / (double) num_active);
        //if (!is_active(index23)) func(index23, 1.0 / (double) num_active);
        bool is_active_13 = is_active(index13);
        bool is_active_23 = is_active(index23);
        if (is_active_13 && is_active_23) return;
        if (is_active_13 || is_active_23) {
          Y key12 = (std::min(v1, v2) & mask);
          key12 = key12 << shift_factor;
          key12 |= (std::max(v1, v2) & mask);
          auto index12 = table.find_index(key12);
          if (is_active_13 && index12 < index13) func(index23, 1.0);
          else if (is_active_23 && index12 < index23) func(index13, 1.0);
          return;
        }
        func(index23, 1.0);
        func(index13, 1.0);
      }

      // Assume v1, v2, v3 is active
      template<class HH, class HG, class I>
      void extract_indices_threefour(uintE v1, uintE v2, uintE v3, uintE v4, 
        HH is_active, HG is_inactive, I func, int r, int k) {
        unsigned __int128 mask = (1ULL << (shift_factor)) - 1;
        size_t num_active = 1;
        // Assume v1, v2 is the active edge
        // Need to get indices for v1,v2,v4..v1,v3,v4..v2,v3,v4
        Y key124 = (std::min({v1, v2, v4}) & mask);
        key124 = key124 << shift_factor;
        key124 |= (middleOfThree(v1, v2, v4) & mask);
        key124 = key124 << shift_factor;
        key124 |= (std::max({v1, v2, v4}) & mask);
        auto index124 = table.find_index(key124);
        if (is_inactive(index124)) return;

        Y key134 = (std::min({v1, v3, v4}) & mask);
        key134 = key134 << shift_factor;
        key134 |= (middleOfThree(v1, v3, v4) & mask);
        key134 = key134 << shift_factor;
        key134 |= (std::max({v1, v3, v4}) & mask);
        auto index134 = table.find_index(key134);
        if (is_inactive(index134)) return;

        Y key234 = (std::min({v2, v3, v4}) & mask);
        key234 = key234 << shift_factor;
        key234 |= (middleOfThree(v2, v3, v4) & mask);
        key234 = key234 << shift_factor;
        key234 |= (std::max({v2, v3, v4}) & mask);
        auto index234 = table.find_index(key234);
        if (is_inactive(index234)) return;

        bool is_active_124 = is_active(index124);
        bool is_active_134 = is_active(index134);
        bool is_active_234 = is_active(index234);
        if (is_active_124 && is_active_134 && is_active_234) return;
        if (is_active_124 || is_active_134 || is_active_234) {
          Y key123 = (std::min({v1, v2, v3}) & mask);
          key123 = key123 << shift_factor;
          key123 |= (middleOfThree(v1, v2, v3) & mask);
          key123 = key123 << shift_factor;
          key123 |= (std::max({v1, v2, v3}) & mask);
          auto index123 = table.find_index(key123);
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
      }
    
    template<class S, class Graph>
    void extract_clique(S index, uintE* base, Graph& G, int k) {
      auto vert = std::get<0>(table.table[index]);
      for (int j = 0; j < rr; ++j) {
        unsigned __int128 mask = (1ULL << (shift_factor)) - 1;
        uintE extract = (uintE) (vert & mask); // vert & mask
        /*if (static_cast<uintE>(extract) >= G.n) {
          std::cout << "Vert: " << static_cast<uintE>(extract) << ", n: " << G.n << std::endl;
        }*/
        //assert(static_cast<uintE>(extract) < G.n);
        if (j == rr - 1) base[0] = static_cast<uintE>(extract);
        else base[k - j] = static_cast<uintE>(extract);
        vert = vert >> shift_factor;
      }
    }

    template<class S>
    std::tuple<uintE, uintE> extract_clique_two(S index, int k) {
      auto vert = std::get<0>(table.table[index]);
      uintE v1 = UINT_E_MAX;
      uintE v2 = UINT_E_MAX;
      for (int j = 0; j < rr; ++j) {
        unsigned __int128 mask = (1ULL << (shift_factor)) - 1;
        uintE extract = (uintE) (vert & mask); // vert & mask
        if (j == rr - 1) v1 = static_cast<uintE>(extract);
        else v2 = static_cast<uintE>(extract);
        vert = vert >> shift_factor;
      }
      assert(v1 != UINT_E_MAX); assert(v2 != UINT_E_MAX);
      return std::tuple(v1, v2);
    }

    template<class S>
    std::tuple<uintE, uintE, uintE> extract_clique_three(S index, int k) {
      auto vert = std::get<0>(table.table[index]);
      uintE v1 = UINT_E_MAX;
      uintE v2 = UINT_E_MAX;
      uintE v3 = UINT_E_MAX;
      for (int j = 0; j < rr; ++j) {
        unsigned __int128 mask = (1ULL << (shift_factor)) - 1;
        uintE extract = (uintE) (vert & mask); // vert & mask
        if (j == rr - 1) v1 = static_cast<uintE>(extract);
        else if (j == 0) v2 = static_cast<uintE>(extract);
        else v3 = static_cast<uintE>(extract);
        vert = vert >> shift_factor;
      }
      assert(v1 != UINT_E_MAX); assert(v2 != UINT_E_MAX); assert(v3 != UINT_E_MAX);
      return std::tuple(v1, v2, v3);
    }
  };

} // end namespace onetable

} //end namespace gbbs