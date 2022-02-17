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

namespace twotable_nosearch {

  template <class Y, class H, class C>
  struct EndTable {
    gbbs::sparse_table<Y, C, H> table;
    uintE vtx;
    //MidTable* up_table;
  };

  template <class Y, class H, class C>
  struct MidTable {
    using EndTableY = EndTable<Y, H, C>;
    gbbs::sparse_table<uintE, EndTableY*, std::hash<uintE>> table;
    sequence<EndTableY*> arr;
  };

  template<class Y>
  bool is_max_val(Y max_val) {
    std::size_t max_bit = sizeof(Y) * 8;
    Y check_bit = (max_val >> (max_bit - 1)) & 1U;
    return (check_bit != 0);
  }

  template<class Y, class C, class S, class EndSpace>
  uintE get_mtable(S index, EndSpace* end_space) {
    using X = std::tuple<Y, C>;
    while (true) {
      auto max_val = std::get<0>(static_cast<X>(end_space[index]));
      std::size_t max_bit = sizeof(Y) * 8;
      Y one = 1;
      Y check_bit = (max_val >> (max_bit - 1)) & 1U;
      if (check_bit != 0) {
        max_val ^= (one << (max_bit - 1));
        //max_val &= ~(1UL << (max_bit - 1));
        return static_cast<uintE>(max_val);
      }
      index++;
      //index = index % mtable.m;
    }
  }
  
  template <class Y, class H, class C>
  class TwolevelHash {
    public:
      using T = gbbs::sparse_table<Y, C, H>;
      using X = std::tuple<Y, C>;
      using EndTableY = EndTable<Y, H, C>;
      using MidTableY = MidTable<Y, H, C>;
      MidTableY top_table;
      sequence<C> top_table_sizes;
      int rr;
      std::size_t total = 0;
      X* space = nullptr;
      int shift_factor;
      size_t nx;
  
      template <class Graph>
      TwolevelHash(int r, Graph& DG, size_t max_deg, bool relabel, int _shift_factor) {
        using W = typename Graph::weight_type;
        shift_factor = _shift_factor;
        nx = DG.n;
        rr = r;
        //top_table.up_table = nullptr;
        // How many vert in top level?
        // For each vert in top level, how many pairs of vert finish it?
        auto tmp_table = gbbs::sparse_additive_map<uintE, C>(
          DG.n, std::make_tuple(UINT_E_MAX, C{0}));
        auto base_f = [&](sequence<uintE>& base){
          auto min_vert = relabel ? base[0] : parlay::reduce_min(base);
          auto tmp = std::make_tuple<uintE, C>(static_cast<uintE>(min_vert), C{1});
          tmp_table.insert(tmp);
        };
        if (r == 2) {
          parallel_for(0, DG.n, [&](std::size_t i){
            if (DG.get_vertex(i).out_degree() != 0) {
              auto map_f = [&](const uintE& src, const uintE& ngh, const W& wgh) {
                auto base = sequence<uintE>(r);
                base[0] = i;
                base[1] = ngh;
                base_f(base);
              };
              DG.get_vertex(i).out_neighbors().map(map_f, false);
            }
          });
        } else {
          auto init_induced = [&](HybridSpace_lw* induced) { induced->alloc(max_deg, r-1, DG.n, true, true); };
          auto finish_induced = [&](HybridSpace_lw* induced) { if (induced != nullptr) { delete induced; } };
          parallel_for_alloc<HybridSpace_lw>(init_induced, finish_induced, 0, DG.n, [&](size_t i, HybridSpace_lw* induced) {
            if (DG.get_vertex(i).out_degree() != 0) {
              induced->setup(DG, r-1, i);
              auto base = sequence<uintE>(r);
              base[0] = i;
              NKCliqueDir_fast_hybrid_rec(DG, 1, r-1, induced, base_f, base);
            }
          }, 1, false);
        }
        auto top_table_sizes2 = tmp_table.entries();
        // sort by key
        parlay::sample_sort_inplace (make_slice(top_table_sizes2), [&](const std::tuple<uintE, C>& u, const std::tuple<uintE, long>&  v) {
          return std::get<0>(u) < std::get<0>(v);
        });
        sequence<long> actual_sizes(top_table_sizes2.size() + 1);
        // Modify top_table_sizes2 to be appropriately oversized
        parallel_for(0, top_table_sizes2.size(), [&](std::size_t i) {
          auto m = 1 + ((size_t)1 << parlay::log2_up((size_t)(1.1 * std::get<1>(top_table_sizes2[i])) + 1));
          actual_sizes[i] = m;
        });
        actual_sizes[top_table_sizes2.size()] = 0;
        // Do a scan inplace
        //auto add_tuple_monoid = pbbs::make_monoid([](std::tuple<uintE, long> a, std::tuple<uintE, long> b){
        //  return std::make_tuple(std::get<0>(b), std::get<1>(a) + std::get<1>(b));
        //}, std::make_tuple(UINT_E_MAX, 0));
        long total_top_table_sizes2 = parlay::scan_inplace(make_slice(actual_sizes));
        // Allocate space for the second level tables
        space = gbbs::new_array_no_init<X>(total_top_table_sizes2);
        tmp_table.del();
  
        //*** for arr
        top_table.arr = sequence<EndTableY*>::from_function(DG.n, [](std::size_t i){return nullptr;});
        top_table_sizes = sequence<C>(DG.n + 1, C{0});
        /*top_table.table = pbbslib::sparse_table<uintE, EndTable*, std::hash<uintE>>(
          top_table_sizes2.size(),
          std::make_tuple<uintE, EndTable*>(UINT_E_MAX, static_cast<EndTable*>(nullptr)),
          std::hash<uintE>());
        top_table_sizes = sequence<long>(top_table.table.m + 1, long{0});*/
  
        parallel_for(0, top_table_sizes2.size(), [&](std::size_t i){
          auto vtx = std::get<0>(top_table_sizes2[i]);
          //if (i != top_table_sizes2.size() - 1) assert(vtx < std::get<0>(top_table_sizes2[i + 1]));
          auto upper_size = actual_sizes[i + 1];
          auto size = upper_size - actual_sizes[i];
          EndTableY* end_table = new EndTableY();

          Y max_val = static_cast<Y>(vtx); 
          assert((max_val >> (max_bit - 1)) & 1U == 0U);
          std::size_t max_bit = sizeof(Y) * 8;
          Y one = 1;
          max_val |= (one << (max_bit - 1));

          end_table->vtx = vtx;
          end_table->table = gbbs::sparse_table<Y, C, H>(
            size - 1, 
            std::make_tuple<Y, C>(static_cast<Y>(max_val), static_cast<C>(0)),
            H{},
            space + actual_sizes[i]
            );
          space[actual_sizes[i] + size - 1] = std::make_tuple<Y, C>(static_cast<Y>(max_val), static_cast<C>(0));

          //uintE vtest = get_mtable<Y>(actual_sizes[i], space);
          //vtest = get_mtable<Y>(actual_sizes[i] + size - 1, space);
          //assert(vtest == vtx);
          /*top_table.table.insert(std::make_tuple(vtx, end_table));
          std::size_t l = top_table.table.find_index(vtx);
          top_table_sizes[l] = end_table->table.m;*/
          //***for arr
          top_table.arr[vtx] = end_table;
          //assert(size == 1 + end_table->table.m);
          top_table_sizes[vtx] = size; //1 + end_table->table.m;
          //assert((end_table->table).table == space + actual_sizes[i]);
        });
        total = parlay::scan_inplace(make_slice(top_table_sizes));
        //assert(total == total_top_table_sizes2);
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
        size_t data_structure_size = sizeof(*this) + sizeof(MidTableY) + total * sizeof(X) +  (DG.n) * sizeof(EndTableY*) + (DG.n + 1) * sizeof(C);
        std::cout << "Data Structure Size: " << data_structure_size << std::endl;
      }

      void insert_twothree(uintE v1, uintE v2, uintE v3, int r, int k) {
        auto add_f = [&] (C* ct, const std::tuple<Y, C>& tup) {
          gbbs::fetch_and_add(ct, (C)1);
        };
        EndTableY* end_table12 = top_table.arr[std::min(v1, v2)];
        (end_table12->table).insert_f(std::make_tuple(Y{std::max(v1, v2)}, (C) 1), add_f);

        EndTableY* end_table13 = top_table.arr[std::min(v1, v3)];
        (end_table13->table).insert_f(std::make_tuple(Y{std::max(v1, v3)}, (C) 1), add_f);

        EndTableY* end_table23 = top_table.arr[std::min(v2, v3)];
        (end_table23->table).insert_f(std::make_tuple(Y{std::max(v2, v3)}, (C) 1), add_f);
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
        unsigned __int128 mask = (1ULL << (shift_factor)) - 1;

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
                key |= (base[i] & mask);
              }
            }
          }
          //EndTable* end_table = top_table.table.find(vtx, nullptr);
          //***for arr
          EndTableY* end_table = top_table.arr[vtx];
          //assert(end_table != nullptr);
          (end_table->table).insert_f(std::make_tuple(key, (C) 1), add_f);

          /*auto index2 = (end_table->table).find_index(key);
          uintE vtest1 = get_mtable<Y>(index2, (end_table->table).table);
          assert(vtest1 == vtx);*/

          //EndTableY* end_table2 = top_table.arr[vtx];
          //auto prefix = top_table_sizes[vtx];
          /*if (end_table != nullptr) {
            assert(end_table->vtx == vtx);
            assert((end_table->table).m + 1 == top_table_sizes[vtx + 1] - top_table_sizes[vtx]);
            assert((end_table->table).table == space + prefix);
          }*/
          //auto index = (end_table->table).find_index(key);
          /*uintE vtest = get_mtable<Y>(index + prefix, space);
          assert(vtest == vtx);*/


        } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
      }

      std::size_t return_total() {return total;}

      size_t get_top_index(std::size_t index) {
        // This gives the first i such that top_table_sizes[i] >= index
        auto idx = parlay::binary_search(top_table_sizes, C{index}, std::less<C>());
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

      C get_count(std::size_t index) {
        if (is_max_val<Y>(std::get<0>(space[index]))) return 0;
        return std::get<1>(space[index]);
      }

      C update_count(std::size_t index, C update){
        if (is_max_val<Y>(std::get<0>(space[index]))) return 0;
        if (get_count(index) < update) {
          std::cout << "i: " << index << ", count: " << get_count(index) << ", update twotablenosearch: " << update << std::endl;
          fflush(stdout);
          exit(0);
        }

        auto val = std::get<1>(space[index]) - update;
        space[index] =
          std::make_tuple(std::get<0>(space[index]), val);
        return val;
      }

      C update_count_atomic(std::size_t index, C update){
        if (is_max_val<Y>(std::get<0>(space[index]))) return 0;
        if (get_count(index) < update) {
          std::cout << "i: " << index << ", count: " << get_count(index) << ", update twotablenosearch: " << update << std::endl;
          fflush(stdout);
          exit(0);
        }

        gbbs::write_add(&std::get<1>(space[index]), -1 * update);
        return std::get<1>(space[index]);
      }

      void clear_count(std::size_t index) {
        space[index] = std::make_tuple(std::get<0>(space[index]), 0);
      }

      void set_count(std::size_t index, C update) {
        space[index] = std::make_tuple(std::get<0>(space[index]), update);
      }

      Y extract_indices_check(uintE* base2, int r) {
        // Size of base2 should be r + 1
        uintE base[10];
        assert(10 > r + 1);
        for(std::size_t i = 0; i < r + 1; i++) {
          base[i] = base2[i];
        }
        std::sort(base, base + r + 1,std::less<uintE>());

        bool use_vtx = false;
        uintE vtx = 0;
        Y key = 0;
        unsigned __int128 mask = (1ULL << (shift_factor)) - 1;
        for (int i = 0; i < static_cast<int>(r)+1; ++i) {
          if (!use_vtx) {
            use_vtx = true;
            vtx = base[i];
          } else {
            key = key << shift_factor;
            key |= (base[i] & mask);
          }
        }
        EndTableY* end_table = top_table.arr[vtx];
        auto prefix = top_table_sizes[vtx];
        auto index = (end_table->table).find_index(key);
        return prefix + index;
      }

      Y extract_indices_two(uintE v1, uintE v3) {
        unsigned __int128 mask = (1ULL << (shift_factor)) - 1;
        // Assume v1, v2 is the active edge
        // Need to get indices for v1, v3 and v2, v3
        auto vtx13 = std::min(v1, v3);
        Y key13 = (std::max(v1, v3) & mask);
  
        EndTableY* end_table13 = top_table.arr[vtx13];
        return top_table_sizes[vtx13] + (end_table13->table).find_index(key13);
      }

      template<class HH, class HG, class I>
      void extract_indices_twothree(uintE v1, uintE v2, uintE v3, HH is_active, 
        HG is_inactive, I func, int r, int k) {
        unsigned __int128 mask = (1ULL << (shift_factor)) - 1;
        size_t num_active = 1;
        // Assume v1, v2 is the active edge
        // Need to get indices for v1, v3 and v2, v3
        auto vtx13 = std::min(v1, v3);
        auto vtx23 = std::min(v2, v3);
        Y key13 = (std::max(v1, v3) & mask);
        Y key23 = (std::max(v2, v3) & mask);
  
        EndTableY* end_table13 = top_table.arr[vtx13];
        auto index13 = top_table_sizes[vtx13] + (end_table13->table).find_index(key13);
        if (is_inactive(index13)) return;

        EndTableY* end_table23 = top_table.arr[vtx23];
        auto index23 = top_table_sizes[vtx23] + (end_table23->table).find_index(key23);
        if (is_inactive(index23)) return;

        bool is_active_13 = is_active(index13);
        bool is_active_23 = is_active(index23);
        if (is_active_13 && is_active_23) return;
        if (is_active_13 || is_active_23) {
          auto vtx12 = std::min(v1, v2);
          Y key12 = (std::max(v1, v2) & mask);
          EndTableY* end_table12 = top_table.arr[vtx12];
          auto index12 = top_table_sizes[vtx12] + (end_table12->table).find_index(key12);
          if (is_active_13 && index12 < index13) func(index23, 1.0);
          else if (is_active_23 && index12 < index23) func(index13, 1.0);
          return;
        }
        func(index23, 1.0);
        func(index13, 1.0);
      }

      template<class HH, class HG, class I>
      void extract_indices_threefour(uintE v1, uintE v2, uintE v3, uintE v4, HH is_active, 
        HG is_inactive, I func, int r, int k) {
        unsigned __int128 mask = (1ULL << (shift_factor)) - 1;
        size_t num_active = 1;
        // Assume v1, v2 is the active edge
        // Need to get indices for v1, v3 and v2, v3
        auto vtx124 = std::min({v1, v2, v4});
        Y key124 = middleOfThree(v1, v2, v4) & mask;
        key124 = key124 << shift_factor;
        key124 |= (std::max({v1, v2, v4}) & mask);

        EndTableY* end_table124 = top_table.arr[vtx124];
        auto index124 = top_table_sizes[vtx124] + (end_table124->table).find_index(key124);
        if (is_inactive(index124)) return;

        auto vtx134 = std::min({v1, v3, v4});
        Y key134 = middleOfThree(v1, v3, v4) & mask;
        key134 = key134 << shift_factor;
        key134 |= (std::max({v1, v3, v4}) & mask);

        EndTableY* end_table134 = top_table.arr[vtx134];
        auto index134 = top_table_sizes[vtx134] + (end_table134->table).find_index(key134);
        if (is_inactive(index134)) return;

        auto vtx234 = std::min({v2, v3, v4});
        Y key234 = middleOfThree(v2, v3, v4) & mask;
        key234 = key234 << shift_factor;
        key234 |= (std::max({v2, v3, v4}) & mask);

        EndTableY* end_table234 = top_table.arr[vtx234];
        auto index234 = top_table_sizes[vtx234] + (end_table234->table).find_index(key234);
        if (is_inactive(index234)) return;

        bool is_active_124 = is_active(index124);
        bool is_active_134 = is_active(index134);
        bool is_active_234 = is_active(index234);
        if (is_active_124 && is_active_134 && is_active_234) return;
        if (is_active_124 || is_active_134 || is_active_234) {
          auto vtx123 = std::min({v1, v2, v3});
          Y key123 = middleOfThree(v1, v2, v3) & mask;
          key123 = key123 << shift_factor;
          key123 |= (std::max({v1, v2, v3}) & mask);
          EndTableY* end_table123 = top_table.arr[vtx123];
          auto index123 = top_table_sizes[vtx123] + (end_table123->table).find_index(key123);
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

      template<class HH, class HG, class I>
      void extract_indices(uintE* base2, HH is_active, HG is_inactive, I func, int r, int k, Y xxx = 0) {
        /*if (std::numeric_limits<Y>::max() == 0) {
          std::cout << "why is this 0" << std::endl; fflush(stdout);
          exit(0);
        }*/
        assert(xxx == 0 || is_active(xxx));
        /*if (xxx != 0) {
          if (!is_active(xxx)) {
            std::cout << "xxx should be active" << std::endl; fflush(stdout); exit(0);
          }
        }*/
        // Sort base
        uintE base[10];
        assert(10 > k);
        for(std::size_t i = 0; i < k + 1; i++) {
          base[i] = base2[i];
        }
        std::sort(base, base + k + 1,std::less<uintE>());

        std::vector<size_t> indices;
        size_t num_active = 0;
        __uint128_t min_active = __uint128_t(__int128_t(-1L));
        bool use_func = true;
        bool one_should_be_xxx = false;

        std::string bitmask(r+1, 1); // K leading 1's
        bitmask.resize(k+1, 0); // N-K trailing 0's
        unsigned __int128 mask = (1ULL << (shift_factor)) - 1;
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
                key |= (base[i] & mask);
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

          indices.push_back(prefix + index);
          if (is_active(prefix + index)) {
            num_active++;
            if (prefix + index < min_active) min_active = prefix + index;
            if (prefix + index == xxx) one_should_be_xxx = true;
            assert(prefix + index <= min_active && min_active != __uint128_t(__int128_t(-1L)));
            /*if (prefix + index > min_active && min_active == __uint128_t(__int128_t(-1L))){
              std::cout << "we have a min problem" << std::endl; fflush(stdout);
              exit(0);
            }*/
            /*if (prefix + index > std::numeric_limits<Y>::max()) {
              std::cout << "greater than max??" << std::endl; fflush(stdout);
              exit(0);
            }*/
          }
          if (is_inactive(prefix + index)) return;
          //func(prefix + index);
        } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
        
        assert(xxx == 0 || one_should_be_xxx);
        /*if (xxx != 0 && !one_should_be_xxx) {
          std::cout << "one is not xxx" << std::endl; fflush(stdout);
          exit(0);
        }*/
        assert(num_active != 1 || xxx == 0 || xxx == min_active);
        /*
        if (num_active == 1 && xxx != 0){
          if (xxx != min_active) {
            std::cout << "only one active so xxx should be min" << std::endl;
            fflush(stdout); exit(0);
          }
        }*/

        //assert(num_active != 0);
        if (use_func && (xxx == 0 || num_active == 1)){// && matches_base == min_active) {
          for (std::size_t i = 0; i < indices.size(); i++) {
            if (!is_active(indices[i]) && !is_inactive(indices[i]))
              func(indices[i], 1.0 / (double) num_active);
              //func(indices[i], 1);
          }
        } else if (use_func && xxx == min_active) {
          for (std::size_t i = 0; i < indices.size(); i++) {
            if (!is_active(indices[i]) && !is_inactive(indices[i]))
              //func(indices[i], 1.0 / (double) num_active);
              func(indices[i], 1);
          }
        }
      }

      // Given an index, get the clique
      template<class S, class Graph>
      void extract_clique(S index, uintE* base, Graph& G, int k) {
        Y vert;
        uintE v = get_mtable<Y, C>(index, space);

        assert(v < nx);

        /*uintE v2 = get_top_index(index); 
        if (v != v2) {
          std::cout << "v: " << v << ", v2: " << v2 << std::endl; fflush(stdout);
        }
        assert(v == v2);*/
        base[0] = v;
        vert = std::get<0>(space[index]);
        for (int j = 0; j < rr - 1; ++j) {
          unsigned __int128 mask = (1ULL << (shift_factor)) - 1;
          uintE extract = (uintE) (vert & mask); // vert & mask
          /*if (static_cast<uintE>(extract) >= G.n) {
            std::cout << "Vert: " << static_cast<uintE>(extract) << ", n: " << G.n << std::endl;
          }*/
          //assert(static_cast<uintE>(extract) < G.n);
          base[k - j] = static_cast<uintE>(extract);
          vert = vert >> shift_factor;
        }
      }

    template<class S>
    std::tuple<uintE, uintE> extract_clique_two(S index, int k) {
      uintE v = get_mtable<Y, C>(index, space);
      Y vert = std::get<0>(space[index]);
      unsigned __int128 mask = (1ULL << (shift_factor)) - 1;
      uintE extract = (uintE) (vert & mask); // vert & mask
      uintE v2 = static_cast<uintE>(extract);
      return std::make_tuple(v, v2);
    }

    template<class S>
    std::tuple<uintE, uintE, uintE> extract_clique_three(S index, int k) {
      Y vert;
      uintE v = get_mtable<Y, C>(index, space);
      vert = std::get<0>(space[index]);
      uintE v2 = UINT_E_MAX; uintE v3 = UINT_E_MAX;
        for (int j = 0; j < rr - 1; ++j) {
          unsigned __int128 mask = (1ULL << (shift_factor)) - 1;
          uintE extract = (uintE) (vert & mask); // vert & mask
          if (j == 0) v2 = static_cast<uintE>(extract);
          else v3 = static_cast<uintE>(extract);
          vert = vert >> shift_factor;
        }
        return std::make_tuple(v, v2, v3);
    }
  };

} // end namespace twotable

} //end namespace gbbs