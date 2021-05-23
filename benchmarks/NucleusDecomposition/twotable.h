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

namespace twotable {

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
      bool contiguous_space;
      int shift_factor;
  
      template <class Graph>
      TwolevelHash(int r, Graph& DG, size_t max_deg, bool _contiguous_space, bool relabel, int _shift_factor) {
        using W = typename Graph::weight_type;
        shift_factor = _shift_factor;
        contiguous_space = _contiguous_space;
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

        if (r == 2) {
          parallel_for(0, DG.n, [&](std::size_t i){
            if (DG.get_vertex(i).getOutDegree() != 0) {
              auto map_f = [&](const uintE& src, const uintE& ngh, const W& wgh) {
                auto base = sequence<uintE>(r);
                base[0] = i;
                base[1] = ngh;
                base_f(base);
              };
              DG.get_vertex(i).mapOutNgh(i, map_f, false);
            }
          });
        } else {
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
        }

        auto top_table_sizes2 = tmp_table.entries();
        // sort by key
        pbbslib::sample_sort_inplace (top_table_sizes2.slice(), [&](const std::tuple<uintE, long>& u, const std::tuple<uintE, long>&  v) {
          return std::get<0>(u) < std::get<0>(v);
        });

        sequence<long> actual_sizes(top_table_sizes2.size() + 1);
        // Modify top_table_sizes2 to be appropriately oversized
        parallel_for(0, top_table_sizes2.size(), [&](std::size_t i) {
          auto m = (size_t)1 << pbbslib::log2_up((size_t)(1.1 * std::get<1>(top_table_sizes2[i])) + 1);
          actual_sizes[i] = m;
        });
        actual_sizes[top_table_sizes2.size()] = 0;
        // Do a scan inplace
        //auto add_tuple_monoid = pbbs::make_monoid([](std::tuple<uintE, long> a, std::tuple<uintE, long> b){
        //  return std::make_tuple(std::get<0>(b), std::get<1>(a) + std::get<1>(b));
        //}, std::make_tuple(UINT_E_MAX, 0));
        long total_top_table_sizes2 = scan_inplace(actual_sizes.slice(), pbbs::addm<long>());
        // Allocate space for the second level tables
        if (contiguous_space) space = pbbslib::new_array_no_init<X>(total_top_table_sizes2);
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
          end_table->vtx = vtx;
          end_table->table = contiguous_space ? pbbslib::sparse_table<Y, long, H>(
            size, 
            std::make_tuple<Y, long>(std::numeric_limits<Y>::max(), static_cast<long>(0)),
            H{},
            space + actual_sizes[i]
            ) :
            pbbslib::sparse_table<Y, long, H>(
            size, 
            std::make_tuple<Y, long>(std::numeric_limits<Y>::max(), static_cast<long>(0)),
            H{}, 1, true);
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
        size_t data_structure_size = sizeof(*this) + sizeof(MidTableY) + total * sizeof(X) +  (DG.n) * sizeof(EndTableY*) + (DG.n + 1) * sizeof(long);
        std::cout << "Data Structure Size: " << data_structure_size << std::endl;
      }

      void insert_twothree(uintE v1, uintE v2, uintE v3, int r, int k) {
        auto add_f = [&] (long* ct, const std::tuple<Y, long>& tup) {
          pbbs::fetch_and_add(ct, (long)1);
        };
        EndTableY* end_table12 = top_table.arr[std::min(v1, v2)];
        (end_table12->table).insert_f(std::make_tuple(Y{std::max(v1, v2)}, (long) 1), add_f);

        EndTableY* end_table13 = top_table.arr[std::min(v1, v3)];
        (end_table13->table).insert_f(std::make_tuple(Y{std::max(v1, v3)}, (long) 1), add_f);

        EndTableY* end_table23 = top_table.arr[std::min(v2, v3)];
        (end_table23->table).insert_f(std::make_tuple(Y{std::max(v2, v3)}, (long) 1), add_f);
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
        if (contiguous_space) {
          if (std::get<0>(space[index]) == std::numeric_limits<Y>::max()) return 0;
          return std::get<1>(space[index]);
        }
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
        EndTableY* end_table = top_table.arr[top_index];
        if (end_table == nullptr) return 0;
        size_t bottom_index = index - top_table_sizes[top_index];
        /*if (bottom_index >= (end_table->table).m) {
          std::cout << "Bottom: " << bottom_index << std::endl;
          std::cout << "M: " << (end_table->table).m << std::endl;
          fflush(stdout);
        }
        assert(bottom_index < (end_table->table).m);*/
        if (std::get<0>((end_table->table).table[bottom_index]) == std::numeric_limits<Y>::max()) return UINT_E_MAX;
        return std::get<1>((end_table->table).table[bottom_index]);
      }

      size_t update_count(std::size_t index, size_t update){
        if (contiguous_space) {
          auto val = std::get<1>(space[index]) - update;
          space[index] =
            std::make_tuple(std::get<0>(space[index]), val);
          return val;
        }
        size_t top_index = get_top_index(index);
        //EndTable* end_table = std::get<1>(top_table.table.table[top_index]);
        //***for arr
        EndTableY* end_table = top_table.arr[top_index];
        size_t bottom_index = index - top_table_sizes[top_index];
        auto val = std::get<1>((end_table->table).table[bottom_index]) - update;
        (end_table->table).table[bottom_index] = std::make_tuple(
          std::get<0>((end_table->table).table[bottom_index]), val
        );
        return val;
      }

      void clear_count(std::size_t index) {
        if (contiguous_space) {
          space[index] = std::make_tuple(std::get<0>(space[index]), 0);
          return;
        }
        size_t top_index = get_top_index(index);
        //***for arr
        //EndTable* end_table = std::get<1>(top_table.table.table[top_index]);
        EndTableY* end_table = top_table.arr[top_index];
        size_t bottom_index = index - top_table_sizes[top_index];
        (end_table->table).table[bottom_index] = std::make_tuple(
          std::get<0>((end_table->table).table[bottom_index]), 0
        );
      }

      void set_count(std::size_t index, size_t update) {
        if (contiguous_space) {
          space[index] = std::make_tuple(std::get<0>(space[index]), update);
          return;
        }
        size_t top_index = get_top_index(index);
        //***for arr
        //EndTable* end_table = std::get<1>(top_table.table.table[top_index]);
        EndTableY* end_table = top_table.arr[top_index];
        size_t bottom_index = index - top_table_sizes[top_index];
        (end_table->table).table[bottom_index] = std::make_tuple(
          std::get<0>((end_table->table).table[bottom_index]), update
        );
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
      void extract_indices(uintE* base2, HH is_active, HG is_inactive, I func, int r, int k) {
        // Sort base
        uintE base[10];
        assert(10 > k);
        for(std::size_t i = 0; i < k + 1; i++) {
          base[i] = base2[i];
        }
        std::sort(base, base + k + 1,std::less<uintE>());

        std::vector<size_t> indices;
        size_t num_active = 0;
        bool use_func = true;
        unsigned __int128 mask = (1ULL << (shift_factor)) - 1;

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
          if (is_active(prefix + index)) num_active++;
          if (is_inactive(prefix + index)) return;
          //func(prefix + index);
        } while (std::prev_permutation(bitmask.begin(), bitmask.end()));

        assert(num_active != 0);
        if (use_func) {
          for (std::size_t i = 0; i < indices.size(); i++) {
            if (!is_active(indices[i]) && !is_inactive(indices[i]))
              func(indices[i], 1.0 / (double) num_active);
          }
        }
      }

      // Given an index, get the clique
      template<class S, class Graph>
      void extract_clique(S index, uintE* base, Graph& G, int k) {
        Y vert;
        // First, do a binary search for index in prefix
        size_t top_index = get_top_index(index);
        /*base[0] = std::get<0>(top_table.table.table[top_index]);
        EndTable* end_table = std::get<1>(top_table.table.table[top_index]);*/
        //***for arr
        base[0] = top_index;
        if (!contiguous_space) {
          EndTableY* end_table = top_table.arr[top_index];
          size_t bottom_index = index - top_table_sizes[top_index];
          vert = std::get<0>((end_table->table).table[bottom_index]);
        } else {
          vert = std::get<0>(space[index]);
        }
        for (int j = 0; j < rr - 1; ++j) {
          unsigned __int128 mask = (1ULL << (shift_factor)) - 1;
          uintE extract = (uintE) (vert & mask); // vert & mask
          /*if (static_cast<uintE>(extract) >= G.n) {
            std::cout << "Vert: " << static_cast<uintE>(extract) << ", n: " << G.n << std::endl;
          }*/
          assert(static_cast<uintE>(extract) < G.n);
          base[k - j] = static_cast<uintE>(extract);
          vert = vert >> shift_factor;
        }
      }

    template<class S>
    std::tuple<uintE, uintE> extract_clique_two(S index, int k) {
      Y vert;
        // First, do a binary search for index in prefix
        size_t top_index = get_top_index(index);
        if (!contiguous_space) {
          EndTableY* end_table = top_table.arr[top_index];
          size_t bottom_index = index - top_table_sizes[top_index];
          vert = std::get<0>((end_table->table).table[bottom_index]);
        } else {
          vert = std::get<0>(space[index]);
        }
        unsigned __int128 mask = (1ULL << (shift_factor)) - 1;
        uintE extract = (uintE) (vert & mask); // vert & mask
        uintE v2 = static_cast<uintE>(extract);
        return std::make_tuple(static_cast<uintE>(top_index), v2);
    }

    template<class S>
    std::tuple<uintE, uintE, uintE> extract_clique_three(S index, int k) {
      Y vert;
        // First, do a binary search for index in prefix
        size_t top_index = get_top_index(index);
        if (!contiguous_space) {
          EndTableY* end_table = top_table.arr[top_index];
          size_t bottom_index = index - top_table_sizes[top_index];
          vert = std::get<0>((end_table->table).table[bottom_index]);
        } else {
          vert = std::get<0>(space[index]);
        }
        uintE v2 = UINT_E_MAX; uintE v3 = UINT_E_MAX;
        for (int j = 0; j < rr - 1; ++j) {
          unsigned __int128 mask = (1ULL << (shift_factor)) - 1;
          uintE extract = (uintE) (vert & mask); // vert & mask
          if (j == 0) v2 = static_cast<uintE>(extract);
          else v3 = static_cast<uintE>(extract);
          vert = vert >> shift_factor;
        }
        return std::make_tuple(static_cast<uintE>(top_index), v2, v3);
    }
  };

} // end namespace twotable

} //end namespace gbbs