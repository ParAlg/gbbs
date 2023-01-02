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

namespace gbbs {

unsigned nChoosek( unsigned n, unsigned k )
{
    if (k > n) return 0;
    if (k * 2 > n) k = n-k;
    if (k == 0) return 1;

    int result = n;
    for( int i = 2; i <= k; ++i ) {
        result *= (n-i+1);
        result /= i;
    }
    return result;
}

class list_buffer {
  public:
    int buffer;
    sequence<uintE> list;
    sequence<size_t> starts;
    sequence<bool> to_pack;
    size_t next;
    size_t num_workers2;
    size_t ss;
    size_t efficient = 1;

    // for dynamic
    using ListType = sequence<uintE>;
    sequence<ListType> dyn_lists;
    sequence<size_t> dyn_list_starts;
    sequence<size_t> nexts;
    sequence<char> dyn_list_init;
    gbbs::dyn_arr<bool> dyn_to_pack;
    size_t init_size;
    size_t next_dyn_list;

    // for easier dynamic?
    std::vector<std::vector<uintE>> ddyn_lists;

    void print_size() {
      size_t space = 0;
      if (efficient == 1) {
        space += (ss + 1024 * num_workers2) * sizeof(uintE); // for list
        space += (num_workers2) * sizeof(size_t); // for starts
        space += sizeof(bool) * (ss + 1024 * num_workers2); // for topack
      } else if (efficient == 0) {
        space += sizeof(uintE) * list.size();
      } else if (efficient == 5) {
        space += sizeof(size_t) * (num_workers2 + 1);
        for (size_t i = 0; i < ddyn_lists.size(); i++) {
          space += sizeof(uintE) * ddyn_lists[i].size();
        }
      }
      std::cout << "List buffer space: " << space << std::endl;
    }

// option to have a thread local vector that's resizable per thread to hold
// the new r cliques; you need to get your thread id
    list_buffer(size_t s, size_t _efficient = 1){
      efficient = _efficient;
      if (efficient == 1) {
        ss = s;
        num_workers2 = num_workers();
        buffer = 1024;
        int buffer2 = 1024;
        list = sequence<uintE>(s + buffer2 * num_workers2);
        starts = sequence<size_t>::from_function(num_workers2, [&](size_t i){return i * buffer2;});
        std::cout << "list size: " << sizeof(size_t) * num_workers2 + sizeof(uintE) * list.size() << std::endl;
        next = num_workers2 * buffer2;
        to_pack = sequence<bool>(s + buffer2 * num_workers2, true);
      } else if (efficient == 0) {
        list = sequence<uintE>(s, static_cast<uintE>(UINT_E_MAX));
        std::cout << "list size: " << sizeof(uintE) * list.size() << std::endl;
        next = 0;
      } else if (efficient == 4) {
        // version of list buffer that's dynamically sized
        ss = s;
        num_workers2 = num_workers();
        buffer = 1024;
        int buffer2 = 1024;
        int multiplier = 10000;
        init_size = (1 + (s / multiplier) / buffer2) * buffer2  + buffer2 * num_workers2;
        dyn_lists = sequence<ListType>::from_function(multiplier, [](size_t i){return ListType();});
        dyn_list_init = sequence<char>::from_function(multiplier, [](size_t i){return false;});
        dyn_lists[0] = ListType(init_size);
        dyn_list_init[0] = true;
        next_dyn_list = 0;
        //dyn_list = dyn_arr<uintE>(init_size);
        dyn_list_starts = sequence<size_t>::from_function(num_workers2, [&](size_t i){return 0;});
        starts = sequence<size_t>::from_function(num_workers2, [&](size_t i){return i * buffer2;});
        std::cout << "list size: " << sizeof(size_t) * num_workers2 + sizeof(uintE) * init_size << std::endl;
        nexts = sequence<size_t>::from_function(multiplier, [](size_t i){return 0;});
        nexts[0] = num_workers2 * buffer2;
        dyn_to_pack = gbbs::dyn_arr<bool>(init_size);
        dyn_to_pack.copyInF([](size_t i){return true;}, init_size);
      } else if (efficient == 5) {
        // easier dynamically sized list buffer
        ss = s;
        num_workers2 = num_workers();
        buffer = 1024;
        int buffer2 = 1024;
        ddyn_lists = std::vector<std::vector<uintE>>(num_workers2, std::vector<uintE>(100, 0));
        starts = sequence<size_t>::from_function(num_workers2 + 1, [&](size_t i){return 0;});
      }
    }

    void resize(size_t num_active, size_t k, size_t r, size_t cur_bkt) {
    }

    void add(size_t index) {
      if (efficient == 1) {
        size_t worker = worker_id();
        list[starts[worker]] = index;
        starts[worker]++;
        if (starts[worker] % buffer == 0) {
          size_t use_next = gbbs::fetch_and_add(&next, buffer);
          starts[worker] = use_next;
        }
      } else if (efficient == 0) {
        size_t use_next = gbbs::fetch_and_add(&next, 1);
        list[use_next] = index;
      } else if (efficient == 5) {
        size_t worker = worker_id();
        if (ddyn_lists[worker].size() < starts[worker] + 1)
          ddyn_lists[worker].resize(2 * (starts[worker] + 1));
        ddyn_lists[worker][starts[worker]] = index;
        starts[worker]++;
      } else if (efficient == 4) {
        size_t worker = worker_id();
        dyn_lists[dyn_list_starts[worker]][starts[worker]] = index;
        starts[worker]++;
        if (starts[worker] % buffer == 0) {
          size_t use_next = gbbs::fetch_and_add(&(nexts[dyn_list_starts[worker]]), buffer);
          while (use_next >= init_size) {
            //while (nexts[dyn_list_starts[worker]] >= init_size) {
            dyn_list_starts[worker]++;
            //}
            use_next = gbbs::fetch_and_add(&(nexts[dyn_list_starts[worker]]), buffer);
          }
          starts[worker] = use_next;
          // But now we need to make sure that dyn_lists[dyn_list_starts[worker]] actually has space....
          // TODO check this is ok esp for contention maybe take a lock instead
          while(dyn_lists[dyn_list_starts[worker]].size() == 0) {
            if (gbbs::atomic_compare_and_swap(&dyn_list_init[dyn_list_starts[worker]], static_cast<char>(false), static_cast<char>(true))) {
              dyn_lists[dyn_list_starts[worker]] = ListType(init_size); 
              break;
            }
          }
          // now need to make sure that dyn_list has space for the buffer we've just claimed
          // this needs to take a lock around dyn_list to work? or use another array to maintain
          // dyn_lists....
          //dyn_list.resize(use_next + buffer);
        }
      }
    }

    size_t num_entries() {
      if (efficient == 1) return next;
      if (efficient == 0) return next;
      if (efficient == 5) {
        starts[num_workers2] = 0;
        size_t total = parlay::scan_inplace(make_slice(starts));
        return total;
      }
      if (efficient == 4){
        auto sizes = sequence<size_t>::from_function(num_workers2, [&](size_t i){return dyn_list_starts[i] * init_size + starts[i];});
        auto max_size = parlay::reduce_max(sizes);
        return max_size;
      }
      return 0;
    }

    void void_v(size_t i, uintE actual_v) {
      if (efficient == 5) {
        // ***TODO this should be a binary search
        for (size_t worker = 0; worker < num_workers2; worker++) {
          auto beginning = starts[worker];
          auto ending = starts[worker + 1];
          if (i >= beginning && i < ending) {
            size_t idx = i - beginning;
            ddyn_lists[worker][idx] = UINT_E_MAX;
            return;
          }
        }
      } else if (efficient == 0) {
        list[i] = UINT_E_MAX;
      }
      else {
        std::cout << "unsupported" << std::endl; fflush(stdout);
        exit(0);
      }
    }
    uintE get_v(size_t i) {
      if (efficient == 5) {
        // ***TODO this should be a binary search
        for (size_t worker = 0; worker < num_workers2; worker++) {
          auto beginning = starts[worker];
          auto ending = starts[worker + 1];
          if (i >= beginning && i < ending) {
            size_t idx = i - beginning;
            return ddyn_lists[worker][idx];
          }
        }
      } else if (efficient == 0) {
        return list[i];
      }
      else {
        std::cout << "unsupported" << std::endl; fflush(stdout);
        exit(0);
      }
      std::cout << "error" << std::endl; fflush(stdout);
      return 0;
    }

    template <class I>
    size_t filter(I& update_changed, sequence<double>& per_processor_counts) {
      if (efficient == 1) {
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
      } else if (efficient == 0) {
        parallel_for(0, next, [&](size_t worker) {
          //assert(list[worker] != UINT_E_MAX);
          assert(per_processor_counts[list[worker]] != 0);
          update_changed(per_processor_counts, worker, list[worker]);
        });
        return next;
      } else if (efficient == 5) {
        starts[num_workers2] = 0;
        size_t total = parlay::scan_inplace(make_slice(starts));
        parallel_for(0, num_workers2, [&](size_t worker){
          parallel_for(starts[worker], starts[worker + 1], [&](size_t j){
            size_t i = j - starts[worker];
            update_changed(per_processor_counts, j, ddyn_lists[worker][i]);
          });
        });
        return total;
      }else if (efficient == 4) {
        //std::cout << "FILTER" << std::endl;
        //fflush(stdout);
        // First ensure that dyn_to_pack is the right size
        // To do this, we need the max of dyn_list_starts[worker] * init_size + starts[worker] for 0 to num_workers2
        auto sizes = sequence<size_t>::from_function(num_workers2, [&](size_t i){return dyn_list_starts[i] * init_size + starts[i];});
        auto max_size = parlay::reduce_max(sizes);
        auto absolute_max = 1 + 10000 * ((1 + (ss / 10000) / 1024) * 1024  + 1024* num_workers());
        assert(max_size <= absolute_max);
        if (max_size > dyn_to_pack.size) dyn_to_pack.copyInF([](size_t i){return true;}, max_size - dyn_to_pack.size);
        assert(dyn_to_pack.size >= max_size);
        parallel_for(0, num_workers2, [&](size_t worker) {
          size_t divide = starts[worker] / buffer;
          size_t offset = dyn_list_starts[worker] * init_size;
          for (size_t j = offset + starts[worker]; j < offset + (divide + 1) * buffer; j++) {
            if (j < dyn_to_pack.size) dyn_to_pack.A[j] = false;
            //assert(dyn_lists[dyn_list_starts[worker]][j - offset] == UINT_E_MAX);
          }
        });
        // Pack out 0 to next of list into pack
        parallel_for(0, max_size, [&] (size_t i) {
          if (dyn_to_pack.A[i]) {
            auto val = dyn_lists[i / init_size][i % init_size];
            //assert(val != UINT_E_MAX);
            update_changed(per_processor_counts, i, val);
          } else {
            //auto val = dyn_lists[i / init_size][i % init_size];
            //assert(val == UINT_E_MAX);
            update_changed(per_processor_counts, i, UINT_E_MAX);
          }
        });
        parallel_for(0, num_workers2, [&](size_t worker) {
          size_t divide = starts[worker] / buffer;
          size_t offset = dyn_list_starts[worker] * init_size;
          for (size_t j = offset + starts[worker]; j < offset + (divide + 1) * buffer; j++) {
            if (j < dyn_to_pack.size) dyn_to_pack.A[j] = true;
          }
        });
        //std::cout << "END FILTER" << std::endl;
        //fflush(stdout);
        return max_size;
      }
      return 0;
    }

    void reset() {
      if (efficient == 1) {
      parallel_for (0, num_workers2, [&] (size_t j) {
        starts[j] = j * buffer;
      });
      /*parallel_for (0, ss + buffer * num_workers2, [&] (size_t j) {
        list[j] = UINT_E_MAX;
      });*/
      next = num_workers2 * buffer;
      } else if (efficient == 0) {
        next = 0;
      } else if(efficient == 5) {
        parallel_for (0, num_workers2, [&] (size_t j) {
          starts[j] = 0;
        });
      }else if (efficient == 4) {
        parallel_for (0, num_workers2, [&] (size_t j) {
          starts[j] = j * buffer;
          dyn_list_starts[j] = 0;
        });
        parallel_for(0, nexts.size(), [](size_t j){return 0;});
        nexts[0] = num_workers2 * buffer;
      }
    }
};


template<class CountType>
class list_count_buffer {
  public:
    using PairType = std::tuple<uintE, CountType>;
    int buffer;
    sequence<PairType> list;
    sequence<size_t> starts;
    sequence<bool> to_pack;
    size_t next;
    size_t num_workers2;
    size_t ss;
    size_t efficient = 1;

    // for easier dynamic?
    std::vector<std::vector<PairType>> ddyn_lists;

// option to have a thread local vector that's resizable per thread to hold
// the new r cliques; you need to get your thread id
    list_count_buffer(size_t s, size_t _efficient = 1){
      efficient = _efficient;
      if (efficient == 1) {
        ss = s;
        num_workers2 = num_workers();
        buffer = 1024;
        int buffer2 = 1024;
        list = sequence<PairType>(s + buffer2 * num_workers2);
        starts = sequence<size_t>::from_function(num_workers2, [&](size_t i){return i * buffer2;});
        std::cout << "list size: " << sizeof(size_t) * num_workers2 + sizeof(PairType) * list.size() << std::endl;
        next = num_workers2 * buffer2;
        to_pack = sequence<bool>(s + buffer2 * num_workers2, true);
      } else if (efficient == 0) {
        list = sequence<PairType>(s, std::make_tuple(UINT_E_MAX, 0));
        std::cout << "list size: " << sizeof(uintE) * list.size() << std::endl;
        next = 0;
      } else if (efficient == 2 || efficient == 4) {
        std::cout << "UNIMPLEMENTED FOR LIST COUNT BUFFER" << std::endl; fflush(stdout);
        exit(0);
      } else if (efficient == 5) {
        // easier dynamically sized list buffer
        ss = s;
        num_workers2 = num_workers();
        buffer = 1024;
        int buffer2 = 1024;
        ddyn_lists = std::vector<std::vector<PairType>>(num_workers2, std::vector<PairType>(100));
        starts = sequence<size_t>::from_function(num_workers2 + 1, [&](size_t i){return 0;});
      }
    }

    void resize(size_t num_active, size_t k, size_t r, size_t cur_bkt) {
    }

    void add(size_t index, PairType val) {
      if (efficient == 1) {
        size_t worker = worker_id();
        list[starts[worker]] = std::make_tuple(index, val);
        starts[worker]++;
        if (starts[worker] % buffer == 0) {
          size_t use_next = gbbs::fetch_and_add(&next, buffer);
          starts[worker] = use_next;
        }
      } else if (efficient == 0) {
        size_t use_next = gbbs::fetch_and_add(&next, 1);
        list[use_next] = std::make_tuple(index, val);
      } else if (efficient == 5) {
        size_t worker = worker_id();
        if (ddyn_lists[worker].size() < starts[worker] + 1)
          ddyn_lists[worker].resize(2 * (starts[worker] + 1));
        ddyn_lists[worker][starts[worker]] = std::make_tuple(index, val);
        starts[worker]++;
      }
    }

    size_t num_entries() {
      if (efficient == 1) return next;
      if (efficient == 0) return next;
      if (efficient == 5) {
        starts[num_workers2] = 0;
        size_t total = parlay::scan_inplace(make_slice(starts));
        return total;
      }
      return 0;
    }

    void void_v(size_t i, uintE actual_v) {
      if (efficient == 5) {
        // ***TODO this should be a binary search
        for (size_t worker = 0; worker < num_workers2; worker++) {
          auto beginning = starts[worker];
          auto ending = starts[worker + 1];
          if (i >= beginning && i < ending) {
            size_t idx = i - beginning;
            ddyn_lists[worker][idx] = std::make_tuple(UINT_E_MAX, 0);
            return;
          }
        }
      } else if (efficient == 0) {
        list[i] = std::make_tuple(UINT_E_MAX, 0);
      }
      else {
        std::cout << "unsupported" << std::endl; fflush(stdout);
        exit(0);
      }
    }
    PairType get_v(size_t i) {
      if (efficient == 5) {
        // ***TODO this should be a binary search
        for (size_t worker = 0; worker < num_workers2; worker++) {
          auto beginning = starts[worker];
          auto ending = starts[worker + 1];
          if (i >= beginning && i < ending) {
            size_t idx = i - beginning;
            return ddyn_lists[worker][idx];
          }
        }
      } else if (efficient == 0) {
        return list[i];
      }
      std::cout << "error" << std::endl; fflush(stdout);
      return std::make_tuple(0, 0);
    }

    template <class I>
    size_t filter(I& update_changed) {
      if (efficient == 1) {
      parallel_for(0, num_workers2, [&](size_t worker) {
        size_t divide = starts[worker] / buffer;
        for (size_t j = starts[worker]; j < (divide + 1) * buffer; j++) {
          to_pack[j] = false;
        }
      });
      // Pack out 0 to next of list into pack
      parallel_for(0, next, [&] (size_t i) {
        if (to_pack[i])
          update_changed(i, list[i]);
        else
          update_changed(i, std::make_tuple(UINT_E_MAX,0));
      });
      parallel_for(0, num_workers2, [&](size_t worker) {
        size_t divide = starts[worker] / buffer;
        for (size_t j = starts[worker]; j < (divide + 1) * buffer; j++) {
          to_pack[j] = true;
        }
      });
      return next;
      } else if (efficient == 0) {
        parallel_for(0, next, [&](size_t worker) {
          update_changed(worker, list[worker]);
        });
        return next;
      } else if (efficient == 5) {
        starts[num_workers2] = 0;
        size_t total = parlay::scan_inplace(make_slice(starts));
        parallel_for(0, num_workers2, [&](size_t worker){
          parallel_for(starts[worker], starts[worker + 1], [&](size_t j){
            size_t i = j - starts[worker];
            update_changed(j, ddyn_lists[worker][i]);
          });
        });
        return total;
      }
      return 0;
    }

    void reset() {
      if (efficient == 1) {
      parallel_for (0, num_workers2, [&] (size_t j) {
        starts[j] = j * buffer;
      });
      /*parallel_for (0, ss + buffer * num_workers2, [&] (size_t j) {
        list[j] = UINT_E_MAX;
      });*/
      next = num_workers2 * buffer;
      } else if (efficient == 0) {
        next = 0;
      } else if(efficient == 5) {
        parallel_for (0, num_workers2, [&] (size_t j) {
          starts[j] = 0;
        });
      }
    }
};

} // namespace gbbs