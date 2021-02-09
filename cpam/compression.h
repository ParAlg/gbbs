#pragma once

#include <functional>
#include <optional>
#include <tuple>
#include <type_traits>

namespace cpam {

struct diffencoded_entry_encoder {

  struct data {};

  template <class Entry, bool is_aug = false>
  struct encoder {
    using ET = typename Entry::entry_t;
    using K = typename Entry::key_t;
    using V = typename Entry::val_t;  // possibly empty (should ensure that default_val in set is empty)
    static constexpr bool is_trivial = false;  // to test

    static inline void print_info(const ET& et) {
    }

    static inline size_t encoded_size(ET* data, size_t size) {
      uint8_t stk[2*sizeof(K)];
      assert(size > 0);
      K prev_key = Entry::get_key(data[0]);
      size_t key_bytes = sizeof(K);  // first key is uncompressed
      for (size_t i=0; i<size; i++) {
        K cur_key = Entry::get_key(data[i]);
        K next_diff = cur_key - prev_key;
        auto bytes_used = encodeUnsigned<K>((uint8_t*)stk, 0, next_diff);  // compressed difference
        key_bytes += bytes_used;
        prev_key = cur_key;
      }
      size_t val_bytes = size * sizeof(V);
      return key_bytes + val_bytes;
    }

//    static inline void encode(ET* data, size_t size, uint8_t* bytes) {
//      V* vals = (V*)bytes;
//      uint8_t* key_bytes = (bytes + size*sizeof(V));
//
//      K prev_key = Entry::get_key(data[0]);
//      *((K*)key_bytes) = prev_key;  // store first key
//      vals[0] = Entry::get_val(data[0]);
//      size_t offset = sizeof(K);  // first key is uncompressed
//      for (size_t i=1; i<size; i++) {
//        vals[i] = Entry::get_val(data[i]);
//        K cur_key = Entry::get_key(data[i]);
//        K next_diff = cur_key - prev_key;
//        offset = encodeUnsigned<K>(key_bytes, offset, next_diff);  // compressed difference
//        prev_key = cur_key;
//      }
//    }

    static inline auto encode(ET* data, size_t size, uint8_t* bytes) {
      V* vals = (V*)bytes;
      uint8_t* key_bytes = (bytes + size*sizeof(V));

      if constexpr (is_aug) {
        using AT = typename Entry::aug_t;
        AT av = Entry::from_entry(data[0]);
        K prev_key = Entry::get_key(data[0]);
        *((K*)key_bytes) = prev_key;  // store first key
        vals[0] = Entry::get_val(data[0]);
        size_t offset = sizeof(K);  // first key is uncompressed
        for (size_t i=1; i<size; i++) {
          vals[i] = Entry::get_val(data[i]);
          K cur_key = Entry::get_key(data[i]);
          K next_diff = cur_key - prev_key;
          offset = encodeUnsigned<K>(key_bytes, offset, next_diff);  // compressed difference
          prev_key = cur_key;
          av = Entry::combine(av, Entry::from_entry(data[i]));
        }
        return av;
      } else {
        K prev_key = Entry::get_key(data[0]);
        *((K*)key_bytes) = prev_key;  // store first key
        vals[0] = Entry::get_val(data[0]);
        size_t offset = sizeof(K);  // first key is uncompressed
        for (size_t i=1; i<size; i++) {
          vals[i] = Entry::get_val(data[i]);
          K cur_key = Entry::get_key(data[i]);
          K next_diff = cur_key - prev_key;
          offset = encodeUnsigned<K>(key_bytes, offset, next_diff);  // compressed difference
          prev_key = cur_key;
        }
      }
    }

    template <class F>
    static inline void decode(uint8_t* bytes, size_t size, const F& f) {
      V* vals = (V*)bytes;
      uint8_t* key_bytes = (bytes + size*sizeof(V));

      K prev_key = *((K*)key_bytes);
      f(Entry::to_entry(prev_key, vals[0]));
      key_bytes += sizeof(K);
      for (size_t i=1; i<size; i++) {
        prev_key += decodeUnsigned<K>(key_bytes);
        f(Entry::to_entry(prev_key, vals[i]));
      }
    }

    template <class F>
    static inline bool decode_cond(uint8_t* bytes, size_t size, const F& f) {
      V* vals = (V*)bytes;
      uint8_t* key_bytes = (bytes + size*sizeof(V));

      K prev_key = *((K*)key_bytes);
      f(Entry::to_entry(prev_key, vals[0]));
      key_bytes += sizeof(K);
      for (size_t i=1; i<size; i++) {
        prev_key += decodeUnsigned<K>(key_bytes);
        if (!f(Entry::to_entry(prev_key, vals[i]))) return false;
      }
      return true;
    }

    static inline void destroy(uint8_t* bytes, size_t size) {
      V* vals = (V*)bytes;
      for (size_t i=0; i<size; i++) {
        vals[i].~V();
      }
    }
  };
};


// Only difference between this and the default entry encoder is that encode
// returns an AT.
//struct aug_entry_encoder {
//  template <class entry>
//  using ET = typename entry::entry_t;
//  using AT = typename entry::aug_t;
//  using Entry = entry;
//
//  static constexpr bool is_trivial = false;  // to test
//
//  static inline void print_info(const ET& et) {
////    std::cout << et.first;
//  }
//
//  static inline size_t encoded_size(ET* data, size_t size) {
//    return sizeof(ET)*size;
//  }
//
//  static inline AT encode(ET* data, size_t size, uint8_t* bytes) {
//    ET* ets = (ET*)bytes;
//    AT av = Entry::from_entry(ets[0]);
//    ets[0] = std::move(data[0]);
//    for (size_t i=1; i<size; i++) {
//      ets[i] = std::move(data[i]);
//      av = Entry::combine(av, Entry::from_entry(ets[i]));
//    }
//    return av;
//  }
//
//  template <class F>
//  static inline void decode(uint8_t* bytes, size_t size, const F& f) {
//    ET* ets = (ET*)bytes;
//    for (size_t i=0; i<size; i++) {
//      f(ets[i]);
//    }
//  }
//
//  static inline void destroy(uint8_t* bytes, size_t size) {
//    ET* ets = (ET*)bytes;
//    for (size_t i=0; i<size; i++) {
//      ets[i].~ET();
//    }
//  }
//};

struct null_encoder {
  template <class ET>
  static inline void print_info(const ET& et) {
    assert(false);
    exit(-1);}
  template <class ET>
  static inline size_t encoded_size(ET* data, size_t size) {
    assert(false);
    exit(-1); return sizeof(ET)*size;}
  template <class ET>
  static inline void encode(ET* data, size_t size, uint8_t* bytes) {
    assert(false);
    exit(-1);}
  template <class F>
  static inline void decode(uint8_t* bytes, size_t size, const F& f) {
    assert(false);
    exit(-1);}
  template <class F>
  static inline bool decode_cond(uint8_t* bytes, size_t size, const F& f) {
    assert(false);
    exit(-1);}
  static inline void destroy(uint8_t* bytes, size_t size) {
    assert(false);
    exit(-1);}
};


struct default_entry_encoder {

  struct data {};

  template <class Entry, bool is_aug=false>
  struct encoder {
    using ET = typename Entry::entry_t;
    static constexpr bool is_trivial = false;  // to test

    static inline void print_info(const ET& et) {
  //    std::cout << et.first;
    }

    // TODO remove ET* data
    static inline size_t encoded_size(ET* data, size_t size) {
      return sizeof(ET)*size;
    }

    static inline auto encode(ET* data, size_t size, uint8_t* bytes) {
      if constexpr (is_aug) {
        using AT = typename Entry::aug_t;
        ET* ets = (ET*)bytes;
        parlay::move_uninitialized(ets[0], data[0]);
        AT av = Entry::from_entry(ets[0]);
        //size_t total_size = 0;
        //total_size += std::get<1>(ets[0]).size();
        for (size_t i=1; i<size; i++) {
          parlay::move_uninitialized(ets[i], data[i]);
          av = Entry::combine(av, Entry::from_entry(ets[i]));
          //total_size += std::get<1>(ets[i]).size();
        }
        //std::cout << "Total size in encoded block = " << total_size << std::endl;
        return av;
      } else {
        ET* ets = (ET*)bytes;
        for (size_t i=0; i<size; i++) {
          parlay::move_uninitialized(ets[i], data[i]);
        }
      }
    }

    template <class F>
    static inline void decode(uint8_t* bytes, size_t size, const F& f) {
      ET* ets = (ET*)bytes;
      for (size_t i=0; i<size; i++) {
        f(ets[i]);
      }
    }

    template <class F>
    static inline bool decode_cond(uint8_t* bytes, size_t size, const F& f) {
      ET* ets = (ET*)bytes;
      for (size_t i=0; i<size; i++) {
        if (!f(ets[i])) return false;
      }
      return true;
    }

    static inline void destroy(uint8_t* bytes, size_t size) {
      ET* ets = (ET*)bytes;
      for (size_t i=0; i<size; i++) {
        ets[i].~ET();
      }
    }
  };
};

}  // namespace cpam
