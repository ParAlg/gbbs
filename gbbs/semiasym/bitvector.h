#pragma once

namespace gbbs {

struct bitvector {
  size_t n;
  uint8_t* data;
  bitvector(size_t n) : n(n) {
    size_t n_bytes = (n + 8 - 1) / 8;  // ceil(n/8);a
    data = pbbs::new_array_no_init<uint8_t>(n_bytes);
    uint8_t zero = 0;
    parallel_for(0, n_bytes, [&](size_t i) { data[i] = zero; });
  }

  inline bool ith_bit(uint8_t byte, uint8_t offset) {
    return byte & (1 << offset);
  }

  // returns true on success
  bool atomic_set_bit(size_t k) {
    size_t byte_id = k >> 3;
    uint8_t offset_in_byte = k & 0x7;

    bool bit_set = false;
    uint8_t val = data[byte_id];
    while (!bit_set && !ith_bit(val, offset_in_byte)) {
      val = data[byte_id];
      uint8_t new_val = val | (1 << offset_in_byte);
      bit_set = pbbs::atomic_compare_and_swap(&data[byte_id], val, new_val);
    }
    return bit_set;
  }

  // Used when application can guarantee a single writer per byte.
  void set_bit(size_t k) {
    size_t byte_id = k >> 3;
    uint8_t offset_in_byte = k & 0x7;

    uint8_t val = data[byte_id];
    uint8_t new_val = val | (1 << offset_in_byte);
    data[byte_id] = new_val;
  }

  bool is_set(size_t k) {
    size_t byte_id = k >> 3;
    uint8_t offset_in_byte = k & 0x7;
    return ith_bit(data[byte_id], offset_in_byte);
  }

  void del() {
    if (data) {
      pbbs::free_array(data);
      data = nullptr;
    }
  }
};

}  // namespace gbbs
