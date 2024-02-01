#ifndef HALO2_PROOFS_SRC_RUST_VEC_H_
#define HALO2_PROOFS_SRC_RUST_VEC_H_

#include <stddef.h>
#include <stdint.h>
#include <string.h>

#include <sstream>
#include <type_traits>

#include "rust/cxx.h"

namespace tachyon::halo2_api {

struct RustVec {
  uintptr_t ptr;
  size_t capacity;
  size_t length;

  void Read(const uint8_t* data) {
    memcpy(&capacity, data, sizeof(size_t));
    data += sizeof(size_t);
    memcpy(&ptr, data, sizeof(uintptr_t));
    data += sizeof(uintptr_t);
    memcpy(&length, data, sizeof(size_t));
    data += sizeof(size_t);
  }

  std::string ToString() const {
    std::stringstream ss;
    ss << std::hex << "ptr: " << std::dec << ptr << " capacity: " << capacity
       << " length: " << length;
    return ss.str();
  }
};

template <typename Container, typename T = typename Container::value_type>
rust::Vec<T> ConvertCppContainerToRustVec(const Container& container) {
  rust::Vec<T> ret;
  ret.reserve(std::size(container));
  for (const T& elem : container) {
    ret.push_back(elem);
  }
  return ret;
}

}  // namespace tachyon::halo2_api

#endif  // HALO2_PROOFS_SRC_RUST_VEC_H_
