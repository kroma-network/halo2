#ifndef HALO2_PROOFS_INCLUDE_XOR_SHIFT_RNG_H_
#define HALO2_PROOFS_INCLUDE_XOR_SHIFT_RNG_H_

#include <stddef.h>
#include <stdint.h>

#include <array>
#include <memory>

#include <tachyon/c/crypto/random/rng.h>

#include "rust/cxx.h"

namespace tachyon::halo2_api {

class XORShiftRng {
 public:
  constexpr static size_t kSeedSize = 16;
  constexpr static size_t kStateSize = 16;

  explicit XORShiftRng(tachyon_rng* rng) : rng_(rng) {}
  explicit XORShiftRng(std::array<uint8_t, kSeedSize> seed);
  XORShiftRng(const XORShiftRng& other) = delete;
  XORShiftRng& operator=(const XORShiftRng& other) = delete;
  ~XORShiftRng();

  uint32_t next_u32();
  std::unique_ptr<XORShiftRng> clone() const;
  rust::Vec<uint8_t> state() const;

 private:
  tachyon_rng* rng_;
};

std::unique_ptr<XORShiftRng> new_xor_shift_rng(
    std::array<uint8_t, XORShiftRng::kSeedSize> seed);

}  // namespace tachyon::halo2_api

#endif  // HALO2_PROOFS_INCLUDE_XOR_SHIFT_RNG_H_
