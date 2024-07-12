#ifndef HALO2_PROOFS_INCLUDE_CHA_CHA20_RNG_H_
#define HALO2_PROOFS_INCLUDE_CHA_CHA20_RNG_H_

#include <stddef.h>
#include <stdint.h>

#include <array>
#include <memory>

#include <tachyon/c/crypto/random/rng.h>

#include "rust/cxx.h"

namespace tachyon::halo2_api {

class ChaCha20Rng {
 public:
  constexpr static size_t kSeedSize = 32;
  constexpr static size_t kStateSize = sizeof(size_t) + 128;

  explicit ChaCha20Rng(tachyon_rng* rng) : rng_(rng) {}
  explicit ChaCha20Rng(std::array<uint8_t, kSeedSize> seed);
  ChaCha20Rng(const ChaCha20Rng& other) = delete;
  ChaCha20Rng& operator=(const ChaCha20Rng& other) = delete;
  ~ChaCha20Rng();

  uint32_t next_u32();
  std::unique_ptr<ChaCha20Rng> clone() const;
  rust::Vec<uint8_t> state() const;

 private:
  tachyon_rng* rng_;
};

std::unique_ptr<ChaCha20Rng> new_cha_cha20_rng(
    std::array<uint8_t, ChaCha20Rng::kSeedSize> seed);

}  // namespace tachyon::halo2_api

#endif  // HALO2_PROOFS_INCLUDE_CHA_CHA20_RNG_H_
