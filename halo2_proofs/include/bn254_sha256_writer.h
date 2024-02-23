#ifndef HALO2_PROOFS_INCLUDE_BN254_SHA256_WRITER_H_
#define HALO2_PROOFS_INCLUDE_BN254_SHA256_WRITER_H_

#include <stddef.h>
#include <stdint.h>

#include <array>
#include <memory>

#include <tachyon/c/zk/plonk/halo2/bn254_transcript.h>

#include "rust/cxx.h"

namespace tachyon::halo2_api::bn254 {

constexpr size_t kSha256DigestLength = 32;
constexpr size_t kSha256StateLength = 112;

class Sha256Writer {
 public:
  Sha256Writer();
  Sha256Writer(const Sha256Writer& other) = delete;
  Sha256Writer& operator=(const Sha256Writer& other) = delete;
  ~Sha256Writer();

  void update(rust::Slice<const uint8_t> data);
  void finalize(std::array<uint8_t, kSha256DigestLength>& result);
  rust::Vec<uint8_t> state() const;

 private:
  tachyon_halo2_bn254_transcript_writer* writer_;
};

std::unique_ptr<Sha256Writer> new_sha256_writer();

}  // namespace tachyon::halo2_api::bn254

#endif  // HALO2_PROOFS_INCLUDE_BN254_SHA256_WRITER_H_
