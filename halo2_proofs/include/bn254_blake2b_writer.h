#ifndef HALO2_PROOFS_INCLUDE_BN254_BLAKE2B_WRITER_H_
#define HALO2_PROOFS_INCLUDE_BN254_BLAKE2B_WRITER_H_

#include <stddef.h>
#include <stdint.h>

#include <array>
#include <memory>

#include <tachyon/c/zk/plonk/halo2/bn254_transcript.h>

#include "rust/cxx.h"

namespace tachyon::halo2_api::bn254 {

constexpr size_t kBlake2bDigestLength = 64;
constexpr size_t kBlake2bStateLength = 216;

class Blake2bWriter {
 public:
  Blake2bWriter();
  Blake2bWriter(const Blake2bWriter& other) = delete;
  Blake2bWriter& operator=(const Blake2bWriter& other) = delete;
  ~Blake2bWriter();

  void update(rust::Slice<const uint8_t> data);
  void finalize(std::array<uint8_t, kBlake2bDigestLength>& result);
  rust::Vec<uint8_t> state() const;

 private:
  tachyon_halo2_bn254_transcript_writer* writer_;
};

std::unique_ptr<Blake2bWriter> new_blake2b_writer();

}  // namespace tachyon::halo2_api::bn254

#endif  // HALO2_PROOFS_INCLUDE_BN254_BLAKE2B_WRITER_H_
