#ifndef HALO2_PROOFS_INCLUDE_BN254_POSEIDON_WRITER_H_
#define HALO2_PROOFS_INCLUDE_BN254_POSEIDON_WRITER_H_

#include <stdint.h>

#include <memory>

#include <tachyon/c/zk/plonk/halo2/bn254_transcript.h>

#include "rust/cxx.h"

namespace tachyon::halo2_api::bn254 {

struct Fr;

class PoseidonWriter {
 public:
  PoseidonWriter();
  PoseidonWriter(const PoseidonWriter& other) = delete;
  PoseidonWriter& operator=(const PoseidonWriter& other) = delete;
  ~PoseidonWriter();

  void update(rust::Slice<const uint8_t> data);
  rust::Box<Fr> squeeze();
  rust::Vec<uint8_t> state() const;

 private:
  tachyon_halo2_bn254_transcript_writer* writer_;
};

std::unique_ptr<PoseidonWriter> new_poseidon_writer();

}  // namespace tachyon::halo2_api::bn254

#endif  // HALO2_PROOFS_INCLUDE_BN254_POSEIDON_WRITER_H_
