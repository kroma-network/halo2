#include "halo2_proofs/include/bn254_blake2b_writer.h"

#include <string.h>

namespace tachyon::halo2_api::bn254 {

Blake2bWriter::Blake2bWriter()
    : writer_(tachyon_halo2_bn254_transcript_writer_create(
          TACHYON_HALO2_BLAKE2B_TRANSCRIPT)) {}

Blake2bWriter::~Blake2bWriter() {
  tachyon_halo2_bn254_transcript_writer_destroy(writer_);
}

void Blake2bWriter::update(rust::Slice<const uint8_t> data) {
  tachyon_halo2_bn254_transcript_writer_update(writer_, data.data(),
                                               data.size());
}

void Blake2bWriter::finalize(
    std::array<uint8_t, kBlake2bDigestLength>& result) {
  uint8_t data[kBlake2bDigestLength];
  size_t data_size;
  tachyon_halo2_bn254_transcript_writer_finalize(writer_, data, &data_size);
  memcpy(result.data(), data, data_size);
}

rust::Vec<uint8_t> Blake2bWriter::state() const {
  rust::Vec<uint8_t> ret;
  // NOTE(chokobole): |rust::Vec<uint8_t>| doesn't have |resize()|.
  ret.reserve(kBlake2bStateLength);
  for (size_t i = 0; i < kBlake2bStateLength; ++i) {
    ret.push_back(0);
  }
  size_t state_size;
  tachyon_halo2_bn254_transcript_writer_get_state(writer_, ret.data(),
                                                  &state_size);
  return ret;
}

std::unique_ptr<Blake2bWriter> new_blake2b_writer() {
  return std::make_unique<Blake2bWriter>();
}

}  // namespace tachyon::halo2_api::bn254
