#include "halo2_proofs/include/bn254_sha256_writer.h"

#include <string.h>

namespace tachyon::halo2_api::bn254 {

Sha256Writer::Sha256Writer()
    : writer_(tachyon_halo2_bn254_transcript_writer_create(
          TACHYON_HALO2_SHA256_TRANSCRIPT)) {}

Sha256Writer::~Sha256Writer() {
  tachyon_halo2_bn254_transcript_writer_destroy(writer_);
}

void Sha256Writer::update(rust::Slice<const uint8_t> data) {
  tachyon_halo2_bn254_transcript_writer_update(writer_, data.data(),
                                               data.size());
}

void Sha256Writer::finalize(std::array<uint8_t, kSha256DigestLength>& result) {
  uint8_t data[kSha256DigestLength];
  size_t data_size;
  tachyon_halo2_bn254_transcript_writer_finalize(writer_, data, &data_size);
  memcpy(result.data(), data, data_size);
}

rust::Vec<uint8_t> Sha256Writer::state() const {
  rust::Vec<uint8_t> ret;
  // NOTE(chokobole): |rust::Vec<uint8_t>| doesn't have |resize()|.
  ret.reserve(kSha256StateLength);
  for (size_t i = 0; i < kSha256StateLength; ++i) {
    ret.push_back(0);
  }
  size_t state_size;
  tachyon_halo2_bn254_transcript_writer_get_state(writer_, ret.data(),
                                                  &state_size);
  return ret;
}

std::unique_ptr<Sha256Writer> new_sha256_writer() {
  return std::make_unique<Sha256Writer>();
}

}  // namespace tachyon::halo2_api::bn254
