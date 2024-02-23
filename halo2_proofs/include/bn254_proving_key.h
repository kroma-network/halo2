#ifndef HALO2_PROOFS_INCLUDE_BN254_PROVING_KEY_H_
#define HALO2_PROOFS_INCLUDE_BN254_PROVING_KEY_H_

#include <stddef.h>
#include <stdint.h>

#include <memory>

#include <tachyon/c/zk/plonk/keys/bn254_plonk_proving_key.h>

#include "rust/cxx.h"

namespace tachyon::halo2_api::bn254 {

struct Fr;
class GWCProver;
class SHPlonkProver;

class ProvingKey {
 public:
  explicit ProvingKey(rust::Slice<const uint8_t> pk_bytes);
  ProvingKey(const ProvingKey& other) = delete;
  ProvingKey& operator=(const ProvingKey& other) = delete;
  ~ProvingKey();

  const tachyon_bn254_plonk_proving_key* pk() const { return pk_; }
  tachyon_bn254_plonk_proving_key* pk() { return pk_; }

  rust::Vec<uint8_t> advice_column_phases() const;
  uint32_t blinding_factors() const;
  rust::Vec<uint8_t> challenge_phases() const;
  rust::Vec<size_t> constants() const;
  size_t num_advice_columns() const;
  size_t num_challenges() const;
  size_t num_instance_columns() const;
  rust::Vec<uint8_t> phases() const;
  rust::Box<Fr> transcript_repr_gwc(const GWCProver& prover);
  rust::Box<Fr> transcript_repr_shplonk(const SHPlonkProver& prover);

 private:
  const tachyon_bn254_plonk_verifying_key* GetVerifyingKey() const;
  const tachyon_bn254_plonk_constraint_system* GetConstraintSystem() const;

  tachyon_bn254_plonk_proving_key* pk_;
};

std::unique_ptr<ProvingKey> new_proving_key(
    rust::Slice<const uint8_t> pk_bytes);

}  // namespace tachyon::halo2_api::bn254

#endif  // HALO2_PROOFS_INCLUDE_BN254_PROVING_KEY_H_
