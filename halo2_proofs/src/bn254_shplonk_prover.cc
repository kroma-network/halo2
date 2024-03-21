#include "halo2_proofs/include/bn254_shplonk_prover.h"

#include <tachyon/c/math/polynomials/univariate/bn254_univariate_evaluation_domain.h>

#include "halo2_proofs/src/bn254.rs.h"
#include "halo2_proofs/src/rust_vec.h"

namespace tachyon::halo2_api::bn254 {

SHPlonkProver::SHPlonkProver(uint8_t transcript_type, uint32_t k, const Fr& s)
    : prover_(tachyon_halo2_bn254_shplonk_prover_create_from_unsafe_setup(
          transcript_type, k, reinterpret_cast<const tachyon_bn254_fr*>(&s))) {}

SHPlonkProver::SHPlonkProver(uint8_t transcript_type, uint32_t k,
                             const uint8_t* params, size_t params_len)
    : prover_(tachyon_halo2_bn254_shplonk_prover_create_from_params(
          transcript_type, k, params, params_len)) {}

SHPlonkProver::~SHPlonkProver() {
  tachyon_halo2_bn254_shplonk_prover_destroy(prover_);
}

uint32_t SHPlonkProver::k() const {
  return tachyon_halo2_bn254_shplonk_prover_get_k(prover_);
}

uint64_t SHPlonkProver::n() const {
  return static_cast<uint64_t>(
      tachyon_halo2_bn254_shplonk_prover_get_n(prover_));
}

rust::Box<G2AffinePoint> SHPlonkProver::s_g2() const {
  return rust::Box<G2AffinePoint>::from_raw(
      reinterpret_cast<G2AffinePoint*>(new tachyon_bn254_g2_affine(
          *tachyon_halo2_bn254_shplonk_prover_get_s_g2(prover_))));
}

rust::Box<G1JacobianPoint> SHPlonkProver::commit(const Poly& poly) const {
  return rust::Box<G1JacobianPoint>::from_raw(
      reinterpret_cast<G1JacobianPoint*>(
          tachyon_halo2_bn254_shplonk_prover_commit(prover_, poly.poly())));
}

rust::Box<G1JacobianPoint> SHPlonkProver::commit_lagrange(
    const Evals& evals) const {
  return rust::Box<G1JacobianPoint>::from_raw(
      reinterpret_cast<G1JacobianPoint*>(
          tachyon_halo2_bn254_shplonk_prover_commit_lagrange(prover_,
                                                             evals.evals())));
}

std::unique_ptr<Evals> SHPlonkProver::empty_evals() const {
  return std::make_unique<Evals>(
      tachyon_bn254_univariate_evaluation_domain_empty_evals(
          tachyon_halo2_bn254_shplonk_prover_get_domain(prover_)));
}

std::unique_ptr<RationalEvals> SHPlonkProver::empty_rational_evals() const {
  return std::make_unique<RationalEvals>(
      tachyon_bn254_univariate_evaluation_domain_empty_rational_evals(
          tachyon_halo2_bn254_shplonk_prover_get_domain(prover_)));
}

std::unique_ptr<Poly> SHPlonkProver::ifft(const Evals& evals) const {
  // NOTE(chokobole): The zero degrees might be removed. This might cause an
  // unexpected error if you use this carelessly. Since this is only used to
  // compute instance polynomial and this is used only in Tachyon side, so it's
  // fine.
  return std::make_unique<Poly>(tachyon_bn254_univariate_evaluation_domain_ifft(
      tachyon_halo2_bn254_shplonk_prover_get_domain(prover_), evals.evals()));
}

void SHPlonkProver::batch_evaluate(
    rust::Slice<const std::unique_ptr<RationalEvals>> rational_evals,
    rust::Slice<std::unique_ptr<Evals>> evals) const {
  for (size_t i = 0; i < rational_evals.size(); ++i) {
    evals[i] = std::make_unique<Evals>(
        tachyon_bn254_univariate_rational_evaluations_batch_evaluate(
            rational_evals[i]->evals()));
  }
}

void SHPlonkProver::set_rng(rust::Slice<const uint8_t> state) {
  tachyon_halo2_bn254_shplonk_prover_set_rng_state(prover_, state.data(),
                                                   state.size());
}

void SHPlonkProver::set_transcript(rust::Slice<const uint8_t> state) {
  tachyon_halo2_bn254_shplonk_prover_set_transcript_state(prover_, state.data(),
                                                          state.size());
}

void SHPlonkProver::set_extended_domain(const ProvingKey& pk) {
  tachyon_halo2_bn254_shplonk_prover_set_extended_domain(prover_, pk.pk());
}

void SHPlonkProver::create_proof(ProvingKey& key,
                                 rust::Slice<InstanceSingle> instance_singles,
                                 rust::Slice<AdviceSingle> advice_singles,
                                 rust::Slice<const Fr> challenges) {
  tachyon_bn254_blinder* blinder =
      tachyon_halo2_bn254_shplonk_prover_get_blinder(prover_);
  const tachyon_bn254_plonk_verifying_key* vk =
      tachyon_bn254_plonk_proving_key_get_verifying_key(key.pk());
  const tachyon_bn254_plonk_constraint_system* cs =
      tachyon_bn254_plonk_verifying_key_get_constraint_system(vk);
  uint32_t blinding_factors =
      tachyon_bn254_plonk_constraint_system_compute_blinding_factors(cs);
  tachyon_halo2_bn254_blinder_set_blinding_factors(blinder, blinding_factors);

  size_t num_circuits = instance_singles.size();

  tachyon_halo2_bn254_argument_data* data =
      tachyon_halo2_bn254_argument_data_create(num_circuits);

  tachyon_halo2_bn254_argument_data_reserve_challenges(data, challenges.size());
  for (size_t i = 0; i < challenges.size(); ++i) {
    tachyon_halo2_bn254_argument_data_add_challenge(
        data, reinterpret_cast<const tachyon_bn254_fr*>(&challenges[i]));
  }

  size_t num_bytes = sizeof(RustVec);
  uint8_t* advice_single_data =
      reinterpret_cast<uint8_t*>(advice_singles.data());
  uint8_t* instance_single_data =
      reinterpret_cast<uint8_t*>(instance_singles.data());
  for (size_t i = 0; i < num_circuits; ++i) {
    RustVec vec;
    vec.Read(advice_single_data);
    size_t num_advice_columns = vec.length;
    uintptr_t* advice_columns_ptr = reinterpret_cast<uintptr_t*>(vec.ptr);
    tachyon_halo2_bn254_argument_data_reserve_advice_columns(
        data, i, num_advice_columns);
    for (size_t j = 0; j < num_advice_columns; ++j) {
      tachyon_halo2_bn254_argument_data_add_advice_column(
          data, i, reinterpret_cast<Evals*>(advice_columns_ptr[j])->release());
    }
    advice_single_data += num_bytes;

    vec.Read(&advice_single_data[0]);
    size_t num_blinds = vec.length;
    const tachyon_bn254_fr* blinds_ptr =
        reinterpret_cast<const tachyon_bn254_fr*>(vec.ptr);
    tachyon_halo2_bn254_argument_data_reserve_advice_blinds(data, i,
                                                            num_blinds);
    for (size_t j = 0; j < num_blinds; ++j) {
      tachyon_halo2_bn254_argument_data_add_advice_blind(data, i,
                                                         &blinds_ptr[j]);
    }
    advice_single_data += num_bytes;

    vec.Read(&instance_single_data[0]);
    size_t num_instance_columns = vec.length;
    uintptr_t* instance_columns_ptr = reinterpret_cast<uintptr_t*>(vec.ptr);
    tachyon_halo2_bn254_argument_data_reserve_instance_columns(
        data, i, num_instance_columns);
    for (size_t j = 0; j < num_instance_columns; ++j) {
      tachyon_halo2_bn254_argument_data_add_instance_column(
          data, i,
          reinterpret_cast<Evals*>(instance_columns_ptr[j])->release());
    }
    instance_single_data += num_bytes;

    vec.Read(&instance_single_data[0]);
    uintptr_t* instance_poly_ptr = reinterpret_cast<uintptr_t*>(vec.ptr);
    tachyon_halo2_bn254_argument_data_reserve_instance_polys(
        data, i, num_instance_columns);
    for (size_t j = 0; j < num_instance_columns; ++j) {
      tachyon_halo2_bn254_argument_data_add_instance_poly(
          data, i, reinterpret_cast<Poly*>(instance_poly_ptr[j])->release());
    }
    instance_single_data += num_bytes;
  }

  tachyon_halo2_bn254_shplonk_prover_create_proof(prover_, key.pk(), data);
  tachyon_halo2_bn254_argument_data_destroy(data);
}

rust::Vec<uint8_t> SHPlonkProver::get_proof() const {
  size_t proof_len;
  tachyon_halo2_bn254_shplonk_prover_get_proof(prover_, nullptr, &proof_len);
  rust::Vec<uint8_t> proof;
  // NOTE(chokobole): |rust::Vec<uint8_t>| doesn't have |resize()|.
  proof.reserve(proof_len);
  for (size_t i = 0; i < proof_len; ++i) {
    proof.push_back(0);
  }
  tachyon_halo2_bn254_shplonk_prover_get_proof(prover_, proof.data(),
                                               &proof_len);
  return proof;
}

std::unique_ptr<SHPlonkProver> new_shplonk_prover(uint8_t transcript_type,
                                                  uint32_t k, const Fr& s) {
  return std::make_unique<SHPlonkProver>(transcript_type, k, s);
}

std::unique_ptr<SHPlonkProver> new_shplonk_prover_from_params(
    uint8_t transcript_type, uint32_t k, rust::Slice<const uint8_t> params) {
  return std::make_unique<SHPlonkProver>(transcript_type, k, params.data(),
                                         params.size());
}

rust::Box<Fr> ProvingKey::transcript_repr_shplonk(const SHPlonkProver& prover) {
  tachyon_halo2_bn254_shplonk_prover_set_transcript_repr(prover.prover(), pk_);
  tachyon_bn254_fr* ret = new tachyon_bn254_fr;
  tachyon_bn254_fr repr = tachyon_bn254_plonk_verifying_key_get_transcript_repr(
      tachyon_bn254_plonk_proving_key_get_verifying_key(pk_));
  memcpy(ret->limbs, repr.limbs, sizeof(uint64_t) * 4);
  return rust::Box<Fr>::from_raw(reinterpret_cast<Fr*>(ret));
}

}  // namespace tachyon::halo2_api::bn254
