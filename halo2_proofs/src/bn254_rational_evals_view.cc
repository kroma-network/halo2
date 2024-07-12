#include <halo2_proofs/include/bn254_rational_evals.h>

#include "halo2_proofs/src/bn254.rs.h"

namespace tachyon::halo2_api::bn254 {

RationalEvalsView::RationalEvalsView(
    tachyon_bn254_univariate_rational_evaluations* evals, size_t start,
    size_t len)
    : evals_(evals), start_(start), len_(len) {}

void RationalEvalsView::set_zero(size_t idx) {
  tachyon_bn254_univariate_rational_evaluations_set_zero(evals_, start_ + idx);
}

void RationalEvalsView::set_trivial(size_t idx, const Fr& numerator) {
  tachyon_bn254_univariate_rational_evaluations_set_trivial(
      evals_, start_ + idx,
      reinterpret_cast<const tachyon_bn254_fr*>(&numerator));
}

void RationalEvalsView::set_rational(size_t idx, const Fr& numerator,
                                     const Fr& denominator) {
  tachyon_bn254_univariate_rational_evaluations_set_rational(
      evals_, start_ + idx,
      reinterpret_cast<const tachyon_bn254_fr*>(&numerator),
      reinterpret_cast<const tachyon_bn254_fr*>(&denominator));
}

void RationalEvalsView::evaluate(size_t idx, Fr& value) const {
  tachyon_bn254_univariate_rational_evaluations_evaluate(
      evals_, start_ + idx, reinterpret_cast<tachyon_bn254_fr*>(&value));
}

}  // namespace tachyon::halo2_api::bn254
