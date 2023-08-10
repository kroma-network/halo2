// clang-format off
#include <tachyon/c/math/msm/msm_gpu.h>
// clang-format on

#include "halo2_proofs/src/lib.rs.h"
#include "halo2_proofs/include/msm.h"

namespace tachyon {

namespace halo2 {

void init_msm_gpu(uint8_t degree) { tachyon_init_msm_gpu(degree); }

void release_msm_gpu() { tachyon_release_msm_gpu(); }

rust::Box<CppG1Jacobian> msm_gpu(rust::Slice<const CppG1Affine> bases,
                                 rust::Slice<const CppFr> scalars) {
  auto ret = tachyon_msm_g1_point2_gpu(
      reinterpret_cast<const tachyon_bn254_point2*>(bases.data()),
      bases.length(), reinterpret_cast<const tachyon_bn254_fr*>(scalars.data()),
      scalars.length());
  return rust::Box<CppG1Jacobian>::from_raw(reinterpret_cast<CppG1Jacobian*>(ret));
}

}  // namespace halo2
}  // namespace tachyon
