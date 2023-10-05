#ifndef HALO2_PROOFS_INCLUDE_MSM_GPU_H_
#define HALO2_PROOFS_INCLUDE_MSM_GPU_H_

#include "rust/cxx.h"

namespace tachyon {
namespace halo2 {

struct G1MSMGpu;

struct G1Point2;
struct G1JacobianPoint;
struct Fr;

rust::Box<G1MSMGpu> create_g1_msm_gpu(uint8_t degree, int algorithm);

void destroy_g1_msm_gpu(rust::Box<G1MSMGpu> msm);

rust::Box<G1JacobianPoint> g1_msm_gpu(G1MSMGpu* msm,
                                      rust::Slice<const G1Point2> bases,
                                      rust::Slice<const Fr> scalars);

}  // namespace halo2
}  // namespace tachyon

#endif  // HALO2_PROOFS_INCLUDE_MSM_GPU_H_
