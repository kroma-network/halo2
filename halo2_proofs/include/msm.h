#ifndef HALO2_PROOFS_INCLUDE_MSM_H_
#define HALO2_PROOFS_INCLUDE_MSM_H_

#include "rust/cxx.h"

namespace tachyon {
namespace halo2 {

struct G1MSM;

struct G1Point2;
struct G1JacobianPoint;
struct Fr;

rust::Box<G1MSM> create_g1_msm(uint8_t degree);

void destroy_g1_msm(rust::Box<G1MSM> msm);

rust::Box<G1JacobianPoint> g1_msm(G1MSM* msm, rust::Slice<const G1Point2> bases,
                                  rust::Slice<const Fr> scalars);

}  // namespace halo2
}  // namespace tachyon

#endif  // HALO2_PROOFS_INCLUDE_MSM_H_
