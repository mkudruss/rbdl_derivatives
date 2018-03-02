#ifndef RBDL_JOINT_ED_h
#define RBDL_JOINT_ED_h

#include <rbdl/Model.h>

#include "ModelED.h"
#include "ModelAD.h"

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------
namespace ED {
// -----------------------------------------------------------------------------

RBDL_DLLAPI void jcalc (
  RigidBodyDynamics::Model &model,
  EDModel &ad_model,
  unsigned int joint_id,
  const Math::VectorNd &q,
  const Math::MatrixNd &q_dirs,
  const Math::VectorNd &qdot,
  const Math::MatrixNd &qdot_dirs
);

RBDL_DLLAPI Math::SpatialTransform jcalc_XJ (
  RigidBodyDynamics::Model &model,
  ADModel &ad_model,
  unsigned int joint_id,
  unsigned int idir,
  const Math::VectorNd &q,
  const Math::MatrixNd &q_dirs);

RBDL_DLLAPI void jcalc_X_lambda_S (
  RigidBodyDynamics::Model &model,
  ADModel &ad_model,
  unsigned int joint_id,
  const Math::VectorNd & q,
  const Math::MatrixNd & q_dirs
);

// -----------------------------------------------------------------------------
} // namespace ED
// -----------------------------------------------------------------------------
} // namespace RigidBodyDynamics
// -----------------------------------------------------------------------------

#endif
