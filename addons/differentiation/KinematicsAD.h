#ifndef RBDL_KINEMATICS_AD_h
#define RBDL_KINEMATICS_AD_h

#include "ModelAD.h"

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------
namespace AD {
// -----------------------------------------------------------------------------

RBDL_DLLAPI void UpdateKinematicsCustom (
        Model &model,
        ADModel & ad_model,
        Math::VectorNd const * q,
        Math::MatrixNd const * q_dirs,
        Math::VectorNd const * qdot,
        Math::MatrixNd const * qdot_dirs,
        Math::VectorNd const * qddot,
        Math::MatrixNd const * qddot_dirs);

RBDL_DLLAPI void UpdateKinematics (
        Model & model,
        ADModel & ad_model,
        Math::VectorNd const & q,
        Math::MatrixNd const & q_dirs,
        Math::VectorNd const & qdot,
        Math::MatrixNd const & qdot_dirs,
        Math::VectorNd const & qddot,
        Math::MatrixNd const & qddot_dirs);

RBDL_DLLAPI Math::Vector3d CalcPointAcceleration (
        Model & model,
        ADModel & ad_model,
        Math::VectorNd const & q,
        Math::MatrixNd const & q_dirs,
        Math::VectorNd const & qdot,
        Math::MatrixNd const & qdot_dirs,
        Math::VectorNd const & qddot,
        Math::MatrixNd const & qddot_dirs,
        unsigned int body_id,
        Math::Vector3d const & point_position,
        Math::MatrixNd const & fd_derivative,
        bool update_kinematics = true);

// -----------------------------------------------------------------------------
} // namespace AD
// -----------------------------------------------------------------------------
} // namespace RigidBodyDynamics
// -----------------------------------------------------------------------------

#endif
