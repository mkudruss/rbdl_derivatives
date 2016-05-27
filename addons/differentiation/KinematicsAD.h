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
        VectorNd const * q,
        MatrixNd const * q_dirs,
        VectorNd const * qdot,
        MatrixNd const * qdot_dirs,
        VectorNd const * qddot,
        MatrixNd const * qddot_dirs);

RBDL_DLLAPI void UpdateKinematics (
        Model & model,
        ADModel & ad_model,
        VectorNd const & q,
        MatrixNd const & q_dirs,
        VectorNd const & qdot,
        MatrixNd const & qdot_dirs,
        VectorNd const & qddot,
        MatrixNd const & qddot_dirs);

RBDL_DLLAPI Vector3d CalcPointAcceleration (
        Model & model,
        ADModel & ad_model,
        VectorNd const & q,
        MatrixNd const & q_dirs,
        VectorNd const & qdot,
        MatrixNd const & qdot_dirs,
        VectorNd const & qddot,
        MatrixNd const & qddot_dirs,
        unsigned int body_id,
        Vector3d const & point_position,
        MatrixNd const & fd_derivative,
        bool update_kinematics = true);

// -----------------------------------------------------------------------------
} // namespace AD
// -----------------------------------------------------------------------------
} // namespace RigidBodyDynamics
// -----------------------------------------------------------------------------

#endif
