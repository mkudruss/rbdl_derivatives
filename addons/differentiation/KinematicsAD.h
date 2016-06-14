#ifndef RBDL_KINEMATICS_AD_h
#define RBDL_KINEMATICS_AD_h

#include "ModelAD.h"

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------

RBDL_DLLAPI Math::Vector3d CalcBodyToBaseCoordinatesSingleFunc (
        Model &model,
        const Math::VectorNd &Q,
        unsigned int body_id,
        const Math::Vector3d &point_body_coordinates);

// -----------------------------------------------------------------------------
namespace AD {
// -----------------------------------------------------------------------------

RBDL_DLLAPI Math::Vector3d CalcBodyToBaseCoordinatesSingleFunc (
		Model &model,
        ADModel &ad_model,
        const Math::VectorNd &q,
        const Math::MatrixNd &q_dirs,
        unsigned int body_id,
        const Math::Vector3d &point_body_coordinates,
        Math::MatrixNd & out);

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

RBDL_DLLAPI Math::Matrix3d CalcBodyWorldOrientation (
        Model & model,
        ADModel & ad_model,
        Math::VectorNd const & q,
        Math::MatrixNd const & q_dirs,
        const unsigned int body_id,
        std::vector<Math::Matrix3d> & ad_derivative,
        bool update_kinematics = true);

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
        Math::MatrixNd & ad_derivative,
        bool update_kinematics = true);

// -----------------------------------------------------------------------------
} // namespace AD
// -----------------------------------------------------------------------------
} // namespace RigidBodyDynamics
// -----------------------------------------------------------------------------

#endif
