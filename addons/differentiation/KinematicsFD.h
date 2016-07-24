#ifndef RBDL_KINEMATICS_FD_h
#define RBDL_KINEMATICS_FD_h

#include <rbdl/Model.h>
#include <rbdl/rbdl_math.h>

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------
namespace FD {
// -----------------------------------------------------------------------------

RBDL_DLLAPI Math::Vector3d CalcBodyToBaseCoordinatesSingleFunc (
		Model &model,
		const Math::VectorNd &q,
		const Math::MatrixNd &q_dirs,
		unsigned int body_id,
		const Math::Vector3d &point_body_coordinates,
		Math::MatrixNd &out);

RBDL_DLLAPI Math::Vector3d CalcBodyToBaseCoordinates (
        Model &model,
        ADModel &ad_model,
        const Math::VectorNd &q,
        const Math::MatrixNd &q_dirs,
        unsigned int body_id,
        const Math::Vector3d &point_body_coordinates,
        std::vector<Math::Vector3d> *fd_body_to_base_coordinates,
        bool update_kinematics=true);

  RBDL_DLLAPI Math::Matrix3d CalcBodyWorldOrientation (
        Model & model,
        Math::VectorNd const & q,
        Math::MatrixNd const & q_dirs,
        const unsigned int body_id,
        std::vector<Math::Matrix3d> & fd_derivative);

RBDL_DLLAPI Math::Vector3d CalcPointAcceleration (
		Model & model,
		Math::VectorNd const & q,
		Math::MatrixNd const & q_dirs,
		Math::VectorNd const & qdot,
		Math::MatrixNd const & qdot_dirs,
		Math::VectorNd const & qddot,
		Math::MatrixNd const & qddot_dirs,
		unsigned int body_id,
		Math::Vector3d const & point_position,
		Math::MatrixNd & fd_derivative);

RBDL_DLLAPI
void CalcPointJacobian (
        Model &model,
        ADModel &ad_model,
        Math::VectorNd const &Q,
        Math::MatrixNd const &Q_dirs,
        unsigned int body_id,
        Math::Vector3d const &point_position,
        Math::MatrixNd &G,
        std::vector<Math::MatrixNd> &G_dirs,
        bool update_kinematics = true
        );

// -----------------------------------------------------------------------------
} // namespace FD
// -----------------------------------------------------------------------------
} // namespace RigidBodyDynamics
// -----------------------------------------------------------------------------

#endif
