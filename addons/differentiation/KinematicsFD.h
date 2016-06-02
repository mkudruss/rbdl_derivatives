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

// -----------------------------------------------------------------------------
} // namespace FD
// -----------------------------------------------------------------------------
} // namespace RigidBodyDynamics
// -----------------------------------------------------------------------------

#endif
