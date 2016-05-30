#ifndef RBDL_UTILS_FD_h
#define RBDL_UTILS_FD_h

#include <rbdl/rbdl_mathutils.h>

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------
namespace Utils {
// -----------------------------------------------------------------------------
namespace FD {
// -----------------------------------------------------------------------------

//extern int fd_counter;
//extern std::vector<std::vector<Math::SpatialVector> > fd_v;
//extern std::vector<std::vector<Math::SpatialMatrix> > fd_I;

RBDL_DLLAPI void CalcCenterOfMass (
        Model & model,
        Math::VectorNd const & q,
        Math::MatrixNd const & q_dirs,
        Math::VectorNd const & qdot,
        Math::MatrixNd const & qdot_dirs,
        double & mass,
        Math::Vector3d & com,
        Math::MatrixNd & fd_com,
        Math::Vector3d * com_velocity,
        Math::MatrixNd * fd_com_velocity,
        Math::Vector3d * angular_momentum,
        Math::MatrixNd * fd_angular_momentum);

RBDL_DLLAPI double CalcPotentialEnergy (
        Model & model,
        Math::VectorNd const & q,
        Math::MatrixNd const & q_dirs,
        Math::MatrixNd & fd_pote);

RBDL_DLLAPI double CalcKineticEnergy (
        Model & model,
        Math::VectorNd const & q,
        Math::MatrixNd const & q_dirs,
        Math::VectorNd const & qdot,
        Math::MatrixNd const & qdot_dirs,
        Math::MatrixNd & fd_kine);

// -----------------------------------------------------------------------------
} // namespace FD
// -----------------------------------------------------------------------------
} // namespace Utils
// -----------------------------------------------------------------------------
} // namespace RigidBodyDynamics
// -----------------------------------------------------------------------------

#endif
