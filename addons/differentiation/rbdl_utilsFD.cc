#include "rbdl_utilsFD.h"

#include <rbdl/rbdl_utils.h>
#include <rbdl/Model.h>

using namespace RigidBodyDynamics::Math;

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------
namespace Utils {
// -----------------------------------------------------------------------------
namespace FD {
// -----------------------------------------------------------------------------

RBDL_DLLAPI void CalcCenterOfMass (
        Model & model,
        VectorNd const & q,
        MatrixNd const & q_dirs,
        VectorNd const & qdot,
        MatrixNd const & qdot_dirs,
        double & mass,
        Vector3d & com,
        MatrixNd & fd_com,
        Vector3d * com_velocity,
        MatrixNd * fd_com_velocity,
        Vector3d * angular_momentum,
        MatrixNd * fd_angular_momentum) {
    assert(q_dirs.cols() == qdot_dirs.cols());
    assert(q_dirs.cols() == fd_com.cols());
    if (com_velocity) {
        assert(q_dirs.cols() == fd_com_velocity->cols());
    }
    if (angular_momentum) {
        assert(q_dirs.cols() == fd_angular_momentum->cols());
    }

    int ndirs = q_dirs.cols();

    RigidBodyDynamics::Utils::CalcCenterOfMass(model, q, qdot, mass, com,
            com_velocity, angular_momentum, true);

    double h = sqrt(1e-16);

    VectorNd q_dir(q);
    VectorNd qdot_dir(qdot);

    for (int i = 0; i < ndirs; i++ ) {
        q_dir = q_dirs.block(0, i, model.q_size, 1);
        qdot_dir = qdot_dirs.block(0, i, model.qdot_size, 1);

        double hd_mass;
        Vector3d hd_com;
        Vector3d * hd_com_velocity = 0;
        Vector3d * hd_angular_momentum = 0;

        if (com_velocity) {
            hd_com_velocity = new Vector3d;
        }

        if (angular_momentum) {
            hd_angular_momentum = new Vector3d;
        }

        RigidBodyDynamics::Utils::CalcCenterOfMass(model, q + h * q_dir,
                qdot + h * qdot_dir, hd_mass, hd_com, hd_com_velocity,
                hd_angular_momentum, true);

        fd_com.block<3,1>(0, i) = (hd_com - com) / h;

        if (com_velocity) {
            fd_com_velocity->block<3,1>(0, i) = (*hd_com_velocity - *com_velocity) / h;
            delete hd_com_velocity;
        }

        if (angular_momentum) {
            fd_angular_momentum->block<3,1>(0, i) = (*hd_angular_momentum - *angular_momentum) / h;
            delete hd_angular_momentum;
        }
    }
}

// -----------------------------------------------------------------------------
} // namespace FD
// -----------------------------------------------------------------------------
} // namespace Utils
// -----------------------------------------------------------------------------
} // namespace RigidBodyDynamics
// -----------------------------------------------------------------------------
