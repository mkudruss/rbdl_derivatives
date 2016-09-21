#include "ContactsFD.h"

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------
namespace FD {
// -----------------------------------------------------------------------------

using namespace RigidBodyDynamics::Math;
using namespace std;

RBDL_DLLAPI
void CalcContactJacobian(
        Model &model,
        ADModel &ad_model,
        const Math::VectorNd &Q,
        const Math::VectorNd &Q_dirs,
        const ConstraintSet &CS,
        ADConstraintSet &ad_CS,
        Math::MatrixNd &G,
        std::vector<Math::MatrixNd> &G_dirs,
        bool update_kinematics
) {
    unsigned int ndirs = Q_dirs.cols();
    ad_model.resize_directions(ndirs);

    // temporary evaluation at current point
    CalcContactJacobian(model, Q, CS, G, update_kinematics);

    double h = 1e-8;
    MatrixNd G_temp = MatrixNd::Zero (3, model.dof_count);
    for (unsigned int idir = 0; idir < ndirs; idir++) {
        VectorNd Q_dir = Q_dirs.block(0, idir, model.q_size, 1);

        // temporary evaluation at perturbed point
        CalcContactJacobian(
            model, Q + h * Q_dir, CS, G_temp, update_kinematics
        );

        // calculate finite difference
        G_dirs[idir] = (G_temp - G) / h;
    }
}

// -----------------------------------------------------------------------------
} // namespace FD
// -----------------------------------------------------------------------------
} // namespace RigidBodyDynamics
// -----------------------------------------------------------------------------
