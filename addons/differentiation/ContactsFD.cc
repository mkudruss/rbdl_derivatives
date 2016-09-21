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
        const Math::MatrixNd &Q_dirs,
        const ConstraintSet &CS,
        ADConstraintSet &ad_CS,
        Math::MatrixNd &G,
        std::vector<Math::MatrixNd> &G_dirs
) {
    // NOTE we provide new Qs every call, therefore update kinematics!
    bool update_kinematics = true;

    unsigned int ndirs = Q_dirs.cols();

    // temporary evaluation at current point
    CalcContactJacobian(model, Q, CS, G, update_kinematics);
    // std::cout << "In function '" << __func__ << "'!" << std::endl;
    // std::cout << "G = \n" << G << std::endl;
    // std::cout << "Leaving function '" << __func__ << "'!" << std::endl;
    // std::cout << std::endl;

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

RBDL_DLLAPI
void ComputeContactImpulsesDirect (
    Model & model,
    const VectorNd & q,
    const MatrixNd & q_dirs,
    const VectorNd & qdot_minus,
    const MatrixNd & qdot_minus_dirs,
    ConstraintSet & CS,
    VectorNd & qdot_plus,
    MatrixNd & fd_qdot_plus
) {
  int ndirs = q_dirs.cols();
  assert(ndirs == qdot_minus_dirs.cols());
  // assert(ndirs == ad_qdot_plus.cols());

  double h = 1e-8;

  ComputeContactImpulsesDirect(model, q, qdot_minus, CS, qdot_plus);
  for (int i = 0; i < ndirs; i++) {
    VectorNd qh  = q + h * q_dirs.col(i);
    VectorNd qdh = qdot_minus + h * qdot_minus_dirs.col(i);
    VectorNd qdot_plush(model.dof_count);
    ComputeContactImpulsesDirect(model, qh, qdh, CS, qdot_plush);
    fd_qdot_plus.col(i) = (qdot_plush - qdot_plus) / h;
  }
}


// -----------------------------------------------------------------------------
} // namespace FD
// -----------------------------------------------------------------------------
} // namespace RigidBodyDynamics
// -----------------------------------------------------------------------------
