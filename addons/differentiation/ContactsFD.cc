#include "ContactsFD.h"
#include "FdModelEntry.h"

using namespace std;

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------
namespace FD {
// -----------------------------------------------------------------------------

using namespace RigidBodyDynamics::Math;

RBDL_DLLAPI
void ForwardDynamicsContactsDirect (Model   & model,
    ADModel * ad_model,
    const VectorNd & q,
    const MatrixNd & q_dirs,
    const VectorNd & qdot,
    const MatrixNd & qdot_dirs,
    const VectorNd & tau,
    const MatrixNd & tau_dirs,
    ConstraintSet   & cs,
    ADConstraintSet & ad_cs,
    Math::VectorNd  & qddot,
    Math::MatrixNd  & ad_qddot) {
  unsigned ndirs = q_dirs.cols();
  assert(ndirs == qdot_dirs.cols());
  assert(ndirs == tau_dirs.cols());
  assert(ndirs == ad_qddot.cols());

  ConstraintSet cs_in = cs;
  ForwardDynamicsConstraintsDirect(model, q, qdot, tau, cs, qddot);

  for (unsigned idir = 0; idir < ndirs; idir++) {

  }

}

RBDL_DLLAPI
void CalcConstraintsJacobian(
    Model & model,
    ADModel * fd_model,
    const VectorNd & q,
    const MatrixNd & q_dirs,
    ConstraintSet & cs,
    ADConstraintSet & fd_CS,
    MatrixNd & G,
    vector<MatrixNd> & G_dirs) {
  unsigned ndirs = q_dirs.cols();
  assert(ndirs == G_dirs.size());
  bool update_kinematics = true;
  ConstraintSet cs_in = cs;
  CalcConstraintsJacobian(model, q, cs, G, update_kinematics);
  double h = 1e-8;
  for (unsigned idir = 0; idir < ndirs; idir++) {
    Model * modelh;
    if (fd_model) {
      modelh = new Model(model);
    } else {
      modelh = &model;
    }
    ConstraintSet csh = cs_in;
    MatrixNd Gh = MatrixNd::Zero (3, model.dof_count);
    VectorNd q_dir = q_dirs.block(0, idir, model.q_size, 1);
    CalcConstraintsJacobian(
          *modelh, q + h * q_dir, csh, Gh, update_kinematics);
    G_dirs[idir] = (Gh - G) / h;
    computeFDEntry(cs, csh, h, idir, fd_CS);
    if (fd_model) {
      computeFDEntry(model, *modelh, h, idir, *fd_model);
      delete modelh;
    }
  }
}

RBDL_DLLAPI
void ComputeConstraintImpulsesDirect (
    Model & model,
    ADModel * fd_model,
    const VectorNd & q,
    const MatrixNd & q_dirs,
    const VectorNd & qdot_minus,
    const MatrixNd & qdot_minus_dirs,
    ConstraintSet & cs,
    ADConstraintSet * fd_cs,
    VectorNd & qdot_plus,
    MatrixNd & fd_qdot_plus
) {
  unsigned ndirs = q_dirs.cols();
  assert(ndirs == qdot_minus_dirs.cols());
  assert(ndirs == fd_qdot_plus.cols());

  double h = 1e-8;
  ConstraintSet cs_in = cs;
  ComputeConstraintImpulsesDirect(model, q, qdot_minus, cs, qdot_plus);
  for (unsigned idir = 0; idir < ndirs; idir++) {
    Model * modelh;
    if (fd_model) {
      modelh = new Model(model);
    } else {
      modelh = &model;
    }
    ConstraintSet csh = cs_in;

    VectorNd qh  = q + h * q_dirs.col(idir);
    VectorNd qdh = qdot_minus + h * qdot_minus_dirs.col(idir);
    VectorNd qdot_plush(model.dof_count);
    ComputeConstraintImpulsesDirect(*modelh, qh, qdh, csh, qdot_plush);
    fd_qdot_plus.col(idir) = (qdot_plush - qdot_plus) / h;

    if (fd_cs) {
      computeFDEntry(cs, csh, h, idir, *fd_cs);
    }
    if (fd_model) {
      computeFDEntry(model, *modelh, h, idir, *fd_model);
      delete modelh;
    }
  }
}

RBDL_DLLAPI
void ComputeConstraintImpulsesDirect (
    Model & model,
    ADModel * fd_model,
    const VectorNd & q,
    const MatrixNd & q_dirs,
    const VectorNd & qdot_minus,
    const MatrixNd & qdot_minus_dirs,
    ConstraintSet & cs,
    ADConstraintSet * fd_cs,
    MatrixNd & fd_b,
    vector<MatrixNd> & fd_A,
    VectorNd & qdot_plus,
    MatrixNd & fd_qdot_plus
) {
  unsigned ndirs = q_dirs.cols();
  assert(ndirs == qdot_minus_dirs.cols());
  assert(ndirs == fd_qdot_plus.cols());

  double h = 1e-8;
  ConstraintSet cs_in = cs;
  ComputeConstraintImpulsesDirect(model, q, qdot_minus, cs, qdot_plus);
  VectorNd b_ref = cs.b;
  MatrixNd A_ref = cs.A;
  for (unsigned idir = 0; idir < ndirs; idir++) {
    Model * modelh;
    if (fd_model) {
      modelh = new Model(model);
    } else {
      modelh = &model;
    }
    ConstraintSet csh = cs_in;

    VectorNd qh  = q + h * q_dirs.col(idir);
    VectorNd qdh = qdot_minus + h * qdot_minus_dirs.col(idir);
    VectorNd qdot_plush(model.dof_count);
    ComputeConstraintImpulsesDirect(*modelh, qh, qdh, csh, qdot_plush);
    fd_qdot_plus.col(idir) = (qdot_plush - qdot_plus) / h;
    fd_b.col(idir) = (cs.b - b_ref) / h;
    fd_A[idir] = (cs.A - A_ref) / h;

    if (fd_cs) {
      computeFDEntry(cs, csh, h, idir, *fd_cs);
    }
    if (fd_model) {
      computeFDEntry(model, *modelh, h, idir, *fd_model);
      delete modelh;
    }
  }
}

RBDL_DLLAPI
void SolveContactSystemDirect (
    const MatrixNd &H,
    const vector<MatrixNd> & H_dirs,
    const MatrixNd &G,
    const vector<MatrixNd> & G_dirs,
    const VectorNd & c,
    const MatrixNd & c_dirs,
    const VectorNd & gamma,
    const MatrixNd & gamma_dirs,
    MatrixNd & A,
    vector<MatrixNd> & A_dirs,
    VectorNd & b,
    MatrixNd & b_dirs,
    VectorNd & x,
    MatrixNd & x_fd,
    LinearSolver & linear_solver,
    int ndirs
) {
  double h = 1e-8;

  VectorNd qddot(H.rows());
  VectorNd lambda(H.rows());

  RigidBodyDynamics::SolveConstrainedSystemDirect(
      H, G, c, gamma, qddot, lambda,
      A, b, x, linear_solver);
  for (int i = 0; i < ndirs; i++) {
    MatrixNd Ah(A);
    VectorNd bh(b);
    VectorNd xh(x);
    RigidBodyDynamics::SolveConstrainedSystemDirect (H + h * H_dirs[i],
        G + h * G_dirs[i], c + h * c_dirs.col(i), gamma + h * gamma_dirs.col(i),
        qddot, lambda, Ah, bh, xh, linear_solver);
    A_dirs[i] = (Ah - A) / h;
    b_dirs.col(i) = (bh - b) / h;
    x_fd.col(i) = (xh - x) / h;
  }
}


// -----------------------------------------------------------------------------
} // namespace FD
// -----------------------------------------------------------------------------
} // namespace RigidBodyDynamics
// -----------------------------------------------------------------------------
