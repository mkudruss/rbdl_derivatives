#include "ConstraintsFDC.h"
#include "ModelEntryFDC.h"

#include <cfloat>
#include <cmath>

using namespace std;

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------
namespace FDC {
// -----------------------------------------------------------------------------

using namespace RigidBodyDynamics::Math;

// define constant difference perturbation according to theory,
// i.e. EPS = eps^(1/3) = cubic_root(eps).
// Here, we us cmath cubic root implementation, i.e. cubic_root = cbrt.
// NOTE assumes that directions are normalized to one
const double EPS = cbrt(DBL_EPSILON);
const double EPSx2 = 2*EPS; // for convenience, we directly compute the denominator

RBDL_DLLAPI void CalcConstrainedSystemVariables (
    Model &model,
    ADModel *fd_model, // NULL means execution without fd_model update
    const VectorNd &q,
    const MatrixNd &q_dirs,
    const VectorNd &qdot,
    const MatrixNd &qdot_dirs,
    const VectorNd &tau,
    const MatrixNd &tau_dirs,
    ConstraintSet   &cs,
    ADConstraintSet &fd_cs) {
  unsigned const ndirs = q_dirs.cols();
  assert(ndirs == qdot_dirs.cols());
  assert(ndirs == tau_dirs.cols());

  ConstraintSet const cs_in = cs;

  CalcConstrainedSystemVariables(model, q, qdot, tau, cs);

  for (unsigned idir = 0; idir < ndirs; idir++) {
    Model *modelh = &model;
    VectorNd qh    = q + EPS * q_dirs.col(idir);
    VectorNd qdoth = qdot + EPS * qdot_dirs.col(idir);
    VectorNd tauh  = tau + EPS * tau_dirs.col(idir);
    ConstraintSet csh = cs_in;

    if (fd_model) {
      modelh = new Model(model);
    }

    CalcConstrainedSystemVariables(*modelh, qh, qdoth, tauh, csh);

    // TODO make it analogue to ADMOdel with pointer arithmetic
    // computeFDEntry(cs, csh, EPS, idir, fd_cs);

    if (fd_model) {
      computeFDEntry(model, *modelh, EPS, idir, *fd_model);
      delete modelh;
    }
  }
}


RBDL_DLLAPI void ForwardDynamicsConstraintsDirect (
  Model &model,
  ADModel *fd_model,
  const VectorNd &q,
  const MatrixNd &q_dirs,
  const VectorNd &qdot,
  const MatrixNd &qdot_dirs,
  const VectorNd &tau,
  const MatrixNd &tau_dirs,
  ConstraintSet   &cs,
  ADConstraintSet &fd_cs,
  VectorNd  &qddot,
  MatrixNd  &fd_qddot
) {
  unsigned const ndirs = q_dirs.cols();
  assert(ndirs == qdot_dirs.cols());
  assert(ndirs == tau_dirs.cols());
  assert(ndirs == fd_qddot.cols());

  ConstraintSet const cs_in = cs;

  ForwardDynamicsConstraintsDirect(model, q, qdot, tau, cs, qddot);

  // temporary variables
  VectorNd qddotph (qddot);
  VectorNd qddotmh (qddot);

  for (unsigned idir = 0; idir < ndirs; idir++) {
    Model * modelh;
    if (fd_model) {
      modelh = new Model(model);
    } else {
      modelh = &model;
    }

    ConstraintSet csh = cs_in;

    ForwardDynamicsConstraintsDirect(
      *modelh,
      q + EPS * q_dirs.col(idir),
      qdot + EPS * qdot_dirs.col(idir),
      tau + EPS * tau_dirs.col(idir),
      csh,
      qddotph
    );

    ForwardDynamicsConstraintsDirect(
      *modelh,
      q - EPS * q_dirs.col(idir),
      qdot - EPS * qdot_dirs.col(idir),
      tau - EPS * tau_dirs.col(idir),
      csh,
      qddotmh
    );

    fd_qddot.col(idir) = (qddotph - qddotmh) / EPSx2;

    // TODO add pointer arithmetic for fd_cs
    // computeFDEntry(cs, csh, EPS, idir, fd_cs);

    if (fd_model) {
      computeFDEntry(model, *modelh, EPS, idir, *fd_model);
      delete modelh;
    }
  }
}

RBDL_DLLAPI void ForwardDynamicsContactsKokkevis (
  Model &model,
  ADModel *fd_model,
  const VectorNd &q,
  const MatrixNd &q_dirs,
  const VectorNd &qdot,
  const MatrixNd &qdot_dirs,
  const VectorNd &tau,
  const MatrixNd &tau_dirs,
  ConstraintSet   &cs,
  ADConstraintSet &fd_cs,
  VectorNd  &qddot,
  MatrixNd  &fd_qddot
) {
  unsigned const ndirs = q_dirs.cols();
  assert(ndirs == qdot_dirs.cols());
  assert(ndirs == tau_dirs.cols());
  assert(ndirs == fd_qddot.cols());

  ConstraintSet const cs_in = cs;

  ForwardDynamicsContactsKokkevis(model, q, qdot, tau, cs, qddot);

  // temporary variables
  VectorNd qddotph (qddot);
  VectorNd qddotmh (qddot);

  for (unsigned idir = 0; idir < ndirs; idir++) {
    Model * modelhp;
    Model * modelhm;
    ConstraintSet * cshp;
    ConstraintSet * cshm;
    if (fd_model) {
      modelhp = new Model(model);
      modelhm = new Model(model);
      cshp = new ConstraintSet(cs);
      cshm = new ConstraintSet(cs);
    } else {
      modelhp = &model;
      modelhm = &model;
      cshp = &cs;
      cshm = &cs;
    }

    ConstraintSet csh = cs_in;

    ForwardDynamicsContactsKokkevis(
      *modelhp,
      q + EPS * q_dirs.col(idir),
      qdot + EPS * qdot_dirs.col(idir),
      tau + EPS * tau_dirs.col(idir),
      csh,
      qddotph
    );

    ForwardDynamicsContactsKokkevis(
      *modelhm,
      q - EPS * q_dirs.col(idir),
      qdot - EPS * qdot_dirs.col(idir),
      tau - EPS * tau_dirs.col(idir),
      csh,
      qddotmh
    );

    fd_qddot.col(idir) = (qddotph - qddotmh) / EPSx2;

    // TODO add pointer arithmetic for fd_cs
    // computeFDEntry(cs, csh, EPS, idir, fd_cs);

    if (fd_model) {
      computeFDEntry(*modelhp, *modelhm, EPS, idir, *fd_model);
      computeFDEntry(*cshp, *cshm, EPS, idir, fd_cs);
      delete modelhp;
      delete modelhm;
      delete cshp;
      delete cshm;
    }
  }
}

RBDL_DLLAPI void CalcConstraintsJacobian(
    Model & model,
    ADModel * fd_model,
    const VectorNd & q,
    const MatrixNd & q_dirs,
    ConstraintSet & cs,
    ADConstraintSet & fd_CS,
    MatrixNd & G,
    vector<MatrixNd> & G_dirs) {
  unsigned const ndirs = q_dirs.cols();
  assert(ndirs == G_dirs.size());

  bool const update_kinematics = true;
  ConstraintSet const cs_in = cs;

  CalcConstraintsJacobian(model, q, cs, G, update_kinematics);

  for (unsigned idir = 0; idir < ndirs; idir++) {
    VectorNd qh = q + EPS * q_dirs.col(idir);
    ConstraintSet csh = cs_in;
    MatrixNd Gh = MatrixNd::Zero (G.rows(), G.cols());

    Model *modelh;
    if (fd_model) {
      modelh = new Model(model);
    } else {
      modelh = &model;
    }

    CalcConstraintsJacobian(*modelh, qh, csh, Gh, update_kinematics);
    G_dirs[idir] = (Gh - G) / EPS;

    if (fd_model) {
      computeFDEntry(model, *modelh, EPS, idir, *fd_model);
      delete modelh;
    }
  }
}

RBDL_DLLAPI void ComputeConstraintImpulsesDirect (
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

    VectorNd qh  = q + EPS * q_dirs.col(idir);
    VectorNd qdh = qdot_minus + EPS * qdot_minus_dirs.col(idir);
    VectorNd qdot_plush(model.dof_count);
    ComputeConstraintImpulsesDirect(*modelh, qh, qdh, csh, qdot_plush);
    fd_qdot_plus.col(idir) = (qdot_plush - qdot_plus) / EPS;

    if (fd_cs) {
      computeFDEntry(cs, csh, EPS, idir, *fd_cs);
    }
    if (fd_model) {
      computeFDEntry(model, *modelh, EPS, idir, *fd_model);
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

    VectorNd qh  = q + EPS * q_dirs.col(idir);
    VectorNd qdh = qdot_minus + EPS * qdot_minus_dirs.col(idir);
    VectorNd qdot_plush(model.dof_count);
    ComputeConstraintImpulsesDirect(*modelh, qh, qdh, csh, qdot_plush);
    fd_qdot_plus.col(idir) = (qdot_plush - qdot_plus) / EPS;
    fd_b.col(idir) = (cs.b - b_ref) / EPS;
    fd_A[idir] = (cs.A - A_ref) / EPS;

    if (fd_cs) {
      computeFDEntry(cs, csh, EPS, idir, *fd_cs);
    }
    if (fd_model) {
      computeFDEntry(model, *modelh, EPS, idir, *fd_model);
      delete modelh;
    }
  }
}

RBDL_DLLAPI
void SolveConstrainedSystemDirect (
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
  VectorNd qddot(H.rows());
  VectorNd lambda(H.rows());

  RigidBodyDynamics::SolveConstrainedSystemDirect(
      H, G, c, gamma, qddot, lambda,
      A, b, x, linear_solver);
  for (int i = 0; i < ndirs; i++) {
    MatrixNd Ah(A);
    VectorNd bh(b);
    VectorNd xh(x);
    RigidBodyDynamics::SolveConstrainedSystemDirect (H + EPS * H_dirs[i],
        G + EPS * G_dirs[i], c + EPS * c_dirs.col(i), gamma + EPS * gamma_dirs.col(i),
        qddot, lambda, Ah, bh, xh, linear_solver);
    A_dirs[i] = (Ah - A) / EPS;
    b_dirs.col(i) = (bh - b) / EPS;
    x_fd.col(i) = (xh - x) / EPS;
  }
}


// -----------------------------------------------------------------------------
} // namespace FD
// -----------------------------------------------------------------------------
} // namespace RigidBodyDynamics
// -----------------------------------------------------------------------------
