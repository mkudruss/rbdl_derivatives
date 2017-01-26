#include <UnitTest++.h>

#include <iostream>
#include <iomanip>

#include "rbdl/Model.h"
#include "rbdl/Kinematics.h"
#include "rbdl/rbdl_mathutils.h"

#include "ConstraintsAD.h"
#include "ConstraintsFD.h"

#include "ModelCheckADvsFD.h"

#include "Fixtures.h"

using namespace std;
using namespace RigidBodyDynamics;
using namespace RigidBodyDynamics::Math;

double const TEST_PREC = 1.0e-8;

// -----------------------------------------------------------------------------

template <typename T>
void CalcContactJacobianTemplate(
    T & obj,
    unsigned trial_count,
    double array_close_prec) {
  Model   ad_model   = obj.model;
  Model   fd_model   = obj.model;
  ADModel ad_d_model = obj.ad_model;
  ADModel fd_d_model = obj.ad_model;

  ConstraintSet ad_cs = obj.constraint_set;
  ConstraintSet fd_cs = obj.constraint_set;
  ADConstraintSet ad_d_cs = obj.ad_constraint_set;
  ADConstraintSet fd_d_cs = obj.ad_constraint_set;

  // set up input quantities
  int const nq = ad_model.dof_count;
  unsigned const ndirs = nq;
  bool update_kinematics = true;

  VectorNd q = VectorNd::Zero(nq);
  MatrixNd q_dirs = MatrixNd::Identity(nq, nq);
  for (unsigned trial = 0; trial < trial_count; trial++) {
    // set up no output quantities
    MatrixNd G = MatrixNd::Zero (3, nq);
    MatrixNd G_ad = MatrixNd::Zero (3, nq);
    MatrixNd G_fd = MatrixNd::Zero (3, nq);

    // set up derivative output quantities
    vector<MatrixNd> ad_d_G (ndirs, G_ad);
    vector<MatrixNd> fd_d_G (ndirs, G_fd);

    // call nominal version
    CalcConstraintsJacobian(ad_model, q, ad_cs, G, update_kinematics);

    AD::CalcConstraintsJacobian(
          ad_model, ad_d_model,
          q, q_dirs,
          ad_cs, ad_d_cs,
          G_ad, ad_d_G,
          update_kinematics
          );

    FD::CalcConstraintsJacobian(
          fd_model, &fd_d_model,
          q, q_dirs,
          fd_cs, fd_d_cs,
          G_fd, fd_d_G
          );

    checkModelsADvsFD(ndirs, ad_model, ad_d_model, fd_model, fd_d_model);
    checkConstraintSetsADvsFD(ndirs, ad_cs, ad_d_cs, fd_cs, fd_d_cs);

    CHECK_ARRAY_CLOSE(G.data(), G_ad.data(), G.size(), array_close_prec);
    CHECK_ARRAY_CLOSE(G.data(), G_fd.data(), G.size(), array_close_prec);
    CHECK_ARRAY_CLOSE(G_fd.data(), G_ad.data(), G.size(), array_close_prec);

    for (unsigned idir = 0; idir < ndirs; idir++) {
      CHECK_ARRAY_CLOSE(
            fd_d_G[idir].data(),
            ad_d_G[idir].data(),
            ad_d_G[idir].size(),
            array_close_prec);
    }
    q.setRandom();
    q_dirs.setRandom();
  }
}

TEST_FIXTURE (FixedBase6DoF, FixedBase6DoFCalcContactJacobian) {
  // add contacts and bind them to constraint set
  constraint_set.AddContactConstraint (contact_body_id, Vector3d (1., 0., 0.), contact_normal);
  constraint_set.AddContactConstraint (contact_body_id, Vector3d (0., 1., 0.), contact_normal);
  constraint_set.Bind (model);
  ad_constraint_set = ADConstraintSet(constraint_set, model.dof_count);
  CalcContactJacobianTemplate(*this, 10, 1e-5);
}

// -----------------------------------------------------------------------------

template <typename T>
void CalcContactSystemVariablesTemplate(
    T & obj,
    unsigned trial_count) {
  Model   model      = obj.model;
  Model   ad_model   = obj.model;
  Model   fd_model   = obj.model;
  ADModel ad_d_model = obj.ad_model;
  ADModel fd_d_model = obj.ad_model;

  ConstraintSet cs = obj.constraint_set;
  ConstraintSet ad_cs = obj.constraint_set;
  ConstraintSet fd_cs = obj.constraint_set;
  ADConstraintSet ad_d_cs = obj.ad_constraint_set;
  ADConstraintSet fd_d_cs = obj.ad_constraint_set;

  // set up input quantities
  int const nq = ad_model.dof_count;
  unsigned const ndirs = nq;

  VectorNd q = VectorNd::Zero(nq);
  MatrixNd q_dirs = MatrixNd::Identity(nq, nq);
  VectorNd qd = VectorNd::Zero(nq);
  MatrixNd qd_dirs = MatrixNd::Identity(nq, nq);
  VectorNd tau = VectorNd::Zero(nq);
  MatrixNd tau_dirs = MatrixNd::Identity(nq, nq);
  for (unsigned trial = 0; trial < trial_count; trial++) {
    CalcConstrainedSystemVariables(model, q, qd, tau, cs);
    AD::CalcConstrainedSystemVariables(ad_model, ad_d_model,
                                   q, q_dirs, qd, qd_dirs,
                                   tau, tau_dirs, ad_cs, ad_d_cs);
    FD::CalcConstrainedSystemVariables(fd_model, &fd_d_model,
                                       q, q_dirs, qd, qd_dirs,
                                       tau, tau_dirs, fd_cs, fd_d_cs);

    checkModelsADvsFD(ndirs, ad_model, ad_d_model, fd_model, fd_d_model);
    checkConstraintSetsADvsFD(ndirs, ad_cs, ad_d_cs, fd_cs, fd_d_cs);

    q.setRandom();
    qd.setRandom();
    tau.setRandom();
    q_dirs.setRandom();
    qd_dirs.setRandom();
    tau_dirs.setRandom();
  }
}

TEST_FIXTURE (FixedBase6DoF, FixedBase6DoFCalcContactSystemVariables) {
  // add contacts and bind them to constraint set
  constraint_set.AddContactConstraint (contact_body_id, Vector3d (1., 0., 0.), contact_normal);
  constraint_set.AddContactConstraint (contact_body_id, Vector3d (0., 1., 0.), contact_normal);
  constraint_set.Bind (model);
  ad_constraint_set = ADConstraintSet(constraint_set, model.dof_count);
  CalcContactSystemVariablesTemplate(*this, 10);
}

// -----------------------------------------------------------------------------

template <typename T>
void ForwardDynamicsConstraintsDirectTemplate(
    T & obj,
    unsigned trial_count,
    double array_close_prec) {
  Model   model    = obj.model;
  Model   ad_model = obj.model;
  Model   fd_model = obj.model;
  ADModel ad_d_model = obj.ad_model;
  ADModel fd_d_model = obj.ad_model;

  ConstraintSet cs = obj.constraint_set;
  ConstraintSet ad_cs = obj.constraint_set;
  ConstraintSet fd_cs = obj.constraint_set;
  ADConstraintSet ad_d_cs = obj.ad_constraint_set;
  ADConstraintSet fd_d_cs = obj.ad_constraint_set;

  // set up input quantities
  int const nq = ad_model.dof_count;
  unsigned const ndirs = nq;

  VectorNd qdd    = VectorNd::Zero(nq);
  VectorNd ad_qdd = VectorNd::Zero(nq);
  VectorNd fd_qdd = VectorNd::Zero(nq);
  MatrixNd ad_qdd_dirs = MatrixNd::Zero(nq, nq);
  MatrixNd fd_qdd_dirs = MatrixNd::Zero(nq, nq);

  MatrixNd q_dirs = MatrixNd::Identity(nq, nq);
  MatrixNd qd_dirs = MatrixNd::Identity(nq, nq);
  MatrixNd tau_dirs = MatrixNd::Identity(nq, nq);

  for (unsigned trial = 0; trial < trial_count; trial++) {
    VectorNd q = VectorNd::Random(nq);
    VectorNd qd = VectorNd::Random(nq);
    VectorNd tau = VectorNd::Random(nq);

    ForwardDynamicsConstraintsDirect(model, q, qd, tau, cs, qdd);

    AD::ForwardDynamicsConstraintsDirect(
          ad_model, ad_d_model,
          q, q_dirs,
          qd, qd_dirs,
          tau, tau_dirs,
          ad_cs, ad_d_cs,
          ad_qdd, ad_qdd_dirs);

    FD::ForwardDynamicsContactsDirect(
          fd_model, &fd_d_model,
          q, q_dirs,
          qd, qd_dirs,
          tau, tau_dirs,
          fd_cs, fd_d_cs,
          fd_qdd, fd_qdd_dirs);

    checkModelsADvsFD(ndirs, ad_model, ad_d_model, fd_model, fd_d_model);
    checkConstraintSetsADvsFD(ndirs, ad_cs, ad_d_cs, fd_cs, fd_d_cs);

    CHECK_ARRAY_CLOSE(qdd.data(), ad_qdd.data(), nq, array_close_prec);
    CHECK_ARRAY_CLOSE(qdd.data(), fd_qdd.data(), nq, array_close_prec);
    CHECK_ARRAY_CLOSE(ad_qdd_dirs.data(), fd_qdd_dirs.data(), nq,
                      array_close_prec);

    q_dirs.setRandom();
    qd_dirs.setRandom();
    tau_dirs.setRandom();
  }
}

TEST_FIXTURE (FixedBase6DoF, FixedBase6DoFForwardDynamicsContactsDirect) {
  // add contacts and bind them to constraint set
  constraint_set.AddContactConstraint (contact_body_id, Vector3d (1., 0., 0.), contact_normal);
  constraint_set.AddContactConstraint (contact_body_id, Vector3d (0., 1., 0.), contact_normal);
  constraint_set.Bind (model);
  ad_constraint_set = ADConstraintSet(constraint_set, model.dof_count);
  ForwardDynamicsConstraintsDirectTemplate(*this, 10, 1e-4);
}

// -----------------------------------------------------------------------------

TEST (SolveContactSystemDirectTest) {
  vector<LinearSolver> lss(6);
  lss[0] = LinearSolverPartialPivLU;
  lss[1] = LinearSolverColPivHouseholderQR;
  lss[2] = LinearSolverHouseholderQR;
  lss[3] = LinearSolverPartialPivLU;
  lss[4] = LinearSolverColPivHouseholderQR;
  lss[5] = LinearSolverHouseholderQR;

  for (unsigned i = 0; i < lss.size() * 3; i++) {
    int nc = 3;
    int nm = 11;
    int ndirs = 20;
    LinearSolver ls = lss[i % lss.size()];

    MatrixNd H      = MatrixNd::Identity(nm, nm);
    MatrixNd G      = MatrixNd::Random(nc, nm);
    G.setRandom(nc, nm);
    VectorNd c      = VectorNd::Ones(nm, 1);
    VectorNd gamma  = VectorNd::Ones(nc, 1);
    vector<MatrixNd> H_dirs(ndirs, H);
    vector<MatrixNd> G_dirs(ndirs, G);
    MatrixNd c_dirs = MatrixNd::Random(nm, ndirs);
    MatrixNd gamma_dirs = MatrixNd::Random(nc, ndirs);

    for (int ir = 0; ir < nm; ir++) {
      for (int ic = ir + 1; ic < nm; ic++) {
        H(ir, ic) = MatrixNd::Random(1, 1)(0, 0);
      }
    }

    for (int j = 0; j < nc; j++) {
      H_dirs[j] = MatrixNd::Random(nm, nm);
      G_dirs[j] = MatrixNd::Random(nc, nm);
    }

    MatrixNd A(nm + nc, nm + nc);
    VectorNd b(nm + nc);
    VectorNd x(nm + nc);
    VectorNd qddot(nm);
    VectorNd lambda(nm);
    A.setZero();
    b.setZero();
    SolveConstrainedSystemDirect(H, G, c, gamma, qddot, lambda, A, b, x, ls);

    MatrixNd ad_A(nm + nc, nm + nc);
    VectorNd ad_b(nm + nc);
    VectorNd ad_x(nm + nc);
    ad_A.setZero();
    ad_b.setZero();
    vector<MatrixNd> ad_A_dirs(ndirs, ad_A);
    MatrixNd ad_b_dirs(nm + nc, ndirs);
    MatrixNd ad_x_dirs(nm + nc, ndirs);
    ad_b_dirs.setZero();
    ad_x_dirs.setZero();
    AD::SolveConstrainedSystemDirect(H, H_dirs, G, G_dirs, c, c_dirs, gamma,
                                 gamma_dirs, ad_A, ad_A_dirs, ad_b, ad_b_dirs,
                                 ad_x, ad_x_dirs, ls, ndirs);

    MatrixNd fd_A(nm + nc, nm + nc);
    VectorNd fd_b(nm + nc);
    VectorNd fd_x(nm + nc);
    fd_A.setZero();
    fd_b.setZero();
    vector<MatrixNd> fd_A_dirs(ndirs, ad_A);
    MatrixNd fd_b_dirs(nm + nc, ndirs);
    MatrixNd fd_x_dirs(nm + nc, ndirs);
    fd_b_dirs.setZero();
    fd_x_dirs.setZero();
    FD::SolveContactSystemDirect(H, H_dirs, G, G_dirs, c, c_dirs, gamma,
                                 gamma_dirs, fd_A, fd_A_dirs, fd_b, fd_b_dirs,
                                 fd_x, fd_x_dirs, ls, ndirs);

    CHECK_ARRAY_CLOSE(x.data(), ad_x.data(), nm + nc, 1e-8);
    CHECK_ARRAY_CLOSE(x.data(), fd_x.data(), nm + nc, 1e-8);
    CHECK_ARRAY_CLOSE(b.data(), ad_b.data(), nm + nc, 1e-8);
    CHECK_ARRAY_CLOSE(b.data(), fd_b.data(), nm + nc, 1e-8);
    CHECK_ARRAY_CLOSE(A.data(), fd_A.data(), (nm + nc) * (nm + nc), 1e-8);
    CHECK_ARRAY_CLOSE(A.data(), ad_A.data(), (nm + nc) * (nm + nc), 1e-8);
    for (int idir = 0; idir < ndirs; idir++) {
      CHECK_ARRAY_CLOSE(ad_A_dirs[idir].data(), fd_A_dirs[idir].data(),
                        (nm + nc) * (nm + nc), 1e-4);
    }

    MatrixNd ad_cwise_normalized = ad_x_dirs.cwiseQuotient(ad_x_dirs.cwiseAbs());
    MatrixNd fd_cwise_normalized = fd_x_dirs.cwiseQuotient(fd_x_dirs.cwiseAbs());

    CHECK_ARRAY_CLOSE(ad_cwise_normalized.data(), fd_cwise_normalized.data(),
                      (nm + nc) * ndirs,
                      1e-5);
    CHECK_ARRAY_CLOSE(ad_b_dirs.data(), fd_b_dirs.data(), (nm + nc) * ndirs,
                      1e-4);
  }
}

// -----------------------------------------------------------------------------

template <typename T>
void ComputeContactImpulsesDirectTestTemplate(
    T & obj,
    unsigned trial_count,
    double const array_close_prec) {
  Model   model      = obj.model;
  Model   ad_model   = obj.model;
  Model   fd_model   = obj.model;
  ADModel ad_d_model = obj.ad_model;
  ADModel fd_d_model = obj.ad_model;

  ConstraintSet   cs      = obj.constraint_set;
  ConstraintSet   ad_cs   = obj.constraint_set;
  ConstraintSet   fd_cs   = obj.constraint_set;
  ADConstraintSet ad_d_cs = obj.ad_constraint_set;
  ADConstraintSet fd_d_cs = obj.ad_constraint_set;

  int const nq       = model.dof_count;
  int const ndirs    = 2 * nq;

  for (unsigned trial = 0; trial < trial_count; trial++) {
    VectorNd q         = VectorNd::Random(nq);
    MatrixNd q_dirs    = MatrixNd::Zero(nq, ndirs);
    q_dirs.block(0, 0, nq, nq) = MatrixNd::Identity(nq, nq);

    VectorNd qd        = VectorNd::Random(nq);
    MatrixNd qd_dirs   = MatrixNd::Zero(nq, ndirs);
    qd_dirs.block(0, nq, nq, nq) = MatrixNd::Identity(nq, nq);

    VectorNd qd_plus(nq);
    VectorNd ad_qd_plus(nq);
    MatrixNd ad_qd_plus_dirs(nq, ndirs);
    VectorNd fd_qd_plus(nq);
    MatrixNd fd_qd_plus_dirs(nq, ndirs);

    ComputeConstraintImpulsesDirect (model, q, qd, cs, qd_plus);
    AD::ComputeConstraintImpulsesDirect(ad_model, ad_d_model,
                                     q, q_dirs, qd, qd_dirs,
                                     ad_cs, ad_d_cs,
                                     ad_qd_plus, ad_qd_plus_dirs);
    FD::ComputeConstraintImpulsesDirect (fd_model, &fd_d_model,
                                         q, q_dirs, qd, qd_dirs,
                                         fd_cs, &fd_d_cs,
                                         fd_qd_plus, fd_qd_plus_dirs);

    checkModelsADvsFD(ndirs, ad_model, ad_d_model, fd_model, fd_d_model);
    checkConstraintSetsADvsFD(ndirs, ad_cs, ad_d_cs, fd_cs, fd_d_cs);

    CHECK_ARRAY_CLOSE(qd_plus.data(), ad_qd_plus.data(), nq, array_close_prec);
    CHECK_ARRAY_CLOSE(qd_plus.data(), fd_qd_plus.data(), nq, array_close_prec);
    CHECK_ARRAY_CLOSE(ad_qd_plus_dirs.data(), fd_qd_plus_dirs.data(),
                      nq * ndirs, array_close_prec);
  }
}

TEST_FIXTURE (FixedBase6DoF, FixedBase6DoFComputeContactImpulsesDirectTest) {
  constraint_set.AddContactConstraint (
        contact_body_id, Vector3d (1., 0., 0.), contact_normal);
  constraint_set.AddContactConstraint (
        contact_body_id, Vector3d (0., 1., 0.), contact_normal);
  constraint_set.Bind (model);
  ad_constraint_set = ADConstraintSet(constraint_set, model.dof_count);
  ComputeContactImpulsesDirectTestTemplate(*this, 5, 1e-5);
}

// -----------------------------------------------------------------------------
