#include <UnitTest++.h>

#include <iostream>

#include "rbdl/Model.h"
#include "rbdl/Kinematics.h"
#include "rbdl/rbdl_mathutils.h"

#include "ContactsAD.h"
#include "ContactsFD.h"

#include "Fixtures.h"

using namespace std;
using namespace RigidBodyDynamics;
using namespace RigidBodyDynamics::Math;

double const TEST_PREC = 1.0e-8;

// -----------------------------------------------------------------------------

template <typename T>
void CalcContactJacobianTemplate(T & obj) {
    // rename for convenience
    Model& model = *(obj.model);
    ADModel& ad_model = *(obj.ad_model);

    ConstraintSet& cs = obj.constraint_set;
    ADConstraintSet& ad_cs = obj.ad_constraint_set;

    // set up input quantities
    int const nq = model.dof_count;
    int const ndirs = nq;

    bool update_kinematics = true;
    VectorNd q = VectorNd::Zero(nq);
    MatrixNd q_dirs = MatrixNd::Identity(nq, nq);

    // set up no output quantities
    MatrixNd G = MatrixNd::Zero (3, nq);
    MatrixNd G_ad = MatrixNd::Zero (3, nq);
    MatrixNd G_fd = MatrixNd::Zero (3, nq);

    // set up derivative output quantities
    vector<MatrixNd> derivative_ad (ndirs, G_ad);
    vector<MatrixNd> derivative_fd (ndirs, G_fd);

    // call nominal version
    RigidBodyDynamics::CalcContactJacobian(model, q, cs, G, update_kinematics);

    // call AD version
    RigidBodyDynamics::AD::CalcContactJacobian(
        model, ad_model,
        q, q_dirs,
        cs, ad_cs,
        G, derivative_ad,
        update_kinematics
    );

    // call FD version
    RigidBodyDynamics::FD::CalcContactJacobian(
        model,
        q, q_dirs,
        cs, ad_cs,
        G, derivative_fd,
        update_kinematics
    );

    // compare nominal results
    CHECK_ARRAY_CLOSE(G_ad.data(), G.data(), G.size(), TEST_PREC);
    CHECK_ARRAY_CLOSE(G_fd.data(), G.data(), G.size(), TEST_PREC);
    CHECK_ARRAY_CLOSE(G_ad.data(), G_fd.data(), G.size(), TEST_PREC);

    for (unsigned idir = 0; idir < ndirs; idir++) {
        CHECK_ARRAY_CLOSE(
            derivative_ad[idir].data(),
            derivative_fd[idir].data(),
            derivative_ad[idir].size(),
            TEST_PREC
        );
    }
}

TEST_FIXTURE (FixedBase6DoF, FixedBase6DoFCalcContactJacobian) {

    // add contacts and bind them to constraint set
    constraint_set.AddConstraint (contact_body_id, Vector3d (1., 0., 0.), contact_normal);
    constraint_set.AddConstraint (contact_body_id, Vector3d (0., 1., 0.), contact_normal);
    constraint_set.Bind (*model);

    ad_constraint_set = ADConstraintSet(constraint_set, model->dof_count);

/*
    VectorNd QDDot_lagrangian = VectorNd::Constant (model->mBodies.size() - 1, 0.);

    ClearLogOutput();

    ForwardDynamicsContactsDirect (*model, Q, QDot, Tau, constraint_set_lagrangian, QDDot_lagrangian);
    ForwardDynamicsContactsKokkevis (*model, Q, QDot, Tau, constraint_set, QDDot_contacts);

    point_accel_lagrangian = CalcPointAcceleration (*model, Q, QDot, QDDot_lagrangian, contact_body_id, contact_point, true);
    point_accel_contacts = CalcPointAcceleration (*model, Q, QDot, QDDot_contacts, contact_body_id, contact_point, true);

    // check whether FDContactsLagrangian and FDContacts match
    CHECK_ARRAY_CLOSE (
            constraint_set_lagrangian.force.data(),
            constraint_set.force.data(),
            constraint_set.size(), TEST_PREC
            );

    // check whether the point accelerations match
    CHECK_ARRAY_CLOSE (point_accel_lagrangian.data(), point_accel_contacts.data(), 3, TEST_PREC);

    // check whether the generalized accelerations match
    CHECK_ARRAY_CLOSE (QDDot_lagrangian.data(), QDDot_contacts.data(), QDDot_lagrangian.size(), TEST_PREC);
    */
    CalcContactJacobianTemplate(*this);
}


/*
TEST_FIXTURE ( CartPendulum, CartPendulumCalcBodyWorldOrientation) {
    CalcBodyWorldOrientationTemplate(*this);
}

TEST_FIXTURE ( Arm2DofX, Arm2DofXCalcBodyWorldOrientation) {
    CalcBodyWorldOrientationTemplate(*this);
}

TEST_FIXTURE ( Arm2DofZ, Arm2DofZCalcBodyWorldOrientation) {
    CalcBodyWorldOrientationTemplate(*this);
}

TEST_FIXTURE ( Arm3DofXZYp, Arm3DofXZYpCalcBodyWorldOrientation) {
    CalcBodyWorldOrientationTemplate(*this);
}

TEST_FIXTURE ( Arm3DofXZZp, Arm3DofXZZpCalcBodyWorldOrientation) {
    CalcBodyWorldOrientationTemplate(*this);
}
*/


// -----------------------------------------------------------------------------

TEST (SolveContactSystemDirectTest) {

  vector<LinearSolver> lss(3);
  lss[0] = LinearSolverPartialPivLU;
  lss[1] = LinearSolverColPivHouseholderQR;
  lss[2] = LinearSolverHouseholderQR;


  for (int i = 0; i < lss.size() * 3; i++) {
    int nc = 3;
    int nm = 11;
    int ndirs = 20;
    LinearSolver ls = lss[i % lss.size()]; //LinearSolverColPivHouseholderQR;

    MatrixNd H      = MatrixNd::Identity(nm, nm);
    MatrixNd G      = MatrixNd::Random(nc, nm);
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
    SolveContactSystemDirect(H, G, c, gamma, qddot, lambda, A, b, x, ls);

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
    AD::SolveContactSystemDirect(H, H_dirs, G, G_dirs, c, c_dirs, gamma,
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
    CHECK_ARRAY_CLOSE(ad_x_dirs.data(), fd_x_dirs.data(), (nm + nc) * ndirs,
                      1e-4);
    CHECK_ARRAY_CLOSE(ad_b_dirs.data(), fd_b_dirs.data(), (nm + nc) * ndirs,
                      1e-4);
  }

}

// -----------------------------------------------------------------------------

template <typename T>
void ComputeContactImpulsesDirectTestTemplate(T & obj,
    double const TEST_PREC = 1e-6) {
  Model   & model    = *(obj.model);
  ADModel & ad_model = *(obj.ad_model);

  ConstraintSet   & cs    = obj.constraint_set;
  ADConstraintSet & ad_cs = obj.ad_constraint_set;

  int const nq       = model.dof_count;
  int const ndirs    = 2 * nq;

  for (int i = 0; i < 5; i++) {
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

    ComputeContactImpulsesDirect(model, q, qd, cs, qd_plus);
    AD::ComputeContactImpulsesDirect(model, ad_model, q, q_dirs, qd, qd_dirs,
                                     cs, ad_cs, ad_qd_plus, ad_qd_plus_dirs);
    FD::ComputeContactImpulsesDirect(model, q, q_dirs, qd, qd_dirs, cs,
                                     fd_qd_plus, fd_qd_plus_dirs);

    CHECK_ARRAY_CLOSE(qd_plus.data(), ad_qd_plus.data(), nq, TEST_PREC);
    CHECK_ARRAY_CLOSE(qd_plus.data(), fd_qd_plus.data(), nq, TEST_PREC);
    CHECK_ARRAY_CLOSE(ad_qd_plus_dirs.data(), fd_qd_plus_dirs.data(),
                      nq * ndirs, TEST_PREC);
  }
}

TEST_FIXTURE (FixedBase6DoF, FixedBase6DoFComputeContactImpulsesDirectTest) {
    constraint_set.AddConstraint (contact_body_id, Vector3d (1., 0., 0.),
                                  contact_normal);
    constraint_set.AddConstraint (contact_body_id, Vector3d (0., 1., 0.),
                                  contact_normal);
    constraint_set.Bind (*model);
    ad_constraint_set = ADConstraintSet(constraint_set, model->dof_count);
    ComputeContactImpulsesDirectTestTemplate(*this);
}

// -----------------------------------------------------------------------------

