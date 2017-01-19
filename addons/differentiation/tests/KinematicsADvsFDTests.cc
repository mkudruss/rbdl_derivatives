#include <UnitTest++.h>

#include <iostream>

#include "rbdl/Model.h"
#include "rbdl/Kinematics.h"
#include "rbdl/rbdl_mathutils.h"

#include "KinematicsAD.h"
#include "KinematicsFD.h"

#include "Fixtures.h"

#include "ModelCheckADvsFD.h"

using namespace std;
using namespace RigidBodyDynamics;
using namespace RigidBodyDynamics::Math;

const double TEST_PREC = 1.0e-8;

// -----------------------------------------------------------------------------

template<typename T>
void UpdateKinematicsCustomTemplate(T & obj) {
  Model   ad_model = obj.model;
  Model   fd_model = obj.model;
  ADModel ad_d_model = obj.ad_model;
  ADModel fd_d_model = obj.ad_model;

  int const nq = obj.model.dof_count;
  unsigned ndirs = 3 * nq;

  VectorNd q    = VectorNd::Random(nq);
  VectorNd qd   = VectorNd::Random(nq);
  VectorNd qdd  = VectorNd::Random(nq);

  MatrixNd q_dirs(nq, ndirs);
  MatrixNd qd_dirs(nq, ndirs);
  MatrixNd qdd_dirs(nq, ndirs);

  q_dirs.setZero();
  qd_dirs.setZero();
  qdd_dirs.setZero();

  for (int i = 0; i < nq; i++) {
    q_dirs(i, i) = 1.0;
    qd_dirs(i, nq + i) = 1.0;
    qdd_dirs(i, 2 * nq + i) = 1.0;
  }

  AD::UpdateKinematicsCustom(ad_model, ad_d_model, &q, &q_dirs, &qd, &qd_dirs, &qdd, &qdd_dirs);

  FD::UpdateKinematicsCustom(fd_model, fd_d_model, q, q_dirs, qd, qd_dirs, qdd, qdd_dirs);

  checkModelsADvsFD(
        ndirs,
        ad_model, ad_d_model,
        fd_model, fd_d_model);
}

TEST_FIXTURE ( CartPendulum, CartPendulumUpdateKinematicsCustom) {
  UpdateKinematicsCustomTemplate(*this);
}

TEST_FIXTURE ( Arm2DofX, Arm2DofXUpdateKinematicsCustom) {
  UpdateKinematicsCustomTemplate(*this);
}

TEST_FIXTURE ( Arm2DofZ, Arm2DofZUpdateKinematicsCustom) {
  UpdateKinematicsCustomTemplate(*this);
}

TEST_FIXTURE ( Arm3DofXZYp, Arm3DofXZYpUpdateKinematicsCustom) {
  UpdateKinematicsCustomTemplate(*this);
}

TEST_FIXTURE ( Arm3DofXZZp, Arm3DofXZZpUpdateKinematicsCustom) {
  UpdateKinematicsCustomTemplate(*this);
}

// -----------------------------------------------------------------------------

template <typename T>
void CalcPointVelocityTemplate(T & obj) {
    Model   & model = obj.model;
    ADModel & ad_model = obj.ad_model;

    // set up input quantities
    int const nq = model.dof_count;
    int const ndirs = 2 * nq;

    Vector3d point_position = Vector3d::Random();

    VectorNd q    = VectorNd::Random(nq);
    VectorNd qdot = VectorNd::Random(nq);

    MatrixNd q_dirs(nq, 2 * nq);
    q_dirs.setZero();
    q_dirs.block(0, 0, nq, nq) = MatrixNd::Identity(nq, nq);

    MatrixNd qdot_dirs(nq, 2 * nq);
    qdot_dirs.setZero();
    qdot_dirs.block(0, nq, nq, nq) = MatrixNd::Identity(nq, nq);

    // set up no output quantities
    Vector3d G    = Vector3d::Zero ();
    Vector3d G_ad = Vector3d::Zero ();
    Vector3d G_fd = Vector3d::Zero ();

    // set up derivative output quantities
    MatrixNd derivative_ad (3, ndirs);
    MatrixNd derivative_fd (3, ndirs);

    for (unsigned i = 1; i < model.mBodies.size(); i++) {
        unsigned int body_id = model.mBodyNameMap[model.GetBodyName(i)];
        // call nominal version
        G = CalcPointVelocity (
            model, q, qdot, body_id, point_position, true);

        // call FD version
        G_fd = FD::CalcPointVelocity (model, q, q_dirs, qdot, qdot_dirs,
            body_id, point_position, derivative_fd);

        // call AD version
        G_ad = AD::CalcPointVelocity (model, ad_model, q, q_dirs, qdot,
            qdot_dirs, body_id, point_position, derivative_ad, true);

        // compare nominal results
        CHECK_ARRAY_CLOSE(G_ad.data(), G.data(), G.size(), TEST_PREC);
        CHECK_ARRAY_CLOSE(G_fd.data(), G.data(), G.size(), TEST_PREC);
        CHECK_ARRAY_CLOSE(G_ad.data(), G_fd.data(), G.size(), TEST_PREC);
        CHECK_ARRAY_CLOSE(derivative_ad.data(), derivative_fd.data(),
                          3 * ndirs, TEST_PREC * 1e1);
    }
}

TEST_FIXTURE ( CartPendulum, CartPendulumCalcPointVelocity) {
    CalcPointVelocityTemplate(*this);
}

TEST_FIXTURE ( Arm2DofX, Arm2DofXCalcPointVelocity) {
    CalcPointVelocityTemplate(*this);
}

TEST_FIXTURE ( Arm2DofZ, Arm2DofZCalcPointVelocity) {
    CalcPointVelocityTemplate(*this);
}

TEST_FIXTURE ( Arm3DofXZYp, Arm3DofXZYpCalcPointVelocity) {
    CalcPointVelocityTemplate(*this);
}

TEST_FIXTURE ( Arm3DofXZZp, Arm3DofXZZpCalcPointVelocity) {
    CalcPointVelocityTemplate(*this);
}

// -----------------------------------------------------------------------------

template <typename T>
void CalcBaseToBodyCoordinatesTemplate(T & obj) {
    Model   & model = obj.model;
    ADModel & ad_model = obj.ad_model;

    // set up input quantities
    int const nq = model.dof_count;
    int const ndirs = nq;

    bool update_kinematics = true;
    Vector3d point_position = Vector3d::Zero();
    MatrixNd point_position_dirs = MatrixNd::Random(3, ndirs);

    VectorNd q = VectorNd::Zero(nq);
    MatrixNd q_dirs = MatrixNd::Identity(nq, nq);


    // set up no output quantities
    Vector3d G    = Vector3d::Zero ();
    Vector3d G_ad = Vector3d::Zero ();
    Vector3d G_fd = Vector3d::Zero ();

    // set up derivative output quantities
    MatrixNd derivative_ad (3, ndirs);
    MatrixNd derivative_fd (3, ndirs);

    for (unsigned i = 1; i < model.mBodies.size(); i++) {
        unsigned int body_id = model.mBodyNameMap[model.GetBodyName(i)];
        // call nominal version
        G = CalcBodyToBaseCoordinates (
            model, q, body_id, point_position, update_kinematics
        );

        // call FD version
        G_fd = FD::CalcBaseToBodyCoordinates (
            model, q, q_dirs, body_id, point_position, point_position_dirs,
            derivative_fd);

        // call AD version
        G_ad = AD::CalcBaseToBodyCoordinates (
            model, ad_model, q, q_dirs, body_id, point_position,
            point_position_dirs, derivative_ad, update_kinematics);

        // compare nominal results
        CHECK_ARRAY_CLOSE(G_ad.data(), G.data(), G.size(), TEST_PREC);
        CHECK_ARRAY_CLOSE(G_fd.data(), G.data(), G.size(), TEST_PREC);
        CHECK_ARRAY_CLOSE(G_ad.data(), G_fd.data(), G.size(), TEST_PREC);
        CHECK_ARRAY_CLOSE(derivative_ad.data(), derivative_fd.data(),
                          3 * ndirs, TEST_PREC);
    }
}

TEST_FIXTURE ( CartPendulum, CartPendulumCalcBaseToBodyCoordinates) {
    CalcBaseToBodyCoordinatesTemplate(*this);
}

TEST_FIXTURE ( Arm2DofX, Arm2DofXCalcBaseToBodyCoordinates) {
    CalcBaseToBodyCoordinatesTemplate(*this);
}

TEST_FIXTURE ( Arm2DofZ, Arm2DofZCalcBaseToBodyCoordinates) {
    CalcBaseToBodyCoordinatesTemplate(*this);
}

TEST_FIXTURE ( Arm3DofXZYp, Arm3DofXZYpCalcBaseToBodyCoordinates) {
    CalcBaseToBodyCoordinatesTemplate(*this);
}

TEST_FIXTURE ( Arm3DofXZZp, Arm3DofXZZpCalcBaseToBodyCoordinates) {
    CalcBaseToBodyCoordinatesTemplate(*this);
}

// -----------------------------------------------------------------------------

template <typename T>
void CalcBodyToBaseCoordinatesTemplate(T & obj) {
    Model &model = obj.model;
    ADModel &ad_model = obj.ad_model;

    // set up input quantities
    int const nq = model.dof_count;
    int const ndirs = nq;

    bool update_kinematics = true;
    Vector3d point_position = Vector3d::Zero();

    VectorNd q = VectorNd::Zero(nq);
    MatrixNd q_dirs = MatrixNd::Identity(nq, nq);

    // set up no output quantities
    Vector3d G = Vector3d::Zero ();
    Vector3d G_ad = Vector3d::Zero ();
    Vector3d G_fd = Vector3d::Zero ();

    // set up derivative output quantities
    MatrixNd derivative_ad (3, ndirs);
    MatrixNd derivative_fd (3, ndirs);

    for (unsigned i = 1; i < model.mBodies.size(); i++) {
        unsigned int body_id = model.mBodyNameMap[model.GetBodyName(i)];
        // call nominal version
        G = CalcBodyToBaseCoordinates (
            model, q, body_id, point_position, update_kinematics
        );

        // call FD version
        G_fd = FD::CalcBodyToBaseCoordinates (
            model, ad_model,
            q, q_dirs,
            body_id,
            point_position,
            derivative_fd
        );

        // call AD version
        G_ad = AD::CalcBodyToBaseCoordinates (
            model, ad_model,
            q, q_dirs,
            body_id,
            point_position,
            derivative_ad,
            update_kinematics
        );

        // compare nominal results
        CHECK_ARRAY_CLOSE(G_ad.data(), G.data(), G.size(), TEST_PREC);
        CHECK_ARRAY_CLOSE(G_fd.data(), G.data(), G.size(), TEST_PREC);
        CHECK_ARRAY_CLOSE(G_ad.data(), G_fd.data(), G.size(), TEST_PREC);
        CHECK_ARRAY_CLOSE(derivative_ad.data(), derivative_fd.data(),
                          3 * ndirs, TEST_PREC);
    }
}

TEST_FIXTURE ( CartPendulum, CartPendulumCalcBodyToBaseCoordinates) {
    CalcBodyToBaseCoordinatesTemplate(*this);
}

TEST_FIXTURE ( Arm2DofX, Arm2DofXCalcBodyToBaseCoordinates) {
    CalcBodyToBaseCoordinatesTemplate(*this);
}

TEST_FIXTURE ( Arm2DofZ, Arm2DofZCalcBodyToBaseCoordinates) {
    CalcBodyToBaseCoordinatesTemplate(*this);
}

TEST_FIXTURE ( Arm3DofXZYp, Arm3DofXZYpCalcBodyToBaseCoordinates) {
    CalcBodyToBaseCoordinatesTemplate(*this);
}

TEST_FIXTURE ( Arm3DofXZZp, Arm3DofXZZpCalcBodyToBaseCoordinates) {
    CalcBodyToBaseCoordinatesTemplate(*this);
}

// -----------------------------------------------------------------------------

template <typename T>
void CalcBodyWorldOrientationTemplate(T & obj) {
    Model & model = obj.model;
    ADModel & ad_model = obj.ad_model;
    int const nq = model.dof_count;

    VectorNd q = VectorNd::Zero(nq);
    MatrixNd q_dirs = MatrixNd::Identity(nq, nq);

    int ndirs = nq;
    vector<Matrix3d> ad_derivative(ndirs);
    vector<Matrix3d> fd_derivative(ndirs);

    int nTrials = 0;
    do {
        for (unsigned i = 1; i < model.mBodies.size(); i++) {
            unsigned id = model.mBodyNameMap[model.GetBodyName(i)];
            Matrix3d ad_E = AD::CalcBodyWorldOrientation(model, ad_model, q, q_dirs, id, ad_derivative);
            Matrix3d fd_E = FD::CalcBodyWorldOrientation(model, q, q_dirs, id, fd_derivative);
            CHECK_ARRAY_CLOSE(ad_E.data(), fd_E.data(), 9, 1e-6);
            for (int idir = 0; idir < ndirs; idir++) {
                CHECK_ARRAY_CLOSE(fd_derivative[idir].data(), ad_derivative[idir].data(), 9, 1e-6);
            }
        }
        q = VectorNd::Random(nq);
    } while(nTrials++ < 10);
}

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

// -----------------------------------------------------------------------------

template <typename T>
void CalcPointAccelerationTemplate(T & obj) {
  Model   & model    = obj.model;
  ADModel & ad_model = obj.ad_model;
  int const nq       = model.dof_count;

  VectorNd q    = VectorNd::Zero(nq);
  VectorNd qd   = VectorNd::Zero(nq);
  VectorNd qdd  = VectorNd::Zero(nq);
  MatrixNd q_dirs   = MatrixNd::Zero(nq, 3 * nq);
  MatrixNd qd_dirs  = MatrixNd::Zero(nq, 3 * nq);
  MatrixNd qdd_dirs = MatrixNd::Zero(nq, 3 * nq);

  q_dirs.block(0, 0, nq, nq)        = MatrixNd::Identity(nq, nq);
  qd_dirs.block(0, nq, nq, nq)      = MatrixNd::Identity(nq, nq);
  qdd_dirs.block(0, 2 * nq, nq, nq) = MatrixNd::Identity(nq, nq);

  int ndirs = 3 * nq;
  MatrixNd ad_derivative(3, ndirs);
  MatrixNd fd_derivative(3, ndirs);

  Vector3d reference_point = Vector3dZero;
  int nTrials = 0;
  do {
    for (unsigned i = 1; i < model.mBodies.size(); i++) {
      unsigned id = model.mBodyNameMap[model.GetBodyName(i)];
      Vector3d ad_a = AD::CalcPointAcceleration(
            model, ad_model,
            q, q_dirs, qd, qd_dirs, qdd, qdd_dirs,
            id, reference_point,
            ad_derivative);
      Vector3d fd_a = FD::CalcPointAcceleration(
            model,
            q, q_dirs, qd, qd_dirs, qdd, qdd_dirs,
            id, reference_point, fd_derivative);
      CHECK_ARRAY_CLOSE(ad_a.data(), fd_a.data(), 3, 1e-12);
      CHECK_ARRAY_CLOSE(fd_derivative.data(), ad_derivative.data(),
                        3 * ndirs, 1e-7);
    }
    q   = VectorNd::Random(nq);
    qd  = VectorNd::Random(nq);
    qdd = VectorNd::Random(nq);
  } while(nTrials++ < 10);
}

TEST_FIXTURE ( CartPendulum, CartPendulumCalcPointAcceleration) {
  CalcPointAccelerationTemplate(*this);
}

TEST_FIXTURE ( Arm2DofX, Arm2DofXCalcPointAcceleration) {
  CalcPointAccelerationTemplate(*this);
}

TEST_FIXTURE ( Arm2DofZ, Arm2DofZCalcPointAcceleration) {
  CalcPointAccelerationTemplate(*this);
}

TEST_FIXTURE ( Arm3DofXZYp, Arm3DofXZYpCalcPointAcceleration) {
  CalcPointAccelerationTemplate(*this);
}

TEST_FIXTURE ( Arm3DofXZZp, Arm3DofXZZpCalcPointAcceleration) {
  CalcPointAccelerationTemplate(*this);
}

// -----------------------------------------------------------------------------

template <typename T>
void CalcPointJacobianTemplate(T & obj) {
  // rename for convenience
  Model& model = obj.model;
  ADModel& ad_model = obj.ad_model;

  // set up input quantities
  int const nq = model.dof_count;
  int const ndirs = nq;

  bool update_kinematics = true;
  Vector3d point_position = Vector3d::Zero();

  VectorNd q = VectorNd::Zero(nq);
  MatrixNd q_dirs = MatrixNd::Identity(nq, nq);

  // set up no output quantities
  MatrixNd G = MatrixNd::Zero (3, nq);
  MatrixNd G_ad = MatrixNd::Zero (3, nq);
  MatrixNd G_fd = MatrixNd::Zero (3, nq);

  // set up derivative output quantities
  vector<MatrixNd> derivative_ad (ndirs, G_ad);
  vector<MatrixNd> derivative_fd (ndirs, G_fd);

  for (unsigned i = 1; i < model.mBodies.size(); i++) {
    unsigned int body_id = model.mBodyNameMap[model.GetBodyName(i)];
    // call nominal version
    CalcPointJacobian (
          model, q, body_id, point_position, G, update_kinematics
          );

    // call AD version
    AD::CalcPointJacobian (
          model, ad_model,
          q, q_dirs,
          body_id,
          point_position,
          G_ad, derivative_ad,
          update_kinematics
          );

    // call FD version
    FD::CalcPointJacobian (
          model, ad_model,
          q, q_dirs,
          body_id,
          point_position,
          G_fd, derivative_fd,
          update_kinematics
          );

    // compare nominal results
    CHECK_ARRAY_CLOSE(G_ad.data(), G.data(), G.size(), TEST_PREC);
    CHECK_ARRAY_CLOSE(G_fd.data(), G.data(), G.size(), TEST_PREC);
    CHECK_ARRAY_CLOSE(G_ad.data(), G_fd.data(), G.size(), TEST_PREC);

    for (int idir = 0; idir < ndirs; idir++) {
      CHECK_ARRAY_CLOSE(
            derivative_ad[idir].data(),
            derivative_fd[idir].data(),
            derivative_ad[idir].size(),
            TEST_PREC
            );
    }
  }
}

TEST_FIXTURE ( CartPendulum, CartPendulumCalcPointJacobian) {
  CalcPointJacobianTemplate(*this);
}

TEST_FIXTURE ( Arm2DofX, Arm2DofXCalcPointJacobian) {
  CalcPointJacobianTemplate(*this);
}

TEST_FIXTURE ( Arm2DofZ, Arm2DofZCalcPointJacobian) {
  CalcPointJacobianTemplate(*this);
}

TEST_FIXTURE ( Arm3DofXZYp, Arm3DofXZYpCalcPointJacobian) {
  CalcPointJacobianTemplate(*this);
}

TEST_FIXTURE ( Arm3DofXZZp, Arm3DofXZZpCalcPointJacobian) {
  CalcPointJacobianTemplate(*this);
}

// -----------------------------------------------------------------------------
