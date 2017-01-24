#include <UnitTest++.h>

#include <iomanip>
#include <iostream>
// #include <random>

#include "Fixtures.h"
#include "ModelAD.h"
#include "KinematicsAD.h"
#include "KinematicsFD.h"

using namespace std;
using namespace RigidBodyDynamics;
using namespace RigidBodyDynamics::Math;
using namespace RigidBodyDynamics::Utils;

// -----------------------------------------------------------------------------

template <typename T>
void CalcBodyToBaseCoordinatesSingleFuncTemplate(
    T & obj,
    unsigned trial_count,
    double array_close_prec,
    unsigned id_ee) {
  Model model = obj.model;
  VectorNd q = obj.q;
  Vector3d point_body_coordinates;

  for (unsigned trial = 0; trial < trial_count; trial++) {
    q.setRandom();
    point_body_coordinates.setRandom();

    Vector3d point_single_func =
        CalcBodyToBaseCoordinatesSingleFunc (model, q, id_ee, point_body_coordinates);
    Vector3d point_default =
        CalcBodyToBaseCoordinates (model, q, id_ee, point_body_coordinates);

    CHECK_ARRAY_CLOSE (point_default.data(), point_single_func.data(), 3, array_close_prec);
  }
}

TEST_FIXTURE ( CartPendulum, CartPendulumCalcBodyToBaseCoordinatesSingleFunc) {
  CalcBodyToBaseCoordinatesSingleFuncTemplate(*this, 10, 1e-6, id_pendulum);
}

TEST_FIXTURE ( Arm2DofX, Arm2DofXCalcBodyToBaseCoordinatesSingleFunc) {
  CalcBodyToBaseCoordinatesSingleFuncTemplate(*this, 10, 1e-6, id_proximal);
}

TEST_FIXTURE ( Arm2DofZ, Arm2DofZCalcBodyToBaseCoordinatesSingleFunc) {
  CalcBodyToBaseCoordinatesSingleFuncTemplate(*this, 10, 1e-6, id_proximal);
}

TEST_FIXTURE ( Arm3DofXZYp, Arm3DofXZYpCalcBodyToBaseCoordinatesSingleFunc) {
  CalcBodyToBaseCoordinatesSingleFuncTemplate(*this, 10, 1e-6, id_proximal);
}

TEST_FIXTURE ( Arm3DofXZZp, Arm3DofXZZpCalcBodyToBaseCoordinatesSingleFunc) {
  CalcBodyToBaseCoordinatesSingleFuncTemplate(*this, 10, 1e-6, id_proximal);
}

// -----------------------------------------------------------------------------

template<typename T>
void JacobianADSimpleTemplate(
    T & obj,
    unsigned trial_count,
    double array_close_prec,
    unsigned int id_ee) {
  Model & model = obj.model;
  ADModel & ad_model = obj.ad_model;

  MatrixNd ad_base_pt_dirs = MatrixNd::Zero(3, model.qdot_size);
  MatrixNd d_base_pt_dirs  = MatrixNd::Zero(3, model.qdot_size);
  MatrixNd fd_base_pt_dirs = MatrixNd::Zero(3, model.qdot_size);

  VectorNd q = VectorNd::Random(model.dof_count);
  Vector3d body_point = Vector3d (1.0, 2.0, 3.0);
  unsigned trial = 0;
  do {
    CalcPointJacobian (model, q, id_ee, body_point, d_base_pt_dirs);

    MatrixNd q_dirs = MatrixNd::Identity (model.qdot_size, model.qdot_size);
    Vector3d base_pt = CalcBodyToBaseCoordinates (
          model, q, id_ee, body_point);
    Vector3d ad_base_pt = AD::CalcBodyToBaseCoordinatesSingleFunc (
          model, ad_model, q, q_dirs, id_ee, body_point, ad_base_pt_dirs);
    Vector3d fd_base_pt = FD::CalcBodyToBaseCoordinatesSingleFunc (
          model, q, q_dirs, id_ee, body_point, fd_base_pt_dirs);

    CHECK_ARRAY_CLOSE (d_base_pt_dirs.data(), ad_base_pt_dirs.data(),
                       3 * model.qdot_size, array_close_prec);
    CHECK_ARRAY_CLOSE (base_pt.data(), ad_base_pt.data(),
                       3, array_close_prec);
    CHECK_ARRAY_CLOSE (base_pt.data(), fd_base_pt.data(),
                       3, array_close_prec);
    CHECK_ARRAY_CLOSE (d_base_pt_dirs.data(), fd_base_pt_dirs.data(),
                       3 * model.qdot_size, array_close_prec);

    q.setRandom();
    body_point.setRandom();
  } while (trial++ < trial_count);
}

TEST_FIXTURE ( CartPendulum, CartPendulumJacobianADSimple ) {
  JacobianADSimpleTemplate(*this, 10, 1e-6, id_pendulum);
}

TEST_FIXTURE ( Arm2DofX, Arm2DofXJacobianADSimple) {
  JacobianADSimpleTemplate(*this, 10, 1e-6, id_proximal);
}

TEST_FIXTURE ( Arm2DofZ, Arm2DofZJacobianADSimple) {
  JacobianADSimpleTemplate(*this, 10, 1e-6, id_proximal);
}

TEST_FIXTURE ( Arm3DofXZYp, Arm3DofXZYpJacobianADSimple) {
  JacobianADSimpleTemplate(*this, 10, 1e-6, id_slider);
}

TEST_FIXTURE ( Arm3DofXZZp, Arm3DofXZZpJacobianADSimple) {
  JacobianADSimpleTemplate(*this, 10, 1e-6, id_slider);
}

// -----------------------------------------------------------------------------

