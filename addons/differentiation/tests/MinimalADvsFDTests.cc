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

  MatrixNd jacobian_ad = MatrixNd::Zero(3, model.qdot_size);
  MatrixNd jacobian_ref = MatrixNd::Zero(3, model.qdot_size);
  MatrixNd jacobian_fd = MatrixNd::Zero(3, model.qdot_size);

  VectorNd q = VectorNd::Random(model.dof_count);
  Vector3d body_point = Vector3d (1.0, 2.0, 3.0);
  unsigned trial = 0;
  do {
    CalcPointJacobian (model, q, id_ee, body_point, jacobian_ref);

    MatrixNd q_dirs = MatrixNd::Identity (model.qdot_size, model.qdot_size);
    Vector3d base_point_standard = CalcBodyToBaseCoordinates (model, q, id_ee, body_point);
    Vector3d base_point_ad = AD::CalcBodyToBaseCoordinatesSingleFunc (model, ad_model, q, q_dirs, id_ee, body_point, jacobian_ad);
    Vector3d base_point_fd = FD::CalcBodyToBaseCoordinatesSingleFunc (model, q, q_dirs, id_ee, body_point, jacobian_fd);

    CHECK_ARRAY_CLOSE (jacobian_ref.data(), jacobian_ad.data(), 3 * model.qdot_size, array_close_prec);
    CHECK_ARRAY_CLOSE (base_point_standard.data(), base_point_ad.data(), 3, array_close_prec);
    CHECK_ARRAY_CLOSE (base_point_standard.data(), base_point_fd.data(), 3, array_close_prec);
    CHECK_ARRAY_CLOSE (jacobian_ref.data(), jacobian_fd.data(), 3 * model.qdot_size, array_close_prec);

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

//TEST_FIXTURE ( Arm2DofZ, Arm2DofZJacobianADSimple) {
//    JacobianADSimpleTemplate(*this, id_proximal);
//}

//TEST_FIXTURE ( Arm3DofXZYp, Arm3DofXZYpJacobianADSimple) {
//    JacobianADSimpleTemplate(*this, id_slider);
//}

//TEST_FIXTURE ( Arm3DofXZZp, Arm3DofXZZpJacobianADSimple) {
//    JacobianADSimpleTemplate(*this, id_slider);
//}

// -----------------------------------------------------------------------------


// -----------------------------------------------------------------------------
// begin deprecated
// -----------------------------------------------------------------------------

// TEST_FIXTURE (CartPendulum, jcalcNominalSolutionTest) {
// 	// set nominal values
// 	q.setZero();
// 	q[0] = 0.3;
// 	q[1] = -0.2;

// 	qdot.setZero();
// 	qdot[0] = 0.3;
// 	qdot[1] = -0.2;

// 	// set directions
// 	unsigned int ndirs = model.q_size + model.qdot_size;
// 	MatrixNd x = MatrixNd::Identity(ndirs, ndirs);
// 	MatrixNd q_dirs = x.block(0, 0, model.q_size, ndirs);
// 	MatrixNd qdot_dirs = x.block(model.q_size, 0, model.qdot_size, ndirs);

// 	// set derivative outputs
// 	std::vector<SpatialMatrix> ad_X_Ji (ndirs, SpatialMatrix::Zero());
// 	std::vector<std::vector<SpatialMatrix> > ad_X_J (model.mBodies.size(), ad_X_Ji);

// 	std::vector<SpatialVector> ad_V (ndirs, SpatialVector::Zero (6));
// 	std::vector<std::vector<SpatialVector> > ad_S (model.mBodies.size(), ad_V);
// 	std::vector<std::vector<SpatialVector> > ad_v_J (model.mBodies.size(), ad_V);
// 	std::vector<std::vector<SpatialVector> > ad_c_J (model.mBodies.size(), ad_V);

// 	for (unsigned int joint_id = 1; joint_id < model.mBodies.size(); joint_id++) {
// 		// evaluate nominal solution
// 		jcalc (model, joint_id, q, qdot);
// 		MatrixNd ref_X_Ji = model.X_J[joint_id].toMatrix();
// 		SpatialVector ref_S_i = model.S[joint_id];
// 		SpatialVector ref_v_Ji = model.v_J[joint_id];
// 		SpatialVector ref_c_Ji = model.c_J[joint_id];

// 		// evaluate AD nominal solution
// 		ad_jcalc (model, ad_model, joint_id, q, q_dirs, qdot, qdot_dirs);
// 		MatrixNd test_X_Ji = model.X_J[joint_id].toMatrix();
// 		SpatialVector test_S_i = model.S[joint_id];
// 		SpatialVector test_v_Ji = model.v_J[joint_id];
// 		SpatialVector test_c_Ji = model.c_J[joint_id];

// 		CHECK_ARRAY_CLOSE (ref_X_Ji.data(), test_X_Ji.data(), 36, TEST_PREC);
// 		CHECK_ARRAY_CLOSE (ref_S_i.data(),  test_S_i.data(),   6, TEST_PREC);
// 		CHECK_ARRAY_CLOSE (ref_v_Ji.data(), test_v_Ji.data(),  6, TEST_PREC);
// 		CHECK_ARRAY_CLOSE (ref_c_Ji.data(), test_c_Ji.data(),  6, TEST_PREC);
// 		/* DEBUG OUTPUT
// 		cout << "===== joint_id: " << joint_id << " =====" << endl;
// 		cout << "ref_X_Ji: " << endl << ref_X_Ji << endl;
// 		cout << "test_X_Ji: " << endl << test_X_Ji << endl;
// 		cout << "error_X_Ji: " << endl << ref_X_Ji - test_X_Ji << endl;
// 		cout << endl;

// 		cout << "ref_S_i: " << endl << ref_S_i << endl;
// 		cout << "test_S_i: " << endl << test_S_i << endl;
// 		cout << "error_S_i: " << endl << ref_S_i - test_S_i << endl;
// 		cout << endl;

// 		cout << "ref_v_Ji: " << endl << ref_v_Ji << endl;
// 		cout << "test_v_Ji: " << endl << test_v_Ji << endl;
// 		cout << "error_v_Ji: " << endl << ref_v_Ji - test_v_Ji << endl;
// 		cout << endl;

// 		cout << "ref_c_Ji: " << endl << ref_c_Ji << endl;
// 		cout << "test_c_Ji: " << endl << test_c_Ji << endl;
// 		cout << "error_c_Ji: " << endl << ref_c_Ji - test_c_Ji << endl;
// 		cout << endl;
// 		cout << endl;
// 		*/
// 	}
// }

// TEST_FIXTURE (CartPendulum, jcalcFDvsADTest) {
// 	double TEST_PREC = 1.0e-08;

// 	// set nominal values
// 	q.setZero();
// 	q[0] = 0.3;
// 	q[1] = -0.2;

// 	qdot.setZero();
// 	qdot[0] = 0.3;
// 	qdot[1] = -0.2;

// 	// set directions
// 	unsigned int ndirs = model.q_size + model.qdot_size;
// 	MatrixNd x = MatrixNd::Identity(ndirs, ndirs);
// 	MatrixNd q_dirs = x.block(0, 0, model.q_size, ndirs);
// 	MatrixNd qdot_dirs = x.block(model.q_size, 0, model.qdot_size, ndirs);

// 	// set derivative outputs
// 	std::vector<SpatialMatrix> ad_X_Ji (ndirs, SpatialMatrix::Zero());
// 	std::vector<std::vector<SpatialMatrix> > fd_X_J (model.mBodies.size(), ad_X_Ji);
// 	std::vector<std::vector<SpatialMatrix> > ad_X_J (model.mBodies.size(), ad_X_Ji);

// 	std::vector<SpatialVector> ad_V (ndirs, SpatialVector::Zero (6));
// 	std::vector<std::vector<SpatialVector> > fd_S (model.mBodies.size(), ad_V);
// 	std::vector<std::vector<SpatialVector> > fd_v_J (model.mBodies.size(), ad_V);
// 	std::vector<std::vector<SpatialVector> > fd_c_J (model.mBodies.size(), ad_V);

// 	for (unsigned int joint_id = 1; joint_id < model.mBodies.size(); joint_id++) {
// 		// evaluate nominal solution
// 		fd_jcalc (model, joint_id, q, q_dirs, qdot, qdot_dirs,
// 			fd_X_J[joint_id], fd_S[joint_id], fd_v_J[joint_id], fd_c_J[joint_id]);

// 		// evaluate AD nominal solution
// 		ad_jcalc (model, ad_model, joint_id, q, q_dirs, qdot, qdot_dirs);

// 		for (int idir = 0; idir < ndirs; ++idir) {
// 			CHECK_ARRAY_CLOSE (fd_X_J[joint_id][idir].data(), ad_model.X_J[joint_id][idir].data(), 36, TEST_PREC);
// 			CHECK_ARRAY_CLOSE (fd_S[joint_id][idir].data(),   ad_model.S[joint_id][idir].data(),    6, TEST_PREC);
// 			CHECK_ARRAY_CLOSE (fd_v_J[joint_id][idir].data(), ad_model.v_J[joint_id][idir].data(),  6, TEST_PREC);
// 			CHECK_ARRAY_CLOSE (fd_c_J[joint_id][idir].data(), ad_model.c_J[joint_id][idir].data(),  6, TEST_PREC);
// 			/* DEBUG OUTPUT
// 			cout << "===== joint_id: " << joint_id << ", idir: " << idir << " =====" << endl;
// 			cout << "fd_X_J[" << joint_id << "][" << idir << "]: " << endl << fd_X_J[joint_id][idir] << endl;
// 			cout << "ad_X_J[" << joint_id << "][" << idir << "]: " << endl << ad_X_J[joint_id][idir] << endl;
// 			cout << "error_X_J: " << endl << fd_X_J[joint_id][idir] - ad_X_J[joint_id][idir] << endl;

// 			cout << "fd_S[" << joint_id << "][" << idir << "]: " << endl << fd_S[joint_id][idir] << endl;
// 			cout << "ad_S[" << joint_id << "][" << idir << "]: " << endl << ad_S[joint_id][idir] << endl;
// 			cout << "error_S: " << endl << fd_S[joint_id][idir] - ad_S[joint_id][idir] << endl;

// 			cout << "fd_v_J[" << joint_id << "][" << idir << "]: " << endl << fd_v_J[joint_id][idir] << endl;
// 			cout << "ad_v_J[" << joint_id << "][" << idir << "]: " << endl << ad_v_J[joint_id][idir] << endl;
// 			cout << "error_v_J: " << endl << fd_v_J[joint_id][idir] - ad_v_J[joint_id][idir] << endl;

// 			cout << "fd_c_J[" << joint_id << "][" << idir << "]: " << endl << fd_c_J[joint_id][idir] << endl;
// 			cout << "ad_c_J[" << joint_id << "][" << idir << "]: " << endl << ad_c_J[joint_id][idir] << endl;
// 			cout << "error_c_J: " << endl << fd_c_J[joint_id][idir] - ad_c_J[joint_id][idir] << endl;
// 			*/
// 		}
// 	}
// }

// TEST_FIXTURE (CartPendulum, ad_UpdateKinematicsNominalTest) {
// 	// set nominal values
// 	q.setZero();
// 	q[0] = 0.3;
// 	q[1] = -0.2;

// 	qdot.setZero();
// 	qdot[0] = 0.3;
// 	qdot[1] = -0.2;

// 	qddot.setZero();
// 	qddot[0] = 0.3;
// 	qddot[1] = -0.2;

// 	// set directions
// 	unsigned int ndirs = 3*model.q_size;
// 	MatrixNd x = MatrixNd::Identity(ndirs, ndirs);
// 	MatrixNd q_dirs = x.block(0, 0, model.q_size, ndirs);
// 	MatrixNd qdot_dirs = x.block(model.q_size, 0, model.q_size, ndirs);
// 	MatrixNd qddot_dirs = x.block(2*model.q_size, 0, model.q_size, ndirs);

// 	// evaluate nominal solution
// 	UpdateKinematics (model, q, qdot, qddot);
// 	Model res = model;

// 	// evaluate AD nominal solution
// 	ad_UpdateKinematics (
// 		model, ad_model, q, q_dirs, qdot, qdot_dirs, qddot, qddot_dirs
// 	);
// 	Model ad_res = model;

// 	for (unsigned int joint_id = 1; joint_id < model.mBodies.size(); joint_id++) {
// 		CHECK_ARRAY_CLOSE (res.X_lambda[joint_id].toMatrix().data(), ad_res.X_lambda[joint_id].toMatrix().data(), 36, TEST_PREC);
// 		CHECK_ARRAY_CLOSE (res.X_base[joint_id].toMatrix().data(), ad_res.X_base[joint_id].toMatrix().data(), 36, TEST_PREC);

// 		CHECK_ARRAY_CLOSE (res.a[joint_id].data(), ad_res.a[joint_id].data(), 6, TEST_PREC);
// 		CHECK_ARRAY_CLOSE (res.v[joint_id].data(), ad_res.v[joint_id].data(), 6, TEST_PREC);
// 		CHECK_ARRAY_CLOSE (res.c[joint_id].data(), ad_res.c[joint_id].data(), 6, TEST_PREC);
// 	}
// }

// TEST_FIXTURE (CartPendulum, ad_UpdateKinematicsADvsFDTest) {
// 	// set nominal values
// 	q.setZero();
// 	q[0] = 0.3;
// 	q[1] = -0.2;

// 	qdot.setZero();
// 	qdot[0] = 0.3;
// 	qdot[1] = -0.2;

// 	qddot.setZero();
// 	qddot[0] = 0.3;
// 	qddot[1] = -0.2;

// 	// set directions
// 	unsigned int ndirs = 3*model.q_size;
// 	MatrixNd x = MatrixNd::Identity(ndirs, ndirs);
// 	MatrixNd q_dirs = x.block(0, 0, model.q_size, ndirs);
// 	MatrixNd qdot_dirs = x.block(model.q_size, 0, model.q_size, ndirs);
// 	MatrixNd qddot_dirs = x.block(2*model.q_size, 0, model.q_size, ndirs);

// 	// set derivative outputs
// 	std::vector<SpatialMatrix> ad_X (ndirs, SpatialMatrix::Zero());
// 	std::vector<std::vector<SpatialMatrix> > fd_X_lambda (model.mBodies.size(), ad_X);
// 	std::vector<std::vector<SpatialMatrix> > fd_X_base (model.mBodies.size(), ad_X);

// 	std::vector<SpatialVector> ad_V (ndirs, SpatialVector::Zero (6));
// 	std::vector<std::vector<SpatialVector> > fd_a (model.mBodies.size(), ad_V);
// 	std::vector<std::vector<SpatialVector> > fd_v (model.mBodies.size(), ad_V);
// 	std::vector<std::vector<SpatialVector> > fd_c (model.mBodies.size(), ad_V);

// 	// evaluate nominal solution
// 	fd_UpdateKinematics (
// 		model, q, q_dirs, qdot, qdot_dirs, qddot, qddot_dirs,
// 		fd_X_lambda, fd_X_base, fd_a, fd_v, fd_c
// 	);

// 	ad_UpdateKinematics (
// 		model, ad_model, q, q_dirs, qdot, qdot_dirs, qddot, qddot_dirs
// 	);

// 	for (unsigned int joint_id = 1; joint_id < model.mBodies.size(); joint_id++) {
// 		for (int idirs = 0; idirs < ndirs; ++idirs) {
// 			CHECK_ARRAY_CLOSE (fd_X_lambda[joint_id][idirs].data(), ad_model.X_lambda[joint_id][idirs].data(), 36, TEST_PREC);
// 			CHECK_ARRAY_CLOSE (fd_X_base[joint_id][idirs].data(), ad_model.X_base[joint_id][idirs].data(), 36, TEST_PREC);

// 			CHECK_ARRAY_CLOSE (fd_a[joint_id][idirs].data(), ad_model.a[joint_id][idirs].data(), 6, TEST_PREC);
// 			CHECK_ARRAY_CLOSE (fd_v[joint_id][idirs].data(), ad_model.v[joint_id][idirs].data(), 6, TEST_PREC);
// 			CHECK_ARRAY_CLOSE (fd_c[joint_id][idirs].data(), ad_model.c[joint_id][idirs].data(), 6, TEST_PREC);
// 			/* DEBUG OUTPUT
// 			// cout << "fd_X_lambda[" << joint_id << "][" << idirs << "]: " << endl << fd_X_lambda[joint_id][idirs] << endl;
// 			//cout << "ad_X_lambda[" << joint_id << "][" << idirs << "]: " << endl << ad_model.X_lambda[joint_id][idirs] << endl;
// 			//cout << "err" << endl << fd_X_lambda[joint_id][idirs] - ad_model.X_lambda[joint_id][idirs] << endl;
// 			// cout << endl;

// 			// cout << "fd_X_base[" << joint_id << "][" << idirs << "]: " << endl << fd_X_base[joint_id][idirs] << endl;
// 			//cout << "ad_X_base[" << joint_id << "][" << idirs << "]: " << endl << ad_model.X_base[joint_id][idirs] << endl;
// 			//cout << "err" << endl << fd_X_base[joint_id][idirs] - ad_model.X_base[joint_id][idirs] << endl;
// 			// cout << endl;

// 			// cout << "fd_a[" << joint_id << "][" << idirs << "]: " << endl << fd_a[joint_id][idirs] << endl;
// 			//cout << "ad_a[" << joint_id << "][" << idirs << "]: " << endl << ad_model.a[joint_id][idirs] << endl;
// 			//cout << "err" << endl << fd_a[joint_id][idirs] - ad_model.a[joint_id][idirs] << endl;
// 			// cout << endl;

// 			// cout << "fd_v[" << joint_id << "][" << idirs << "]: " << endl << fd_v[joint_id][idirs] << endl;
// 			//cout << "ad_v[" << joint_id << "][" << idirs << "]: " << endl << ad_model.v[joint_id][idirs] << endl;
// 			//cout << "err" << endl << fd_v[joint_id][idirs] - ad_model.v[joint_id][idirs] << endl;
// 			// cout << endl;

// 			cout << "fd_c[" << joint_id << "][" << idirs << "]: " << endl << fd_c[joint_id][idirs] << endl;
// 			//cout << "ad_c[" << joint_id << "][" << idirs << "]: " << endl << ad_model.c[joint_id][idirs] << endl;
// 			//cout << "err" << endl << fd_c[joint_id][idirs] - ad_model.c[joint_id][idirs] << endl;
// 			cout << endl;
// 			*/
// 		}
// 	}
// }

// TEST_FIXTURE (CartPendulum, CalcPointAccelerationNominalTest) {
// 	// set nominal values
// 	q.setZero();
// 	q[0] = 0.3;
// 	q[1] = -0.2;

// 	qdot.setZero();
// 	qdot[0] = 0.3;
// 	qdot[1] = -0.2;

// 	qddot.setZero();
// 	qddot[0] = 0.3;
// 	qddot[1] = -0.2;

// 	Vector3d point_body_coordinates (0.1, 3.2, 4.2);

// 	// set directions
// 	unsigned int ndirs = 3*model.q_size;
// 	MatrixNd x = MatrixNd::Identity(ndirs, ndirs);
// 	MatrixNd q_dirs = x.block(0, 0, model.q_size, ndirs);
// 	MatrixNd qdot_dirs = x.block(model.q_size, 0, model.q_size, ndirs);
// 	MatrixNd qddot_dirs = x.block(2*model.q_size, 0, model.q_size, ndirs);

// 	// set derivative output
// 	MatrixNd ad_jacobian = MatrixNd::Zero(3, ndirs);

// 	// evaluate nominal solution
// 	Vector3d ref_acc = CalcPointAcceleration (
// 		model, q, qdot, qddot, id_pendulum, point_body_coordinates
// 	);

// 	Vector3d ad_acc = ad_CalcPointAcceleration (
// 		model, ad_model,
// 		q, q_dirs, qdot, qdot_dirs, qddot, qddot_dirs,
// 		id_pendulum, point_body_coordinates,
// 		ad_jacobian
// 	);

// 	CHECK_ARRAY_CLOSE (ref_acc.data(), ad_acc.data(), 3, TEST_PREC);
// }


// TEST_FIXTURE (CartPendulum, CalcPointAccelerationFDvsADTest) {
// 	// set nominal values
// 	q.setZero();
// 	q[0] = 0.3;
// 	q[1] = -0.2;

// 	qdot.setZero();
// 	qdot[0] = 0.3;
// 	qdot[1] = -0.2;

// 	qddot.setZero();
// 	qddot[0] = 0.3;
// 	qddot[1] = -0.2;

// 	Vector3d point_body_coordinates (0.1, 3.2, 4.2);

// 	// set directions
// 	unsigned int ndirs = model.q_size + model.qdot_size + model.q_size;
// 	MatrixNd x = MatrixNd::Identity(ndirs, ndirs);
// 	MatrixNd q_dirs = x.block(0, 0, model.q_size, ndirs);
// 	MatrixNd qdot_dirs = x.block(model.q_size, 0, model.qdot_size, ndirs);
// 	MatrixNd qddot_dirs = x.block(model.q_size + model.qdot_size, 0, model.q_size, ndirs);


// 	// set derivative output
// 	MatrixNd fd_jacobian = MatrixNd::Zero(3, ndirs);
// 	MatrixNd ad_jacobian = MatrixNd::Zero(3, ndirs);

// 	// evaluate nominal solution
// 	Vector3d ref_acc = fd_CalcPointAcceleration (
// 		model,
// 		q, q_dirs, qdot, qdot_dirs, qddot, qddot_dirs,
// 		id_pendulum, point_body_coordinates,
// 		ad_jacobian
// 	);

// 	Vector3d ad_acc = ad_CalcPointAcceleration (
// 		model, ad_model,
// 		q, q_dirs, qdot, qdot_dirs, qddot, qddot_dirs,
// 		id_pendulum, point_body_coordinates,
// 		ad_jacobian
// 	);

// 	CHECK_ARRAY_CLOSE (ref_acc.data(), ad_acc.data(), 3, TEST_PREC);
// }

// TEST_FIXTURE (CartPendulum, ForwardDynamicsCholesky) {
//   // set nominal values

//   q = VectorNd::Random(model.q_size);
//   qdot = VectorNd::Random(model.q_size);
//   qddot = VectorNd::Random(model.q_size);
//   tau = VectorNd::Random(model.q_size);


//   std::vector<SpatialVector> f_ext (model.mBodies.size(),SpatialVector::Zero(model.q_size));
//   for (int i = 0; i < model.mBodies.size(); ++i) {
// 	f_ext[i]=SpatialVector::Random(model.q_size);
//   }


//   VectorNd qddot_ref (VectorNd::Zero (model.q_size));

//   ForwardDynamicsCholesky(model,q,qdot,tau,qddot,&f_ext);

//   ForwardDynamics(model,q,qdot,tau,qddot_ref,&f_ext);

//   CHECK_ARRAY_CLOSE (qddot, qddot_ref, model.q_size, TEST_PREC);
//   /*
// 	cout << "qddot_ref: " << endl << qddot_ref << endl;
// 	cout << "qddot_test: " << endl << qddot << endl;
// 	cout << "error qddot: " << endl << qddot_ref - qddot << endl;
//   */
// }


// TEST_FIXTURE(CartPendulum, ForwardDynamicsCholeskyADTest){

//   VectorNd q = VectorNd::Random(model.q_size);
//   VectorNd qdot = VectorNd::Random(model.q_size);
//   VectorNd tau = VectorNd::Random(model.q_size);
//   VectorNd qddot = VectorNd::Zero(model.q_size);
//   VectorNd qddot_ref = VectorNd::Zero(model.q_size);

//   unsigned int ndirs = 3 * model.q_size;
//   MatrixNd x = MatrixNd::Identity(ndirs, ndirs);
//   MatrixNd q_dirs = x.block(0, 0, model.q_size, ndirs);
//   MatrixNd qdot_dirs = x.block(model.q_size, 0, model.q_size, ndirs);
//   MatrixNd tau_dirs = x.block(2*model.q_size, 0, model.q_size, ndirs);


//   std::vector<SpatialVector> f_ext (model.mBodies.size(),SpatialVector::Zero(model.q_size));
//   for (int i = 0; i < model.mBodies.size(); ++i) {
// 	f_ext[i]=SpatialVector::Zero(model.q_size);
//   }

//   MatrixNd ad_qddot  = MatrixNd::Zero(model.q_size, ndirs);
//   MatrixNd fd_qddot  = MatrixNd::Zero(model.q_size, ndirs);

//   fd_ForwardDynamics(model,
// 			 q,
// 			 q_dirs,
// 			 qdot,
// 			 qdot_dirs,
// 			 tau,
// 			 tau_dirs,
// 			 qddot_ref,
// 			 fd_qddot,
// 			 &f_ext);

//   ad_ForwardDynamicsCholesky(model,
// 				 ad_model,
// 				 q,
// 				 q_dirs,
// 				 qdot,
// 				 qdot_dirs,
// 				 tau,
// 				 tau_dirs,
// 				 qddot,
// 				 ad_qddot,
// 				 &f_ext);

//   CHECK_ARRAY_CLOSE (qddot_ref.data(), qddot.data(), qddot.size(), 1e-7);
//   CHECK_ARRAY_CLOSE (fd_qddot.data(), ad_qddot.data(), fd_qddot.cols()*fd_qddot.rows(), 1e-7);
// }




// TEST_FIXTURE(CartPendulum, InverseDynamicsADTest_with_external_forces){
// 	for(unsigned int i = 0; i < model.qdot_size; i++){
// 		q[i] = (i+1)*(0.897878435);
// 		qdot[i] = (i+1)*(0.27563682);
// 		qddot[i] = (i+1)*(0.06565644564455);
// 	}

// 	MatrixNd q_dirs 	= MatrixNd::Identity (model.qdot_size, model.qdot_size);
// 	MatrixNd qdot_dirs 	= MatrixNd::Identity (model.qdot_size, model.qdot_size);
// 	MatrixNd qddot_dirs = MatrixNd::Identity (model.qdot_size, model.qdot_size);

// 	std::vector<SpatialVector> f_ext (model.mBodies.size(),SpatialVector::Zero(model.q_size));
// 	for (int i = 0; i < model.mBodies.size(); ++i) {
// 	  f_ext[i]=SpatialVector::Random(model.q_size);
// 	}


// 	VectorNd tau_ref (tau);

// 	MatrixNd ad_tau  = MatrixNd::Zero(model.qdot_size, model.qdot_size);
// 	MatrixNd fd_tau  = MatrixNd::Zero(model.qdot_size, model.qdot_size);

// 	ad_InverseDynamics(model,
// 			ad_model,
// 			q,
// 			q_dirs,
// 			qdot,
// 			qdot_dirs,
// 			qddot,
// 			qddot_dirs,
// 			tau,
// 			ad_tau,
// 			&f_ext);


// 	fd_InverseDynamics(model,
// 			q,
// 			q_dirs,
// 			qdot,
// 			qdot_dirs,
// 			qddot,
// 			qddot_dirs,
// 			tau_ref,
// 			fd_tau,
// 			&f_ext);

// 	CHECK_ARRAY_CLOSE (fd_tau.data(), ad_tau.data(), fd_tau.cols()*fd_tau.rows(), 1e-7);
// 	CHECK_ARRAY_CLOSE (tau_ref.data(), tau.data(), tau_ref.rows(), 1e-7);

// }

// TEST_FIXTURE(CartPendulum, ForwardDynamicsADTest_with_external_forces){
//   srand((unsigned int) time(0));

//   for(unsigned int trial = 0; trial < 10; trial++) {
// 	VectorNd q = VectorNd::Random(model.q_size);
// 	VectorNd qdot = VectorNd::Random(model.q_size);
// 	VectorNd tau = VectorNd::Random(model.q_size);

// 	unsigned int ndirs = 3 * model.q_size;
// 	MatrixNd x = MatrixNd::Identity(ndirs, ndirs);
// 	MatrixNd q_dirs = x.block(0, 0, model.q_size, ndirs);
// 	MatrixNd qdot_dirs = x.block(model.q_size, 0, model.q_size, ndirs);
// 	MatrixNd tau_dirs = x.block(2*model.q_size, 0, model.q_size, ndirs);

// 	std::vector<SpatialVector> f_ext (model.mBodies.size(),SpatialVector::Zero(model.q_size));
// 	for (int i = 0; i < model.mBodies.size(); ++i) {
// 	  f_ext[i]=SpatialVector::Random(model.q_size);
// 	}


// 	VectorNd qddot (VectorNd::Zero(model.q_size));
// 	MatrixNd ad_qddot  = MatrixNd::Zero(model.qdot_size, ndirs);
// 	MatrixNd fd_qddot  = MatrixNd::Zero(model.qdot_size, ndirs);


// 	fd_ForwardDynamics(model, q, q_dirs, qdot, qdot_dirs, tau, tau_dirs, qddot, fd_qddot,&f_ext);

// 	ad_ForwardDynamics(model, ad_model, q, q_dirs, qdot, qdot_dirs, tau, tau_dirs, qddot, ad_qddot, &f_ext);

// 	CHECK_ARRAY_CLOSE (fd_qddot.data(), ad_qddot.data(), fd_qddot.cols()*fd_qddot.rows(), 1e-7);
//   }
//   }



// TEST_FIXTURE(CartPendulum, ForwardDynamicsCholeskyADTest_with_external_forces){

//   VectorNd q = VectorNd::Random(model.q_size);
//   VectorNd qdot = VectorNd::Random(model.q_size);
//   VectorNd tau = VectorNd::Random(model.q_size);
//   VectorNd qddot = VectorNd::Zero(model.q_size);
//   VectorNd qddot_ref = VectorNd::Zero(model.q_size);

//   unsigned int ndirs = 3 * model.q_size;
//   MatrixNd x = MatrixNd::Identity(ndirs, ndirs);
//   MatrixNd q_dirs = x.block(0, 0, model.q_size, ndirs);
//   MatrixNd qdot_dirs = x.block(model.q_size, 0, model.q_size, ndirs);
//   MatrixNd tau_dirs = x.block(2*model.q_size, 0, model.q_size, ndirs);


//   std::vector<SpatialVector> f_ext (model.mBodies.size(),SpatialVector::Zero(model.q_size));
//   for (int i = 0; i < model.mBodies.size(); ++i) {
// 	f_ext[i]=SpatialVector::Random(model.q_size);
//   }

//   MatrixNd ad_qddot  = MatrixNd::Zero(model.q_size, ndirs);
//   MatrixNd fd_qddot  = MatrixNd::Zero(model.q_size, ndirs);

//   fd_ForwardDynamics(model,
// 			 q,
// 			 q_dirs,
// 			 qdot,
// 			 qdot_dirs,
// 			 tau,
// 			 tau_dirs,
// 			 qddot_ref,
// 			 fd_qddot,
// 			 &f_ext);

//   ad_ForwardDynamicsCholesky(model,
// 				 ad_model,
// 				 q,
// 				 q_dirs,
// 				 qdot,
// 				 qdot_dirs,
// 				 tau,
// 				 tau_dirs,
// 				 qddot,
// 				 ad_qddot,
// 				 &f_ext);

//   CHECK_ARRAY_CLOSE (qddot_ref.data(), qddot.data(), qddot.size(), 1e-7);
//   CHECK_ARRAY_CLOSE (fd_qddot.data(), ad_qddot.data(), fd_qddot.cols()*fd_qddot.rows(), 1e-7);
// }
