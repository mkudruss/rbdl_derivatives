#include <UnitTest++.h>

#include <iostream>

#include "rbdl/rbdl_mathutils.h"
#include "rbdl/Logging.h"

#include "rbdl/Model.h"
#include "rbdl/Kinematics.h"

using namespace std;
using namespace RigidBodyDynamics;
using namespace RigidBodyDynamics::Math;

const double TEST_PREC = 1.0e-12;

struct CartPendulum {
	CartPendulum () {
		using namespace RigidBodyDynamics;
		using namespace RigidBodyDynamics::Math;

		ClearLogOutput();
		model = new Model;

		double cart_w = 0.5;
		double cart_d = 0.2;
		double cart_h = 0.2;
		double cart_m = 10.0;
		Vector3d cart_com (0., 0., 0.);

		double pend_l = 0.5;
		double pend_r = 0.1;
		double pend_m = 1.0;
		Vector3d pend_com (0., 0., pend_l);

		body_cart = Body (cart_m, cart_com, Matrix3d (
					1. / 12. * cart_m * (cart_h*cart_h + cart_d * cart_d), 0., 0.,
					0., 1. / 12. * cart_m * (cart_w*cart_w + cart_h * cart_h), 0.,
					0., 0., 1. / 12. * cart_m * (cart_w*cart_w + cart_d * cart_d)
					));
		body_pendulum = Body (pend_m, pend_com, Matrix3d(
					2. / 5. * pend_m * pend_r, 0., 0.,
					0., 2. / 5. * pend_m * pend_r, 0.,
					0., 0., 2. / 5. * pend_m * pend_r
				));

		joint_cart = Joint (SpatialVector (0., 0., 0., 1., 0., 0.));
		joint_pendulum = Joint (SpatialVector (0., 1., 0., 0., 0., 0.));

		id_cart = model->AddBody (0, Xtrans(Vector3d (0., 0., 0.)), joint_cart, body_cart, "cart");
		id_pendulum = model->AddBody (id_cart, Xtrans(Vector3d (0., 0., 0.)), joint_pendulum, body_pendulum, "pendulum");

		q = VectorNd::Constant ((size_t) model->dof_count, 0.);
		qdot = VectorNd::Constant ((size_t) model->dof_count, 0.);
		qddot = VectorNd::Constant ((size_t) model->dof_count, 0.);
		tau = VectorNd::Constant ((size_t) model->dof_count, 0.);

		body_point = Vector3d (0., 0., pend_l);

		ClearLogOutput();
	}
	~CartPendulum () {
		delete model;
	}

	RigidBodyDynamics::Model *model;

	unsigned int id_cart, id_pendulum;
	RigidBodyDynamics::Body body_cart, body_pendulum;
	RigidBodyDynamics::Joint joint_cart, joint_pendulum;

	RigidBodyDynamics::Math::VectorNd q;
	RigidBodyDynamics::Math::VectorNd qdot;
	RigidBodyDynamics::Math::VectorNd qddot;
	RigidBodyDynamics::Math::VectorNd tau;


	Vector3d body_point;
};

RBDL_DLLAPI
Vector3d CalcBodyToBaseCoordinatesSingleFunc (
		Model &model,
		const VectorNd &Q,
		unsigned int body_id,
		const Vector3d &point_body_coordinates) {
	if (body_id >= model.fixed_body_discriminator) {
		std::cerr << "Fixed bodies not yet supported!" << std::endl;
		abort();
	}

	// Update the kinematics
	VectorNd QDot_zero (VectorNd::Zero (model.q_size));

	for (unsigned int i = 1; i < model.mBodies.size(); i++) {
		unsigned int lambda = model.lambda[i];

		// Calculate joint dependent variables
		if (model.mJoints[i].mJointType == JointTypeRevoluteY) {
			model.X_J[i] = Xroty (Q[model.mJoints[i].q_index]);
		} else if (model.S[i] == SpatialVector (0., 0., 0., 1., 0., 0.)) {
			model.X_J[i] = Xtrans (Vector3d (1., 0., 0.) * Q[model.mJoints[i].q_index]);
		} else {
			std::cerr << "Unsupported joint! Only RotY and TransX supported!" << std::endl;
			abort();
		}

		model.X_lambda[i] = model.X_J[i] * model.X_T[i];
		model.X_base[i] = model.X_lambda[i] * model.X_base[lambda];
	}

	Matrix3d body_rotation = model.X_base[body_id].E.transpose();
	Vector3d body_position = model.X_base[body_id].r;

	return body_position + body_rotation * point_body_coordinates;
}

inline SpatialMatrix dq_Xroty (const double &yrot, const double &dot_yrot) {
	SpatialMatrix result (SpatialMatrix::Zero(6,6));

	double s, c;
	s = sin (yrot) * dot_yrot;
	c = cos (yrot) * dot_yrot;
	Matrix3d E(
				-s, 0., -c,
				0., 0., 0.,
				c, 0., -s
				);

	result.block<3,3>(0,0) = E;
	result.block<3,3>(3,3) = E;

	return result;
}

inline SpatialMatrix dq_Xtrans (const Vector3d &dq_trans) {
	SpatialMatrix result (SpatialMatrix::Zero(6,6));

	result(3,1) =  dq_trans[2];
	result(3,2) = -dq_trans[1];

	result(4,0) = -dq_trans[2];
	result(4,2) =  dq_trans[0];

	result(5,0) =  dq_trans[1];
	result(5,1) = -dq_trans[0];

	return result;
}

RBDL_DLLAPI
Vector3d dq_CalcBodyToBaseCoordinatesSingleFunc (
		Model &model,
		const VectorNd &q,
		const MatrixNd &q_dirs,
		unsigned int body_id,
		const Vector3d &point_body_coordinates,
		MatrixNd &out
		) {
	if (body_id >= model.fixed_body_discriminator) {
		std::cerr << "Fixed bodies not yet supported!" << std::endl;
		abort();
	}

	assert (out.rows() == 3 && out.cols() == model.qdot_size);
	unsigned int ndirs = q_dirs.cols();

	cout << "body_id = " << body_id << endl;
	cout << "q_dirs = " << endl << q_dirs.transpose() << endl;

	// Update the kinematics
	VectorNd QDot_zero (VectorNd::Zero (model.q_size));
	VectorNd fd_out (MatrixNd::Zero (3, model.q_size));

	std::vector<MatrixNd> dq_X_J_i (ndirs, MatrixNd::Zero (6,6));
	std::vector<std::vector<MatrixNd> > dq_X_J (model.mBodies.size(), dq_X_J_i);
	std::vector<std::vector<MatrixNd> > fd_X_J (model.mBodies.size(), dq_X_J_i);
	// dq_X_J[3][5] gives for body 3 the 5th direction

	std::vector<MatrixNd> dq_X_lambda_i (ndirs, MatrixNd::Zero (6,6));
	std::vector<std::vector<MatrixNd> > dq_X_lambda (model.mBodies.size(), dq_X_lambda_i);
	std::vector<std::vector<MatrixNd> > fd_X_lambda (model.mBodies.size(), dq_X_lambda_i);
	// dq_X_lambda[3][5] gives for body 3 the 5th direction

	std::vector<MatrixNd> dq_X_base_i (ndirs, MatrixNd::Zero (6,6));
	std::vector<std::vector<MatrixNd> > dq_X_base (model.mBodies.size(), dq_X_base_i);
	std::vector<std::vector<MatrixNd> > fd_X_base (model.mBodies.size(), dq_X_base_i);
	// dq_X_base[3][5] gives for body 3 the 5th direction

	for (unsigned int i = 1; i < model.mBodies.size(); i++) {
		unsigned int lambda = model.lambda[i];
		// Calculate joint dependent variables
		cout << "== BEGIN == " << endl;
		cout << "== X_J == " << endl;
		//if (model.mJoints[i].mJointType == JointTypeRevoluteY) {
		if (model.S[i] == SpatialVector (0., 1., 0., 0., 0., 0.)) {
			cout << "Revolute" << endl;
			for (unsigned int j = 0; j < ndirs; j++) {
				cout << "q[" << i << "] = " << q[model.mJoints[i].q_index] << endl;
				cout << "q_dirs[" << i << ", " << j << "] = " << q_dirs(i-1, j) << endl;
				// FD TEST
				fd_X_J[i][j] = (
					Xroty (q[model.mJoints[i].q_index] + 1e-08 * q_dirs(i-1, j)).toMatrix()
					- Xroty (q[model.mJoints[i].q_index]).toMatrix()
					) / 1e-08;
				//
				dq_X_J[i][j] = dq_Xroty (q[model.mJoints[i].q_index], q_dirs(i-1,j));
				cout << "fd_X_J[" << i << "][" << j << "] = " << endl;
				cout << fd_X_J[i][j] << endl;
				cout << "dq_X_J[" << i << "][" << j << "] = " << endl;
				cout << dq_X_J[i][j] << endl;
			}
			model.X_J[i] = Xroty (q[model.mJoints[i].q_index]);
		} else if (model.S[i] == SpatialVector (0., 0., 0., 1., 0., 0.)) {
			cout << "Trans" << endl;
			for (unsigned int j = 0; j < ndirs; j++) {
				cout << "q[" << i << "] = " << q[model.mJoints[i].q_index] << endl;
				cout << "q_dirs[" << i << ", " << j << "] = " << q_dirs(i-1, j) << endl;
				// FD TEST
				fd_X_J[i][j] = (
					Xtrans (Vector3d(q[model.mJoints[i].q_index] + 1e-08 * q_dirs(i-1, j), 0. , 0.)).toMatrix()
					- Xtrans (Vector3d (q[model.mJoints[i].q_index], 0., 0.)).toMatrix()
					) / 1e-08;
				//
				dq_X_J[i][j] = dq_Xtrans (Vector3d (q_dirs(i-1, j), 0., 0.));
				cout << "fd_X_J[" << i << "][" << j << "] = " << endl;
				cout << fd_X_J[i][j] << endl;
				cout << "dq_X_J[" << i << "][" << j << "] = " << endl;
				cout << dq_X_J[i][j] << endl;
			}
			model.X_J[i] = Xtrans (Vector3d (1., 0., 0.) * q[model.mJoints[i].q_index]);
		} else {
			std::cerr << "Unsupported joint! Only RotY and TransX supported!" << std::endl;
			abort();
		}
		cout << "== END == " << endl;

		cout << "== BEGIN == " << endl;
		cout << "== X_lambda == " << endl;
		for (unsigned int j = 0; j < ndirs; j++) {
			dq_X_lambda[i][j] = dq_X_J[i][j] * model.X_T[i].toMatrix();
			cout << "dq_X_lambda[" << i << "][" << j << "] = " << endl;
			cout << dq_X_lambda[i][j] << endl;
			// BEGIN FD TEST
			fd_X_lambda[i][j] = (
				(model.X_J[i].toMatrix() + 1e-08 * fd_X_J[i][j]) * model.X_T[i].toMatrix()
				- model.X_J[i].toMatrix() * model.X_T[i].toMatrix()
				) / 1e-08;
			//
			cout << "fd_X_lambda[" << i << "][" << j << "] = " << endl;
			cout << fd_X_lambda[i][j] << endl;
			// END FD TEST
		}
		model.X_lambda[i] = model.X_J[i] * model.X_T[i];
		cout << "== END == " << endl;

		cout << "== BEGIN == " << endl;
		cout << "== X_base == " << endl;
		for (unsigned int j = 0; j < ndirs; j++) {
			dq_X_base[i][j] = dq_X_lambda[i][j] * model.X_base[lambda].toMatrix()
				+ model.X_lambda[i].toMatrix() * dq_X_base[lambda][j];
			cout << "dq_X_base[" << i << "][" << j << "] = " << endl;
			cout << dq_X_base[i][j] << endl;
			// BEGIN FD TEST
			fd_X_base[i][j] = (
					(
						model.X_lambda[i].toMatrix() + 1e-08 * dq_X_lambda[i][j] - model.X_lambda[i].toMatrix()
					) * model.X_base[lambda].toMatrix()
					+
					model.X_lambda[i].toMatrix() * (
						model.X_base[lambda].toMatrix() + 1e-08 * dq_X_base[lambda][j] - model.X_base[lambda].toMatrix()
					)
				) / 1e-08;
			//
			cout << "fd_X_base[" << i << "][" << j << "] = " << endl;
			cout << fd_X_base[i][j] << endl;
			// END FD TEST
		}
		model.X_base[i] = model.X_lambda[i] * model.X_base[lambda];
		cout << "== END == " << endl;
	}

	Matrix3d body_rotation = model.X_base[body_id].E.transpose();
	Vector3d body_position = model.X_base[body_id].r;

	cout << "== last loop ==" << endl;

	for (unsigned int j = 0; j < ndirs; j++) {
		SpatialTransform dq_base = SpatialTransform::fromMatrix (dq_X_base[body_id][j]);
		Vector3d vec3 = dq_base.r + dq_base.E.transpose() * point_body_coordinates;
		cout << "col = " <<  vec3.transpose() << endl;
		out.block<3, 1>(0,j) = vec3;
	}

	return model.X_base[body_id].r + model.X_base[body_id].E.transpose() * point_body_coordinates;
}

RBDL_DLLAPI
Vector3d fd_dq_CalcBodyToBaseCoordinatesSingleFunc (
		Model &model,
		const VectorNd &q,
		const MatrixNd &q_dirs,
		unsigned int body_id,
		const Vector3d &point_body_coordinates,
		MatrixNd &out
		) {
	Vector3d ref = CalcBodyToBaseCoordinatesSingleFunc (model, q, body_id, point_body_coordinates);

	unsigned int ndirs = q_dirs.cols();
	double h = 1.0e-8;
	for (unsigned int j = 0; j < ndirs; j++) {
		VectorNd q_dir = q_dirs.block(0,j, model.qdot_size, 1);
		Vector3d res_hd = CalcBodyToBaseCoordinatesSingleFunc (model, q + h * q_dir, body_id, point_body_coordinates);
		Vector3d res_hd_rbdl = CalcBodyToBaseCoordinates (model, q + h * q_dir, body_id, point_body_coordinates);

		cout << "calc_body err = " << (res_hd - res_hd_rbdl).transpose() << endl;

		out.block<3,1>(0,j) = (res_hd - ref) / h;
	}

	return ref;
}

TEST_FIXTURE ( CartPendulum, CartPendulumCalcBodyToBaseSimple) {
	q[0] = 0.3;
	q[1] = -0.2;
	Vector3d point_body_coordinates (0.1, 3.2, 4.2);

	Vector3d point_single_func = CalcBodyToBaseCoordinatesSingleFunc (*model, q, id_pendulum, point_body_coordinates);
	Vector3d point_default = CalcBodyToBaseCoordinates (*model, q, id_pendulum, point_body_coordinates);

	CHECK_ARRAY_CLOSE (point_default.data(), point_single_func.data(), 3, TEST_PREC);
}

TEST_FIXTURE ( CartPendulum, CartPendulumJacobianADSimple ) {
	MatrixNd jacobian_ad = MatrixNd::Zero(3, model->qdot_size);
	MatrixNd jacobian_ref = MatrixNd::Zero(3, model->qdot_size);
	MatrixNd jacobian_fd = MatrixNd::Zero(3, model->qdot_size);

	q.setZero();
	q[0] = 2.0;
	q[1] = 0.3;
	body_point = Vector3d (0., 0., 0.6);

	CalcPointJacobian (*model, q, id_pendulum, body_point, jacobian_ref);

	MatrixNd q_dirs = MatrixNd::Identity (model->qdot_size, model->qdot_size);
	Vector3d base_point_standard = CalcBodyToBaseCoordinates (*model, q, id_pendulum, body_point);
	Vector3d base_point_ad = dq_CalcBodyToBaseCoordinatesSingleFunc (*model, q, q_dirs, id_pendulum, body_point, jacobian_ad);
	Vector3d base_point_fd = fd_dq_CalcBodyToBaseCoordinatesSingleFunc (*model, q, q_dirs, id_pendulum, body_point, jacobian_fd);

	cout << "point err = " << (base_point_standard - base_point_ad).transpose() << endl;

	cout << "jacobian_ref: " << endl << jacobian_ref << endl;
	cout << "jacobian_ad: " << endl << jacobian_ad << endl;
	cout << "Jacobian error (AD, ref):" << endl << (jacobian_ad - jacobian_ref) << endl;
	cout << "Jacobian error (FD, ref):" << endl << (jacobian_fd - jacobian_ref) << endl;
//	cout << "Jacobian error (AD, FD):" << endl << (jacobian_ad - jacobian_fd) << endl;

	CHECK_ARRAY_CLOSE (jacobian_ref.data(), jacobian_ad.data(), 3 * model->qdot_size, TEST_PREC);
//	CHECK_ARRAY_CLOSE (v_fixed_body.data(), v_body.data(), 6, TEST_PREC);
}
