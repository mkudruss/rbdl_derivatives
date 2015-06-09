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

		body_point.setZero();

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
			model.X_J[i] = Xtrans (Vector3d (1., 0., 0.));
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
	

	// Update the kinematics
	VectorNd QDot_zero (VectorNd::Zero (model.q_size));

	for (unsigned int i = 1; i < model.mBodies.size(); i++) {
		unsigned int lambda = model.lambda[i];

		// Calculate joint dependent variables
		if (model.mJoints[i].mJointType == JointTypeRevoluteY) {
			model.X_J[i] = Xroty (q[model.mJoints[i].q_index]);
		} else if (model.S[i] == SpatialVector (0., 0., 0., 1., 0., 0.)) {
			model.X_J[i] = Xtrans (Vector3d (1., 0., 0.) * q[model.mJoints[i].q_index]);
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

TEST_FIXTURE ( CartPendulum, CartPendulumJacobianSimple ) {
	MatrixNd jacobian_ad = MatrixNd::Zero(3, model->qdot_size);
	MatrixNd jacobian_ref = MatrixNd::Zero(3, model->qdot_size);

	q.setZero();

	CalcPointJacobian (*model, q, id_pendulum, body_point, jacobian_ref);

	MatrixNd q_dirs = MatrixNd::Identity (model->qdot_size, model->qdot_size);
	Vector3d base_point = dq_CalcBodyToBaseCoordinatesSingleFunc (*model, q, q_dirs, id_pendulum, body_point, jacobian_ad);

	cout << "Jacobian error:" << endl << (jacobian_ref - jacobian_ad).transpose() << endl;

	CHECK_ARRAY_CLOSE (jacobian_ref.data(), jacobian_ad.data(), 3 * model->qdot_size, TEST_PREC);
//	CHECK_ARRAY_CLOSE (v_fixed_body.data(), v_body.data(), 6, TEST_PREC);
}
