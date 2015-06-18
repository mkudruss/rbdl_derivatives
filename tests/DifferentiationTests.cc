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

struct ADModel {
	ADModel () {};
	void fromModel () {};

	// TODO: Add getter methods with ndirs to check if dimensions fit for functions that write into data structure

	// NOTE: old initialization values
	//std::vector<MatrixNd> ad_X_lambda_i (ndirs, MatrixNd::Zero (6,6));
	//std::vector<MatrixNd> ad_X_Ji (ndirs, MatrixNd::Zero (6,6));
	// std::vector<SpatialVector> ad_V (ndirs, SpatialVector::Zero (6));

	// derivative values
	// TODO: remove ad_ prefix?
	std::vector<std::vector<MatrixNd> > ad_X_lambda; // (model.mBodies.size(), ad_X_lambda_i);
	std::vector<std::vector<MatrixNd> > ad_X_J; // (model->mBodies.size(), ad_X_Ji);

	std::vector<std::vector<SpatialVector> > ad_S;   //(model->mBodies.size(), ad_V);
	std::vector<std::vector<SpatialVector> > ad_v_J; //(model->mBodies.size(), ad_V);
	std::vector<std::vector<SpatialVector> > ad_c_J; //(model->mBodies.size(), ad_V);

};

Matrix3d E_from_Matrix(const SpatialMatrix X) {
	Matrix3d E = Matrix3d::Zero();
	E = X.block<3, 3>(0,0);
	return E;
}

Vector3d r_from_Matrix(const SpatialMatrix X) {
	Matrix3d E = E_from_Matrix(X);
	Matrix3d Erx = Matrix3d::Zero();
	Erx = X.block<3,3>(3,0);
	Matrix3d rx = E.transpose() * Erx;
	Vector3d r = Vector3d::Zero();
	r(0) = -rx(2,1);
	r(1) =  rx(2,0);
	r(2) = -rx(1,0);
	return r;
}

Vector3d ad_r_from_Matrix(const SpatialMatrix X, const SpatialMatrix X_dot ) {
	Matrix3d E_dot = E_from_Matrix(X_dot);
	Matrix3d E = E_from_Matrix(X);

	Matrix3d Erx_dot = Matrix3d::Zero();
	Matrix3d Erx = Matrix3d::Zero();
	Erx_dot = X_dot.block<3,3>(3,0);
	Erx = X.block<3,3>(3,0);

	Matrix3d rx_dot = E_dot.transpose() * Erx + E.transpose() * Erx_dot;
	Matrix3d rx = E.transpose() * Erx;

	Vector3d r_dot = Vector3d::Zero();
	r_dot(0) = -rx_dot(2,1);
	r_dot(1) =  rx_dot(2,0);
	r_dot(2) = -rx_dot(1,0);

	return r_dot;
}

inline SpatialMatrix ad_Xroty (const double &yrot, const double &dot_yrot) {
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

inline SpatialMatrix ad_Xtrans (const Vector3d &ad_trans) {
	SpatialMatrix result (SpatialMatrix::Zero(6,6));

	result(3,1) =  ad_trans[2];
	result(3,2) = -ad_trans[1];

	result(4,0) = -ad_trans[2];
	result(4,2) =  ad_trans[0];

	result(5,0) =  ad_trans[1];
	result(5,1) = -ad_trans[0];

	return result;
}

RBDL_DLLAPI
void ad_jcalc (
	Model &model,
	//ADModel &ad_model,
	unsigned int joint_id,
	const VectorNd &q,
	const MatrixNd &q_dirs,
	const VectorNd &qdot,
	const MatrixNd &qdot_dirs,
	std::vector<MatrixNd> &ad_X_Ji,
	std::vector<SpatialVector> &ad_S_i,
	std::vector<SpatialVector> &ad_v_Ji,
	std::vector<SpatialVector> &ad_c_Ji
) {
	unsigned int ndirs = q_dirs.cols();

	// derivative intermediate values
	std::vector<MatrixNd> ad_X_lambda_i (ndirs, MatrixNd::Zero (6,6));
	std::vector<std::vector<MatrixNd> > ad_X_lambda (model.mBodies.size(), ad_X_lambda_i);

	// check input dimensions
	if (q_dirs.cols() != qdot_dirs.cols()) {
		std::cerr << "directions have different dimensions: " << "#q_dirs = " << q_dirs.cols() << " != " << qdot_dirs.cols() << " = #qdot_dirs." << std::endl;
		std::cerr << "In: " << __func__ << endl;
		abort();
	}
	// check output dimensions
	if (ad_X_Ji.size() != ndirs) {
		std::cerr << "derivative does not have proper dimensions " << "#ad_X_Ji.size() = " << ad_X_Ji.size() << " != " << ndirs << " = ndirs." << std::endl;
		std::cerr << "In: " << __func__ << endl;
		abort();
	}
	if (ad_S_i.size() != ndirs) {
		std::cerr << "derivative does not have proper dimensions " << "#ad_S_i.size() = " << ad_S_i.size() << " != " << ndirs << " = ndirs." << std::endl;
		std::cerr << "In: " << __func__ << endl;
		abort();
	}
	if (ad_c_Ji.size() != ndirs) {
		std::cerr << "derivative does not have proper dimensions " << "#ad_c_Ji.size() = " << ad_c_Ji.size() << " != " << ndirs << " = ndirs." << std::endl;
		std::cerr << "In: " << __func__ << endl;
		abort();
	}

	// exception if we calculate it for the root body
	assert (joint_id > 0);

	if (model.mJoints[joint_id].mJointType == JointTypeRevoluteX) {
		std::cerr << "In: " << __func__ << "this joint type is not supported." << std::endl;
		abort();
	} else if (model.mJoints[joint_id].mJointType == JointTypeRevoluteY) {
		// derivative code
		for (int idir = 0; idir < ndirs; ++idir) {
			ad_X_Ji[idir] = ad_Xroty (q[model.mJoints[joint_id].q_index], q_dirs(model.mJoints[joint_id].q_index, idir));
			ad_S_i[idir]  = SpatialVector::Zero(); // S = [0., 1., 0., 0., 0., 0.]
			ad_v_Ji[idir][1] = qdot_dirs(model.mJoints[joint_id].q_index, idir); // v_J = S*qdot
			ad_c_Ji[idir] = SpatialVector::Zero(); // c_J = Sdot*qdot
		}
		// nominal code
		model.X_J[joint_id] = Xroty (q[model.mJoints[joint_id].q_index]);
		model.v_J[joint_id][1] = qdot[model.mJoints[joint_id].q_index];
	} else if (model.mJoints[joint_id].mJointType == JointTypeRevoluteZ) {
		std::cerr << "In: " << __func__ << "this joint type is not supported." << std::endl;
		abort();
	} else if (model.S[joint_id] == SpatialVector (0., 0., 0., 1., 0., 0.)) {
		// derivative code
		for (int idir = 0; idir < ndirs; ++idir) {
			ad_X_Ji[idir] = ad_Xtrans (Vector3d (q_dirs(model.mJoints[joint_id].q_index, idir), 0.0, 0.0));
			ad_S_i[idir]  = SpatialVector::Zero(); // S = [0., 0., 0., 1., 0., 0.]
			ad_v_Ji[idir][3] = qdot_dirs(model.mJoints[joint_id].q_index, idir); // v_J = S*qdot
			ad_c_Ji[idir] = SpatialVector::Zero(); // v_J = Sdot*qdot
		}
		// nominal code
		model.X_J[joint_id] = Xtrans (Vector3d (q[model.mJoints[joint_id].q_index], 0., 0.));
		model.v_J[joint_id][3] = qdot[model.mJoints[joint_id].q_index];
	} else if (model.mJoints[joint_id].mDoFCount == 1) {
		std::cerr << "In: " << __func__ << "this joint type is not supported." << std::endl;
		abort();
	} else if (model.mJoints[joint_id].mJointType == JointTypeSpherical) {
		std::cerr << "In: " << __func__ << "this joint type is not supported." << std::endl;
		abort();
	} else if (model.mJoints[joint_id].mJointType == JointTypeEulerZYX) {
		std::cerr << "In: " << __func__ << "this joint type is not supported." << std::endl;
		abort();
	} else if (model.mJoints[joint_id].mJointType == JointTypeEulerXYZ) {
		std::cerr << "In: " << __func__ << "this joint type is not supported." << std::endl;
		abort();
	} else if (model.mJoints[joint_id].mJointType == JointTypeEulerYXZ) {
		std::cerr << "In: " << __func__ << "this joint type is not supported." << std::endl;
		abort();
	} else if (model.mJoints[joint_id].mJointType == JointTypeTranslationXYZ) {
		std::cerr << "In: " << __func__ << "this joint type is not supported." << std::endl;
		abort();
	} else {
		std::cerr << "Error: invalid joint type " << model.mJoints[joint_id].mJointType << " at id " << joint_id << std::endl;
		abort();
	}
	// derivative code
	for (int idir = 0; idir < ndirs; ++idir) {
		// TODO: we need ad_X_lambda in ad data structures
		ad_X_lambda[joint_id][idir] = ad_X_Ji[idir] * model.X_T[joint_id].toMatrix();
	}
	// nominal code
	model.X_lambda[joint_id] = model.X_J[joint_id] * model.X_T[joint_id];
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

	Matrix3d body_rotation = E_from_Matrix(model.X_base[body_id].toMatrix());
	Vector3d body_position = r_from_Matrix(model.X_base[body_id].toMatrix());

	return body_position + body_rotation.transpose() * point_body_coordinates;
}

RBDL_DLLAPI
Vector3d ad_CalcBodyToBaseCoordinatesSingleFunc (
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

	// Update the kinematics
	VectorNd QDot_zero (VectorNd::Zero (model.q_size));
	VectorNd fd_out (MatrixNd::Zero (3, model.q_size));

	std::vector<MatrixNd> ad_X_J_i (ndirs, MatrixNd::Zero (6,6));
	std::vector<std::vector<MatrixNd> > ad_X_J (model.mBodies.size(), ad_X_J_i);
	std::vector<std::vector<MatrixNd> > fd_X_J (model.mBodies.size(), ad_X_J_i);
	// ad_X_J[3][5] gives for body 3 the 5th direction

	std::vector<MatrixNd> ad_X_lambda_i (ndirs, MatrixNd::Zero (6,6));
	std::vector<std::vector<MatrixNd> > ad_X_lambda (model.mBodies.size(), ad_X_lambda_i);
	std::vector<std::vector<MatrixNd> > fd_X_lambda (model.mBodies.size(), ad_X_lambda_i);
	// ad_X_lambda[3][5] gives for body 3 the 5th direction

	std::vector<MatrixNd> ad_X_base_i (ndirs, MatrixNd::Zero (6,6));
	std::vector<std::vector<MatrixNd> > ad_X_base (model.mBodies.size(), ad_X_base_i);
	std::vector<std::vector<MatrixNd> > fd_X_base (model.mBodies.size(), ad_X_base_i);
	// ad_X_base[3][5] gives for body 3 the 5th direction

	for (unsigned int i = 1; i < model.mBodies.size(); i++) {
		unsigned int lambda = model.lambda[i];
		// Calculate joint dependent variables
		if (model.mJoints[i].mJointType == JointTypeRevoluteY) {
			for (unsigned int j = 0; j < ndirs; j++) {
				ad_X_J[i][j] = ad_Xroty (q[model.mJoints[i].q_index], q_dirs(i-1,j));
			}
			model.X_J[i] = Xroty (q[model.mJoints[i].q_index]);
		} else if (model.S[i] == SpatialVector (0., 0., 0., 1., 0., 0.)) {
			for (unsigned int j = 0; j < ndirs; j++) {
				ad_X_J[i][j] = ad_Xtrans (Vector3d (q_dirs(i-1, j), 0., 0.));
			}
			model.X_J[i] = Xtrans (Vector3d (1., 0., 0.) * q[model.mJoints[i].q_index]);
		} else {
			std::cerr << "Unsupported joint! Only RotY and TransX supported!" << std::endl;
			abort();
		}

		for (unsigned int j = 0; j < ndirs; j++) {
			ad_X_lambda[i][j] = ad_X_J[i][j] * model.X_T[i].toMatrix();
		}
		model.X_lambda[i] = model.X_J[i] * model.X_T[i];

		for (unsigned int j = 0; j < ndirs; j++) {
			ad_X_base[i][j] = ad_X_lambda[i][j] * model.X_base[lambda].toMatrix() + model.X_lambda[i].toMatrix() * ad_X_base[lambda][j];
		}
		model.X_base[i] = model.X_lambda[i] * model.X_base[lambda];
	}

	for (unsigned int j = 0; j < ndirs; j++) {
		SpatialMatrix X_base_ib = model.X_base[body_id].toMatrix();
		Matrix3d ad_E = E_from_Matrix(ad_X_base[body_id][j]);
		Vector3d ad_r = ad_r_from_Matrix(X_base_ib, ad_X_base[body_id][j]);

		out.block<3,1>(0,j) = ad_r + ad_E.transpose() * point_body_coordinates;
	}

	Matrix3d body_rotation = E_from_Matrix(model.X_base[body_id].toMatrix());
	Vector3d body_position = r_from_Matrix(model.X_base[body_id].toMatrix());

	return body_position + body_rotation.transpose() * point_body_coordinates;
}

RBDL_DLLAPI
Vector3d fd_jcalc(
	Model &model,
	unsigned int joint_id,
	const VectorNd &q,
	const MatrixNd &q_dirs,
	const VectorNd &qdot,
	const MatrixNd &qdot_dirs,
	std::vector<MatrixNd> &ad_X_Ji,
	std::vector<SpatialVector> &ad_S_i,
	std::vector<SpatialVector> &ad_v_Ji,
	std::vector<SpatialVector> &ad_c_Ji
) {
	double h = 1.0e-8;
	unsigned int ndirs = q_dirs.cols();

	// check input dimensions
	if (q_dirs.cols() != qdot_dirs.cols()) {
		std::cerr << "directions have different dimensions: " << "#q_dirs = " << q_dirs.cols() << " != " << qdot_dirs.cols() << " = #qdot_dirs." << std::endl;
		std::cerr << "In: " << __func__ << endl;
		abort();
	}
	// check output dimensions
	if (ad_X_Ji.size() != ndirs) {
		std::cerr << "derivative does not have proper dimensions " << "#ad_X_Ji.size() = " << ad_X_Ji.size() << " != " << ndirs << " = ndirs." << std::endl;
		std::cerr << "In: " << __func__ << endl;
		abort();
	}
	if (ad_S_i.size() != ndirs) {
		std::cerr << "derivative does not have proper dimensions " << "#ad_S_i.size() = " << ad_S_i.size() << " != " << ndirs << " = ndirs." << std::endl;
		std::cerr << "In: " << __func__ << endl;
		abort();
	}
	if (ad_c_Ji.size() != ndirs) {
		std::cerr << "derivative does not have proper dimensions " << "#ad_c_Ji.size() = " << ad_c_Ji.size() << " != " << ndirs << " = ndirs." << std::endl;
		std::cerr << "In: " << __func__ << endl;
		abort();
	}

	// calculate y(t)
	jcalc (model, joint_id, q, qdot);
	MatrixNd ref_X_Ji = model.X_J[joint_id].toMatrix();
	SpatialVector ref_S_i = model.S[joint_id];
	SpatialVector ref_v_Ji = model.v_J[joint_id];
	SpatialVector ref_c_Ji = model.c_J[joint_id];

	for (unsigned int j = 0; j < ndirs; j++) {
		VectorNd q_dir = q_dirs.block(0,j, model.qdot_size, 1);
		VectorNd qdot_dir = qdot_dirs.block(0,j, model.qdot_size, 1);
		jcalc (model, joint_id, q + h * q_dir, qdot + h * qdot_dir);

		MatrixNd hd_X_Ji = model.X_J[joint_id].toMatrix();
		SpatialVector hd_S_i = model.S[joint_id];
		SpatialVector hd_v_Ji = model.v_J[joint_id];
		SpatialVector hd_c_Ji = model.c_J[joint_id];

		ad_X_Ji[j] = (hd_X_Ji - ref_X_Ji) / h;
		ad_S_i[j]  = (hd_S_i  - ref_S_i)  / h;
		ad_v_Ji[j] = (hd_v_Ji - ref_v_Ji) / h;
		ad_c_Ji[j] = (hd_c_Ji - ref_c_Ji) / h;
	}
};

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

		//cout << "calc_body err = " << (res_hd - res_hd_rbdl).transpose() << endl;

		out.block<3,1>(0,j) = (res_hd - ref) / h;
	}

	return ref;
}

TEST_FIXTURE (CartPendulum, jcalcNominalSolutionTest) {
	cout << "Nominal Evaluation Test" << endl;

	// set nominal values
	q.setZero();
	q[0] = 0.3;
	q[1] = -0.2;

	qdot.setZero();
	qdot[0] = 0.3;
	qdot[1] = -0.2;

	// set directions
	unsigned int ndirs = model->q_size + model->qdot_size;
	MatrixNd x = MatrixNd::Identity(ndirs, ndirs);
	MatrixNd q_dirs = x.block(0, 0, model->q_size, ndirs);
	MatrixNd qdot_dirs = x.block(model->q_size, 0, model->qdot_size, ndirs);

	// set derivative outputs
	std::vector<MatrixNd> ad_X_Ji (ndirs, MatrixNd::Zero (6,6));
	std::vector<std::vector<MatrixNd> > ad_X_J (model->mBodies.size(), ad_X_Ji);

	std::vector<SpatialVector> ad_V (ndirs, SpatialVector::Zero (6));
	std::vector<std::vector<SpatialVector> > ad_S (model->mBodies.size(), ad_V);
	std::vector<std::vector<SpatialVector> > ad_v_J (model->mBodies.size(), ad_V);
	std::vector<std::vector<SpatialVector> > ad_c_J (model->mBodies.size(), ad_V);

	for (unsigned int joint_id = 1; joint_id < model->mBodies.size(); joint_id++) {
		// evaluate nominal solution
		jcalc (*model, joint_id, q, qdot);
		MatrixNd ref_X_Ji = model->X_J[joint_id].toMatrix();
		SpatialVector ref_S_i = model->S[joint_id];
		SpatialVector ref_v_Ji = model->v_J[joint_id];
		SpatialVector ref_c_Ji = model->c_J[joint_id];

		// evaluate AD nominal solution
		ad_jcalc (*model, joint_id, q, q_dirs, qdot, qdot_dirs,
			ad_X_J[joint_id], ad_S[joint_id], ad_v_J[joint_id], ad_c_J[joint_id]);
		MatrixNd test_X_Ji = model->X_J[joint_id].toMatrix();
		SpatialVector test_S_i = model->S[joint_id];
		SpatialVector test_v_Ji = model->v_J[joint_id];
		SpatialVector test_c_Ji = model->c_J[joint_id];

		CHECK_ARRAY_CLOSE (ref_X_Ji.data(), test_X_Ji.data(), 36, TEST_PREC);
		CHECK_ARRAY_CLOSE (ref_S_i.data(),  test_S_i.data(),   6, TEST_PREC);
		CHECK_ARRAY_CLOSE (ref_v_Ji.data(), test_v_Ji.data(),  6, TEST_PREC);
		CHECK_ARRAY_CLOSE (ref_c_Ji.data(), test_c_Ji.data(),  6, TEST_PREC);
		/* DEBUG OUTPUT
		cout << "===== joint_id: " << joint_id << " =====" << endl;
		cout << "ref_X_Ji: " << endl << ref_X_Ji << endl;
		cout << "test_X_Ji: " << endl << test_X_Ji << endl;
		cout << "error_X_Ji: " << endl << ref_X_Ji - test_X_Ji << endl;
		cout << endl;

		cout << "ref_S_i: " << endl << ref_S_i << endl;
		cout << "test_S_i: " << endl << test_S_i << endl;
		cout << "error_S_i: " << endl << ref_S_i - test_S_i << endl;
		cout << endl;

		cout << "ref_v_Ji: " << endl << ref_v_Ji << endl;
		cout << "test_v_Ji: " << endl << test_v_Ji << endl;
		cout << "error_v_Ji: " << endl << ref_v_Ji - test_v_Ji << endl;
		cout << endl;

		cout << "ref_c_Ji: " << endl << ref_c_Ji << endl;
		cout << "test_c_Ji: " << endl << test_c_Ji << endl;
		cout << "error_c_Ji: " << endl << ref_c_Ji - test_c_Ji << endl;
		cout << endl;
		cout << endl;
		*/
	}
}

TEST_FIXTURE (CartPendulum, jcalcFDvsADTest) {
	cout << "Derivative Evaluation Test" << endl;
	double TEST_PREC = 1.0e-08;

	// set nominal values
	q.setZero();
	q[0] = 0.3;
	q[1] = -0.2;

	qdot.setZero();
	qdot[0] = 0.3;
	qdot[1] = -0.2;

	// set directions
	unsigned int ndirs = model->q_size + model->qdot_size;
	MatrixNd x = MatrixNd::Identity(ndirs, ndirs);
	MatrixNd q_dirs = x.block(0, 0, model->q_size, ndirs);
	MatrixNd qdot_dirs = x.block(model->q_size, 0, model->qdot_size, ndirs);

	// set derivative outputs
	std::vector<MatrixNd> ad_X_Ji (ndirs, MatrixNd::Zero (6,6));
	std::vector<std::vector<MatrixNd> > fd_X_J (model->mBodies.size(), ad_X_Ji);
	std::vector<std::vector<MatrixNd> > ad_X_J (model->mBodies.size(), ad_X_Ji);

	std::vector<SpatialVector> ad_V (ndirs, SpatialVector::Zero (6));
	std::vector<std::vector<SpatialVector> > fd_S (model->mBodies.size(), ad_V);
	std::vector<std::vector<SpatialVector> > fd_v_J (model->mBodies.size(), ad_V);
	std::vector<std::vector<SpatialVector> > fd_c_J (model->mBodies.size(), ad_V);

	std::vector<std::vector<SpatialVector> > ad_S (model->mBodies.size(), ad_V);
	std::vector<std::vector<SpatialVector> > ad_v_J (model->mBodies.size(), ad_V);
	std::vector<std::vector<SpatialVector> > ad_c_J (model->mBodies.size(), ad_V);

	for (unsigned int joint_id = 1; joint_id < model->mBodies.size(); joint_id++) {
		// evaluate nominal solution
		fd_jcalc (*model, joint_id, q, q_dirs, qdot, qdot_dirs,
			fd_X_J[joint_id], fd_S[joint_id], fd_v_J[joint_id], fd_c_J[joint_id]);

		// evaluate AD nominal solution
		ad_jcalc (*model, joint_id, q, q_dirs, qdot, qdot_dirs,
			ad_X_J[joint_id], ad_S[joint_id], ad_v_J[joint_id], ad_c_J[joint_id]);

		for (int idir = 0; idir < ndirs; ++idir) {
			CHECK_ARRAY_CLOSE (fd_X_J[joint_id][idir].data(), ad_X_J[joint_id][idir].data(), 36, TEST_PREC);
			CHECK_ARRAY_CLOSE (fd_S[joint_id][idir].data(),   ad_S[joint_id][idir].data(),    6, TEST_PREC);
			CHECK_ARRAY_CLOSE (fd_v_J[joint_id][idir].data(), ad_v_J[joint_id][idir].data(),  6, TEST_PREC);
			CHECK_ARRAY_CLOSE (fd_c_J[joint_id][idir].data(), ad_c_J[joint_id][idir].data(),  6, TEST_PREC);
			/* DEBUG OUTPUT
			cout << "===== joint_id: " << joint_id << ", idir: " << idir << " =====" << endl;
			cout << "fd_X_J[" << joint_id << "][" << idir << "]: " << endl << fd_X_J[joint_id][idir] << endl;
			cout << "ad_X_J[" << joint_id << "][" << idir << "]: " << endl << ad_X_J[joint_id][idir] << endl;
			cout << "error_X_J: " << endl << fd_X_J[joint_id][idir] - ad_X_J[joint_id][idir] << endl;

			cout << "fd_S[" << joint_id << "][" << idir << "]: " << endl << fd_S[joint_id][idir] << endl;
			cout << "ad_S[" << joint_id << "][" << idir << "]: " << endl << ad_S[joint_id][idir] << endl;
			cout << "error_S: " << endl << fd_S[joint_id][idir] - ad_S[joint_id][idir] << endl;

			cout << "fd_v_J[" << joint_id << "][" << idir << "]: " << endl << fd_v_J[joint_id][idir] << endl;
			cout << "ad_v_J[" << joint_id << "][" << idir << "]: " << endl << ad_v_J[joint_id][idir] << endl;
			cout << "error_v_J: " << endl << fd_v_J[joint_id][idir] - ad_v_J[joint_id][idir] << endl;

			cout << "fd_c_J[" << joint_id << "][" << idir << "]: " << endl << fd_c_J[joint_id][idir] << endl;
			cout << "ad_c_J[" << joint_id << "][" << idir << "]: " << endl << ad_c_J[joint_id][idir] << endl;
			cout << "error_c_J: " << endl << fd_c_J[joint_id][idir] - ad_c_J[joint_id][idir] << endl;
			*/
		}
	}
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
	q[0] = -1.0;
	q[1] = -0.3;
	body_point = Vector3d (1.0, 2.0, 3.0);

	CalcPointJacobian (*model, q, id_pendulum, body_point, jacobian_ref);

	MatrixNd q_dirs = MatrixNd::Identity (model->qdot_size, model->qdot_size);
	Vector3d base_point_standard = CalcBodyToBaseCoordinates (*model, q, id_pendulum, body_point);
	Vector3d base_point_ad = ad_CalcBodyToBaseCoordinatesSingleFunc (*model, q, q_dirs, id_pendulum, body_point, jacobian_ad);
	Vector3d base_point_fd = fd_dq_CalcBodyToBaseCoordinatesSingleFunc (*model, q, q_dirs, id_pendulum, body_point, jacobian_fd);

	CHECK_ARRAY_CLOSE (jacobian_ref.data(), jacobian_ad.data(), 3 * model->qdot_size, TEST_PREC);
//	CHECK_ARRAY_CLOSE (v_fixed_body.data(), v_body.data(), 6, TEST_PREC);
}
