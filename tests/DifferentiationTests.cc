#include <UnitTest++.h>

#include <iostream>

#include "rbdl/rbdl_mathutils.h"
#include "rbdl/Logging.h"

#include "rbdl/Model.h"
#include "rbdl/Kinematics.h"
#include "rbdl/Dynamics.h"

using namespace std;
using namespace RigidBodyDynamics;
using namespace RigidBodyDynamics::Math;

const double TEST_PREC = 1.0e-12;

struct ADModel {
	ADModel () {
		// TODO: should not be called, I think!
		//std::cerr << "Data structure can only be initialized with argument RigidBodyDynamics::Model!" << std::endl;
		//throw std::exception;
	};

	ADModel (Model& model) {
		ndirs = 4 * model.dof_count;

		// NOTE: old initialization values
		std::vector<SpatialMatrix> ad_X(ndirs, SpatialMatrix::Zero());
		std::vector<SpatialVector> ad_vec(ndirs, SpatialVector::Zero());
		//std::vector<double> ad_double(ndirs, 0.0);

		ad_X_lambda.resize(model.mBodies.size(), ad_X);
		ad_X_base.resize(model.mBodies.size(), ad_X);
		ad_X_J.resize(model.mBodies.size(), ad_X);

		ad_S.resize(model.mBodies.size(), ad_vec);
		ad_v_J.resize(model.mBodies.size(), ad_vec);
		ad_v.resize(model.mBodies.size(), ad_vec);
		ad_c_J.resize(model.mBodies.size(), ad_vec);
		ad_a_J.resize(model.mBodies.size(), ad_vec);
		ad_a.resize(model.mBodies.size(), ad_vec);
		
		ad_I_ci.resize(model.mBodies.size(), ad_X);
		ad_F.resize(model.mBodies.size(), ad_vec);
		ad_c.resize(model.mBodies.size(), ad_vec);
		ad_f.resize(model.mBodies.size(), ad_vec);
		
		ad_pA.resize(model.mBodies.size(), ad_vec);
		ad_U.resize(model.mBodies.size(), ad_vec);
		ad_IA.resize(model.mBodies.size(), ad_X);
		//ad_u.resize(model.mBodies.size(), ad_double);
		//ad_d.resize(model.mBodies.size(), ad_double);
		ad_u.resize(model.q_size, ndirs);
		ad_d.resize(model.q_size, ndirs);
	};

	unsigned int ndirs;

	// derivative values
	// TODO: remove ad_ prefix?
	std::vector<std::vector<SpatialMatrix> > ad_X_lambda;
	std::vector<std::vector<SpatialMatrix> > ad_X_base;
	std::vector<std::vector<SpatialMatrix> > ad_X_J;

	std::vector<std::vector<SpatialVector> > ad_S;
	std::vector<std::vector<SpatialVector> > ad_v_J;
	std::vector<std::vector<SpatialVector> > ad_v;
	std::vector<std::vector<SpatialVector> > ad_c_J;
	std::vector<std::vector<SpatialVector> > ad_a_J;
	std::vector<std::vector<SpatialVector> > ad_a;
	std::vector<std::vector<SpatialVector> > ad_c;
	std::vector<std::vector<SpatialVector> > ad_f;

	std::vector<std::vector<SpatialMatrix> > ad_I_ci;
	std::vector<std::vector<SpatialVector> > ad_F;
	
	std::vector<std::vector<SpatialVector> > ad_pA;
	std::vector<std::vector<SpatialVector> > ad_U;
	std::vector<std::vector<SpatialMatrix> > ad_IA;
	//std::vector<std::vector<double> > ad_u;
	//std::vector<std::vector<double> > ad_d;
	MatrixNd ad_u;
	MatrixNd ad_d;

	void resize_directions (unsigned requested_ndirs){
		if (ndirs < requested_ndirs) {
			ndirs = requested_ndirs;

			for (int i = 0; i < ad_X_lambda.size(); ++i) {
				ad_X_lambda[i].resize(ndirs, SpatialMatrix::Zero());
			}

			for (int i = 0; i < ad_X_base.size(); ++i) {
				ad_X_base[i].resize(ndirs, SpatialMatrix::Zero());
			}

			for (int i = 0; i < ad_X_J.size(); ++i) {
				ad_X_J[i].resize(ndirs, SpatialMatrix::Zero());
			}

			for (int i = 0; i < ad_S.size(); ++i) {
				ad_S[i].resize(ndirs, SpatialVector::Zero());
			}

			for (int i = 0; i < ad_v_J.size(); ++i) {
				ad_v_J[i].resize(ndirs, SpatialVector::Zero());
			}

			for (int i = 0; i < ad_v_J.size(); ++i) {
				ad_v[i].resize(ndirs, SpatialVector::Zero());
			}

			for (int i = 0; i < ad_c_J.size(); ++i) {
				ad_c_J[i].resize(ndirs, SpatialVector::Zero());
			}

			for (int i = 0; i < ad_c_J.size(); ++i) {
				ad_a_J[i].resize(ndirs, SpatialVector::Zero());
			}

			for (int i = 0; i < ad_c_J.size(); ++i) {
				ad_a[i].resize(ndirs, SpatialVector::Zero());
			}

			for (int i = 0; i < ad_c_J.size(); ++i) {
				ad_c[i].resize(ndirs, SpatialVector::Zero());
			}

			for (int i = 0; i < ad_c_J.size(); ++i) {
				ad_f[i].resize(ndirs, SpatialVector::Zero());
			}
			
			for (int i = 0; i < ad_I_ci.size(); ++i) {
				ad_I_ci[i].resize(ndirs, SpatialMatrix::Zero());
			}

			for (int i = 0; i < ad_F.size(); ++i) {
				ad_F[i].resize(ndirs, SpatialVector::Zero());
			}
			
			for(int i = 0; i < ad_pA.size(); ++i) {
				ad_pA[i].resize(ndirs, SpatialVector::Zero());
			}
			
			for(int i = 0; i < ad_U.size(); ++i) {
				ad_U[i].resize(ndirs, SpatialVector::Zero());
			}
			
			for (int i = 0; i < ad_IA.size(); ++i) {
				ad_IA[i].resize(ndirs, SpatialMatrix::Zero());
			}
			
			// for(int i = 0; i < ad_u.size(); ++i) {
			// 	ad_u[i].resize(ndirs, 0.0);
			// }
			
			// for (int i = 0; i < ad_d.size(); ++i) {
			// 	ad_d[i].resize(ndirs, 0.0);
			// }

			ad_u.resize(ad_u.rows(), ndirs);
			ad_d.resize(ad_d.rows(), ndirs);
		}
	}
};

struct CartPendulum {
	CartPendulum () {
		using namespace RigidBodyDynamics;
		using namespace RigidBodyDynamics::Math;

		ClearLogOutput();

		double cart_w = 0.5;
		double cart_d = 0.2;
		double cart_h = 0.2;
		double cart_m = 100.0;
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

		id_cart = model.AddBody (0, Xtrans(Vector3d (0., 0., 0.)), joint_cart, body_cart, "cart");
		id_pendulum = model.AddBody (id_cart, Xtrans(Vector3d (0., 0., 0.)), joint_pendulum, body_pendulum, "pendulum");

		q = VectorNd::Constant ((size_t) model.dof_count, 0.);
		qdot = VectorNd::Constant ((size_t) model.dof_count, 0.);
		qddot = VectorNd::Constant ((size_t) model.dof_count, 0.);
		tau = VectorNd::Constant ((size_t) model.dof_count, 0.);

		ad_model = ADModel(model);

		body_point = Vector3d (0., 0., pend_l);

		ClearLogOutput();
	}

	RigidBodyDynamics::Model model;
	ADModel ad_model;

	unsigned int id_cart, id_pendulum;
	RigidBodyDynamics::Body body_cart, body_pendulum;
	RigidBodyDynamics::Joint joint_cart, joint_pendulum;

	RigidBodyDynamics::Math::VectorNd q;
	RigidBodyDynamics::Math::VectorNd qdot;
	RigidBodyDynamics::Math::VectorNd qddot;
	RigidBodyDynamics::Math::VectorNd tau;

	Vector3d body_point;
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
	ADModel &ad_model,
	unsigned int joint_id,
	const VectorNd &q,
	const MatrixNd &q_dirs,
	const VectorNd &qdot,
	const MatrixNd &qdot_dirs
) {
	unsigned int ndirs = q_dirs.cols();
	ad_model.resize_directions(ndirs);

	// check input dimensions
	if (q_dirs.cols() != qdot_dirs.cols()) {
		std::cerr << "directions have different dimensions: " << "#q_dirs = " << q_dirs.cols() << " != " << qdot_dirs.cols() << " = #qdot_dirs." << std::endl;
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
			ad_model.ad_X_J[joint_id][idir] = ad_Xroty (q[model.mJoints[joint_id].q_index], q_dirs(model.mJoints[joint_id].q_index, idir));
			ad_model.ad_S[joint_id][idir]  = SpatialVector::Zero(); // S = [0., 1., 0., 0., 0., 0.]
			ad_model.ad_v_J[joint_id][idir][1] = qdot_dirs(model.mJoints[joint_id].q_index, idir); // v_J = S*qdot
			ad_model.ad_c_J[joint_id][idir] = SpatialVector::Zero(); // c_J = Sdot*qdot
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
			ad_model.ad_X_J[joint_id][idir] = ad_Xtrans (Vector3d (q_dirs(model.mJoints[joint_id].q_index, idir), 0.0, 0.0));
			ad_model.ad_S[joint_id][idir]  = SpatialVector::Zero(); // S = [0., 0., 0., 1., 0., 0.]
			ad_model.ad_v_J[joint_id][idir][3] = qdot_dirs(model.mJoints[joint_id].q_index, idir); // v_J = S*qdot
			ad_model.ad_c_J[joint_id][idir] = SpatialVector::Zero(); // v_J = Sdot*qdot
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
		ad_model.ad_X_lambda[joint_id][idir] = ad_model.ad_X_J[joint_id][idir] * model.X_T[joint_id].toMatrix();
        //cout << "ad_jcalc:  ad_X_lambda[" << joint_id << "][" << idir << "] =" << endl << ad_X_lambda[joint_id][idir] << endl;
        //cout << "ad_jcalc:  ad_X_Ji[" << joint_id << "][" << idir << "] =" << endl << ad_X_Ji[idir] << endl;
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
		ADModel &ad_model,
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

	ad_model.resize_directions(ndirs);

	std::vector<MatrixNd> ad_X_J_i (ndirs, MatrixNd::Zero (6,6));
	std::vector<std::vector<MatrixNd> > fd_X_J (model.mBodies.size(), ad_X_J_i);
	// ad_X_J[3][5] gives for body 3 the 5th direction

	std::vector<MatrixNd> ad_X_lambda_i (ndirs, MatrixNd::Zero (6,6));
	std::vector<std::vector<MatrixNd> > fd_X_lambda (model.mBodies.size(), ad_X_lambda_i);
	// ad_X_lambda[3][5] gives for body 3 the 5th direction

	std::vector<MatrixNd> ad_X_base_i (ndirs, MatrixNd::Zero (6,6));
	std::vector<std::vector<MatrixNd> > fd_X_base (model.mBodies.size(), ad_X_base_i);
	// ad_X_base[3][5] gives for body 3 the 5th direction

	for (unsigned int i = 1; i < model.mBodies.size(); i++) {
		unsigned int lambda = model.lambda[i];
		// Calculate joint dependent variables
		if (model.mJoints[i].mJointType == JointTypeRevoluteY) {
			for (unsigned int j = 0; j < ndirs; j++) {
				ad_model.ad_X_J[i][j] = ad_Xroty (q[model.mJoints[i].q_index], q_dirs(i-1,j));
			}
			model.X_J[i] = Xroty (q[model.mJoints[i].q_index]);
		} else if (model.S[i] == SpatialVector (0., 0., 0., 1., 0., 0.)) {
			for (unsigned int j = 0; j < ndirs; j++) {
				ad_model.ad_X_J[i][j] = ad_Xtrans (Vector3d (q_dirs(i-1, j), 0., 0.));
			}
			model.X_J[i] = Xtrans (Vector3d (1., 0., 0.) * q[model.mJoints[i].q_index]);
		} else {
			std::cerr << "Unsupported joint! Only RotY and TransX supported!" << std::endl;
			abort();
		}

		for (unsigned int j = 0; j < ndirs; j++) {
			ad_model.ad_X_lambda[i][j] = ad_model.ad_X_J[i][j] * model.X_T[i].toMatrix();
		}
		model.X_lambda[i] = model.X_J[i] * model.X_T[i];

		for (unsigned int j = 0; j < ndirs; j++) {
			ad_model.ad_X_base[i][j] = ad_model.ad_X_lambda[i][j] * model.X_base[lambda].toMatrix() + model.X_lambda[i].toMatrix() * ad_model.ad_X_base[lambda][j];
		}
		model.X_base[i] = model.X_lambda[i] * model.X_base[lambda];
	}

	for (unsigned int j = 0; j < ndirs; j++) {
		SpatialMatrix X_base_ib = model.X_base[body_id].toMatrix();
		Matrix3d ad_E = E_from_Matrix(ad_model.ad_X_base[body_id][j]);
		Vector3d ad_r = ad_r_from_Matrix(X_base_ib, ad_model.ad_X_base[body_id][j]);

		out.block<3,1>(0,j) = ad_r + ad_E.transpose() * point_body_coordinates;
	}

	Matrix3d body_rotation = E_from_Matrix(model.X_base[body_id].toMatrix());
	Vector3d body_position = r_from_Matrix(model.X_base[body_id].toMatrix());

	return body_position + body_rotation.transpose() * point_body_coordinates;
}

RBDL_DLLAPI
void ad_UpdateKinematics (
	Model &model,
	ADModel &ad_model,
	const VectorNd &Q,
	const MatrixNd &q_dirs,
	const VectorNd &QDot,
	const MatrixNd &qdot_dirs,
	const VectorNd &QDDot,
	const MatrixNd &qddot_dirs
) {
	LOG << "-------- " << __func__ << " --------" << std::endl;

	unsigned int i;

	SpatialVector spatial_gravity (0., 0., 0., model.gravity[0], model.gravity[1], model.gravity[2]);

	model.a[0].setZero();
	//model.a[0] = spatial_gravity;

	for (i = 1; i < model.mBodies.size(); i++) {
		unsigned int q_index = model.mJoints[i].q_index;

		Joint joint = model.mJoints[i];
		unsigned int lambda = model.lambda[i];

		jcalc (model, i, Q, QDot);

		model.X_lambda[i] = model.X_J[i] * model.X_T[i];

		if (lambda != 0) {
			model.X_base[i] = model.X_lambda[i] * model.X_base[lambda];
			model.v[i] = model.X_lambda[i].apply(model.v[lambda]) + model.v_J[i];
		}	else {
			model.X_base[i] = model.X_lambda[i];
			model.v[i] = model.v_J[i];
		}

		model.c[i] = model.c_J[i] + crossm(model.v[i],model.v_J[i]);
		model.a[i] = model.X_lambda[i].apply(model.a[lambda]) + model.c[i];

		if (model.mJoints[i].mDoFCount == 3) {
			Vector3d omegadot_temp (QDDot[q_index], QDDot[q_index + 1], QDDot[q_index + 2]);
			model.a[i] = model.a[i] + model.multdof3_S[i] * omegadot_temp;
		} else {
			model.a[i] = model.a[i] + model.S[i] * QDDot[q_index];
		}
	}

	for (i = 1; i < model.mBodies.size(); i++) {
		LOG << "a[" << i << "] = " << model.a[i].transpose() << std::endl;
	}
}

Math::Vector3d ad_CalcPointAcceleration (
		Model &model,
		ADModel &ad_model,
		const Math::VectorNd &q,
		const Math::MatrixNd &q_dirs,
		const Math::VectorNd &qdot,
		const Math::MatrixNd &qdot_dirs,
		const Math::VectorNd &qddot,
		const Math::MatrixNd &qddot_dirs,
		unsigned int body_id,
		const Math::Vector3d &point_position,
		const Math::MatrixNd &fd_derivative,
		bool update_kinematics = true
) {
	LOG << "-------- " << __func__ << " --------" << std::endl;

	unsigned int ndirs = q_dirs.cols();
	ad_model.resize_directions(ndirs);

	// Reset the velocity of the root body
	model.v[0].setZero();
	model.a[0].setZero();

	if (update_kinematics) {
		ad_UpdateKinematics (model, ad_model, q, q_dirs, qdot, qdot_dirs, qddot, qddot_dirs);
		UpdateKinematics (model, q, qdot, qddot);
	}

	LOG << std::endl;

	unsigned int reference_body_id = body_id;
	Vector3d reference_point = point_position;

	if (model.IsFixedBodyId(body_id)) {
		unsigned int fbody_id = body_id - model.fixed_body_discriminator;
		reference_body_id = model.mFixedBodies[fbody_id].mMovableParent;
		Vector3d base_coords = CalcBodyToBaseCoordinates (model, q, body_id, point_position, false);
		reference_point = CalcBaseToBodyCoordinates (model, q, reference_body_id, base_coords, false);
	}

	SpatialTransform p_X_i (CalcBodyWorldOrientation (model, q, reference_body_id, false).transpose(), reference_point);

	SpatialVector p_v_i = p_X_i.apply(model.v[reference_body_id]);
	Vector3d a_dash = Vector3d (p_v_i[0], p_v_i[1], p_v_i[2]).cross(Vector3d (p_v_i[3], p_v_i[4], p_v_i[5]));
	SpatialVector p_a_i = p_X_i.apply(model.a[reference_body_id]);

	return Vector3d (
			p_a_i[3] + a_dash[0],
			p_a_i[4] + a_dash[1],
			p_a_i[5] + a_dash[2]
			);
};

RBDL_DLLAPI
Vector3d fd_jcalc(
	Model &model,
	unsigned int joint_id,
	const VectorNd &q,
	const MatrixNd &q_dirs,
	const VectorNd &qdot,
	const MatrixNd &qdot_dirs,
	std::vector<SpatialMatrix> &ad_X_Ji,
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
Vector3d fd_CalcBodyToBaseCoordinatesSingleFunc (
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
//==============================================================================
void fd_InverseDynamics(Model& model,
						const Math::VectorNd& q,
						const Math::MatrixNd& q_dirs,
						const Math::VectorNd& qdot,
						const Math::MatrixNd& qdot_dirs,
						const Math::VectorNd& qddot,
						const Math::MatrixNd& qddot_dirs,
						Math::VectorNd& tau,
						Math::MatrixNd& fd_tau,
						std::vector<Math::SpatialVector> *f_ext){

	assert(q_dirs.cols() == qdot_dirs.cols() && 
           q_dirs.cols() == qddot_dirs.cols() && 
		   "q_dirs, qdot_dirs, qddot_dirs have different dimensions");

	assert( fd_tau.cols() == q_dirs.cols() && 
			fd_tau.rows() == tau.rows() && 
			"fd_tau and tau have different dimensions");

	double h = sqrt(1e-16);

	unsigned int ndirs   = q_dirs.cols();

	InverseDynamics(model, q, qdot, qddot, tau, f_ext);

	VectorNd hd_tau(tau);
	VectorNd q_dir(q);
	VectorNd qdot_dir(qdot);
	VectorNd qddot_dir(qddot);

	for (int i = 0; i < ndirs; i++ ){
			q_dir = q_dirs.block(0,i, model.q_size, 1);
			qdot_dir = qdot_dirs.block(0,i, model.q_size, 1);
			qddot_dir = qddot_dirs.block(0,i, model.q_size, 1);
			InverseDynamics(model,         q + h*q_dir, 
						                qdot + h*qdot_dir, 
						               qddot + h*qddot_dir, 
						                         hd_tau,  f_ext);

		fd_tau.block(0,i,hd_tau.rows(),1) = (hd_tau - tau) / h;
	}
}
						
//==============================================================================
void ad_InverseDynamics(Model& model,
			ADModel &ad_model,
			const Math::VectorNd& q,
			const Math::MatrixNd& q_dirs,
			const Math::VectorNd& qdot,
			const Math::MatrixNd& qdot_dirs,
			const Math::VectorNd& qddot,
			const Math::MatrixNd& qddot_dirs,
			Math::VectorNd& tau,
			Math::MatrixNd& ad_tau,
			std::vector<Math::SpatialVector> *f_ext){

	model.v[0].setZero();
	model.a[0].set (0., 0., 0., 
			-model.gravity[0], 
			-model.gravity[1], 
			-model.gravity[2]);

	unsigned int ndirs = q_dirs.cols();

	for (unsigned int i = 1; i < model.mBodies.size(); i++) {
		unsigned int lambda = model.lambda[i];
		unsigned int q_index = model.mJoints[i].q_index;

		ad_jcalc (model,
			  ad_model,
			  i,
			  q,
			  q_dirs,
			  qdot,
			  qdot_dirs);

// 		Done in ad_jcalc
// 		for(unsigned int j = 0; j < ndirs; j++) {
// 			ad_model.ad_X_lambda[i][j] = ad_model.ad_X_J[i][j] * model.X_T[i].toMatrix();
// 		}

		if (lambda != 0) {
			for(unsigned int j = 0; j < ndirs; j++) {
				ad_model.ad_X_base[i][j] = ad_model.ad_X_lambda[i][j] * model.X_base[lambda].toMatrix() 
					+ model.X_lambda[i].toMatrix() * ad_model.ad_X_base[i][j];
			}
			model.X_base[i] = model.X_lambda[i] * model.X_base[lambda];
		} else {
			for(unsigned int j = 0; j < ndirs; j++) {
				ad_model.ad_X_base[i][j] = ad_model.ad_X_lambda[i][j];
			}
			model.X_base[i] = model.X_lambda[i];
		}

		for(unsigned int j = 0; j < ndirs; j++) {
			ad_model.ad_v[i][j] = ad_model.ad_X_lambda[i][j] * model.v[lambda] 
				+ model.X_lambda[i].apply (ad_model.ad_v[lambda][j]) 
				+ ad_model.ad_v_J[i][j];
		}
		model.v[i] = model.X_lambda[i].apply(model.v[lambda]) + model.v_J[i];

		for(unsigned int j = 0; j < ndirs; j++) {
			ad_model.ad_c[i][j] = ad_model.ad_c_J[i][j]
				+ crossm(ad_model.ad_v[i][j],model.v_J[i]) 
				+ crossm(model.v[i],ad_model.ad_v_J[i][j]);
		}
		model.c[i] = model.c_J[i] + crossm(model.v[i],model.v_J[i]);

		if (model.mJoints[i].mDoFCount == 3) {
			cerr << "Multi-dof not supported." << endl;
			abort();
			model.a[i] = model.X_lambda[i].apply(model.a[lambda]) 
				+ model.c[i] 
				+ model.multdof3_S[i] * Vector3d (qddot[q_index], qddot[q_index + 1], qddot[q_index + 2]);
		} else {
			for(unsigned int j = 0; j < ndirs; j++) {
				ad_model.ad_a[i][j] = ad_model.ad_X_lambda[i][j] * model.a[lambda] + model.X_lambda[i].apply(ad_model.ad_a[lambda][j]) 
					+ ad_model.ad_c[i][j] 
					+ ad_model.ad_S[i][j] * qddot[q_index] + model.S[i] * qddot_dirs(q_index,j);
			}
			model.a[i] = model.X_lambda[i].apply(model.a[lambda]) + model.c[i] + model.S[i] * qddot[q_index];
		}	

		if (!model.mBodies[i].mIsVirtual) {
			for(unsigned int j = 0; j < ndirs; j++) {
				ad_model.ad_f[i][j] = model.I[i] * ad_model.ad_a[i][j] 
					+ crossf(ad_model.ad_v[i][j],model.I[i] * model.v[i]) 
					+ crossf(model.v[i],model.I[i] * ad_model.ad_v[i][j]);
			}
			model.f[i] = model.I[i] * model.a[i] + crossf(model.v[i],model.I[i] * model.v[i]);
		} else {
			for(unsigned int j = 0; j < ndirs; j++) {
				ad_model.ad_f[i][j].setZero();
			}
			model.f[i].setZero();
		}

		if (f_ext != NULL && (*f_ext)[i] != SpatialVectorZero) {
			for(unsigned int j = 0; j < ndirs; j++) {
				SpatialMatrix ad_X_base_force(ad_model.ad_X_base[i][j]);
				ad_X_base_force.block<3,3>(3,0) = Matrix3d::Zero();
				ad_X_base_force.block<3,3>(0,3) = ad_model.ad_X_base[i][j].block<3,3>(3,0);
				ad_model.ad_f[i][j] -= ad_X_base_force * (*f_ext)[i];
			}
			model.f[i] -= model.X_base[i].toMatrixAdjoint() * (*f_ext)[i];
		}
	}

	for (unsigned int i = model.mBodies.size() - 1; i > 0; i--) {
		if (model.mJoints[i].mDoFCount == 3) {
			cerr << "Multi-dof not supported." << endl;
			abort();
			tau.block<3,1>(model.mJoints[i].q_index, 0) = model.multdof3_S[i].transpose() * model.f[i];
		} else {
			for(unsigned int j = 0; j < ndirs; j++) {
				ad_tau(model.mJoints[i].q_index,j) = ad_model.ad_S[i][j].dot(model.f[i]) 
					+ model.S[i].dot(ad_model.ad_f[i][j]);
			}
			tau[model.mJoints[i].q_index] = model.S[i].dot(model.f[i]);
		}

		if (model.lambda[i] != 0) {
			for(unsigned int j = 0; j < ndirs; j++) {
				ad_model.ad_f[model.lambda[i]][j] = ad_model.ad_f[model.lambda[i]][j]
					+ ad_model.ad_X_lambda[i][j].transpose() * model.f[i] 
					+ model.X_lambda[i].applyTranspose(ad_model.ad_f[i][j]);
			}
			model.f[model.lambda[i]] = model.f[model.lambda[i]] 
				+ model.X_lambda[i].applyTranspose(model.f[i]);
		}
	}
}

void fd_CompositeRigidBodyAlgorithm (
					Model &model,
					const VectorNd &q,
					const MatrixNd &q_dirs,
					std::vector<MatrixNd> &fd_out
					) {

	MatrixNd inertia_ref = MatrixNd::Zero(q.size(),q.size());
	MatrixNd inertia_fd = MatrixNd::Zero(q.size(),q.size());

	CompositeRigidBodyAlgorithm(model, q, inertia_ref, true);
	unsigned int ndirs = q_dirs.cols();
	double h = 1.0e-8;
	for (unsigned int j = 0; j < ndirs; j++) {
		VectorNd q_dir = q_dirs.block(0,j, model.qdot_size, 1);
		CompositeRigidBodyAlgorithm(model, q+h*q_dir, inertia_fd, true);
		fd_out.push_back( (inertia_fd - inertia_ref) / h);
	}
}


void ad_CompositeRigidBodyAlgorithm (
						 Model &model,
						 ADModel &ad_model,
						 const VectorNd &q,
						 const MatrixNd &q_dirs,
						 MatrixNd &H,
						 std::vector<MatrixNd> &out,
						 bool update_kinematics = true
						 ) {
	
// 	double h = 1.0e-8;
	size_t ndirs = q_dirs.cols();
	
	ad_model.resize_directions(ndirs);
	
	for(unsigned i = 1; i < model.mBodies.size(); i++)
		ad_jcalc(model,ad_model,i,q,q_dirs,VectorNd::Zero(model.q_size),MatrixNd::Zero(model.q_size,ndirs));
	
	// --- derivative quantities ----------------------------------------------------------------------
	
	assert (H.rows() == model.dof_count && H.cols() == model.dof_count);

	for (unsigned int i = 1; i < model.mBodies.size(); i++) {
		if (update_kinematics) {
			jcalc_X_lambda_S (model, i, q);
			// this works for single dof joint
			cerr << "No ad_jcalc_X_lambda_S implemented yet" << endl;
		}
		// ad code
		for (size_t j = 0; j < ndirs; j++) {
			ad_model.ad_I_ci[i][j] = SpatialMatrix::Zero();
//			 fd_I_ci[i][j] = SpatialMatrix::Zero();
		}
		// normal quantity code
		model.Ic[i] = model.I[i];
	}
	
	for (unsigned int i = model.mBodies.size() - 1; i > 0; i--) {
		//if (model.lambda[i] != 0) {
			// ad code
			for(size_t j = 0; j < ndirs; j++) {
				ad_model.ad_I_ci[model.lambda[i]][j] = ad_model.ad_I_ci[model.lambda[i]][j] + model.X_lambda[i].toMatrixTranspose()*ad_model.ad_I_ci[i][j]*model.X_lambda[i].toMatrix() + ad_model.ad_X_lambda[i][j].transpose()*model.Ic[i].toMatrix()*model.X_lambda[i].toMatrix() + model.X_lambda[i].toMatrixTranspose()*model.Ic[i].toMatrix()*ad_model.ad_X_lambda[i][j];
// 	fd_I_ci[model.lambda[i]][j] = fd_I_ci[model.lambda[i]][j] + model.X_lambda[i].toMatrixTranspose()*fd_I_ci[i][j]*model.X_lambda[i].toMatrix() + (fd_per_X_lambda[i][j].transpose()*model.Ic[i].toMatrix()*fd_per_X_lambda[i][j] - model.X_lambda[i].toMatrixTranspose()*model.Ic[i].toMatrix()*model.X_lambda[i].toMatrix())/h;

// 	cout << "fd - ad: " << endl << fd_I_ci[model.lambda[i]][j] - ad_I_ci[model.lambda[i]][j] << endl;
			}
			// normal quantity code
			model.Ic[model.lambda[i]] = model.Ic[model.lambda[i]] + model.X_lambda[i].applyTranspose(model.Ic[i]);
			
			
		//}

		unsigned int dof_index_i = model.mJoints[i].q_index;
		assert( model.mJoints[i].mDoFCount == 1);
		// if (model.mJoints[i].mDoFCount == 3) {
		// 	Matrix63 F_63 = model.Ic[i].toMatrix() * model.multdof3_S[i];
		// 	H.block<3,3>(dof_index_i, dof_index_i) = model.multdof3_S[i].toMatrixTranspose() * F_63;

		// 	unsigned int j = i;
		// 	unsigned int dof_index_j = dof_index_i;

		// 	while (model.lambda[j] != 0) {
		// 		F_63 = model.X_lambda[j].toMatrixTranspose() * (F_63);
		// 		j = model.lambda[j];
		// 		dof_index_j = model.mJoints[j].q_index;

		// 		if (model.mJoints[j].mDoFCount == 3) {
		// 			Matrix3d H_temp2 = F_63.transpose() * (model.multdof3_S[j]);

		// 			H.block<3,3>(dof_index_i,dof_index_j) = H_temp2;
		// 			H.block<3,3>(dof_index_j,dof_index_i) = H_temp2.transpose();
		// 		} else {
		// 			Vector3d H_temp2 = F_63.transpose() * (model.S[j]);

		// 			H.block<3,1>(dof_index_i,dof_index_j) = H_temp2;
		// 			H.block<1,3>(dof_index_j,dof_index_i) = H_temp2.transpose();
		// 		}
		// 	}
		// } else {
		
		// ad code
		for (size_t j = 0; j < ndirs; j++) {
			ad_model.ad_F[i][j] = ad_model.ad_I_ci[i][j]*model.S[i]+model.Ic[i].toMatrix()*ad_model.ad_S[i][j];
		}
		// normal quantity code
		SpatialVector F = model.Ic[i] * model.S[i];
		
		// ad code
		for (size_t j = 0; j < ndirs; j++) {
			out[j](dof_index_i, dof_index_i) = ad_model.ad_S[i][j].dot(F)+model.S[i].dot(ad_model.ad_F[i][j]);
		}
		// normal quantity code
		H(dof_index_i, dof_index_i) = model.S[i].dot(F);
		unsigned int j = i;
		unsigned int dof_index_j = dof_index_i;

		while (model.lambda[j] != 0) {
			// ad code
			for (size_t ndir = 0; ndir < ndirs; ndir++) {
				ad_model.ad_F[i][ndir] = ad_model.ad_X_lambda[j][ndir].transpose()*F+model.X_lambda[j].toMatrixTranspose()*ad_model.ad_F[i][ndir];
			}
			// normal quantity code
			F = model.X_lambda[j].applyTranspose(F);
			j = model.lambda[j];
			dof_index_j = model.mJoints[j].q_index;
			assert(model.mJoints[j].mDoFCount == 1);
			// if (model.mJoints[j].mDoFCount == 3) {
			//	 Vector3d H_temp2 = (F.transpose() * model.multdof3_S[j]).transpose();
			//	 H.block<1,3>(dof_index_i,dof_index_j) = H_temp2.transpose();
			//	 H.block<3,1>(dof_index_j,dof_index_i) = H_temp2;
			// } else {
			// normal quantity code
			H(dof_index_i,dof_index_j) = F.dot(model.S[j]);
			H(dof_index_j,dof_index_i) = H(dof_index_i,dof_index_j);
			// ad code
			for (size_t ndir = 0; ndir < ndirs; ndir++) {
				out[ndir](dof_index_i,dof_index_j) = ad_model.ad_F[i][ndir].dot(model.S[j])+F.dot(ad_model.ad_S[j][ndir]);
				out[ndir](dof_index_j,dof_index_i) = out[ndir](dof_index_i,dof_index_j);
			}
		
			//}
		}
		// }
	}
}

/*
void fd_UpdateKinematics (
		Model &model,
		const Math::VectorNd &q,
		const Math::VectorNd &q_dirs,
		const Math::VectorNd &qdot,
		const Math::VectorNd &qdot_dirs,
		const Math::VectorNd &qddot,
		const Math::VectorNd &qddot_dirs,
		unsigned int body_id,
		const Math::Vector3d &point_position,
		const Math::MatrixNd &fd_derivative,
		bool update_kinematics = true
) {
	// evaluate y(t)
	Vector3d ref = CalcPointAcceleration (model, q, qdot, qddot, body_id, point_position);

	unsigned int ndirs = q_dirs.cols();
	double h = 1.0e-8;

	for (unsigned int j = 0; j < ndirs; j++) {
		VectorNd q_dir = q_dirs.block(0, j, model.q_size, 1);
		VectorNd qdot_dir = qdot_dirs.block(0, j, model.q_size, 1);
		VectorNd qddot_dir = qddot_dirs.block(0, j, model.q_size, 1);

		// evaluate y(t+h*d)
		Vector3d res_hd = CalcPointAcceleration (
			model,
			q + h * q_dir,
			qdot + h * qdot_dir,
			qddot + h * qddot_dir,
			body_id, point_position
		);

		//fd_derivative.block<3,1>(0, j) = (res_hd - ref) / h;
	}

	return ref;
};
*/





/*This functin shall compute the Forward Dynamics by solving M qddot = tau - N(q,qdot)*/
void ForwardDynamicsCholesky (
			      Model &model,
			      const VectorNd &q,
			      const VectorNd &qdot,
			      const VectorNd &tau,
			      VectorNd &qddot,
			      std::vector<SpatialVector>* f_ext
			      ){

  VectorNd zero_vector (VectorNd::Zero (tau.size()));

  //Here we geht the zero_vector by calling Inverse dynamics with qddot = 0
  InverseDynamics(model,q,qdot,zero_vector,qddot,f_ext);
  // cout << "qddot: " << qddot << endl;
  // cout << "zero_vector: " << zero_vector << endl;
  // cout << "tau: " << tau << endl;

  MatrixNd M ( MatrixNd::Zero(qddot.size(),qddot.size()));
  CompositeRigidBodyAlgorithm (model, q, M);  //bool update kinematics?? Works without it
  //  cout << "Mass Matrix: " << M << endl;

  SparseFactorizeLTL(model,M);
  
  qddot = tau - qddot;
  SparseSolveLTx(model,M,qddot);
  SparseSolveLx(model,M,qddot);

};

void ad_ForwardDynamicsCholesky (
				 Model& model,
				 ADModel& ad_model,
				 const VectorNd& q,
				 const MatrixNd& q_dirs,
				 const VectorNd& qdot,
				 const MatrixNd& qdot_dirs,
				 const VectorNd& tau,
				 const MatrixNd& tau_dirs,
				 VectorNd& qddot,
				 MatrixNd& ad_qddot,
				 std::vector<SpatialVector>* f_ext) {

  unsigned int ndirs = q_dirs.cols();
  ad_model.resize_directions(ndirs);

	
  VectorNd zero_vector (VectorNd::Zero (q.size()));
  MatrixNd zero_vector_dirs (MatrixNd::Zero(q.size(),ndirs));

  /*Here we get the coriolis term by calling Inverse dynamics with accelerations = 0.
    The result is stored in qddot. */

  ad_InverseDynamics(model,
		     ad_model,
		     q,
		     q_dirs,
		     qdot,
		     qdot_dirs,
		     zero_vector,
		     zero_vector_dirs,
		     qddot,
		     ad_qddot,
		     f_ext);    

  MatrixNd M (MatrixNd::Zero(qddot.size(),qddot.size()));
  std::vector<MatrixNd> ad_M (q_dirs.cols(),MatrixNd::Zero(model.q_size,model.q_size));

  ad_CompositeRigidBodyAlgorithm(model,
				 ad_model,
				 q,
				 q_dirs,
				 M,
				 ad_M);
	
  //	CompositeRigidBodyAlgorithm (model, q, M);  //bool update kinematics?? Works without it

  //  cout << "Mass Matrix: " << M << endl;

  SparseFactorizeLTL(model,M);
  
  qddot = tau - qddot;
  SparseSolveLTx(model,M,qddot);
  SparseSolveLx(model,M,qddot);



  // I would like to do that without the local copy, maybe Martin knows, how to do that and whether it's a good idea.
  VectorNd local_column_copy;
  
  for (unsigned int j = 0; j < ndirs; j++)
    {
      ad_qddot.col(j)=tau_dirs.col(j) - ad_qddot.col(j)-ad_M[j]*qddot;
      local_column_copy=ad_qddot.col(j);
      SparseSolveLTx(model,M,local_column_copy);
      SparseSolveLx(model,M,local_column_copy);
      ad_qddot.col(j)=local_column_copy;
    }
	

};



Math::Vector3d fd_CalcPointAcceleration (
		Model &model,
		const Math::VectorNd &q,
		const Math::VectorNd &q_dirs,
		const Math::VectorNd &qdot,
		const Math::VectorNd &qdot_dirs,
		const Math::VectorNd &qddot,
		const Math::VectorNd &qddot_dirs,
		unsigned int body_id,
		const Math::Vector3d &point_position,
		const Math::MatrixNd &fd_derivative,
		bool update_kinematics = true
) {
	// evaluate y(t)
	Vector3d ref = CalcPointAcceleration (model, q, qdot, qddot, body_id, point_position);

	unsigned int ndirs = q_dirs.cols();
	double h = 1.0e-8;

	for (unsigned int j = 0; j < ndirs; j++) {
		VectorNd q_dir = q_dirs.block(0, j, model.q_size, 1);
		VectorNd qdot_dir = qdot_dirs.block(0, j, model.q_size, 1);
		VectorNd qddot_dir = qddot_dirs.block(0, j, model.q_size, 1);

		// evaluate y(t+h*d)
		Vector3d res_hd = CalcPointAcceleration (
			model,
			q + h * q_dir,
			qdot + h * qdot_dir,
			qddot + h * qddot_dir,
			body_id, point_position
		);

		//fd_derivative.block<3,1>(0, j) = (res_hd - ref) / h;
	}

	return ref;
};

void fd_ForwardDynamics(
	Model& model,
	const VectorNd& q,
	const MatrixNd& q_dirs,
	const VectorNd& qdot,
	const MatrixNd& qdot_dirs,
	const VectorNd& tau,
	const MatrixNd& tau_dirs,
	VectorNd& qddot,
	MatrixNd& fd_qddot,
	std::vector<SpatialVector>* f_ext) {
	assert(q_dirs.cols() == qdot_dirs.cols() 
		&& q_dirs.cols() == fd_qddot.cols() 
		&& tau_dirs.cols() == q_dirs.cols()
		&& "q_dirs, qdot_dirs, tau_dirs, fd_qddot have different dimensions");
	
	double h = sqrt(1e-16);
	unsigned int ndirs = q_dirs.cols();
	
	ForwardDynamics(model, q, qdot, tau, qddot,f_ext);
	
	VectorNd hd_qddot(qddot);
	VectorNd q_dir(q);
	VectorNd qdot_dir(qdot);
	VectorNd tau_dir(tau);
	
	for(unsigned int i = 0; i < ndirs; ++i) {
	  q_dir = q_dirs.col(i);
	  qdot_dir = qdot_dirs.col(i);
	  tau_dir = tau_dirs.col(i);
	  ForwardDynamics(model, q + h*q_dir, qdot + h*qdot_dir, tau + h*tau_dir, hd_qddot,f_ext);
		
	  fd_qddot.col(i) = (hd_qddot - qddot) / h;
	}
}

void ad_ForwardDynamics (
	Model& model,
	ADModel& ad_model,
	const VectorNd& q,
	const MatrixNd& q_dirs,
	const VectorNd& qdot,
	const MatrixNd& qdot_dirs,
	const VectorNd& tau,
	const MatrixNd& tau_dirs,
	VectorNd& qddot,
	MatrixNd& ad_qddot,
	std::vector<SpatialVector>* f_ext) {
	SpatialVector spatial_gravity (0., 0., 0., model.gravity[0], model.gravity[1], model.gravity[2]);

	unsigned int ndirs = q_dirs.cols();
	ad_model.resize_directions(ndirs);
	
	unsigned int i = 0;

	// Reset the velocity of the root body
	model.v[0].setZero();

	for (i = 1; i < model.mBodies.size(); i++) {
		unsigned int lambda = model.lambda[i];
		
		ad_jcalc(model, ad_model, i, q, q_dirs, qdot, qdot_dirs);
		
		if (lambda != 0) {
			for(unsigned int j = 0; j < ndirs; j++) {
				ad_model.ad_X_base[i][j] = ad_model.ad_X_lambda[i][j] * model.X_base[lambda].toMatrix() 
					+ model.X_lambda[i].toMatrix() * ad_model.ad_X_base[i][j];
			}
			model.X_base[i] = model.X_lambda[i] * model.X_base[lambda];
		} else {
			for(unsigned int j = 0; j < ndirs; j++) {
				ad_model.ad_X_base[i][j] = ad_model.ad_X_lambda[i][j];
			}
			model.X_base[i] = model.X_lambda[i];
		}

		for(unsigned int j = 0; j < ndirs; j++) {
			ad_model.ad_v[i][j] = ad_model.ad_X_lambda[i][j] * model.v[lambda] 
				+ model.X_lambda[i].apply(ad_model.ad_v[lambda][j]) 
				+ ad_model.ad_v_J[i][j];
		}
		model.v[i] = model.X_lambda[i].apply(model.v[lambda]) + model.v_J[i];

		for(unsigned int j = 0; j < ndirs; j++) {
			ad_model.ad_c[i][j] = ad_model.ad_c_J[i][j]
				+ crossm(ad_model.ad_v[i][j],model.v_J[i]) 
				+ crossm(model.v[i],ad_model.ad_v_J[i][j]);
		}
		model.c[i] = model.c_J[i] + crossm(model.v[i],model.v_J[i]);
		
		for(unsigned int j = 0; j < ndirs; j++) {
				ad_model.ad_IA[i][j].setZero();
		}
		model.I[i].setSpatialMatrix(model.IA[i]);

		for(unsigned int j = 0; j < ndirs; j++) {
			ad_model.ad_pA[i][j] = crossf(ad_model.ad_v[i][j],model.I[i] * model.v[i])
				+ crossf(model.v[i],model.I[i] * ad_model.ad_v[i][j]);;
		}
		model.pA[i] = crossf(model.v[i],model.I[i] * model.v[i]);

		if (f_ext != NULL && (*f_ext)[i] != SpatialVectorZero) {
			for(unsigned int j = 0; j < ndirs; j++) {
				SpatialMatrix ad_X_base_force(ad_model.ad_X_base[i][j]);
				ad_X_base_force.block<3,3>(3,0) = Matrix3d::Zero();
				ad_X_base_force.block<3,3>(0,3) = ad_model.ad_X_base[i][j].block<3,3>(3,0);
				ad_model.ad_pA[i][j] -= ad_X_base_force * (*f_ext)[i];
			}
			model.pA[i] -= model.X_base[i].toMatrixAdjoint() * (*f_ext)[i];
		}
	}

	for (i = model.mBodies.size() - 1; i > 0; i--) {
		unsigned int q_index = model.mJoints[i].q_index;

		if (model.mJoints[i].mDoFCount == 3) {
			cerr << "Multi-dof not supported." << endl;
			abort();
// 			model.multdof3_U[i] = model.IA[i] * model.multdof3_S[i];
// #ifdef EIGEN_CORE_H
// 			model.multdof3_Dinv[i] = (model.multdof3_S[i].transpose() * model.multdof3_U[i]).inverse().eval();
// #else
// 			model.multdof3_Dinv[i] = (model.multdof3_S[i].transpose() * model.multdof3_U[i]).inverse();
// #endif
// 			Vector3d tau_temp (Tau[q_index], Tau[q_index + 1], Tau[q_index + 2]);
// 
// 			model.multdof3_u[i] = tau_temp - model.multdof3_S[i].transpose() * model.pA[i];
// 
// //			LOG << "multdof3_u[" << i << "] = " << model.multdof3_u[i].transpose() << std::endl;
// 			unsigned int lambda = model.lambda[i];
// 			if (lambda != 0) {
// 				SpatialMatrix Ia = model.IA[i] - model.multdof3_U[i] * model.multdof3_Dinv[i] * model.multdof3_U[i].transpose();
// 				SpatialVector pa = model.pA[i] + Ia * model.c[i] + model.multdof3_U[i] * model.multdof3_Dinv[i] * model.multdof3_u[i];
// #ifdef EIGEN_CORE_H
// 				model.IA[lambda].noalias() += model.X_lambda[i].toMatrixTranspose() * Ia * model.X_lambda[i].toMatrix();
// 				model.pA[lambda].noalias() += model.X_lambda[i].applyTranspose(pa);
// #else
// 				model.IA[lambda] += model.X_lambda[i].toMatrixTranspose() * Ia * model.X_lambda[i].tadoMatrix();
// 				model.pA[lambda] += model.X_lambda[i].applyTranspose(pa);
// #endif
// 				LOG << "pA[" << lambda << "] = " << model.pA[lambda].transpose() << std::endl;
// 			}
		} 
		else {
			for(unsigned int j = 0; j < ndirs; j++) {
				ad_model.ad_U[i][j] = model.IA[i] * ad_model.ad_S[i][j];
				//	+ ad_model.ad_IA[i][j] * model.S[i];
			}
			model.U[i] = model.IA[i] * model.S[i];
			
			for(unsigned int j = 0; j < ndirs; j++) {
				ad_model.ad_d(i,j) = ad_model.ad_S[i][j].dot(model.U[i]) + model.S[i].dot(ad_model.ad_U[i][j]);
			}
			model.d[i] = model.S[i].dot(model.U[i]);
			
			for(unsigned int j = 0; j < ndirs; j++) {
				ad_model.ad_u(i,j) = tau_dirs(q_index,j) 
					- ad_model.ad_S[i][j].dot(model.pA[i]) 
					- model.S[i].dot(ad_model.ad_pA[i][j]);
			}
			model.u[i] = tau[q_index] - model.S[i].dot(model.pA[i]);

			unsigned int lambda = model.lambda[i];
			if (lambda != 0) {
				vector<SpatialMatrix> ad_Ia(ndirs, SpatialMatrix::Zero());
				for(unsigned int j = 0; j < ndirs; j++) {
					ad_Ia[j] = ad_model.ad_IA[i][j] 
						- ad_model.ad_U[i][j] * (model.U[i] / model.d[i]).transpose()
						- model.U[i] * (ad_model.ad_U[i][j] / model.d[i]).transpose()
						- model.U[i] * (model.U[i] * ad_model.ad_d(i,j) / model.d[i]).transpose();
				}
				SpatialMatrix Ia = model.IA[i] - model.U[i] * (model.U[i] / model.d[i]).transpose();
				
				vector<SpatialVector> ad_pa(ndirs, SpatialVector::Zero());
				for(unsigned int j = 0; j < ndirs; j++) {
					ad_pa[j] = ad_model.ad_pA[i][j]
						+ ad_Ia[j] * model.c[i]
						+ Ia * ad_model.ad_c[i][j] 
						+ ad_model.ad_U[i][j] * model.u[i] / model.d[i]
						+ model.U[i] * ad_model.ad_u(i,j) / model.d[i]
						+ model.U[i] * model.u[i] * ad_model.ad_d(i,j) / model.d[i];
				}
				SpatialVector pa = model.pA[i] + Ia * model.c[i] + model.U[i] * model.u[i] / model.d[i];
#ifdef EIGEN_CORE_H
				for(unsigned int j = 0; j < ndirs; j++) {
					ad_model.ad_IA[lambda][j].noalias() += 
						ad_model.ad_X_lambda[i][j].transpose() * Ia * model.X_lambda[i].toMatrix()
						+ model.X_lambda[i].toMatrixTranspose() * ad_Ia[j] * model.X_lambda[i].toMatrix()
						+ model.X_lambda[i].toMatrixTranspose() * Ia * ad_model.ad_X_lambda[i][j];
				}
				model.IA[lambda].noalias() += model.X_lambda[i].toMatrixTranspose() * Ia * model.X_lambda[i].toMatrix();
				
				for(unsigned int j = 0; j < ndirs; j++) {
					ad_model.ad_pA[lambda][j].noalias() += ad_model.ad_X_lambda[i][j].transpose() * pa
						+ model.X_lambda[i].applyTranspose(ad_pa[j]);
				}
				model.pA[lambda].noalias() += model.X_lambda[i].applyTranspose(pa);
#else
				cerr << "Simple math not yet supported." << endl;
				abort();
				//model.IA[lambda] += model.X_lambda[i].toMatrixTranspose() * Ia * model.X_lambda[i].toMatrix();
				//model.pA[lambda] += model.X_lambda[i].applyTranspose(pa);
#endif
			}
		}
	}

	model.a[0] = spatial_gravity * -1.;

	for (i = 1; i < model.mBodies.size(); i++) {
		unsigned int q_index = model.mJoints[i].q_index;
		unsigned int lambda = model.lambda[i];
		SpatialTransform X_lambda = model.X_lambda[i];

		for(unsigned int j = 0; j < ndirs; j++) {
			ad_model.ad_a[i][j] = ad_model.ad_X_lambda[i][j] * model.a[lambda]
				+ X_lambda.apply(ad_model.ad_a[lambda][j]) 
				+ ad_model.ad_c[i][j];
		}
		model.a[i] = X_lambda.apply(model.a[lambda]) + model.c[i];

		if (model.mJoints[i].mDoFCount == 3) {
			cerr << "Multi-dof joints not supported." << endl;
			abort();
// 			Vector3d qdd_temp = model.multdof3_Dinv[i] * (model.multdof3_u[i] - model.multdof3_U[i].transpose() * model.a[i]);
// 			qddot[q_index] = qdd_temp[0];
// 			qddot[q_index + 1] = qdd_temp[1];
// 			qddot[q_index + 2] = qdd_temp[2];
// 			model.a[i] = model.a[i] + model.multdof3_S[i] * qdd_temp;
		} 
		else {
			for(unsigned int j = 0; j < ndirs; j++) {
				ad_qddot(q_index,j) = (ad_model.ad_d(i,j)/model.d[i]) * (model.u[i] - model.U[i].dot(model.a[i]))
					+ (1./model.d[i]) * (ad_model.ad_u(i,j) - ad_model.ad_U[i][j].dot(model.a[i]) - model.U[i].dot(ad_model.ad_a[i][j]));
			}
			qddot[q_index] = (1./model.d[i]) * (model.u[i] - model.U[i].dot(model.a[i]));
			
			for(unsigned int j = 0; j < ndirs; j++) {
				ad_model.ad_a[i][j] = ad_model.ad_a[i][j]
					+ ad_model.ad_S[i][j] * qddot[q_index]
					+ model.S[i] * ad_qddot(q_index,j);
			}
			model.a[i] = model.a[i] + model.S[i] * qddot[q_index];
		}
	}
}

TEST_FIXTURE (CartPendulum, jcalcNominalSolutionTest) {
	// set nominal values
	q.setZero();
	q[0] = 0.3;
	q[1] = -0.2;

	qdot.setZero();
	qdot[0] = 0.3;
	qdot[1] = -0.2;

	// set directions
	unsigned int ndirs = model.q_size + model.qdot_size;
	MatrixNd x = MatrixNd::Identity(ndirs, ndirs);
	MatrixNd q_dirs = x.block(0, 0, model.q_size, ndirs);
	MatrixNd qdot_dirs = x.block(model.q_size, 0, model.qdot_size, ndirs);

	// set derivative outputs
	std::vector<SpatialMatrix> ad_X_Ji (ndirs, SpatialMatrix::Zero());
	std::vector<std::vector<SpatialMatrix> > ad_X_J (model.mBodies.size(), ad_X_Ji);

	std::vector<SpatialVector> ad_V (ndirs, SpatialVector::Zero (6));
	std::vector<std::vector<SpatialVector> > ad_S (model.mBodies.size(), ad_V);
	std::vector<std::vector<SpatialVector> > ad_v_J (model.mBodies.size(), ad_V);
	std::vector<std::vector<SpatialVector> > ad_c_J (model.mBodies.size(), ad_V);

	for (unsigned int joint_id = 1; joint_id < model.mBodies.size(); joint_id++) {
		// evaluate nominal solution
		jcalc (model, joint_id, q, qdot);
		MatrixNd ref_X_Ji = model.X_J[joint_id].toMatrix();
		SpatialVector ref_S_i = model.S[joint_id];
		SpatialVector ref_v_Ji = model.v_J[joint_id];
		SpatialVector ref_c_Ji = model.c_J[joint_id];

		// evaluate AD nominal solution
		ad_jcalc (model, ad_model, joint_id, q, q_dirs, qdot, qdot_dirs);
		MatrixNd test_X_Ji = model.X_J[joint_id].toMatrix();
		SpatialVector test_S_i = model.S[joint_id];
		SpatialVector test_v_Ji = model.v_J[joint_id];
		SpatialVector test_c_Ji = model.c_J[joint_id];

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
	double TEST_PREC = 1.0e-08;

	// set nominal values
	q.setZero();
	q[0] = 0.3;
	q[1] = -0.2;

	qdot.setZero();
	qdot[0] = 0.3;
	qdot[1] = -0.2;

	// set directions
	unsigned int ndirs = model.q_size + model.qdot_size;
	MatrixNd x = MatrixNd::Identity(ndirs, ndirs);
	MatrixNd q_dirs = x.block(0, 0, model.q_size, ndirs);
	MatrixNd qdot_dirs = x.block(model.q_size, 0, model.qdot_size, ndirs);

	// set derivative outputs
	std::vector<SpatialMatrix> ad_X_Ji (ndirs, SpatialMatrix::Zero());
	std::vector<std::vector<SpatialMatrix> > fd_X_J (model.mBodies.size(), ad_X_Ji);
	std::vector<std::vector<SpatialMatrix> > ad_X_J (model.mBodies.size(), ad_X_Ji);

	std::vector<SpatialVector> ad_V (ndirs, SpatialVector::Zero (6));
	std::vector<std::vector<SpatialVector> > fd_S (model.mBodies.size(), ad_V);
	std::vector<std::vector<SpatialVector> > fd_v_J (model.mBodies.size(), ad_V);
	std::vector<std::vector<SpatialVector> > fd_c_J (model.mBodies.size(), ad_V);

	for (unsigned int joint_id = 1; joint_id < model.mBodies.size(); joint_id++) {
		// evaluate nominal solution
		fd_jcalc (model, joint_id, q, q_dirs, qdot, qdot_dirs,
			fd_X_J[joint_id], fd_S[joint_id], fd_v_J[joint_id], fd_c_J[joint_id]);

		// evaluate AD nominal solution
		ad_jcalc (model, ad_model, joint_id, q, q_dirs, qdot, qdot_dirs);

		for (int idir = 0; idir < ndirs; ++idir) {
			CHECK_ARRAY_CLOSE (fd_X_J[joint_id][idir].data(), ad_model.ad_X_J[joint_id][idir].data(), 36, TEST_PREC);
			CHECK_ARRAY_CLOSE (fd_S[joint_id][idir].data(),   ad_model.ad_S[joint_id][idir].data(),    6, TEST_PREC);
			CHECK_ARRAY_CLOSE (fd_v_J[joint_id][idir].data(), ad_model.ad_v_J[joint_id][idir].data(),  6, TEST_PREC);
			CHECK_ARRAY_CLOSE (fd_c_J[joint_id][idir].data(), ad_model.ad_c_J[joint_id][idir].data(),  6, TEST_PREC);
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

TEST_FIXTURE (CartPendulum, CalcPointAccelerationNominalTest) {
	// set nominal values
	q.setZero();
	q[0] = 0.3;
	q[1] = -0.2;

	qdot.setZero();
	qdot[0] = 0.3;
	qdot[1] = -0.2;

	qddot.setZero();
	qddot[0] = 0.3;
	qddot[1] = -0.2;

	Vector3d point_body_coordinates (0.1, 3.2, 4.2);

	// set directions
	unsigned int ndirs = 3*model.q_size;
	MatrixNd x = MatrixNd::Identity(ndirs, ndirs);
	MatrixNd q_dirs = x.block(0, 0, model.q_size, ndirs);
	MatrixNd qdot_dirs = x.block(model.q_size, 0, model.q_size, ndirs);
	MatrixNd qddot_dirs = x.block(2*model.q_size, 0, model.q_size, ndirs);

	// set derivative output
	MatrixNd ad_jacobian = MatrixNd::Zero(3, ndirs);

	// evaluate nominal solution
	Vector3d ref_acc = CalcPointAcceleration (
		model, q, qdot, qddot, id_pendulum, point_body_coordinates
	);

	Vector3d ad_acc = ad_CalcPointAcceleration (
		model, ad_model,
		q, q_dirs, qdot, qdot_dirs, qddot, qddot_dirs,
		id_pendulum, point_body_coordinates,
		ad_jacobian
	);

	CHECK_ARRAY_CLOSE (ref_acc.data(), ad_acc.data(), 3, TEST_PREC);
}


TEST_FIXTURE (CartPendulum, CalcPointAccelerationFDvsADTest) {
	// set nominal values
	q.setZero();
	q[0] = 0.3;
	q[1] = -0.2;

	qdot.setZero();
	qdot[0] = 0.3;
	qdot[1] = -0.2;

	qddot.setZero();
	qddot[0] = 0.3;
	qddot[1] = -0.2;

	Vector3d point_body_coordinates (0.1, 3.2, 4.2);
	
	// set directions
	unsigned int ndirs = model.q_size + model.qdot_size + model.q_size;
	MatrixNd x = MatrixNd::Identity(ndirs, ndirs);
	MatrixNd q_dirs = x.block(0, 0, model.q_size, ndirs);
	MatrixNd qdot_dirs = x.block(model.q_size, 0, model.qdot_size, ndirs);
	MatrixNd qddot_dirs = x.block(model.q_size + model.qdot_size, 0, model.q_size, ndirs);
	
	
	// set derivative output
	MatrixNd fd_jacobian = MatrixNd::Zero(3, ndirs);
	MatrixNd ad_jacobian = MatrixNd::Zero(3, ndirs);

	// evaluate nominal solution
	Vector3d ref_acc = fd_CalcPointAcceleration (
		model,
		q, q_dirs, qdot, qdot_dirs, qddot, qddot_dirs,
		id_pendulum, point_body_coordinates,
		ad_jacobian
	);

	Vector3d ad_acc = ad_CalcPointAcceleration (
		model, ad_model,
		q, q_dirs, qdot, qdot_dirs, qddot, qddot_dirs,
		id_pendulum, point_body_coordinates,
		ad_jacobian
	);

	CHECK_ARRAY_CLOSE (ref_acc.data(), ad_acc.data(), 3, TEST_PREC);
}


TEST_FIXTURE ( CartPendulum, CartPendulumCalcBodyToBaseSimple) {
	q[0] = 0.3;
	q[1] = -0.2;
	Vector3d point_body_coordinates (0.1, 3.2, 4.2);

	Vector3d point_single_func = CalcBodyToBaseCoordinatesSingleFunc (model, q, id_pendulum, point_body_coordinates);
	Vector3d point_default = CalcBodyToBaseCoordinates (model, q, id_pendulum, point_body_coordinates);

	CHECK_ARRAY_CLOSE (point_default.data(), point_single_func.data(), 3, TEST_PREC);
}

TEST_FIXTURE ( CartPendulum, CartPendulumJacobianADSimple ) {
	MatrixNd jacobian_ad = MatrixNd::Zero(3, model.qdot_size);
	MatrixNd jacobian_ref = MatrixNd::Zero(3, model.qdot_size);
	MatrixNd jacobian_fd = MatrixNd::Zero(3, model.qdot_size);

	q.setZero();
	q[0] = -1.0;
	q[1] = -0.3;
	body_point = Vector3d (1.0, 2.0, 3.0);

	CalcPointJacobian (model, q, id_pendulum, body_point, jacobian_ref);

	MatrixNd q_dirs = MatrixNd::Identity (model.qdot_size, model.qdot_size);
	Vector3d base_point_standard = CalcBodyToBaseCoordinates (model, q, id_pendulum, body_point);
	Vector3d base_point_ad = ad_CalcBodyToBaseCoordinatesSingleFunc (model, ad_model, q, q_dirs, id_pendulum, body_point, jacobian_ad);
	Vector3d base_point_fd = fd_CalcBodyToBaseCoordinatesSingleFunc (model, q, q_dirs, id_pendulum, body_point, jacobian_fd);

	CHECK_ARRAY_CLOSE (jacobian_ref.data(), jacobian_ad.data(), 3 * model.qdot_size, TEST_PREC);
//	CHECK_ARRAY_CLOSE (v_fixed_body.data(), v_body.data(), 6, TEST_PREC);
}

TEST_FIXTURE( CartPendulum, CompositeRigidBodyAlgorithmADTest) {
	MatrixNd q_dirs = MatrixNd::Identity (model.qdot_size, model.qdot_size);

	MatrixNd H = MatrixNd::Zero(model.q_size,model.q_size);
	std::vector<MatrixNd> fd_out;
	std::vector<MatrixNd> ad_out(q_dirs.cols(),MatrixNd::Zero(model.q_size,model.q_size));
	ADModel ad_model(model);
	
	VectorNd qrandom = VectorNd::Random(model.qdot_size);
	fd_CompositeRigidBodyAlgorithm(model, qrandom, q_dirs, fd_out);
	ad_CompositeRigidBodyAlgorithm(model, ad_model, qrandom, q_dirs, H, ad_out);
		
	MatrixNd inertia_test = MatrixNd::Zero(qrandom.size(),qrandom.size());
	CompositeRigidBodyAlgorithm(model, qrandom, inertia_test, true);
	for (size_t nIdx = 0; nIdx < fd_out.size(); nIdx++) {
// 		cout << "Dir " << nIdx << endl; 
// 		cout << "fd_CRBA: " << endl << fd_out[nIdx] << std::endl << std::endl;
// 		cout << "ad_CRBA: " << endl << ad_out[nIdx] << std::endl << std::endl;
// 		cout << "ad_CRBA error (fd,ad)" << endl << (fd_out[nIdx] - ad_out[nIdx]) << std::endl;
		CHECK_ARRAY_CLOSE (inertia_test.data(), H.data(), model.q_size*model.q_size, 1e-08);
		CHECK_ARRAY_CLOSE (fd_out[nIdx].data(), ad_out[nIdx].data(), model.q_size*model.q_size, 1e-08);
	}
}


TEST_FIXTURE(CartPendulum, InverseDynamicsADTest){
	for(unsigned int i = 0; i < model.qdot_size; i++){
		q[i] = (i+1)*(0.897878435);
		qdot[i] = (i+1)*(0.27563682);
		qddot[i] = (i+1)*(0.06565644564455);
	}

	MatrixNd q_dirs 	= MatrixNd::Identity (model.qdot_size, model.qdot_size);
	MatrixNd qdot_dirs 	= MatrixNd::Identity (model.qdot_size, model.qdot_size);
	MatrixNd qddot_dirs = MatrixNd::Identity (model.qdot_size, model.qdot_size);

	std::vector<SpatialVector> f_ext (model.mBodies.size(),SpatialVector::Zero(model.q_size));
	for (int i = 0; i < model.mBodies.size(); ++i) {
	  f_ext[i]=SpatialVector::Zero(model.q_size);
	}

	
	double h = sqrt(1.0e-16);
	VectorNd tau_ref (tau);

	MatrixNd ad_tau  = MatrixNd::Zero(model.qdot_size, model.qdot_size);
	MatrixNd fd_tau  = MatrixNd::Zero(model.qdot_size, model.qdot_size);

	ad_InverseDynamics(model,
			ad_model,
			q,
			q_dirs,
			qdot,
			qdot_dirs,
			qddot,
			qddot_dirs,
			tau,
			ad_tau,
			&f_ext);    


	fd_InverseDynamics(model,
			q,
			q_dirs,
			qdot,
			qdot_dirs,
			qddot,
			qddot_dirs,
			tau,
			fd_tau,
			&f_ext);    

	CHECK_ARRAY_CLOSE (fd_tau.data(), ad_tau.data(), fd_tau.cols()*fd_tau.rows(), 1e-7);
}

TEST_FIXTURE(CartPendulum, ForwardDynamicsADTest){
  srand((unsigned int) time(0));

  for(unsigned int trial = 0; trial < 10; trial++) {
    VectorNd q = VectorNd::Random(model.q_size);
    VectorNd qdot = VectorNd::Random(model.q_size);
    VectorNd tau = VectorNd::Random(model.q_size);

    unsigned int ndirs = 3 * model.q_size;
    MatrixNd x = MatrixNd::Identity(ndirs, ndirs);
    MatrixNd q_dirs = x.block(0, 0, model.q_size, ndirs);
    MatrixNd qdot_dirs = x.block(model.q_size, 0, model.q_size, ndirs);
    MatrixNd tau_dirs = x.block(2*model.q_size, 0, model.q_size, ndirs);

    std::vector<SpatialVector> f_ext (model.mBodies.size(),SpatialVector::Zero(model.q_size));
    for (int i = 0; i < model.mBodies.size(); ++i) {
      f_ext[i]=SpatialVector::Zero(model.q_size);
    }
    
    double h = sqrt(1.0e-16);

    MatrixNd ad_qddot  = MatrixNd::Zero(model.qdot_size, ndirs);
    MatrixNd fd_qddot  = MatrixNd::Zero(model.qdot_size, ndirs); 


    fd_ForwardDynamics(model, q, q_dirs, qdot, qdot_dirs, tau, tau_dirs, qddot, fd_qddot,&f_ext);    

    ad_ForwardDynamics(model, ad_model, q, q_dirs, qdot, qdot_dirs, tau, tau_dirs, qddot, ad_qddot, &f_ext);

    CHECK_ARRAY_CLOSE (fd_qddot.data(), ad_qddot.data(), fd_qddot.cols()*fd_qddot.rows(), 1e-7);
  }
  }

TEST_FIXTURE (CartPendulum, ForwardDynamicsCholesky) {
  // set nominal values

  q = VectorNd::Random(model.q_size);
  qdot = VectorNd::Random(model.q_size);
  qddot = VectorNd::Random(model.q_size);
  tau = VectorNd::Random(model.q_size);

  
  std::vector<SpatialVector> f_ext (model.mBodies.size(),SpatialVector::Zero(model.q_size));
  for (int i = 0; i < model.mBodies.size(); ++i) {
    f_ext[i]=SpatialVector::Random(model.q_size);
  }


  VectorNd qddot_ref (VectorNd::Zero (model.q_size));
	
  ForwardDynamicsCholesky(model,q,qdot,tau,qddot,&f_ext);
		
  ForwardDynamics(model,q,qdot,tau,qddot_ref,&f_ext);
	
  CHECK_ARRAY_CLOSE (qddot, qddot_ref, model.q_size, TEST_PREC);
  /*
    cout << "qddot_ref: " << endl << qddot_ref << endl;
    cout << "qddot_test: " << endl << qddot << endl;
    cout << "error qddot: " << endl << qddot_ref - qddot << endl;
  */
}


TEST_FIXTURE(CartPendulum, ForwardDynamicsCholeskyADTest){

  VectorNd q = VectorNd::Random(model.q_size);
  VectorNd qdot = VectorNd::Random(model.q_size);
  VectorNd tau = VectorNd::Random(model.q_size);
  VectorNd qddot = VectorNd::Zero(model.q_size);
  VectorNd qddot_ref = VectorNd::Zero(model.q_size);

  unsigned int ndirs = 3 * model.q_size;
  MatrixNd x = MatrixNd::Identity(ndirs, ndirs);
  MatrixNd q_dirs = x.block(0, 0, model.q_size, ndirs);
  MatrixNd qdot_dirs = x.block(model.q_size, 0, model.q_size, ndirs);
  MatrixNd tau_dirs = x.block(2*model.q_size, 0, model.q_size, ndirs);


  std::vector<SpatialVector> f_ext (model.mBodies.size(),SpatialVector::Zero(model.q_size));
  for (int i = 0; i < model.mBodies.size(); ++i) {
    f_ext[i]=SpatialVector::Zero(model.q_size);
  }

  MatrixNd ad_qddot  = MatrixNd::Zero(model.q_size, ndirs);
  MatrixNd fd_qddot  = MatrixNd::Zero(model.q_size, ndirs); 

  fd_ForwardDynamics(model,
		     q,
		     q_dirs,
		     qdot,
		     qdot_dirs,
		     tau,
		     tau_dirs,
		     qddot_ref,
		     fd_qddot,
		     &f_ext);    

  ad_ForwardDynamicsCholesky(model,
			     ad_model,
			     q,
			     q_dirs,
			     qdot,
			     qdot_dirs,
			     tau,
			     tau_dirs,
			     qddot,
			     ad_qddot,
			     &f_ext);

  CHECK_ARRAY_CLOSE (qddot_ref.data(), qddot.data(), qddot.size(), 1e-7);
  CHECK_ARRAY_CLOSE (fd_qddot.data(), ad_qddot.data(), fd_qddot.cols()*fd_qddot.rows(), 1e-7);
}




TEST_FIXTURE(CartPendulum, InverseDynamicsADTest_with_external_forces){
	for(unsigned int i = 0; i < model.qdot_size; i++){
		q[i] = (i+1)*(0.897878435);
		qdot[i] = (i+1)*(0.27563682);
		qddot[i] = (i+1)*(0.06565644564455);
	}

	MatrixNd q_dirs 	= MatrixNd::Identity (model.qdot_size, model.qdot_size);
	MatrixNd qdot_dirs 	= MatrixNd::Identity (model.qdot_size, model.qdot_size);
	MatrixNd qddot_dirs = MatrixNd::Identity (model.qdot_size, model.qdot_size);

	std::vector<SpatialVector> f_ext (model.mBodies.size(),SpatialVector::Zero(model.q_size));
	for (int i = 0; i < model.mBodies.size(); ++i) {
	  f_ext[i]=SpatialVector::Random(model.q_size);
	}

	
	double h = sqrt(1.0e-16);
	VectorNd tau_ref (tau);

	MatrixNd ad_tau  = MatrixNd::Zero(model.qdot_size, model.qdot_size);
	MatrixNd fd_tau  = MatrixNd::Zero(model.qdot_size, model.qdot_size);

	ad_InverseDynamics(model,
			ad_model,
			q,
			q_dirs,
			qdot,
			qdot_dirs,
			qddot,
			qddot_dirs,
			tau,
			ad_tau,
			&f_ext);    


	fd_InverseDynamics(model,
			q,
			q_dirs,
			qdot,
			qdot_dirs,
			qddot,
			qddot_dirs,
			tau_ref,
			fd_tau,
			&f_ext);    

	CHECK_ARRAY_CLOSE (fd_tau.data(), ad_tau.data(), fd_tau.cols()*fd_tau.rows(), 1e-7);
	CHECK_ARRAY_CLOSE (tau_ref.data(), tau.data(), tau_ref.rows(), 1e-7);

}

TEST_FIXTURE(CartPendulum, ForwardDynamicsADTest_with_external_forces){
  srand((unsigned int) time(0));

  for(unsigned int trial = 0; trial < 10; trial++) {
    VectorNd q = VectorNd::Random(model.q_size);
    VectorNd qdot = VectorNd::Random(model.q_size);
    VectorNd tau = VectorNd::Random(model.q_size);

    unsigned int ndirs = 3 * model.q_size;
    MatrixNd x = MatrixNd::Identity(ndirs, ndirs);
    MatrixNd q_dirs = x.block(0, 0, model.q_size, ndirs);
    MatrixNd qdot_dirs = x.block(model.q_size, 0, model.q_size, ndirs);
    MatrixNd tau_dirs = x.block(2*model.q_size, 0, model.q_size, ndirs);

    std::vector<SpatialVector> f_ext (model.mBodies.size(),SpatialVector::Zero(model.q_size));
    for (int i = 0; i < model.mBodies.size(); ++i) {
      f_ext[i]=SpatialVector::Random(model.q_size);
    }
    
    double h = sqrt(1.0e-16);

    MatrixNd ad_qddot  = MatrixNd::Zero(model.qdot_size, ndirs);
    MatrixNd fd_qddot  = MatrixNd::Zero(model.qdot_size, ndirs); 


    fd_ForwardDynamics(model, q, q_dirs, qdot, qdot_dirs, tau, tau_dirs, qddot, fd_qddot,&f_ext);    

    ad_ForwardDynamics(model, ad_model, q, q_dirs, qdot, qdot_dirs, tau, tau_dirs, qddot, ad_qddot, &f_ext);

    CHECK_ARRAY_CLOSE (fd_qddot.data(), ad_qddot.data(), fd_qddot.cols()*fd_qddot.rows(), 1e-7);
  }
  }



TEST_FIXTURE(CartPendulum, ForwardDynamicsCholeskyADTest_with_external_forces){

  VectorNd q = VectorNd::Random(model.q_size);
  VectorNd qdot = VectorNd::Random(model.q_size);
  VectorNd tau = VectorNd::Random(model.q_size);
  VectorNd qddot = VectorNd::Zero(model.q_size);
  VectorNd qddot_ref = VectorNd::Zero(model.q_size);

  unsigned int ndirs = 3 * model.q_size;
  MatrixNd x = MatrixNd::Identity(ndirs, ndirs);
  MatrixNd q_dirs = x.block(0, 0, model.q_size, ndirs);
  MatrixNd qdot_dirs = x.block(model.q_size, 0, model.q_size, ndirs);
  MatrixNd tau_dirs = x.block(2*model.q_size, 0, model.q_size, ndirs);


  std::vector<SpatialVector> f_ext (model.mBodies.size(),SpatialVector::Zero(model.q_size));
  for (int i = 0; i < model.mBodies.size(); ++i) {
    f_ext[i]=SpatialVector::Random(model.q_size);
  }

  MatrixNd ad_qddot  = MatrixNd::Zero(model.q_size, ndirs);
  MatrixNd fd_qddot  = MatrixNd::Zero(model.q_size, ndirs); 

  fd_ForwardDynamics(model,
		     q,
		     q_dirs,
		     qdot,
		     qdot_dirs,
		     tau,
		     tau_dirs,
		     qddot_ref,
		     fd_qddot,
		     &f_ext);    

  ad_ForwardDynamicsCholesky(model,
			     ad_model,
			     q,
			     q_dirs,
			     qdot,
			     qdot_dirs,
			     tau,
			     tau_dirs,
			     qddot,
			     ad_qddot,
			     &f_ext);

  CHECK_ARRAY_CLOSE (qddot_ref.data(), qddot.data(), qddot.size(), 1e-7);
  CHECK_ARRAY_CLOSE (fd_qddot.data(), ad_qddot.data(), fd_qddot.cols()*fd_qddot.rows(), 1e-7);
}
