/*
 * RBDL - Rigid Body Library
 * Copyright (c) 2011 Martin Felis <martin.felis@iwr.uni-heidelberg.de>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#include <iostream>
#include <limits>
#include <cstring>
#include <assert.h>

#include "mathutils.h"
#include "Logging.h"

#include "Model.h"
#include "Kinematics.h"

using namespace SpatialAlgebra;
using namespace SpatialAlgebra::Operators;

namespace RigidBodyDynamics {

void ForwardKinematics (Model &model,
		const VectorNd &Q,
		const VectorNd &QDot,
		const VectorNd &QDDot
		) {
	LOG << "-------- " << __func__ << " --------" << std::endl;

	unsigned int i;

	assert (model.q.size() == Q.size() + 1);
	assert (model.qdot.size() == QDot.size() + 1);
	assert (model.qddot.size() == QDDot.size() + 1);

	if (model.experimental_floating_base) {
		assert (0 && !"ForwardKinematics not supported yet for experimental floating bases");
	}
	
	// positions
	for (i = 0; i < model.dof_count; i++) {
		model.q[i + 1] = Q[i];
	}

	// velocities
	for (i = 0; i < model.dof_count; i++) {
		model.qdot[i + 1] = QDot[i];
	}

	// accelerations
	for (i = 0; i < model.dof_count; i++) {
		model.qddot[i + 1] = QDDot[i];
	}

	for (i = 1; i < model.mBodies.size(); i++) {
		SpatialMatrix X_J;
		SpatialVector v_J;
		SpatialVector c_J;
		Joint joint = model.mJoints[i];
		unsigned int lambda = model.lambda[i];

		jcalc (model, i, X_J, model.S[i], v_J, c_J, model.q[i], model.qdot[i]);

		model.X_lambda[i] = X_J * model.X_T[i];

		if (lambda != 0) {
			model.X_base[i] = model.X_lambda[i] * model.X_base.at(lambda);
			model.v[i] = model.X_lambda[i] * model.v[lambda] + v_J;
			model.c[i] = c_J + crossm(model.v[i],v_J);
			model.a[i] = model.X_lambda[i] * model.a[lambda] + model.c[i];
		}	else {
			model.X_base[i] = model.X_lambda[i];
			model.v[i] = v_J;
			model.c[i].setZero();
			model.a[i].setZero();
		}

		model.a[i] = model.a[i] + model.S[i] * model.qddot[i];
	}
}

void CalcPointJacobian (
		Model &model,
		const VectorNd &Q,
		unsigned int body_id,
		const Vector3d &point_position,
		MatrixNd &G,
		bool update_kinematics
	) {
	LOG << "-------- " << __func__ << " --------" << std::endl;
	if (model.experimental_floating_base) {
		assert (0 && "CalcPointJacobian() not yet supported for experimental floating base models");
	};

	// update the Kinematics with zero acceleration
	if (update_kinematics) {
		VectorNd QDDot_zero = VectorNd::Zero(Q.size());
		VectorNd QDot_zero = VectorNd::Zero(Q.size());
		
		ForwardKinematics (model, Q, QDot_zero, QDDot_zero);
	}

	Vector3d point_base_pos = model.CalcBodyToBaseCoordinates(body_id, point_position);
	SpatialMatrix point_trans = Xtrans (point_base_pos);

	assert (G.rows() == 3 && G.cols() == model.dof_count );

	G.setZero();

	// we have to make sure that only the joints that contribute to the
	// bodies motion also get non-zero columns in the jacobian.
	//VectorNd e = VectorNd::Zero(Q.size() + 1);

	unsigned int j = body_id;

	char e[(Q.size() + 1)];
	memset (&e[0], 0, Q.size() + 1);

	// e will contain 
	while (j != 0) {
		e[j] = 1;
		j = model.lambda[j];
	}

	for (j = 1; j < model.mBodies.size(); j++) {
		if (e[j] == 1) {
			SpatialVector S_base;
			S_base = point_trans * spatial_inverse(model.X_base[j]) * model.S[j];

			G(0, j - 1) = S_base[3];
			G(1, j - 1) = S_base[4];
			G(2, j - 1) = S_base[5];
		}
	}
}

Vector3d CalcPointVelocity (
		Model &model,
		const VectorNd &Q,
		const VectorNd &QDot,
		unsigned int body_id,
		const Vector3d &point_position,
		bool update_kinematics
	) {
	LOG << "-------- " << __func__ << " --------" << std::endl;
	unsigned int i;

	if (model.experimental_floating_base) {
		// set the transformation for the base body
		model.X_base[0] = XtransRotZYXEuler (Vector3d (Q[0], Q[1], Q[2]), Vector3d (Q[3], Q[4], Q[5]));
		model.X_lambda[0] = model.X_base[0];

		// in this case the appropriate function has to be called, see
		// ForwardDynamicsFloatingBase
		model.v[0].set (QDot[5], QDot[4], QDot[3], QDot[0], QDot[1], QDot[2]);

		model.v[0] = spatial_inverse(model.X_base[0]) * model.v[0];
		
		// global_velocities[0] = model.v[0];

		if (model.mBodies.size() > 1) {
			// Copy state values from the input to the variables in model
			for (i = 1; i < model.mBodies.size(); i++) {
				model.q[i] = Q[6 + i];
				model.qdot[i] = QDot[6 + i];
			}
		}
	} else {
		assert (body_id > 0 && body_id < model.mBodies.size());
		assert (model.q.size() == Q.size() + 1);
		assert (model.qdot.size() == QDot.size() + 1);

		// Reset the velocity of the root body
		model.v[0].setZero();

		// Copy state values from the input to the variables in model
		for (i = 0; i < Q.size(); i++) {
			model.q[i+1] = Q[i];
			model.qdot[i+1] = QDot[i];
		}
	}

	// update the Kinematics with zero acceleration
	if (update_kinematics) {
		VectorNd QDDot_zero = VectorNd::Zero(Q.size());
		
		ForwardKinematics (model, Q, QDot, QDDot_zero);
	}

	Vector3d point_abs_pos = model.CalcBodyToBaseCoordinates(body_id, point_position); 

	LOG << "body_index     = " << body_id << std::endl;
	LOG << "point_pos      = " << point_position << std::endl;
//	LOG << "global_velo    = " << global_velocities.at(body_id) << std::endl;
	LOG << "body_transf    = " << model.X_base[body_id] << std::endl;
	LOG << "point_abs_ps   = " << point_abs_pos << std::endl;
	LOG << "X   = " << Xtrans (point_abs_pos) * spatial_inverse(model.X_base[body_id]) << std::endl;
	LOG << "v   = " << model.v[body_id] << std::endl;

	// Now we can compute the spatial velocity at the given point
//	SpatialVector body_global_velocity (global_velocities.at(body_id));
	SpatialVector point_spatial_velocity = Xtrans (point_abs_pos) * spatial_inverse(model.X_base[body_id]) * model.v[body_id];

	LOG << "point_velocity = " <<	Vector3d (
			point_spatial_velocity[3],
			point_spatial_velocity[4],
			point_spatial_velocity[5]
			) << std::endl;

	return Vector3d (
			point_spatial_velocity[3],
			point_spatial_velocity[4],
			point_spatial_velocity[5]
			);
}

Vector3d CalcPointAcceleration (
		Model &model,
		const VectorNd &Q,
		const VectorNd &QDot,
		const VectorNd &QDDot,
		unsigned int body_id,
		const Vector3d &point_position,
			bool update_kinematics
	)
{
	LOG << "-------- " << __func__ << " --------" << std::endl;
	unsigned int i;

	// Reset the velocity of the root body
	model.v[0].setZero();
	model.a[0].setZero();

	if (model.experimental_floating_base) {
		// set the transformation for the base body
		model.X_base[0] = XtransRotZYXEuler (Vector3d (Q[0], Q[1], Q[2]), Vector3d (Q[3], Q[4], Q[5]));
		model.X_lambda[0] = model.X_base[0];

		// in this case the appropriate function has to be called, see
		// ForwardDynamicsFloatingBase.
		model.v[0].set (QDot[5], QDot[4], QDot[3], QDot[0], QDot[1], QDot[2]);
		model.a[0].set (QDDot[5], QDDot[4], QDDot[3], QDDot[0], QDDot[1], QDDot[2]);
		
		if (model.mBodies.size() > 1) {
			// Copy state values from the input to the variables in model
			for (i = 1; i < model.mBodies.size(); i++) {
				model.q[i] = Q[6 + i];
				model.qdot[i] = QDot[6 + i];
				model.qddot[i] = QDDot[6 + i];
			}
		}
	} else {
		if (update_kinematics)
			ForwardKinematics (model, Q, QDot, QDDot);
	}

	LOG << std::endl;

	// computation of the global position of the point
	Vector3d point_abs_pos = model.CalcBodyToBaseCoordinates (body_id, point_position);
	LOG << "point_abs_ps = " << point_abs_pos << std::endl;

	// The whole computation looks in formulae like the following:
	SpatialVector body_global_velocity (spatial_inverse(model.X_base[body_id]) * model.v[body_id]);
	SpatialVector body_global_acceleration (spatial_inverse(model.X_base[body_id]) * model.a[body_id]);
	SpatialMatrix global_point_transform (Xtrans (point_abs_pos));
	SpatialMatrix local_point_transform (Xtrans (point_position));

	LOG << " orientation " << std::endl << model.GetBodyWorldOrientation(body_id) << std::endl;
	LOG << " orientationT " << std::endl <<  model.GetBodyWorldOrientation(body_id).transpose() << std::endl;

	Matrix3d global_body_orientation_inv = model.GetBodyWorldOrientation (body_id).transpose();
	SpatialMatrix p_X_i = SpatialMatrixZero;

	p_X_i.block<3,3>(0,0) = global_body_orientation_inv;
	p_X_i.block<3,3>(3,3) = global_body_orientation_inv;

	LOG << " p_X_i = " << std::endl << p_X_i << std::endl;
	LOG << " xtrans = " << std::endl << Xtrans (point_position) << std::endl;

	p_X_i *= Xtrans (point_position);

	SpatialVector p_v_i = p_X_i * model.v[body_id];
	SpatialVector p_a_i = p_X_i * model.a[body_id];

	SpatialVector frame_acceleration = 
		crossm( SpatialVector(0., 0., 0., p_v_i[3], p_v_i[4], p_v_i[5]), (body_global_velocity));

	LOG << model.X_base[body_id] << std::endl;
	LOG << "v_i                = " << model.v[body_id] << std::endl;
	LOG << "a_i                = " << model.a[body_id] << std::endl;
	LOG << "p_X_i              = " << p_X_i << std::endl;
	LOG << "p_v_i              = " << p_v_i << std::endl;
	LOG << "p_a_i              = " << p_a_i << std::endl;
	LOG << "body_global_vel    = " << body_global_velocity << std::endl;
	LOG << "frame_acceleration = " << frame_acceleration << std::endl;

	SpatialVector p_a_i_dash = p_a_i - frame_acceleration;

	LOG << "point_acceleration = " <<	Vector3d (
			p_a_i_dash[3],
			p_a_i_dash[4],
			p_a_i_dash[5]
			) << std::endl;

	return Vector3d (
			p_a_i_dash[3],
			p_a_i_dash[4],
			p_a_i_dash[5]
			);
}

}
