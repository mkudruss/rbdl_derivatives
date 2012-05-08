/*
 * RBDL - Rigid Body Dynamics Library
 * Copyright (c) 2011-2012 Martin Felis <martin.felis@iwr.uni-heidelberg.de>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#include <iostream>
#include <limits>
#include <assert.h>

#include "rbdl_mathutils.h"

#include "Logging.h"

#include "Model.h"
#include "Body.h"
#include "Joint.h"

using namespace RigidBodyDynamics;
using namespace RigidBodyDynamics::Math;

void Model::Init() {
	Body root_body;
	Joint root_joint;

	Vector3d zero_position (0., 0., 0.);
	SpatialVector zero_spatial (0., 0., 0., 0., 0., 0.);

	// structural information
	lambda.push_back(0.);
	mu.push_back(std::vector<unsigned int>());
	dof_count = 0;

	gravity = Vector3d (0., -9.81, 0.);

	// state information
	v.push_back(zero_spatial);
	a.push_back(zero_spatial);

	// Joints
	mJoints.push_back(root_joint);
	S.push_back (zero_spatial);
	X_T.push_back(SpatialTransform());
	
	// Dynamic variables
	c.push_back(zero_spatial);
	IA.push_back(SpatialMatrixIdentity);
	pA.push_back(zero_spatial);
	U.push_back(zero_spatial);

	u = VectorNd::Zero(1);
	d = VectorNd::Zero(1);

	f.push_back (zero_spatial);
	Ic.push_back (
			SpatialRigidBodyInertia(
				0.,
				Vector3d (0., 0., 0.),
				Matrix3d::Zero(3,3)
				)
			);

	// Bodies
	X_lambda.push_back(SpatialTransform());
	X_base.push_back(SpatialTransform());

	mBodies.push_back(root_body);
	mBodyNames.push_back("ROOT");
}

unsigned int Model::AddBody (const unsigned int parent_id,
		const SpatialTransform &joint_frame,
		const Joint &joint,
		const Body &body,
		std::string body_name) {
	// Here we emulate multi DoF joints by simply adding nullbodies. This
	// may be changed in the future, but so far it is reasonably fast.
	if (joint.mJointType != JointTypePrismatic 
			&& joint.mJointType != JointTypeRevolute) {
		unsigned int joint_count;
		if (joint.mJointType == JointType1DoF)
			joint_count = 1;
		else if (joint.mJointType == JointType2DoF)
			joint_count = 2;
		else if (joint.mJointType == JointType3DoF)
			joint_count = 3;
		else if (joint.mJointType == JointType4DoF)
			joint_count = 4;
		else if (joint.mJointType == JointType5DoF)
			joint_count = 5;
		else if (joint.mJointType == JointType6DoF)
			joint_count = 6;
		else {
			std::cerr << "Error: Invalid joint type: " << joint.mJointType << std::endl;
			assert (0 && !"Invalid joint type!");
		}

		Body null_body (0., Vector3d (0., 0., 0.), Vector3d (0., 0., 0.));
		unsigned int null_parent = parent_id;
		SpatialTransform joint_frame_transform;

		for (unsigned int j = 0; j < joint_count; j++) {
			Joint single_dof_joint;

			Vector3d rotation (
					joint.mJointAxes[j][0],
					joint.mJointAxes[j][1],
					joint.mJointAxes[j][2]);
			Vector3d translation (
					joint.mJointAxes[j][3],
					joint.mJointAxes[j][4],
					joint.mJointAxes[j][5]);

			if (rotation == Vector3d (0., 0., 0.)) {
				single_dof_joint = Joint (JointTypePrismatic, translation);
			} else if (translation == Vector3d (0., 0., 0.)) {
				single_dof_joint = Joint (JointTypeRevolute, rotation);
			}

			// the first joint has to be transformed by joint_frame, all the
			// others must have a null transformation
			if (j == 0)
				joint_frame_transform = joint_frame;
			else
				joint_frame_transform = SpatialTransform();

			if (j == joint_count - 1)
				// if we are at the last we must add the real body
				return AddBody (null_parent, joint_frame_transform, single_dof_joint, body, body_name);
			else
				// otherwise we just add an intermediate body
				null_parent = AddBody (null_parent, joint_frame_transform, single_dof_joint, null_body);
		}
	}

	assert (lambda.size() > 0);
	assert (joint.mJointType != JointTypeUndefined);

	// structural information
	lambda.push_back(parent_id);
	mu.push_back(std::vector<unsigned int>());
	mu.at(parent_id).push_back(mBodies.size());

	// Bodies
	X_lambda.push_back(SpatialTransform());
	X_base.push_back(SpatialTransform());
	mBodies.push_back(body);
	mBodyNames.push_back(body_name);

	dof_count++;

	// state information
	v.push_back(SpatialVector(0., 0., 0., 0., 0., 0.));
	a.push_back(SpatialVector(0., 0., 0., 0., 0., 0.));

	// Joints
	mJoints.push_back(joint);
	S.push_back (joint.mJointAxes[0]);
	// we have to invert the transformation as it is later always used from the
	// child bodies perspective.
	X_T.push_back(joint_frame);

	// Dynamic variables
	c.push_back(SpatialVector(0., 0., 0., 0., 0., 0.));
	IA.push_back(body.mSpatialInertia);
	pA.push_back(SpatialVector(0., 0., 0., 0., 0., 0.));
	U.push_back(SpatialVector(0., 0., 0., 0., 0., 0.));

	d = VectorNd::Zero (mBodies.size());
	u = VectorNd::Zero (mBodies.size());

	f.push_back (SpatialVector (0., 0., 0., 0., 0., 0.));
	Ic.push_back (
			SpatialRigidBodyInertia(
				body.mMass,
				body.mCenterOfMass,
				body.mInertia
				)
			);

	return mBodies.size() - 1;
}

unsigned int Model::AppendBody (
		const Math::SpatialTransform &joint_frame,
		const Joint &joint,
		const Body &body,
		std::string body_name 
		) {
	unsigned int prev_body_id = mBodies.size() - 1;
	return Model::AddBody (prev_body_id, joint_frame,
			joint, body, body_name);
}

unsigned int Model::SetFloatingBaseBody (const Body &body) {
	assert (lambda.size() >= 0);

	// Add five zero weight bodies and append the given body last. Order of
	// the degrees of freedom is:
	// 		tx ty tz rz ry rx
	//

	Joint floating_base_joint (
			SpatialVector (0., 0., 0., 1., 0., 0.),
			SpatialVector (0., 0., 0., 0., 1., 0.),
			SpatialVector (0., 0., 0., 0., 0., 1.),
			SpatialVector (0., 0., 1., 0., 0., 0.),
			SpatialVector (0., 1., 0., 0., 0., 0.),
			SpatialVector (1., 0., 0., 0., 0., 0.)
			);

	unsigned int body_id = this->AddBody (0, Xtrans (Vector3d (0., 0., 0.)), floating_base_joint, body);

	return body_id;
}

unsigned int Model::GetBodyId (const char *id) const {
	for (unsigned int i = 0; i < mBodyNames.size(); i++) {
		if (mBodyNames[i] == id)
			return i;
	}

	return std::numeric_limits<unsigned int>::max();
}
