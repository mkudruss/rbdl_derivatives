/*
 * RBDL - Rigid Body Dynamics Library
 * Copyright (c) 2011-2015 Martin Felis <martin.felis@iwr.uni-heidelberg.de>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#ifndef RBDL_DYNAMICS_AD_H
#define RBDL_DYNAMICS_AD_H

#include <iostream>
#include <limits>
#include <assert.h>
#include <string.h>

#include "rbdl/rbdl_mathutils.h"
#include "rbdl/Logging.h"

#include "rbdl/Model.h"
#include "rbdl/Joint.h"
#include "rbdl/Body.h"
#include "rbdl/Dynamics.h"
#include "rbdl/Kinematics.h"

#include "ModelAD.h"
#include "JointAD.h"

using namespace RigidBodyDynamics::Math;

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------
namespace AD {
// -----------------------------------------------------------------------------

/*
RBDL_DLLAPI
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
	std::vector<SpatialVector>* f_ext
);

RBDL_DLLAPI
void ad_InverseDynamics(
	Model& model,
	ADModel &ad_model,
	const Math::VectorNd& q,
	const Math::MatrixNd& q_dirs,
	const Math::VectorNd& qdot,
	const Math::MatrixNd& qdot_dirs,
	const Math::VectorNd& qddot,
	const Math::MatrixNd& qddot_dirs,
	Math::VectorNd& tau,
	Math::MatrixNd& ad_tau,
	std::vector<Math::SpatialVector> *f_ext
);
*/

RBDL_DLLAPI
void CompositeRigidBodyAlgorithm (
	Model &model,
	ADModel &ad_model,
	const VectorNd &q,
	const MatrixNd &q_dirs,
	MatrixNd &H,
	std::vector<MatrixNd> &out,
	bool update_kinematics = true
);

// -----------------------------------------------------------------------------
} /* AD */
// -----------------------------------------------------------------------------
} /* RigidBodyDynamics */
// -----------------------------------------------------------------------------

/* RBDL_MATH_AD_H */
#endif
