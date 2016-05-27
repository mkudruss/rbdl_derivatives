/*
 * RBDL - Rigid Body Dynamics Library
 * Copyright (c) 2011-2015 Martin Felis <martin.felis@iwr.uni-heidelberg.de>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

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

using namespace RigidBodyDynamics::Math;

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------
namespace FD {
// -----------------------------------------------------------------------------

RBDL_DLLAPI
void ForwardDynamics(
	Model& model,
	const VectorNd& q,
	const MatrixNd& q_dirs,
	const VectorNd& qdot,
	const MatrixNd& qdot_dirs,
	const VectorNd& tau,
	const MatrixNd& tau_dirs,
	VectorNd& qddot,
	MatrixNd& fd_qddot,
	std::vector<SpatialVector>* f_ext
);

RBDL_DLLAPI
void InverseDynamics(
	Model& model,
	const Math::VectorNd& q,
	const Math::MatrixNd& q_dirs,
	const Math::VectorNd& qdot,
	const Math::MatrixNd& qdot_dirs,
	const Math::VectorNd& qddot,
	const Math::MatrixNd& qddot_dirs,
	Math::VectorNd& tau,
	Math::MatrixNd& fd_tau,
	std::vector<Math::SpatialVector> *f_ext
);

RBDL_DLLAPI
void CompositeRigidBodyAlgorithm (
	Model &model,
	const VectorNd &q,
	const MatrixNd &q_dirs,
	std::vector<MatrixNd> &fd_out
);

// -----------------------------------------------------------------------------
} /* FD */
// -----------------------------------------------------------------------------
} /* RigidBodyDynamics */
// -----------------------------------------------------------------------------
