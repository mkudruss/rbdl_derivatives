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

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------
namespace FD {
// -----------------------------------------------------------------------------

RBDL_DLLAPI
void ForwardDynamics(
    Model& model,
    const Math::VectorNd& q,
    const Math::MatrixNd& q_dirs,
    const Math::VectorNd& qdot,
    const Math::MatrixNd& qdot_dirs,
    const Math::VectorNd& tau,
    const Math::MatrixNd& tau_dirs,
    Math::VectorNd& qddot,
    Math::MatrixNd& fd_qddot,
    std::vector<Math::SpatialVector>* f_ext
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
void NonlinearEffects (
    Model & model,
    const Math::VectorNd & q,
    const Math::MatrixNd & q_dirs,
    const Math::VectorNd & qDot,
    const Math::MatrixNd & qDot_dirs,
    Math::VectorNd & tau,
    Math::MatrixNd & fd_tau
    );

RBDL_DLLAPI
void CompositeRigidBodyAlgorithm (
    Model &model,
    const Math::VectorNd &q,
    const Math::MatrixNd &q_dirs,
    std::vector<Math::MatrixNd> &fd_out
    );

// -----------------------------------------------------------------------------
} /* FD */
// -----------------------------------------------------------------------------
} /* RigidBodyDynamics */
// -----------------------------------------------------------------------------
