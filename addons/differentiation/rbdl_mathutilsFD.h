/*
 * RBDL - Rigid Body Dynamics Library
 * Copyright (c) 2011-2015 Martin Felis <martin.felis@iwr.uni-heidelberg.de>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#include <cmath>
#include <limits>

#include <iostream>
#include <assert.h>

#include <rbdl/Model.h>
#include <rbdl/rbdl_mathutils.h>
#include <rbdl/SpatialAlgebraOperators.h>

#include "rbdl/Logging.h"

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------
namespace Math {
// -----------------------------------------------------------------------------
namespace FD {
// -----------------------------------------------------------------------------

SpatialVector crossm (
        const SpatialVector &v1, const SpatialVector &v1_dirs,
        const SpatialVector &v2, const SpatialVector &v2_dirs
) {
    double h = 1.0e-8;

    // evaluate y(t+h*d)
    SpatialVector res = RigidBodyDynamics::Math::crossm (v1, v2);

    // evaluate y(t+h*d)
    SpatialVector res_hd = RigidBodyDynamics::Math::crossm (v1 + h * v1_dirs, v2 + h * v2_dirs);

    return (res_hd - res) / h;
};

SpatialMatrix crossm (
        const SpatialVector &v, const SpatialVector &v_dirs
) {
    double h = 1.0e-8;

    // evaluate y(t+h*d)
    SpatialMatrix res = RigidBodyDynamics::Math::crossm (v);

    // evaluate y(t+h*d)
    SpatialMatrix res_hd = RigidBodyDynamics::Math::crossm (v + h * v_dirs);

    return (res_hd - res) / h;
};

// -----------------------------------------------------------------------------
} /* FD */
// -----------------------------------------------------------------------------
} /* Math */
// -----------------------------------------------------------------------------
} /* RigidBodyDynamics */
// -----------------------------------------------------------------------------
