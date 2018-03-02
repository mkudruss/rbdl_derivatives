/*
 * RBDL - Rigid Body Dynamics Library
 * Copyright (c) 2011-2015 Martin Felis <martin.felis@iwr.uni-heidelberg.de>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#ifndef RBDL_MODEL_ED_H
#define RBDL_MODEL_ED_H

#include <assert.h>
#include <iostream>
#include <limits>
#include <cstring>

#include "rbdl/rbdl.h"
#include "rbdl/rbdl_math.h"

#include "SpatialAlgebraOperatorsED.h"

//using namespace std;
//using namespace RigidBodyDynamics;
//using namespace RigidBodyDynamics::Math;

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------
// namespace ED {
// -----------------------------------------------------------------------------

struct EDModel {

    unsigned int ndirs;

    // derivative values
    std::vector<std::vector<RigidBodyDynamics::Math::SpatialTransform> > X_lambda;
    std::vector<std::vector<RigidBodyDynamics::Math::SpatialTransform> > X_base;
    std::vector<std::vector<RigidBodyDynamics::Math::ED::SpatialTransformDot> > X_J;

    std::vector<std::vector<RigidBodyDynamics::Math::SpatialVector> > v_J;
    std::vector<std::vector<RigidBodyDynamics::Math::SpatialVector> > v;
    std::vector<std::vector<RigidBodyDynamics::Math::SpatialVector> > c_J;
    std::vector<std::vector<RigidBodyDynamics::Math::SpatialVector> > a;
    std::vector<std::vector<RigidBodyDynamics::Math::SpatialVector> > c;
    std::vector<std::vector<RigidBodyDynamics::Math::SpatialVector> > f;

//    std::vector<std::vector<RigidBodyDynamics::Math::SpatialMatrix> > Ic;

    std::vector<std::vector<RigidBodyDynamics::Math::SpatialRigidBodyInertia> > Ic;

    std::vector<std::vector<RigidBodyDynamics::Math::SpatialVector> > hc;
    std::vector<std::vector<RigidBodyDynamics::Math::SpatialVector> > F;

    std::vector<std::vector<RigidBodyDynamics::Math::SpatialVector> > pA;
    std::vector<std::vector<RigidBodyDynamics::Math::SpatialVector> > U;
    std::vector<std::vector<RigidBodyDynamics::Math::SpatialMatrix> > IA;

    //std::vector<std::vector<double> > u;
    //std::vector<std::vector<double> > d;
    RigidBodyDynamics::Math::MatrixNd u;
    RigidBodyDynamics::Math::MatrixNd d;

    EDModel ();
    EDModel (RigidBodyDynamics::Model& model);

    void resize_directions (unsigned requested_ndirs);

};

// -----------------------------------------------------------------------------
// } // ED
// -----------------------------------------------------------------------------
} // RigidBodyDynamics
// -----------------------------------------------------------------------------

/* RBDL_MODEL_ED_H */
#endif
