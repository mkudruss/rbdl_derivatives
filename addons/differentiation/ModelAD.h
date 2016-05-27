/*
 * RBDL - Rigid Body Dynamics Library
 * Copyright (c) 2011-2015 Martin Felis <martin.felis@iwr.uni-heidelberg.de>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#ifndef RBDL_MODEL_AD_H
#define RBDL_MODEL_AD_H

#include <assert.h>
#include <iostream>
#include <limits>
#include <cstring>

#include "rbdl/rbdl.h"
#include "rbdl/rbdl_math.h"

// #include "rbdl/Logging.h"
// #include "rbdl/Joint.h"
// #include "rbdl/Body.h"

// #include "rbdl/Model.h"

using namespace std;
using namespace RigidBodyDynamics;
//using namespace RigidBodyDynamics::Math;

// // -----------------------------------------------------------------------------
// namespace RigidBodyDynamics {
// // -----------------------------------------------------------------------------
// namespace AD {
// -----------------------------------------------------------------------------

struct ADModel {

    unsigned int ndirs;

    // derivative values
    vector<vector<Math::SpatialMatrix> > X_lambda;
    vector<vector<Math::SpatialMatrix> > X_base;
    vector<vector<Math::SpatialMatrix> > X_J;

    vector<vector<Math::SpatialVector> > S;
    vector<vector<Math::SpatialVector> > v_J;
    vector<vector<Math::SpatialVector> > v;
    vector<vector<Math::SpatialVector> > c_J;
    vector<vector<Math::SpatialVector> > a_J;
    vector<vector<Math::SpatialVector> > a;
    vector<vector<Math::SpatialVector> > c;
    vector<vector<Math::SpatialVector> > f;

    vector<vector<Math::SpatialMatrix> > Ic;
    vector<vector<Math::SpatialVector> > hc;
    vector<vector<Math::SpatialVector> > F;

    vector<vector<Math::SpatialVector> > pA;
    vector<vector<Math::SpatialVector> > U;
    vector<vector<Math::SpatialMatrix> > IA;

    //vector<vector<double> > u;
    //vector<vector<double> > d;
    Math::MatrixNd u;
    Math::MatrixNd d;

    ADModel ();
    ADModel (RigidBodyDynamics::Model& model);

    void resize_directions (unsigned requested_ndirs);

};

// -----------------------------------------------------------------------------
// } /* AD */
// -----------------------------------------------------------------------------
// } /* RigidBodyDynamics */
// -----------------------------------------------------------------------------

/* RBDL_MODEL_AD_H */
#endif
