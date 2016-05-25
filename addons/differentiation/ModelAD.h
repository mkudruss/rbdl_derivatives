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
using namespace RigidBodyDynamics::Math;

// // -----------------------------------------------------------------------------
// namespace RigidBodyDynamics {
// // -----------------------------------------------------------------------------
// namespace AD {
// -----------------------------------------------------------------------------

struct ADModel {

    unsigned int ndirs;

    // derivative values
    vector<vector<SpatialMatrix> > X_lambda;
    vector<vector<SpatialMatrix> > X_base;
    vector<vector<SpatialMatrix> > X_J;

    vector<vector<SpatialVector> > S;
    vector<vector<SpatialVector> > v_J;
    vector<vector<SpatialVector> > v;
    vector<vector<SpatialVector> > c_J;
    vector<vector<SpatialVector> > a_J;
    vector<vector<SpatialVector> > a;
    vector<vector<SpatialVector> > c;
    vector<vector<SpatialVector> > f;

    vector<vector<SpatialMatrix> > Ic;
    vector<vector<SpatialVector> > hc;
    vector<vector<SpatialVector> > F;

    vector<vector<SpatialVector> > pA;
    vector<vector<SpatialVector> > U;
    vector<vector<SpatialMatrix> > IA;

    //vector<vector<double> > u;
    //vector<vector<double> > d;
    MatrixNd u;
    MatrixNd d;

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
