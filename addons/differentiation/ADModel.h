/*
 * RBDL - Rigid Body Dynamics Library
 * Copyright (c) 2011-2015 Martin Felis <martin.felis@iwr.uni-heidelberg.de>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#ifndef RBDL_MODEL_H
#define RBDL_MODEL_H

#include "rbdl/rbdl_math.h"
#include <map>
#include <list>
#include <assert.h>
#include <iostream>
#include <limits>
#include <cstring>

#include "rbdl/Logging.h"
#include "rbdl/Joint.h"
#include "rbdl/Body.h"

// std::vectors containing any objectst that have Eigen matrices or vectors
// as members need to have a special allocater. This can be achieved with
// the following macro.

using namespace std;
using namespace RigidBodyDynamics;
using namespace RigidBodyDynamics::Math;

namespace RigidBodyDynamics {

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

        ad_Ic.resize(model.mBodies.size(), ad_X);
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

    std::vector<std::vector<SpatialMatrix> > ad_Ic;
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

            for (int i = 0; i < ad_Ic.size(); ++i) {
                ad_Ic[i].resize(ndirs, SpatialMatrix::Zero());
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
            //  ad_u[i].resize(ndirs, 0.0);
            // }

            // for (int i = 0; i < ad_d.size(); ++i) {
            //  ad_d[i].resize(ndirs, 0.0);
            // }

            ad_u.resize(ad_u.rows(), ndirs);
            ad_d.resize(ad_d.rows(), ndirs);
        }
    }
};


/** @} */
}

/* _ADMODEL_H */
#endif
