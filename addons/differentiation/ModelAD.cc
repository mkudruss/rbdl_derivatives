/*
 * RBDL - Rigid Body Dynamics Library
 * Copyright (c) 2011-2015 Martin Felis <martin.felis@iwr.uni-heidelberg.de>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#include "ModelAD.h"

// std::vectors containing any objectst that have Eigen matrices or vectors
// as members need to have a special allocater. This can be achieved with
// the following macro.

ADModel::ADModel () {
    // TODO: should not be called, I think!
    //std::cerr << "Data structure can only be initialized with argument RigidBodyDynamics::Model!" << std::endl;
    //throw std::exception;
};

ADModel::ADModel (RigidBodyDynamics::Model& model) {
    ndirs = 4 * model.dof_count;

    // NOTE: old initialization values
    std::vector<SpatialMatrix> X(ndirs, SpatialMatrix::Zero());
    std::vector<SpatialVector> vec(ndirs, SpatialVector::Zero());
    //std::vector<double> double(ndirs, 0.0);

    X_lambda.resize(model.mBodies.size(), X);
    X_base.resize(model.mBodies.size(), X);
    X_J.resize(model.mBodies.size(), X);

    S.resize(model.mBodies.size(), vec);
    v_J.resize(model.mBodies.size(), vec);
    v.resize(model.mBodies.size(), vec);
    c_J.resize(model.mBodies.size(), vec);
    a_J.resize(model.mBodies.size(), vec);
    a.resize(model.mBodies.size(), vec);

    Ic.resize(model.mBodies.size(), X);
    F.resize(model.mBodies.size(), vec);
    c.resize(model.mBodies.size(), vec);
    f.resize(model.mBodies.size(), vec);

    pA.resize(model.mBodies.size(), vec);
    U.resize(model.mBodies.size(), vec);
    IA.resize(model.mBodies.size(), X);
    //u.resize(model.mBodies.size(), double);
    //d.resize(model.mBodies.size(), double);
    u.resize(model.q_size, ndirs);
    d.resize(model.q_size, ndirs);
};

void ADModel::resize_directions (unsigned requested_ndirs){
    if (ndirs < requested_ndirs) {
        ndirs = requested_ndirs;

        for (int i = 0; i < X_lambda.size(); ++i) {
            X_lambda[i].resize(ndirs, SpatialMatrix::Zero());
        }

        for (int i = 0; i < X_base.size(); ++i) {
            X_base[i].resize(ndirs, SpatialMatrix::Zero());
        }

        for (int i = 0; i < X_J.size(); ++i) {
            X_J[i].resize(ndirs, SpatialMatrix::Zero());
        }

        for (int i = 0; i < S.size(); ++i) {
            S[i].resize(ndirs, SpatialVector::Zero());
        }

        for (int i = 0; i < v_J.size(); ++i) {
            v_J[i].resize(ndirs, SpatialVector::Zero());
        }

        for (int i = 0; i < v_J.size(); ++i) {
            v[i].resize(ndirs, SpatialVector::Zero());
        }

        for (int i = 0; i < c_J.size(); ++i) {
            c_J[i].resize(ndirs, SpatialVector::Zero());
        }

        for (int i = 0; i < c_J.size(); ++i) {
            a_J[i].resize(ndirs, SpatialVector::Zero());
        }

        for (int i = 0; i < c_J.size(); ++i) {
            a[i].resize(ndirs, SpatialVector::Zero());
        }

        for (int i = 0; i < c_J.size(); ++i) {
            c[i].resize(ndirs, SpatialVector::Zero());
        }

        for (int i = 0; i < c_J.size(); ++i) {
            f[i].resize(ndirs, SpatialVector::Zero());
        }

        for (int i = 0; i < Ic.size(); ++i) {
            Ic[i].resize(ndirs, SpatialMatrix::Zero());
        }

        for (int i = 0; i < F.size(); ++i) {
            F[i].resize(ndirs, SpatialVector::Zero());
        }

        for(int i = 0; i < pA.size(); ++i) {
            pA[i].resize(ndirs, SpatialVector::Zero());
        }

        for(int i = 0; i < U.size(); ++i) {
            U[i].resize(ndirs, SpatialVector::Zero());
        }

        for (int i = 0; i < IA.size(); ++i) {
            IA[i].resize(ndirs, SpatialMatrix::Zero());
        }

        // for(int i = 0; i < u.size(); ++i) {
        //  u[i].resize(ndirs, 0.0);
        // }

        // for (int i = 0; i < d.size(); ++i) {
        //  d[i].resize(ndirs, 0.0);
        // }

        u.resize(u.rows(), ndirs);
        d.resize(d.rows(), ndirs);
    }
}
