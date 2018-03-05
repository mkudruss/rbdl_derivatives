/*
 * RBDL - Rigid Body Dynamics Library
 * Copyright (c) 2011-2015 Martin Felis <martin.felis@iwr.uni-heidelberg.de>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#include "ModelED.h"

using namespace RigidBodyDynamics::Math;

// std::vectors containing any objectst that have Eigen matrices or vectors
// as members need to have a special allocater. This can be achieved with
// the following macro.

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------
// namespace ED {
// -----------------------------------------------------------------------------

EDModel::EDModel () {
    // TODO: should not be called, I think!
    //std::cerr << "Data structure can only be initialized with argument RigidBodyDynamics::Model!" << std::endl;
    //throw std::exception;
}

EDModel::EDModel (RigidBodyDynamics::Model& model) {
    ndirs = 4 * model.dof_count;

    // NOTE: old initialization values
    std::vector<SpatialRigidBodyInertia> SRBI(ndirs);
    std::vector<SpatialTransform> T(ndirs, SpatialTransform::Zero());
    std::vector<SpatialTransform> Tdot(ndirs, SpatialTransform::Zero());
    std::vector<SpatialMatrix> X(ndirs, SpatialMatrix::Zero());
    std::vector<SpatialVector> vec(ndirs, SpatialVector::Zero());
    //std::vector<double> double(ndirs, 0.0);

    X_lambda.resize(model.mBodies.size(), T);
    X_base.resize(model.mBodies.size(), T);
    X_J.resize(model.mBodies.size(), T);

    v_J.resize(model.mBodies.size(), vec);
    v.resize(model.mBodies.size(), vec);
    c_J.resize(model.mBodies.size(), vec);
    a.resize(model.mBodies.size(), vec);

    Ic.resize(model.mBodies.size(), SRBI);
    hc.resize(model.mBodies.size(), vec);
    F.resize(model.mBodies.size(), vec);
    c.resize(model.mBodies.size(), vec);
    f.resize(model.mBodies.size(), vec);

    pA.resize(model.mBodies.size(), vec);
    U.resize(model.mBodies.size(), vec);
    IA.resize(model.mBodies.size(), X);
    //u.resize(model.mBodies.size(), double);
    //d.resize(model.mBodies.size(), double);
    u.resize(model.u.rows(), ndirs);
    u.setZero();
    d.resize(model.d.rows(), ndirs);
    d.setZero();
}

void EDModel::resize_directions (unsigned int requested_ndirs){
    if (ndirs < requested_ndirs) {
        ndirs = requested_ndirs;

        for (unsigned int i = 0; i < X_lambda.size(); ++i) {
            X_lambda[i].resize(ndirs, SpatialTransform::Zero());
        }

        for (unsigned int i = 0; i < X_base.size(); ++i) {
            X_base[i].resize(ndirs, SpatialTransform::Zero());
        }

        for (unsigned int i = 0; i < X_J.size(); ++i) {
            X_J[i].resize(ndirs, SpatialTransform::Zero());
        }

        for (unsigned int i = 0; i < v_J.size(); ++i) {
            v_J[i].resize(ndirs, SpatialVector::Zero());
        }

        for (unsigned int i = 0; i < v_J.size(); ++i) {
            v[i].resize(ndirs, SpatialVector::Zero());
        }

        for (unsigned int i = 0; i < c_J.size(); ++i) {
            c_J[i].resize(ndirs, SpatialVector::Zero());
        }

        for (unsigned int i = 0; i < c_J.size(); ++i) {
            a[i].resize(ndirs, SpatialVector::Zero());
        }

        for (unsigned int i = 0; i < c_J.size(); ++i) {
            c[i].resize(ndirs, SpatialVector::Zero());
        }

        for (unsigned int i = 0; i < c_J.size(); ++i) {
            f[i].resize(ndirs, SpatialVector::Zero());
        }

        for (unsigned int i = 0; i < Ic.size(); ++i) {
            Ic[i].resize(ndirs);
        }

        for (unsigned int i = 0; i < F.size(); ++i) {
            F[i].resize(ndirs, SpatialVector::Zero());
        }

        for(unsigned int i = 0; i < pA.size(); ++i) {
            pA[i].resize(ndirs, SpatialVector::Zero());
        }

        for(unsigned int i = 0; i < U.size(); ++i) {
            U[i].resize(ndirs, SpatialVector::Zero());
        }

        for (unsigned int i = 0; i < IA.size(); ++i) {
            IA[i].resize(ndirs, SpatialMatrix::Zero());
        }

        for (unsigned int i = 0; i < hc.size(); ++i) {
            hc[i].resize(ndirs, SpatialVector::Zero());
        }

        u.resize(u.rows(), ndirs);
        u.setZero();
        d.resize(d.rows(), ndirs);
        d.setZero();
    }
}

// -----------------------------------------------------------------------------
// } // ED
// -----------------------------------------------------------------------------
} // RigidBodyDynamics
// -----------------------------------------------------------------------------
