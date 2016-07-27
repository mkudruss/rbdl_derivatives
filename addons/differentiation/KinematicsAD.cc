/*
 * RBDL - Rigid Body Dynamics Library
 * Copyright (c) 2011-2015 Martin Felis <martin.felis@iwr.uni-heidelberg.de>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#include <iostream>
#include <limits>
#include <cstring>
#include <assert.h>

#include "rbdl/rbdl_mathutils.h"
#include "rbdl/Logging.h"

#include "rbdl/Model.h"
#include "rbdl/Kinematics.h"

#include "JointAD.h"
#include "ModelAD.h"
#include "KinematicsAD.h"

#include "rbdl_mathutilsAD.h"
#include "rbdl_mathutilsFD.h"

using std::cerr;
using std::endl;
using std::vector;

using namespace RigidBodyDynamics::Math;

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------

RBDL_DLLAPI Vector3d CalcBodyToBaseCoordinatesSingleFunc (
        Model &model,
        const VectorNd &Q,
        unsigned int body_id,
        const Vector3d &point_body_coordinates) {
    if (body_id >= model.fixed_body_discriminator) {
        std::cerr << "Fixed bodies not yet supported!" << std::endl;
        abort();
    }

    // Update the kinematics
    VectorNd QDot_zero (VectorNd::Zero (model.q_size));

    for (unsigned int i = 1; i < model.mBodies.size(); i++) {
        unsigned int lambda = model.lambda[i];

        // Calculate joint dependent variables
        if (model.mJoints[i].mJointType == JointTypeRevoluteX) {
            model.X_J[i] = Xrotx (Q[model.mJoints[i].q_index]);
        } else if (model.mJoints[i].mJointType == JointTypeRevoluteY) {
            model.X_J[i] = Xroty (Q[model.mJoints[i].q_index]);
        } else if (model.mJoints[i].mJointType == JointTypeRevoluteZ) {
            model.X_J[i] = Xrotz (Q[model.mJoints[i].q_index]);
        } else if (model.mJoints[i].mDoFCount == 1) {
            model.X_J[i] = Xtrans (model.S[i].block<3,1>(3,0) * Q[model.mJoints[i].q_index]);
        } else {
            std::cerr << "Unsupported joint! Only RotX, RotY, RotZ and TransX, TransY, TransZ supported!" << std::endl;
            abort();
        }

        model.X_lambda[i] = model.X_J[i] * model.X_T[i];
        model.X_base[i] = model.X_lambda[i] * model.X_base[lambda];
    }

    Matrix3d body_rotation = Math::AD::E_from_Matrix(model.X_base[body_id].toMatrix());
    Vector3d body_position = Math::AD::r_from_Matrix(model.X_base[body_id].toMatrix());

    return body_position + body_rotation.transpose() * point_body_coordinates;
}

// -----------------------------------------------------------------------------
namespace AD {
// -----------------------------------------------------------------------------

RBDL_DLLAPI Vector3d CalcBodyToBaseCoordinatesSingleFunc (
        Model &model,
        ADModel &ad_model,
        const VectorNd &q,
        const MatrixNd &q_dirs,
        unsigned int body_id,
        const Vector3d &point_body_coordinates,
        MatrixNd & out) {
    if (body_id >= model.fixed_body_discriminator) {
        std::cerr << "Fixed bodies not yet supported!" << std::endl;
        abort();
    }

    assert (out.rows() == 3 && out.cols() == model.qdot_size);
    unsigned int ndirs = q_dirs.cols();

    // Update the kinematics
    VectorNd QDot_zero (VectorNd::Zero (model.q_size));
    MatrixNd fd_out (MatrixNd::Zero (3, model.q_size));

    ad_model.resize_directions(ndirs);

    std::vector<MatrixNd> ad_X_J_i (ndirs, MatrixNd::Zero (6,6));
    std::vector<std::vector<MatrixNd> > fd_X_J (model.mBodies.size(), ad_X_J_i);
    // ad_X_J[3][5] gives for body 3 the 5th direction

    std::vector<MatrixNd> ad_X_lambda_i (ndirs, MatrixNd::Zero (6,6));
    std::vector<std::vector<MatrixNd> > fd_X_lambda (model.mBodies.size(), ad_X_lambda_i);
    // ad_X_lambda[3][5] gives for body 3 the 5th direction

    std::vector<MatrixNd> ad_X_base_i (ndirs, MatrixNd::Zero (6,6));
    std::vector<std::vector<MatrixNd> > fd_X_base (model.mBodies.size(), ad_X_base_i);
    // ad_X_base[3][5] gives for body 3 the 5th direction

    for (unsigned int i = 1; i < model.mBodies.size(); i++) {
        unsigned int lambda = model.lambda[i];
        // Calculate joint dependent variables

        if (model.mJoints[i].mJointType == JointTypeRevoluteX) {
            for (unsigned int j = 0; j < ndirs; j++) {
                ad_model.X_J[i][j] = Math::AD::Xrotx (q[model.mJoints[i].q_index], q_dirs(i-1,j));
            }
            model.X_J[i] = Xrotx (q[model.mJoints[i].q_index]);
        } else if (model.mJoints[i].mJointType == JointTypeRevoluteY) {
            for (unsigned int j = 0; j < ndirs; j++) {
                ad_model.X_J[i][j] = Math::AD::Xroty (q[model.mJoints[i].q_index], q_dirs(i-1,j));
            }
            model.X_J[i] = Xroty (q[model.mJoints[i].q_index]);
        } else if (model.mJoints[i].mJointType == JointTypeRevoluteZ) {
            for (unsigned int j = 0; j < ndirs; j++) {
                ad_model.X_J[i][j] = Math::AD::Xrotz (q[model.mJoints[i].q_index], q_dirs(i-1,j));
            }
            model.X_J[i] = Xrotz (q[model.mJoints[i].q_index]);
        } else if (model.mJoints[i].mDoFCount == 1) {
            for (unsigned int j = 0; j < ndirs; j++) {
                ad_model.X_J[i][j] = Math::AD::Xtrans (
                                        model.S[i].block<3,1>(3,0) * q[model.mJoints[i].q_index],
                                        model.S[i].block<3,1>(3,0) * q_dirs(i-1, j));
            }
            model.X_J[i] = Xtrans (model.S[i].block<3,1>(3,0) * q[model.mJoints[i].q_index]);
//            } else if (model.S[i] == SpatialVector (0., 0., 0., 1., 0., 0.)) {
//            for (unsigned int j = 0; j < ndirs; j++) {
//                ad_model.X_J[i][j] = Math::AD::Xtrans (
//                    Vector3d (1., 0., 0.) * q[model.mJoints[i].q_index],
//                    Vector3d (q_dirs(i-1, j), 0., 0.)
//                );
//            }
//            model.X_J[i] = Xtrans (Vector3d (1., 0., 0.) * q[model.mJoints[i].q_index]);
        } else {
            std::cerr << "Unsupported joint! Only RotX, RotY, RotZ and TransX supported!" << std::endl;
            abort();
        }

        for (unsigned int j = 0; j < ndirs; j++) {
            ad_model.X_lambda[i][j] = ad_model.X_J[i][j] * model.X_T[i].toMatrix();
        }
        model.X_lambda[i] = model.X_J[i] * model.X_T[i];

        for (unsigned int j = 0; j < ndirs; j++) {
            ad_model.X_base[i][j] = ad_model.X_lambda[i][j] * model.X_base[lambda].toMatrix() + model.X_lambda[i].toMatrix() * ad_model.X_base[lambda][j];
        }
        model.X_base[i] = model.X_lambda[i] * model.X_base[lambda];
    }

    for (unsigned int j = 0; j < ndirs; j++) {
        SpatialMatrix X_base_ib = model.X_base[body_id].toMatrix();
        Matrix3d ad_E = Math::AD::E_from_Matrix(ad_model.X_base[body_id][j]);
        Vector3d ad_r = Math::AD::r_from_Matrix(X_base_ib, ad_model.X_base[body_id][j]);

        out.block<3,1>(0,j) = ad_r + ad_E.transpose() * point_body_coordinates;
    }

    Matrix3d body_rotation = Math::AD::E_from_Matrix(model.X_base[body_id].toMatrix());
    Vector3d body_position = Math::AD::r_from_Matrix(model.X_base[body_id].toMatrix());

    return body_position + body_rotation.transpose() * point_body_coordinates;
}

RBDL_DLLAPI
Vector3d CalcBodyToBaseCoordinates (
        Model &model,
        ADModel &ad_model,
        const VectorNd &Q,
        const MatrixNd &Q_dirs,
        unsigned int body_id,
        const Vector3d &point_body_coordinates,
        std::vector<Math::Vector3d> *ad_body_to_base_coordinates,
        bool update_kinematics
) {
    unsigned int ndirs = Q_dirs.cols();

    std::vector<Matrix3d> ad_body_rotation (ndirs, Matrix3d::Zero());
    std::vector<Vector3d> ad_body_position (ndirs, Vector3d::Zero());

    // update the Kinematics if necessary
    if (update_kinematics) {
        // derivative evaluation
        AD::UpdateKinematicsCustom (
            model, ad_model,
            &Q, &Q_dirs,
            NULL, NULL,
            NULL, NULL
        );
        // nominal evaluation
        // NOTE: nominal evaluation is already done in AD::UpdateKinematicsCustom
        // UpdateKinematicsCustom (model, &Q, NULL, NULL);
    }

    if (body_id >= model.fixed_body_discriminator) {
        cerr << "fixed bodies not supported yet." << endl;
        abort();

        unsigned int fbody_id = body_id - model.fixed_body_discriminator;
        unsigned int parent_id = model.mFixedBodies[fbody_id].mMovableParent;

        // derivative evaluation
        for (unsigned int idirs = 0; idirs < ndirs; ++idirs) {
        // Matrix3d fixed_rotation = model.mFixedBodies[fbody_id].mParentTransform.E.transpose();
        // Vector3d fixed_position = model.mFixedBodies[fbody_id].mParentTransform.r;
        }
        // nominal evaluation
        Matrix3d fixed_rotation = model.mFixedBodies[fbody_id].mParentTransform.E.transpose();
        Vector3d fixed_position = model.mFixedBodies[fbody_id].mParentTransform.r;

        // derivative evaluation
        for (unsigned int idirs = 0; idirs < ndirs; ++idirs) {
        // Matrix3d parent_body_rotation = model.X_base[parent_id].E.transpose();
        // Vector3d parent_body_position = model.X_base[parent_id].r;
        }
        // nominal evaluation
        Matrix3d parent_body_rotation = model.X_base[parent_id].E.transpose();
        Vector3d parent_body_position = model.X_base[parent_id].r;

        // derivative evaluation
        for (unsigned int idirs = 0; idirs < ndirs; ++idirs) {
            // *(ad_body_to_base_coordinates)[idirs] =
            // parent_body_position + parent_body_rotation * (fixed_position + fixed_rotation * (point_body_coordinates));
        }
        // nominal evaluation
        return parent_body_position + parent_body_rotation * (fixed_position + fixed_rotation * (point_body_coordinates));
    }

    // derivative evaluation
    for (unsigned int idirs = 0; idirs < ndirs; ++idirs) {
        ad_body_rotation[idirs] = Math::AD::E_from_Matrix(ad_model.X_base[body_id][idirs]).transpose();
        ad_body_position[idirs] = Math::AD::r_from_Matrix(
            model.X_base[body_id].toMatrix(), ad_model.X_base[body_id][idirs]
        );

        // NOTE point_body_coordinates is a constant, no derivative needed!
        (*ad_body_to_base_coordinates)[idirs] = ad_body_position[idirs]
             + ad_body_rotation[idirs] * point_body_coordinates;
    }

    // nominal evaluation
    Matrix3d body_rotation = model.X_base[body_id].E.transpose();
    Vector3d body_position = model.X_base[body_id].r;
    return body_position + body_rotation * point_body_coordinates;
}


RBDL_DLLAPI void UpdateKinematicsCustom (
        Model &model,
        ADModel & ad_model,
        const VectorNd * q,
        const MatrixNd * q_dirs,
        const VectorNd * qdot,
        const MatrixNd * qdot_dirs,
        const VectorNd * qddot,
        const MatrixNd * qddot_dirs) {
    LOG << "-------- " << __func__ << " --------" << std::endl;
    unsigned int i;
    unsigned int ndirs = 0;

    if (q) {
        ndirs = q_dirs->cols();
    }

    ad_model.resize_directions(ndirs);

    SpatialVector spatial_gravity (0., 0., 0., model.gravity[0], model.gravity[1], model.gravity[2]);

    // derivative evaluation
    for (unsigned int idirs = 0; idirs < ndirs; ++idirs) {
        ad_model.a[0][idirs].setZero();
    }
    // nominal evaluation
    model.a[0].setZero();

    if (q) {
        for (i = 1; i < model.mBodies.size(); i++) {
            unsigned int lambda = model.lambda[i];

            VectorNd QDot_zero (VectorNd::Zero (model.q_size));
            MatrixNd QDot_zero_dirs (MatrixNd::Zero (model.q_size, ndirs));

            // Derivative evaluation and nominal evaluation
            jcalc (model, ad_model, i, *q, *q_dirs, QDot_zero, QDot_zero_dirs);
            // derivative evaluation
            for (unsigned int idirs = 0; idirs < ndirs; ++idirs) {
                // NOTE: X_T is a constant model dependent transformation
                ad_model.X_lambda[i][idirs] = ad_model.X_J[i][idirs] * model.X_T[i].toMatrix();
            }
            // nominal evaluation
            model.X_lambda[i] = model.X_J[i] * model.X_T[i];

            if (lambda != 0) {
                // derivative evaluation
                for (unsigned int idirs = 0; idirs < ndirs; ++idirs) {
                    ad_model.X_base[i][idirs] =
                        ad_model.X_lambda[i][idirs] * model.X_base[lambda].toMatrix()
                        + model.X_lambda[i].toMatrix() * ad_model.X_base[lambda][idirs];
                }
                // nominal evaluation
                model.X_base[i] = model.X_lambda[i] * model.X_base[lambda];
            } else {
                // derivative evaluation
                for (unsigned int idirs = 0; idirs < ndirs; ++idirs) {
                    ad_model.X_base[i][idirs] = ad_model.X_lambda[i][idirs];
                }
                // nominal evaluation
                model.X_base[i] = model.X_lambda[i];
            }
        }
    }

    if (qdot) {
        for (i = 1; i < model.mBodies.size(); i++) {
            unsigned int lambda = model.lambda[i];

            // Derivative evaluation and nominal evaluation
            jcalc (model, ad_model, i, *q, *q_dirs, *qdot, *qdot_dirs);

            if (lambda != 0) {
                // derivative evaluation
                for (unsigned int idirs = 0; idirs < ndirs; ++idirs) {
                    ad_model.v[i][idirs] =
                        ad_model.X_lambda[i][idirs] * model.v[lambda]
                        + model.X_lambda[i].apply(ad_model.v[lambda][idirs])
                        + ad_model.v_J[i][idirs];
                }
                // nominal evaluation
                model.v[i] = model.X_lambda[i].apply(model.v[lambda]) + model.v_J[i];
            } else {
                // derivative evaluation
                for (unsigned int idirs = 0; idirs < ndirs; ++idirs) {
                    ad_model.v[i][idirs] = ad_model.v_J[i][idirs];
                }
                // nominal evaluation
                model.v[i] = model.v_J[i];
            }
            // derivative evaluation
            for (unsigned int idirs = 0; idirs < ndirs; ++idirs) {
                ad_model.c[i][idirs] = ad_model.c_J[i][idirs]
                    + Math::AD::crossm(
                        model.v[i], ad_model.v[i][idirs],
                        model.v_J[i], ad_model.v_J[i][idirs]);
            }
            // nominal evaluation
            model.c[i] = model.c_J[i] + crossm(model.v[i],model.v_J[i]);
        }
    }

    if (qddot) {
        for (i = 1; i < model.mBodies.size(); i++) {
            unsigned int q_index = model.mJoints[i].q_index;
            unsigned int lambda = model.lambda[i];

            if (lambda != 0) {
                // derivative evaluation
                for (unsigned int idirs = 0; idirs < ndirs; ++idirs) {
                    ad_model.a[i][idirs] = ad_model.X_lambda[i][idirs] * model.a[lambda]
                        + model.X_lambda[i].apply(ad_model.a[lambda][idirs]) + ad_model.c[i][idirs];
                }
                // nominal evaluation
                model.a[i] = model.X_lambda[i].apply(model.a[lambda]) + model.c[i];
            } else {
                // derivative evaluation
                for (unsigned int idirs = 0; idirs < ndirs; ++idirs) {
                    ad_model.a[i][idirs] = ad_model.c[i][idirs];
                }
                // nominal evaluation
                model.a[i] = model.c[i];
            }

            if (model.mJoints[i].mDoFCount == 3) {
                cerr << "Multi-DoF not supported." << endl;
                abort();
                /*
                // nominal evaluation
                Vector3d omegadot_temp ((*qddot)[q_index], (*qddot)[q_index + 1], (*qddot)[q_index + 2]);
                model.a[i] = model.a[i] + model.multdof3_S[i] * omegadot_temp;
                */
            } else {
                // derivative evaluation
                for (unsigned int idirs = 0; idirs < ndirs; ++idirs) {
                    ad_model.a[i][idirs] = ad_model.a[i][idirs]
                        + ad_model.S[i][idirs] * (*qddot)[q_index]
                        + model.S[i] * (*qddot_dirs)(q_index, idirs);
                }
                // nominal evaluation
                model.a[i] = model.a[i] + model.S[i] * (*qddot)[q_index];
            }
        }
    }
}

RBDL_DLLAPI void UpdateKinematics (
    Model &model,
    ADModel &ad_model,
    const VectorNd &q,
    const MatrixNd &q_dirs,
    const VectorNd &qdot,
    const MatrixNd &qdot_dirs,
    const VectorNd &qddot,
    const MatrixNd &qddot_dirs) {
    LOG << "-------- " << __func__ << " --------" << std::endl;

    unsigned int i;
    unsigned int ndirs = q_dirs.cols();

    ad_model.resize_directions(ndirs);

    SpatialVector spatial_gravity (0., 0., 0., model.gravity[0], model.gravity[1], model.gravity[2]);

    // derivative evaluation
    for (unsigned int idirs = 0; idirs < ndirs; ++idirs) {
        ad_model.a[0][idirs].setZero();
    }
    // nominal evaluation
    model.a[0].setZero();
    //model.a[0] = spatial_gravity;

    for (i = 1; i < model.mBodies.size(); i++) {
        unsigned int q_index = model.mJoints[i].q_index;

        Joint joint = model.mJoints[i];
        unsigned int lambda = model.lambda[i];

        // derivative evaluation
        jcalc (model, ad_model, i, q, q_dirs, qdot, qdot_dirs);
        // nominal evaluation
        // NOTE: nominal evaluation is not needed, because already done in ad_jcalc
        //jcalc (model, i, q, qdot);

        // derivative evaluation
        for (unsigned int idirs = 0; idirs < ndirs; ++idirs) {
            // NOTE: X_T is a constant model dependent transformation
            ad_model.X_lambda[i][idirs] = ad_model.X_J[i][idirs] * model.X_T[i].toMatrix();
        }
        // nominal evaluation
        model.X_lambda[i] = model.X_J[i] * model.X_T[i];

        if (lambda != 0) {
            // derivative evaluation
            for (unsigned int idirs = 0; idirs < ndirs; ++idirs) {
                ad_model.X_base[i][idirs] =
                    ad_model.X_lambda[i][idirs] * model.X_base[lambda].toMatrix()
                    + model.X_lambda[i].toMatrix() * ad_model.X_base[lambda][idirs];
                ad_model.v[i][idirs] =
                    ad_model.X_lambda[i][idirs] * model.v[lambda]
                    + model.X_lambda[i].apply(ad_model.v[lambda][idirs])
                    + ad_model.v_J[i][idirs];
            }
            // nominal evaluation
            model.X_base[i] = model.X_lambda[i] * model.X_base[lambda];
            model.v[i] = model.X_lambda[i].apply(model.v[lambda]) + model.v_J[i];
        } else {
            // derivative evaluation
            for (unsigned int idirs = 0; idirs < ndirs; ++idirs) {
                ad_model.X_base[i][idirs] = ad_model.X_lambda[i][idirs];
                ad_model.v[i][idirs] = ad_model.v_J[i][idirs];
            }
            // nominal evaluation
            model.X_base[i] = model.X_lambda[i];
            model.v[i] = model.v_J[i];
        }

        // derivative evaluation
        for (unsigned int idirs = 0; idirs < ndirs; ++idirs) {
            ad_model.c[i][idirs] = ad_model.c_J[i][idirs]
                + Math::AD::crossm(
                    model.v[i], ad_model.v[i][idirs],
                    model.v_J[i], ad_model.v_J[i][idirs]
                );
            ad_model.a[i][idirs] = ad_model.X_lambda[i][idirs] * model.a[lambda]
                + model.X_lambda[i].apply(ad_model.a[lambda][idirs]) + ad_model.c[i][idirs];
        }
        // nominal evaluation
        model.c[i] = model.c_J[i] + crossm(model.v[i],model.v_J[i]);
        model.a[i] = model.X_lambda[i].apply(model.a[lambda]) + model.c[i];

        if (model.mJoints[i].mDoFCount == 3) {
            cerr << "Multi-DoF not supported." << endl;
            abort();
            /*
            // derivative evaluation
            for (int idirs = 0; idirs < ndirs; ++idirs) {
            }
            // nominal evaluation
            Vector3d omegadot_temp (qddot[q_index], qddot[q_index + 1], qddot[q_index + 2]);
            model.a[i] = model.a[i] + model.multdof3_S[i] * omegadot_temp;
            */
        } else {
            // derivative evaluation
            for (unsigned int idirs = 0; idirs < ndirs; ++idirs) {
                ad_model.a[i][idirs] = ad_model.a[i][idirs]
                    + ad_model.S[i][idirs] * qddot[q_index]
                    + model.S[i] * qddot_dirs(q_index, idirs);
            }
            // nominal evaluation
            model.a[i] = model.a[i] + model.S[i] * qddot[q_index];
        }
    }

    for (i = 1; i < model.mBodies.size(); i++) {
        // derivative evaluation
        for (unsigned int idirs = 0; idirs < ndirs; ++idirs) {
            LOG << "ad_a[" << i << "][" << idirs << "] = " << ad_model.a[i][idirs].transpose() << std::endl;
        }
        // nominal evaluation
        LOG << "a[" << i << "] = " << model.a[i].transpose() << std::endl;
    }
}


RBDL_DLLAPI Matrix3d CalcBodyWorldOrientation (
        Model & model,
        ADModel & ad_model,
        VectorNd const & q,
        MatrixNd const & q_dirs,
        const unsigned int body_id,
        vector<Matrix3d> & ad_derivative,
        bool update_kinematics) {
    unsigned int ndirs = q_dirs.cols();
    assert(ad_derivative.size() == ndirs);

    // update the Kinematics if necessary
    if (update_kinematics) {
        UpdateKinematicsCustom (model, ad_model, &q, &q_dirs, 0, 0, 0, 0);
    }

    if (body_id >= model.fixed_body_discriminator) {
        std::cerr << "Fixed bodies not yet supported!" << std::endl;
        abort();
    }

    for (unsigned int idir = 0; idir < ndirs; idir++) {
        ad_derivative[idir] = ad_model.X_base[body_id][idir].block<3,3>(0,0);
    }
    // nominal evaluation
    return model.X_base[body_id].E;
}


RBDL_DLLAPI Vector3d CalcPointAcceleration (
        Model &model,
        ADModel &ad_model,
        const Math::VectorNd &q,
        const Math::MatrixNd &q_dirs,
        const Math::VectorNd &qdot,
        const Math::MatrixNd &qdot_dirs,
        const Math::VectorNd &qddot,
        const Math::MatrixNd &qddot_dirs,
        unsigned int body_id,
        const Math::Vector3d &point_position,
        Math::MatrixNd &ad_derivative,
        bool update_kinematics
) {
    LOG << "-------- " << __func__ << " --------" << std::endl;

    int const ndirs = q_dirs.cols();
    assert(ndirs == ad_derivative.cols());
    assert(3     == ad_derivative.rows());

    ad_model.resize_directions(ndirs);

    // Reset the velocity of the root body
    model.v[0].setZero();
    model.a[0].setZero();

    if (update_kinematics) {
        // derivative evaluation
        UpdateKinematics (model, ad_model, q, q_dirs, qdot, qdot_dirs, qddot, qddot_dirs);
        // nominal evaluation
        // NOTE: Kinematics are already updated in ad_UpdateKinematics
        //UpdateKinematics (model, q, qdot, qddot);
    }

    LOG << std::endl;

    unsigned int reference_body_id = body_id;
    Vector3d reference_point = point_position;

    if (model.IsFixedBodyId(body_id)) {
        std::cerr << "Fixed bodies not yet supported!" << std::endl;
        abort();
    }

    vector<Matrix3d> ad_E(ndirs);
    Matrix3d E = CalcBodyWorldOrientation(model, ad_model, q, q_dirs, reference_body_id, ad_E, false);

    // derivative evaluation
    vector<SpatialMatrix> ad_p_X_i(ndirs);
    for (int idir = 0; idir < ndirs; idir++) {
        Matrix3d ad_ETrx = ad_E[idir].transpose() * Matrix3d(
                              0., -reference_point[2],  reference_point[1],
              reference_point[2],                  0., -reference_point[0],
             -reference_point[1],  reference_point[0],                  0.);
             // + E.transpose() * Matrix3dZero;
        ad_p_X_i[idir].block<3,3>(0, 0) = ad_E[idir].transpose();
        ad_p_X_i[idir].block<3,3>(0, 3) = Matrix3dZero;
        ad_p_X_i[idir].block<3,3>(3, 0) = -ad_ETrx;
        ad_p_X_i[idir].block<3,3>(3, 3) = ad_E[idir].transpose();
    }
    // nominal evaluation
    SpatialTransform p_X_i (E.transpose(), reference_point);

    // derivative evaluation
    vector<SpatialVector> ad_p_v_i(ndirs);
    for (int idir = 0; idir < ndirs; idir++) {
        ad_p_v_i[idir] = ad_p_X_i[idir] * model.v[reference_body_id]
                        + p_X_i.apply(ad_model.v[reference_body_id][idir]);
    }
    // nominal evaluation
    SpatialVector p_v_i = p_X_i.apply(model.v[reference_body_id]);

    Vector3d p_v_i_0t2(p_v_i[0], p_v_i[1], p_v_i[2]);
    Vector3d p_v_i_3t5(p_v_i[3], p_v_i[4], p_v_i[5]);
    // derivative evaluation
    vector<Vector3d> ad_a_dash(ndirs);
    for (int idir = 0; idir < ndirs; idir++) {
        Vector3d ad_p_v_i_0t2_idir = ad_p_v_i[idir].block<3,1>(0, 0);
        Vector3d ad_p_v_i_3t5_idir = ad_p_v_i[idir].block<3,1>(3, 0);
		ad_a_dash[idir] = p_v_i_0t2.cross(ad_p_v_i_3t5_idir)
						  + ad_p_v_i_0t2_idir.cross(p_v_i_3t5);
    }
    // nominal evaluation
    Vector3d a_dash = p_v_i_0t2.cross(p_v_i_3t5);

    // derivative evaluation
    vector<SpatialVector> ad_p_a_i(ndirs);
    for (int idir = 0; idir < ndirs; idir++) {
        ad_p_a_i[idir] = ad_p_X_i[idir] * model.a[reference_body_id]
                + p_X_i.apply(ad_model.a[reference_body_id][idir]);
    }
    // nominal evaluation
    SpatialVector p_a_i = p_X_i.apply(model.a[reference_body_id]);

    // derivative evaluation
    for (int idir = 0; idir < ndirs; idir++) {
        ad_derivative(0, idir) = ad_p_a_i[idir][3] + ad_a_dash[idir][0];
        ad_derivative(1, idir) = ad_p_a_i[idir][4] + ad_a_dash[idir][1];
        ad_derivative(2, idir) = ad_p_a_i[idir][5] + ad_a_dash[idir][2];
    }
    // nominal evaluation
    return Vector3d (
            p_a_i[3] + a_dash[0],
            p_a_i[4] + a_dash[1],
            p_a_i[5] + a_dash[2]);
}

RBDL_DLLAPI
void CalcPointJacobian (
        Model &model,
        ADModel &ad_model,
        Math::VectorNd const &Q,
        Math::MatrixNd const &Q_dirs,
        unsigned int body_id,
        Math::Vector3d const &point_position,
        Math::MatrixNd &G,
        std::vector<Math::MatrixNd> &G_dirs,
        bool update_kinematics
) {
    LOG << "-------- " << __func__ << " --------" << std::endl;

    unsigned int ndirs = Q_dirs.cols();
    ad_model.resize_directions(ndirs);

    // update the Kinematics if necessary
    if (update_kinematics) {
        // derivative evaluation
        UpdateKinematicsCustom (
            model, ad_model,
            &Q, &Q_dirs,
            NULL, NULL,
            NULL, NULL
        );
        // nominal evaluation
        // NOTE kinematics are already updated in the AD version
        // UpdateKinematicsCustom (model, &Q, NULL, NULL);
    }

    std::vector<SpatialMatrix> ad_point_trans (ndirs);
    // derivative evaluation
    for (unsigned int idirs = 0; idirs < ndirs; idirs++) {
        // NOTE point_trans is spatial transform from E, r, i.e,
        //
        //      X = [   E 0], E = I
        //          [-Erx E]
        //
        //      for the derivative the blocks E vanish as they are constant
        //      and for the dot product the product rule applies, i.e.,
        //
        //      X_ad = [                   0 0]
        //             [-E_dot*rx - E*rx_dot 0]
        //
        // SpatialMatrix& ad_point_trans_i = ad_point_trans[idirs];
        // Matrix3d Erx = Matrix3d::Zero();
        // Erx = X.block<3,3>(3,0);
        // = SpatialTransform (
        //     Matrix3d::Identity(),
        //     CalcBodyToBaseCoordinates (
        //         model, Q, body_id, point_position, false
        //     )
        // );
    }
    // nominal evaluation
    SpatialTransform point_trans = SpatialTransform (
            Matrix3d::Identity(),
            CalcBodyToBaseCoordinates (
                model, Q, body_id, point_position, false
            )
    );

    assert (G.rows() == 3 && G.cols() == model.qdot_size );

    unsigned int reference_body_id = body_id;

    if (model.IsFixedBodyId(body_id)) {
        unsigned int fbody_id = body_id - model.fixed_body_discriminator;
        reference_body_id = model.mFixedBodies[fbody_id].mMovableParent;
    }

    unsigned int j = reference_body_id;

    // e[j] is set to 1 if joint j contributes to the jacobian that we are
    // computing. For all other joints the column will be zero.
    while (j != 0) {
        unsigned int q_index = model.mJoints[j].q_index;

        if (model.mJoints[j].mDoFCount == 3) {
            std::cout << "3DoF joints are not yet supported!" << std::endl;
            std::cout << "bailing out ..." << std::endl;
            abort();
            G.block(0, q_index, 3, 3) = ((point_trans * model.X_base[j].inverse()).toMatrix() * model.multdof3_S[j]).block(3,0,3,3);
        } else {
            // derivative evaluation
            // nominal evaluation
            G.block(0, q_index, 3, 1) = point_trans.apply(
                model.X_base[j].inverse().apply(model.S[j])
            ).block(3,0,3,1);
        }

        j = model.lambda[j];
    }
}

// -----------------------------------------------------------------------------
} // namespace AD
// -----------------------------------------------------------------------------
} // namespace RigidBodyDynamics
// -----------------------------------------------------------------------------
