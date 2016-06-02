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
    VectorNd fd_out (MatrixNd::Zero (3, model.q_size));

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
        const Math::MatrixNd &fd_derivative,
        bool update_kinematics = true) {
    LOG << "-------- " << __func__ << " --------" << std::endl;

    unsigned int ndirs = q_dirs.cols();
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
        unsigned int fbody_id = body_id - model.fixed_body_discriminator;
        reference_body_id = model.mFixedBodies[fbody_id].mMovableParent;
        Vector3d base_coords = CalcBodyToBaseCoordinates (model, q, body_id, point_position, false);
        reference_point = CalcBaseToBodyCoordinates (model, q, reference_body_id, base_coords, false);
    }

    SpatialTransform p_X_i (CalcBodyWorldOrientation (model, q, reference_body_id, false).transpose(), reference_point);

    SpatialVector p_v_i = p_X_i.apply(model.v[reference_body_id]);
    Vector3d a_dash = Vector3d (p_v_i[0], p_v_i[1], p_v_i[2]).cross(Vector3d (p_v_i[3], p_v_i[4], p_v_i[5]));
    SpatialVector p_a_i = p_X_i.apply(model.a[reference_body_id]);

    return Vector3d (
            p_a_i[3] + a_dash[0],
            p_a_i[4] + a_dash[1],
            p_a_i[5] + a_dash[2]
            );
}

// -----------------------------------------------------------------------------
} // namespace AD
// -----------------------------------------------------------------------------
} // namespace RigidBodyDynamics
// -----------------------------------------------------------------------------
