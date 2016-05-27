#include <UnitTest++.h>

#include <iostream>

#include "rbdl/Logging.h"
#include "rbdl/rbdl_utils.h"
#include "rbdl/rbdl_mathutils.h"

#include "rbdl/Model.h"
#include "rbdl/Kinematics.h"
#include "rbdl/Dynamics.h"

#include "Fixtures.h"
#include "ModelAD.h"
#include "JointAD.h"
#include "rbdl_mathutilsAD.h"

using namespace std;
using namespace RigidBodyDynamics;
using namespace RigidBodyDynamics::Math;
using namespace RigidBodyDynamics::Utils;

const double TEST_PREC = 1.0e-12;

RBDL_DLLAPI
Vector3d CalcBodyToBaseCoordinatesSingleFunc (
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
        if (model.mJoints[i].mJointType == JointTypeRevoluteY) {
            model.X_J[i] = Xroty (Q[model.mJoints[i].q_index]);
        } else if (model.S[i] == SpatialVector (0., 0., 0., 1., 0., 0.)) {
            model.X_J[i] = Xtrans (Vector3d (1., 0., 0.) * Q[model.mJoints[i].q_index]);
        } else {
            std::cerr << "Unsupported joint! Only RotY and TransX supported!" << std::endl;
            abort();
        }

        model.X_lambda[i] = model.X_J[i] * model.X_T[i];
        model.X_base[i] = model.X_lambda[i] * model.X_base[lambda];
    }

    Matrix3d body_rotation = AD::E_from_Matrix(model.X_base[body_id].toMatrix());
    Vector3d body_position = AD::r_from_Matrix(model.X_base[body_id].toMatrix());

    return body_position + body_rotation.transpose() * point_body_coordinates;
}

RBDL_DLLAPI
Vector3d ad_CalcBodyToBaseCoordinatesSingleFunc (
        Model &model,
        ADModel &ad_model,
        const VectorNd &q,
        const MatrixNd &q_dirs,
        unsigned int body_id,
        const Vector3d &point_body_coordinates,
        MatrixNd &out
        ) {
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
        if (model.mJoints[i].mJointType == JointTypeRevoluteY) {
            for (unsigned int j = 0; j < ndirs; j++) {
                ad_model.X_J[i][j] = AD::Xroty (q[model.mJoints[i].q_index], q_dirs(i-1,j));
            }
            model.X_J[i] = Xroty (q[model.mJoints[i].q_index]);
        } else if (model.S[i] == SpatialVector (0., 0., 0., 1., 0., 0.)) {
            for (unsigned int j = 0; j < ndirs; j++) {
                ad_model.X_J[i][j] = AD::Xtrans (
                    Vector3d (1., 0., 0.) * q[model.mJoints[i].q_index],
                    Vector3d (q_dirs(i-1, j), 0., 0.)
                );
            }
            model.X_J[i] = Xtrans (Vector3d (1., 0., 0.) * q[model.mJoints[i].q_index]);
        } else {
            std::cerr << "Unsupported joint! Only RotY and TransX supported!" << std::endl;
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
        Matrix3d ad_E = AD::E_from_Matrix(ad_model.X_base[body_id][j]);
        Vector3d ad_r = AD::r_from_Matrix(X_base_ib, ad_model.X_base[body_id][j]);

        out.block<3,1>(0,j) = ad_r + ad_E.transpose() * point_body_coordinates;
    }

    Matrix3d body_rotation = AD::E_from_Matrix(model.X_base[body_id].toMatrix());
    Vector3d body_position = AD::r_from_Matrix(model.X_base[body_id].toMatrix());

    return body_position + body_rotation.transpose() * point_body_coordinates;
}


RBDL_DLLAPI
Vector3d fd_CalcBodyToBaseCoordinatesSingleFunc (
		Model &model,
		const VectorNd &q,
		const MatrixNd &q_dirs,
		unsigned int body_id,
		const Vector3d &point_body_coordinates,
		MatrixNd &out
		) {
	Vector3d ref = CalcBodyToBaseCoordinatesSingleFunc (model, q, body_id, point_body_coordinates);

	unsigned int ndirs = q_dirs.cols();
	double h = 1.0e-8;
	for (unsigned int j = 0; j < ndirs; j++) {
		VectorNd q_dir = q_dirs.block(0,j, model.qdot_size, 1);
		Vector3d res_hd = CalcBodyToBaseCoordinatesSingleFunc (model, q + h * q_dir, body_id, point_body_coordinates);
		Vector3d res_hd_rbdl = CalcBodyToBaseCoordinates (model, q + h * q_dir, body_id, point_body_coordinates);

		//cout << "calc_body err = " << (res_hd - res_hd_rbdl).transpose() << endl;

		out.block<3,1>(0,j) = (res_hd - ref) / h;
	}

	return ref;
}

TEST_FIXTURE ( CartPendulum, CartPendulumCalcBodyToBaseCoordinatesSingleFunc) {
	q[0] = 0.3;
	q[1] = -0.2;
	Vector3d point_body_coordinates (0.1, 3.2, 4.2);

	Vector3d point_single_func = CalcBodyToBaseCoordinatesSingleFunc (model, q, id_pendulum, point_body_coordinates);
	Vector3d point_default = CalcBodyToBaseCoordinates (model, q, id_pendulum, point_body_coordinates);

	CHECK_ARRAY_CLOSE (point_default.data(), point_single_func.data(), 3, TEST_PREC);
}

RBDL_DLLAPI void fd_CalcCenterOfMass (
        Model & model,
        const VectorNd & q,
        const MatrixNd & q_dirs,
        const VectorNd & qdot,
        const MatrixNd & qdot_dirs,
        double & mass,
        Vector3d & com,
        MatrixNd & fd_com,
        Vector3d * com_velocity,
        MatrixNd * fd_com_velocity,
        Vector3d * angular_momentum,
        MatrixNd * fd_angular_momentum) {
    assert(q_dirs.cols() == qdot_dirs.cols());
    assert(q_dirs.cols() == fd_com.cols());
    if (com_velocity) {
        assert(q_dirs.cols() == fd_com_velocity->cols());
    }
    if (angular_momentum) {
        assert(q_dirs.cols() == fd_angular_momentum->cols());
    }

    int ndirs = q_dirs.cols();

    CalcCenterOfMass(model, q, qdot, mass, com, com_velocity, angular_momentum);

    double h = sqrt(1e-16);

    VectorNd q_dir(q);
    VectorNd qdot_dir(qdot);

    for (int i = 0; i < ndirs; i++ ) {
        q_dir = q_dirs.block(0, i, model.q_size, 1);
        qdot_dir = qdot_dirs.block(0, i, model.qdot_size, 1);

        double hd_mass;
        Vector3d hd_com;
        Vector3d * hd_com_velocity = 0;
        Vector3d * hd_angular_momentum = 0;

        if (com_velocity) {
            hd_com_velocity = new Vector3d;
        }

        if (angular_momentum) {
            hd_angular_momentum = new Vector3d;
        }

        CalcCenterOfMass(model, q + h * q_dir, qdot + h * qdot_dir, hd_mass, hd_com, hd_com_velocity, hd_angular_momentum);

        fd_com.block<3,1>(0, i) = (hd_com - com) / h;

        if (com_velocity) {
            fd_com_velocity->block<3,1>(0, i) = (*hd_com_velocity - *com_velocity) / h;
            delete hd_com_velocity;
        }

        if (angular_momentum) {
            fd_angular_momentum->block<3,1>(0, i) = (*hd_angular_momentum - *angular_momentum) / h;
            delete hd_angular_momentum;
        }
    }
}

RBDL_DLLAPI void ad_UpdateKinematicsCustom (
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
    unsigned ndirs = 0;

    if (q) {
        ndirs = q_dirs->cols();
    }

    ad_model.resize_directions(ndirs);

    SpatialVector spatial_gravity (0., 0., 0., model.gravity[0], model.gravity[1], model.gravity[2]);

    // derivative evaluation
    for (int idirs = 0; idirs < ndirs; ++idirs) {
        ad_model.a[0][idirs].setZero();
    }
    // nominal evaluation
    model.a[0].setZero();

    if (q) {
        for (i = 1; i < model.mBodies.size(); i++) {
            unsigned int lambda = model.lambda[i];

            VectorNd QDot_zero (VectorNd::Zero (model.q_size));
            MatrixNd QDot_zero_dirs (MatrixNd::Zero (model.q_size, ndirs));
            cout << 1131 << endl;

            // Derivative evaluation and nominal evaluation
            ad_jcalc (model, ad_model, i, *q, *q_dirs, QDot_zero, QDot_zero_dirs);
            cout << 1132 << endl;
            // derivative evaluation
            for (int idirs = 0; idirs < ndirs; ++idirs) {
                // NOTE: X_T is a constant model dependent transformation
                ad_model.X_lambda[i][idirs] = ad_model.X_J[i][idirs] * model.X_T[i].toMatrix();
            }
            // nominal evaluation
            model.X_lambda[i] = model.X_J[i] * model.X_T[i];
            cout << 1133 << endl;

            if (lambda != 0) {
                // derivative evaluation
                for (int idirs = 0; idirs < ndirs; ++idirs) {
                    ad_model.X_base[i][idirs] =
                        ad_model.X_lambda[i][idirs] * model.X_base[lambda].toMatrix()
                        + model.X_lambda[i].toMatrix() * ad_model.X_base[lambda][idirs];
                }
                // nominal evaluation
                model.X_base[i] = model.X_lambda[i] * model.X_base[lambda];
            }	else {
                // derivative evaluation
                for (int idirs = 0; idirs < ndirs; ++idirs) {
                    ad_model.X_base[i][idirs] = ad_model.X_lambda[i][idirs];
                }
                // nominal evaluation
                model.X_base[i] = model.X_lambda[i];
            }
            cout << 1134 << endl;
        }
    }

    if (qdot) {
        for (i = 1; i < model.mBodies.size(); i++) {
            unsigned int lambda = model.lambda[i];

            // Derivative evaluation and nominal evaluation
            ad_jcalc (model, ad_model, i, *q, *q_dirs, *qdot, *qdot_dirs);

            if (lambda != 0) {
                // derivative evaluation
                for (int idirs = 0; idirs < ndirs; ++idirs) {
                    ad_model.v[i][idirs] =
                        ad_model.X_lambda[i][idirs] * model.v[lambda]
                        + model.X_lambda[i].apply(ad_model.v[lambda][idirs])
                        + ad_model.v_J[i][idirs];
                }
                // nominal evaluation
                model.v[i] = model.X_lambda[i].apply(model.v[lambda]) + model.v_J[i];
            } else {
                // derivative evaluation
                for (int idirs = 0; idirs < ndirs; ++idirs) {
                    ad_model.v[i][idirs] = ad_model.v_J[i][idirs];
                }
                // nominal evaluation
                model.v[i] = model.v_J[i];
            }
            // derivative evaluation
            for (int idirs = 0; idirs < ndirs; ++idirs) {
                ad_model.c[i][idirs] = ad_model.c_J[i][idirs]
                    + AD::crossm(
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
                for (int idirs = 0; idirs < ndirs; ++idirs) {
                    ad_model.a[i][idirs] = ad_model.X_lambda[i][idirs] * model.a[lambda]
                        + model.X_lambda[i].apply(ad_model.a[lambda][idirs]) + ad_model.c[i][idirs];
                }
                // nominal evaluation
                model.a[i] = model.X_lambda[i].apply(model.a[lambda]) + model.c[i];
            } else {
                // derivative evaluation
                for (int idirs = 0; idirs < ndirs; ++idirs) {
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
                for (int idirs = 0; idirs < ndirs; ++idirs) {
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

RBDL_DLLAPI void ad_CalcCenterOfMass (
        Model & model,
        ADModel & ad_model,
        const VectorNd & q,
        const MatrixNd & q_dirs,
        const VectorNd & qdot,
        const MatrixNd & qdot_dirs,
        double & mass,
        Vector3d & com,
        MatrixNd & ad_com,
        Vector3d * com_velocity,
        MatrixNd * ad_com_velocity,
        Vector3d * angular_momentum,
        MatrixNd * ad_angular_momentum,
        bool update_kinematics) {

    assert(q_dirs.cols() == qdot_dirs.cols());
    assert(q_dirs.cols() == ad_com.cols());
    if (com_velocity) {
        assert(q_dirs.cols() == ad_com_velocity->cols());
    }
    if (angular_momentum) {
        assert(q_dirs.cols() == ad_angular_momentum->cols());
    }

    if (update_kinematics) {
        ad_UpdateKinematicsCustom (model, ad_model, &q, &q_dirs, &qdot, &qdot_dirs, NULL, NULL);
    }

    int ndirs = q_dirs.cols();

    for (size_t i = 1; i < model.mBodies.size(); i++) {
        // derivative evaluation
        for(size_t idir = 0; idir < ndirs; idir++) {
            ad_model.Ic[i][idir] = SpatialMatrix::Zero();
        }
        // nominal evaluation
        model.Ic[i] = model.I[i];

        // derivative evaluation (can be shortened, as 2nd summand is currently always zero.
        for(size_t idir = 0; idir < ndirs; idir++) {
            ad_model.hc[i][idir] = model.Ic[i].toMatrix() * ad_model.v[i][idir]
                    + ad_model.Ic[i][idir] * model.v[i];
        }
        // nominal evaluation
        model.hc[i] = model.Ic[i].toMatrix() * model.v[i];
    }

    SpatialRigidBodyInertia Itot (0., Vector3d (0., 0., 0.), Matrix3d::Zero(3,3));
    vector<SpatialRigidBodyInertia> ad_Itot(ndirs, Itot);

    SpatialVector htot (SpatialVector::Zero(6));
    vector<SpatialVector> ad_htot(ndirs, htot);

    for (size_t i = model.mBodies.size() - 1; i > 0; i--) {
        unsigned int lambda = model.lambda[i];
        if (lambda != 0) {
            // derivative evaluation
            for (size_t idir = 0; idir < ndirs; idir++) {
                ad_model.Ic[lambda][idir] = ad_model.Ic[lambda][idir]
                        + model.X_lambda[i].toMatrixTranspose() * ad_model.Ic[i][idir] * model.X_lambda[i].toMatrix()
                        + model.X_lambda[i].toMatrixTranspose() * model.Ic[i].toMatrix() * ad_model.X_lambda[i][idir]
                        + ad_model.X_lambda[i][idir].transpose() * model.Ic[i].toMatrix() * model.X_lambda[i].toMatrix();
            }
            // nominal evaluation
            model.Ic[lambda] = model.Ic[lambda] + model.X_lambda[i].applyTranspose (model.Ic[i]);

            // derivative evaluation
            for (size_t idir = 0; idir < ndirs; idir++) {
                ad_model.hc[lambda][idir] = ad_model.hc[lambda][idir]
                        + model.X_lambda[i].toMatrixTranspose() * ad_model.hc[i][idir]
                        + ad_model.X_lambda[i][idir].transpose() * model.hc[i];
            }
            // nominal evaluation
            model.hc[lambda] = model.hc[lambda] + model.X_lambda[i].applyTranspose (model.hc[i]);
        } else {
            // derivative evaluation
            for (size_t idir = 0; idir < ndirs; idir++) {
                SpatialMatrix m =
                        model.X_lambda[i].toMatrixTranspose() * ad_model.Ic[i][idir] * model.X_lambda[i].toMatrix()
                        + model.X_lambda[i].toMatrixTranspose() * model.Ic[i].toMatrix() * ad_model.X_lambda[i][idir]
                        + ad_model.X_lambda[i][idir].transpose() * model.Ic[i].toMatrix() * model.X_lambda[i].toMatrix();

                SpatialRigidBodyInertia m2i;
                m2i.createFromMatrix(model.X_lambda[i].toMatrixTranspose() * ad_model.Ic[i][idir] * model.X_lambda[i].toMatrix());
                ad_Itot[idir] = ad_Itot[idir] + m2i;
            }
            // nominal evaluation
            Itot = Itot + model.X_lambda[i].applyTranspose (model.Ic[i]);

            // derivative evaluation
            for (size_t idir = 0; idir < ndirs; idir++) {
                ad_htot[idir] += model.X_lambda[i].toMatrixTranspose() * ad_model.hc[i][idir]
                        + ad_model.X_lambda[i][idir].transpose() * model.hc[i];
            }
            // nominal evaluation
            htot = htot + model.X_lambda[i].applyTranspose (model.hc[i]);
        }
    }

    mass = Itot.m;
    com = Itot.h / mass;
    for (size_t idir = 0; idir < ndirs; idir++) {
        ad_com.block<3,1>(0, idir) = ad_Itot[idir].h / mass;
    }

    LOG << "mass = " << mass << " com = " << com.transpose() << " htot = " << htot.transpose() << std::endl;

    if (com_velocity) {
        // derivative evaluation
        for (size_t idir = 0; idir < ndirs; idir++) {
            ad_com_velocity->block<3, 1>(0, idir) = Vector3d(ad_htot[idir][3] / mass, ad_htot[idir][4] / mass, ad_htot[idir][5] / mass);
        }
        // nominal evaluation
        *com_velocity = Vector3d (htot[3] / mass, htot[4] / mass, htot[5] / mass);
    }

    if (angular_momentum) {
        // derivative evaluation
        for (size_t idir = 0; idir < ndirs; idir++) {
            ad_htot[idir] = Xtrans(com).toMatrixAdjoint() * ad_htot[idir]
                            + Xtrans(ad_com.block<3,1>(0, idir)).applyAdjoint(ad_htot[idir]);
        }
        // nominal evaluation
        htot = Xtrans (com).applyAdjoint (htot);

        // derivative evaluation
        for (size_t idir = 0; idir < ndirs; idir++) {
            ad_angular_momentum->block<3,1>(0, idir) = Vector3d(ad_htot[idir][0], ad_htot[idir][1], ad_htot[idir][2]);
        }
        // nominal evaluation
        angular_momentum->set (htot[0], htot[1], htot[2]);
    }
}

TEST_FIXTURE ( CartPendulum, CartPendulumCalcCenterOfMass) {
    q[0] = 0.3;
    q[1] = -0.2;
    Vector3d point_body_coordinates (0.1, 3.2, 4.2);
    qdot[0] = .1;
    qdot[1] = -.15;


    MatrixNd q_dirs = MatrixNd::Identity(model.dof_count, model.dof_count);
    MatrixNd qdot_dirs = MatrixNd::Identity(model.dof_count, model.dof_count);

    double   ad_mass;
    Vector3d ad_com;
    MatrixNd ad_d_com(3, q_dirs.cols());
    ad_CalcCenterOfMass(model, ad_model, q, q_dirs, qdot, qdot_dirs, ad_mass,
                        ad_com, ad_d_com,
                        0, 0, 0, 0, true);

//    double   fd_mass;
//    Vector3d fd_com;
//    MatrixNd fd_d_com(3, q_dirs.cols());
//    fd_CalcCenterOfMass(model, q, q_dirs, qdot, qdot_dirs, fd_mass, fd_com,
//                        fd_d_com, 0, 0, 0, 0);

//    CHECK_EQUAL(ad_mass, fd_mass);

}


// TEST_FIXTURE ( CartPendulum, CartPendulumJacobianADSimple ) {
// 	MatrixNd jacobian_ad = MatrixNd::Zero(3, model.qdot_size);
// 	MatrixNd jacobian_ref = MatrixNd::Zero(3, model.qdot_size);
// 	MatrixNd jacobian_fd = MatrixNd::Zero(3, model.qdot_size);

// 	q.setZero();
// 	q[0] = -1.0;
// 	q[1] = -0.3;
// 	body_point = Vector3d (1.0, 2.0, 3.0);

// 	CalcPointJacobian (model, q, id_pendulum, body_point, jacobian_ref);

// 	MatrixNd q_dirs = MatrixNd::Identity (model.qdot_size, model.qdot_size);
// 	Vector3d base_point_standard = CalcBodyToBaseCoordinates (model, q, id_pendulum, body_point);
// 	Vector3d base_point_ad = ad_CalcBodyToBaseCoordinatesSingleFunc (model, ad_model, q, q_dirs, id_pendulum, body_point, jacobian_ad);
// 	Vector3d base_point_fd = fd_CalcBodyToBaseCoordinatesSingleFunc (model, q, q_dirs, id_pendulum, body_point, jacobian_fd);

// 	CHECK_ARRAY_CLOSE (jacobian_ref.data(), jacobian_ad.data(), 3 * model.qdot_size, TEST_PREC);
// //	CHECK_ARRAY_CLOSE (v_fixed_body.data(), v_body.data(), 6, TEST_PREC);
// }





// TEST_FIXTURE (CartPendulum, jcalcNominalSolutionTest) {
// 	// set nominal values
// 	q.setZero();
// 	q[0] = 0.3;
// 	q[1] = -0.2;

// 	qdot.setZero();
// 	qdot[0] = 0.3;
// 	qdot[1] = -0.2;

// 	// set directions
// 	unsigned int ndirs = model.q_size + model.qdot_size;
// 	MatrixNd x = MatrixNd::Identity(ndirs, ndirs);
// 	MatrixNd q_dirs = x.block(0, 0, model.q_size, ndirs);
// 	MatrixNd qdot_dirs = x.block(model.q_size, 0, model.qdot_size, ndirs);

// 	// set derivative outputs
// 	std::vector<SpatialMatrix> ad_X_Ji (ndirs, SpatialMatrix::Zero());
// 	std::vector<std::vector<SpatialMatrix> > ad_X_J (model.mBodies.size(), ad_X_Ji);

// 	std::vector<SpatialVector> ad_V (ndirs, SpatialVector::Zero (6));
// 	std::vector<std::vector<SpatialVector> > ad_S (model.mBodies.size(), ad_V);
// 	std::vector<std::vector<SpatialVector> > ad_v_J (model.mBodies.size(), ad_V);
// 	std::vector<std::vector<SpatialVector> > ad_c_J (model.mBodies.size(), ad_V);

// 	for (unsigned int joint_id = 1; joint_id < model.mBodies.size(); joint_id++) {
// 		// evaluate nominal solution
// 		jcalc (model, joint_id, q, qdot);
// 		MatrixNd ref_X_Ji = model.X_J[joint_id].toMatrix();
// 		SpatialVector ref_S_i = model.S[joint_id];
// 		SpatialVector ref_v_Ji = model.v_J[joint_id];
// 		SpatialVector ref_c_Ji = model.c_J[joint_id];

// 		// evaluate AD nominal solution
// 		ad_jcalc (model, ad_model, joint_id, q, q_dirs, qdot, qdot_dirs);
// 		MatrixNd test_X_Ji = model.X_J[joint_id].toMatrix();
// 		SpatialVector test_S_i = model.S[joint_id];
// 		SpatialVector test_v_Ji = model.v_J[joint_id];
// 		SpatialVector test_c_Ji = model.c_J[joint_id];

// 		CHECK_ARRAY_CLOSE (ref_X_Ji.data(), test_X_Ji.data(), 36, TEST_PREC);
// 		CHECK_ARRAY_CLOSE (ref_S_i.data(),  test_S_i.data(),   6, TEST_PREC);
// 		CHECK_ARRAY_CLOSE (ref_v_Ji.data(), test_v_Ji.data(),  6, TEST_PREC);
// 		CHECK_ARRAY_CLOSE (ref_c_Ji.data(), test_c_Ji.data(),  6, TEST_PREC);
// 		/* DEBUG OUTPUT
// 		cout << "===== joint_id: " << joint_id << " =====" << endl;
// 		cout << "ref_X_Ji: " << endl << ref_X_Ji << endl;
// 		cout << "test_X_Ji: " << endl << test_X_Ji << endl;
// 		cout << "error_X_Ji: " << endl << ref_X_Ji - test_X_Ji << endl;
// 		cout << endl;

// 		cout << "ref_S_i: " << endl << ref_S_i << endl;
// 		cout << "test_S_i: " << endl << test_S_i << endl;
// 		cout << "error_S_i: " << endl << ref_S_i - test_S_i << endl;
// 		cout << endl;

// 		cout << "ref_v_Ji: " << endl << ref_v_Ji << endl;
// 		cout << "test_v_Ji: " << endl << test_v_Ji << endl;
// 		cout << "error_v_Ji: " << endl << ref_v_Ji - test_v_Ji << endl;
// 		cout << endl;

// 		cout << "ref_c_Ji: " << endl << ref_c_Ji << endl;
// 		cout << "test_c_Ji: " << endl << test_c_Ji << endl;
// 		cout << "error_c_Ji: " << endl << ref_c_Ji - test_c_Ji << endl;
// 		cout << endl;
// 		cout << endl;
// 		*/
// 	}
// }

// TEST_FIXTURE (CartPendulum, jcalcFDvsADTest) {
// 	double TEST_PREC = 1.0e-08;

// 	// set nominal values
// 	q.setZero();
// 	q[0] = 0.3;
// 	q[1] = -0.2;

// 	qdot.setZero();
// 	qdot[0] = 0.3;
// 	qdot[1] = -0.2;

// 	// set directions
// 	unsigned int ndirs = model.q_size + model.qdot_size;
// 	MatrixNd x = MatrixNd::Identity(ndirs, ndirs);
// 	MatrixNd q_dirs = x.block(0, 0, model.q_size, ndirs);
// 	MatrixNd qdot_dirs = x.block(model.q_size, 0, model.qdot_size, ndirs);

// 	// set derivative outputs
// 	std::vector<SpatialMatrix> ad_X_Ji (ndirs, SpatialMatrix::Zero());
// 	std::vector<std::vector<SpatialMatrix> > fd_X_J (model.mBodies.size(), ad_X_Ji);
// 	std::vector<std::vector<SpatialMatrix> > ad_X_J (model.mBodies.size(), ad_X_Ji);

// 	std::vector<SpatialVector> ad_V (ndirs, SpatialVector::Zero (6));
// 	std::vector<std::vector<SpatialVector> > fd_S (model.mBodies.size(), ad_V);
// 	std::vector<std::vector<SpatialVector> > fd_v_J (model.mBodies.size(), ad_V);
// 	std::vector<std::vector<SpatialVector> > fd_c_J (model.mBodies.size(), ad_V);

// 	for (unsigned int joint_id = 1; joint_id < model.mBodies.size(); joint_id++) {
// 		// evaluate nominal solution
// 		fd_jcalc (model, joint_id, q, q_dirs, qdot, qdot_dirs,
// 			fd_X_J[joint_id], fd_S[joint_id], fd_v_J[joint_id], fd_c_J[joint_id]);

// 		// evaluate AD nominal solution
// 		ad_jcalc (model, ad_model, joint_id, q, q_dirs, qdot, qdot_dirs);

// 		for (int idir = 0; idir < ndirs; ++idir) {
// 			CHECK_ARRAY_CLOSE (fd_X_J[joint_id][idir].data(), ad_model.X_J[joint_id][idir].data(), 36, TEST_PREC);
// 			CHECK_ARRAY_CLOSE (fd_S[joint_id][idir].data(),   ad_model.S[joint_id][idir].data(),    6, TEST_PREC);
// 			CHECK_ARRAY_CLOSE (fd_v_J[joint_id][idir].data(), ad_model.v_J[joint_id][idir].data(),  6, TEST_PREC);
// 			CHECK_ARRAY_CLOSE (fd_c_J[joint_id][idir].data(), ad_model.c_J[joint_id][idir].data(),  6, TEST_PREC);
// 			/* DEBUG OUTPUT
// 			cout << "===== joint_id: " << joint_id << ", idir: " << idir << " =====" << endl;
// 			cout << "fd_X_J[" << joint_id << "][" << idir << "]: " << endl << fd_X_J[joint_id][idir] << endl;
// 			cout << "ad_X_J[" << joint_id << "][" << idir << "]: " << endl << ad_X_J[joint_id][idir] << endl;
// 			cout << "error_X_J: " << endl << fd_X_J[joint_id][idir] - ad_X_J[joint_id][idir] << endl;

// 			cout << "fd_S[" << joint_id << "][" << idir << "]: " << endl << fd_S[joint_id][idir] << endl;
// 			cout << "ad_S[" << joint_id << "][" << idir << "]: " << endl << ad_S[joint_id][idir] << endl;
// 			cout << "error_S: " << endl << fd_S[joint_id][idir] - ad_S[joint_id][idir] << endl;

// 			cout << "fd_v_J[" << joint_id << "][" << idir << "]: " << endl << fd_v_J[joint_id][idir] << endl;
// 			cout << "ad_v_J[" << joint_id << "][" << idir << "]: " << endl << ad_v_J[joint_id][idir] << endl;
// 			cout << "error_v_J: " << endl << fd_v_J[joint_id][idir] - ad_v_J[joint_id][idir] << endl;

// 			cout << "fd_c_J[" << joint_id << "][" << idir << "]: " << endl << fd_c_J[joint_id][idir] << endl;
// 			cout << "ad_c_J[" << joint_id << "][" << idir << "]: " << endl << ad_c_J[joint_id][idir] << endl;
// 			cout << "error_c_J: " << endl << fd_c_J[joint_id][idir] - ad_c_J[joint_id][idir] << endl;
// 			*/
// 		}
// 	}
// }

// TEST_FIXTURE (CartPendulum, ad_UpdateKinematicsNominalTest) {
// 	// set nominal values
// 	q.setZero();
// 	q[0] = 0.3;
// 	q[1] = -0.2;

// 	qdot.setZero();
// 	qdot[0] = 0.3;
// 	qdot[1] = -0.2;

// 	qddot.setZero();
// 	qddot[0] = 0.3;
// 	qddot[1] = -0.2;

// 	// set directions
// 	unsigned int ndirs = 3*model.q_size;
// 	MatrixNd x = MatrixNd::Identity(ndirs, ndirs);
// 	MatrixNd q_dirs = x.block(0, 0, model.q_size, ndirs);
// 	MatrixNd qdot_dirs = x.block(model.q_size, 0, model.q_size, ndirs);
// 	MatrixNd qddot_dirs = x.block(2*model.q_size, 0, model.q_size, ndirs);

// 	// evaluate nominal solution
// 	UpdateKinematics (model, q, qdot, qddot);
// 	Model res = model;

// 	// evaluate AD nominal solution
// 	ad_UpdateKinematics (
// 		model, ad_model, q, q_dirs, qdot, qdot_dirs, qddot, qddot_dirs
// 	);
// 	Model ad_res = model;

// 	for (unsigned int joint_id = 1; joint_id < model.mBodies.size(); joint_id++) {
// 		CHECK_ARRAY_CLOSE (res.X_lambda[joint_id].toMatrix().data(), ad_res.X_lambda[joint_id].toMatrix().data(), 36, TEST_PREC);
// 		CHECK_ARRAY_CLOSE (res.X_base[joint_id].toMatrix().data(), ad_res.X_base[joint_id].toMatrix().data(), 36, TEST_PREC);

// 		CHECK_ARRAY_CLOSE (res.a[joint_id].data(), ad_res.a[joint_id].data(), 6, TEST_PREC);
// 		CHECK_ARRAY_CLOSE (res.v[joint_id].data(), ad_res.v[joint_id].data(), 6, TEST_PREC);
// 		CHECK_ARRAY_CLOSE (res.c[joint_id].data(), ad_res.c[joint_id].data(), 6, TEST_PREC);
// 	}
// }

// TEST_FIXTURE (CartPendulum, ad_UpdateKinematicsADvsFDTest) {
// 	// set nominal values
// 	q.setZero();
// 	q[0] = 0.3;
// 	q[1] = -0.2;

// 	qdot.setZero();
// 	qdot[0] = 0.3;
// 	qdot[1] = -0.2;

// 	qddot.setZero();
// 	qddot[0] = 0.3;
// 	qddot[1] = -0.2;

// 	// set directions
// 	unsigned int ndirs = 3*model.q_size;
// 	MatrixNd x = MatrixNd::Identity(ndirs, ndirs);
// 	MatrixNd q_dirs = x.block(0, 0, model.q_size, ndirs);
// 	MatrixNd qdot_dirs = x.block(model.q_size, 0, model.q_size, ndirs);
// 	MatrixNd qddot_dirs = x.block(2*model.q_size, 0, model.q_size, ndirs);

// 	// set derivative outputs
// 	std::vector<SpatialMatrix> ad_X (ndirs, SpatialMatrix::Zero());
// 	std::vector<std::vector<SpatialMatrix> > fd_X_lambda (model.mBodies.size(), ad_X);
// 	std::vector<std::vector<SpatialMatrix> > fd_X_base (model.mBodies.size(), ad_X);

// 	std::vector<SpatialVector> ad_V (ndirs, SpatialVector::Zero (6));
// 	std::vector<std::vector<SpatialVector> > fd_a (model.mBodies.size(), ad_V);
// 	std::vector<std::vector<SpatialVector> > fd_v (model.mBodies.size(), ad_V);
// 	std::vector<std::vector<SpatialVector> > fd_c (model.mBodies.size(), ad_V);

// 	// evaluate nominal solution
// 	fd_UpdateKinematics (
// 		model, q, q_dirs, qdot, qdot_dirs, qddot, qddot_dirs,
// 		fd_X_lambda, fd_X_base, fd_a, fd_v, fd_c
// 	);

// 	ad_UpdateKinematics (
// 		model, ad_model, q, q_dirs, qdot, qdot_dirs, qddot, qddot_dirs
// 	);

// 	for (unsigned int joint_id = 1; joint_id < model.mBodies.size(); joint_id++) {
// 		for (int idirs = 0; idirs < ndirs; ++idirs) {
// 			CHECK_ARRAY_CLOSE (fd_X_lambda[joint_id][idirs].data(), ad_model.X_lambda[joint_id][idirs].data(), 36, TEST_PREC);
// 			CHECK_ARRAY_CLOSE (fd_X_base[joint_id][idirs].data(), ad_model.X_base[joint_id][idirs].data(), 36, TEST_PREC);

// 			CHECK_ARRAY_CLOSE (fd_a[joint_id][idirs].data(), ad_model.a[joint_id][idirs].data(), 6, TEST_PREC);
// 			CHECK_ARRAY_CLOSE (fd_v[joint_id][idirs].data(), ad_model.v[joint_id][idirs].data(), 6, TEST_PREC);
// 			CHECK_ARRAY_CLOSE (fd_c[joint_id][idirs].data(), ad_model.c[joint_id][idirs].data(), 6, TEST_PREC);
// 			/* DEBUG OUTPUT
// 			// cout << "fd_X_lambda[" << joint_id << "][" << idirs << "]: " << endl << fd_X_lambda[joint_id][idirs] << endl;
// 			//cout << "ad_X_lambda[" << joint_id << "][" << idirs << "]: " << endl << ad_model.X_lambda[joint_id][idirs] << endl;
// 			//cout << "err" << endl << fd_X_lambda[joint_id][idirs] - ad_model.X_lambda[joint_id][idirs] << endl;
// 			// cout << endl;

// 			// cout << "fd_X_base[" << joint_id << "][" << idirs << "]: " << endl << fd_X_base[joint_id][idirs] << endl;
// 			//cout << "ad_X_base[" << joint_id << "][" << idirs << "]: " << endl << ad_model.X_base[joint_id][idirs] << endl;
// 			//cout << "err" << endl << fd_X_base[joint_id][idirs] - ad_model.X_base[joint_id][idirs] << endl;
// 			// cout << endl;

// 			// cout << "fd_a[" << joint_id << "][" << idirs << "]: " << endl << fd_a[joint_id][idirs] << endl;
// 			//cout << "ad_a[" << joint_id << "][" << idirs << "]: " << endl << ad_model.a[joint_id][idirs] << endl;
// 			//cout << "err" << endl << fd_a[joint_id][idirs] - ad_model.a[joint_id][idirs] << endl;
// 			// cout << endl;

// 			// cout << "fd_v[" << joint_id << "][" << idirs << "]: " << endl << fd_v[joint_id][idirs] << endl;
// 			//cout << "ad_v[" << joint_id << "][" << idirs << "]: " << endl << ad_model.v[joint_id][idirs] << endl;
// 			//cout << "err" << endl << fd_v[joint_id][idirs] - ad_model.v[joint_id][idirs] << endl;
// 			// cout << endl;

// 			cout << "fd_c[" << joint_id << "][" << idirs << "]: " << endl << fd_c[joint_id][idirs] << endl;
// 			//cout << "ad_c[" << joint_id << "][" << idirs << "]: " << endl << ad_model.c[joint_id][idirs] << endl;
// 			//cout << "err" << endl << fd_c[joint_id][idirs] - ad_model.c[joint_id][idirs] << endl;
// 			cout << endl;
// 			*/
// 		}
// 	}
// }

// TEST_FIXTURE (CartPendulum, CalcPointAccelerationNominalTest) {
// 	// set nominal values
// 	q.setZero();
// 	q[0] = 0.3;
// 	q[1] = -0.2;

// 	qdot.setZero();
// 	qdot[0] = 0.3;
// 	qdot[1] = -0.2;

// 	qddot.setZero();
// 	qddot[0] = 0.3;
// 	qddot[1] = -0.2;

// 	Vector3d point_body_coordinates (0.1, 3.2, 4.2);

// 	// set directions
// 	unsigned int ndirs = 3*model.q_size;
// 	MatrixNd x = MatrixNd::Identity(ndirs, ndirs);
// 	MatrixNd q_dirs = x.block(0, 0, model.q_size, ndirs);
// 	MatrixNd qdot_dirs = x.block(model.q_size, 0, model.q_size, ndirs);
// 	MatrixNd qddot_dirs = x.block(2*model.q_size, 0, model.q_size, ndirs);

// 	// set derivative output
// 	MatrixNd ad_jacobian = MatrixNd::Zero(3, ndirs);

// 	// evaluate nominal solution
// 	Vector3d ref_acc = CalcPointAcceleration (
// 		model, q, qdot, qddot, id_pendulum, point_body_coordinates
// 	);

// 	Vector3d ad_acc = ad_CalcPointAcceleration (
// 		model, ad_model,
// 		q, q_dirs, qdot, qdot_dirs, qddot, qddot_dirs,
// 		id_pendulum, point_body_coordinates,
// 		ad_jacobian
// 	);

// 	CHECK_ARRAY_CLOSE (ref_acc.data(), ad_acc.data(), 3, TEST_PREC);
// }


// TEST_FIXTURE (CartPendulum, CalcPointAccelerationFDvsADTest) {
// 	// set nominal values
// 	q.setZero();
// 	q[0] = 0.3;
// 	q[1] = -0.2;

// 	qdot.setZero();
// 	qdot[0] = 0.3;
// 	qdot[1] = -0.2;

// 	qddot.setZero();
// 	qddot[0] = 0.3;
// 	qddot[1] = -0.2;

// 	Vector3d point_body_coordinates (0.1, 3.2, 4.2);

// 	// set directions
// 	unsigned int ndirs = model.q_size + model.qdot_size + model.q_size;
// 	MatrixNd x = MatrixNd::Identity(ndirs, ndirs);
// 	MatrixNd q_dirs = x.block(0, 0, model.q_size, ndirs);
// 	MatrixNd qdot_dirs = x.block(model.q_size, 0, model.qdot_size, ndirs);
// 	MatrixNd qddot_dirs = x.block(model.q_size + model.qdot_size, 0, model.q_size, ndirs);


// 	// set derivative output
// 	MatrixNd fd_jacobian = MatrixNd::Zero(3, ndirs);
// 	MatrixNd ad_jacobian = MatrixNd::Zero(3, ndirs);

// 	// evaluate nominal solution
// 	Vector3d ref_acc = fd_CalcPointAcceleration (
// 		model,
// 		q, q_dirs, qdot, qdot_dirs, qddot, qddot_dirs,
// 		id_pendulum, point_body_coordinates,
// 		ad_jacobian
// 	);

// 	Vector3d ad_acc = ad_CalcPointAcceleration (
// 		model, ad_model,
// 		q, q_dirs, qdot, qdot_dirs, qddot, qddot_dirs,
// 		id_pendulum, point_body_coordinates,
// 		ad_jacobian
// 	);

// 	CHECK_ARRAY_CLOSE (ref_acc.data(), ad_acc.data(), 3, TEST_PREC);
// }



// TEST_FIXTURE(CartPendulum, InverseDynamicsADTest){
// 	for(unsigned int i = 0; i < model.qdot_size; i++){
// 		q[i] = (i+1)*(0.897878435);
// 		qdot[i] = (i+1)*(0.27563682);
// 		qddot[i] = (i+1)*(0.06565644564455);
// 	}

// 	MatrixNd q_dirs 	= MatrixNd::Identity (model.qdot_size, model.qdot_size);
// 	MatrixNd qdot_dirs 	= MatrixNd::Identity (model.qdot_size, model.qdot_size);
// 	MatrixNd qddot_dirs = MatrixNd::Identity (model.qdot_size, model.qdot_size);

// 	std::vector<SpatialVector> f_ext (model.mBodies.size(),SpatialVector::Zero(model.q_size));
// 	for (int i = 0; i < model.mBodies.size(); ++i) {
// 	  f_ext[i]=SpatialVector::Zero(model.q_size);
// 	}



// 	VectorNd tau_ref (tau);

// 	MatrixNd ad_tau  = MatrixNd::Zero(model.qdot_size, model.qdot_size);
// 	MatrixNd fd_tau  = MatrixNd::Zero(model.qdot_size, model.qdot_size);

// 	ad_InverseDynamics(model,
// 			ad_model,
// 			q,
// 			q_dirs,
// 			qdot,
// 			qdot_dirs,
// 			qddot,
// 			qddot_dirs,
// 			tau,
// 			ad_tau,
// 			&f_ext);


// 	fd_InverseDynamics(model,
// 			q,
// 			q_dirs,
// 			qdot,
// 			qdot_dirs,
// 			qddot,
// 			qddot_dirs,
// 			tau_ref,
// 			fd_tau,
// 			&f_ext);

// 	CHECK_ARRAY_CLOSE (tau_ref.data(), tau.data(), tau_ref.rows(), 1e-7);
// 	CHECK_ARRAY_CLOSE (fd_tau.data(), ad_tau.data(), fd_tau.cols()*fd_tau.rows(), 1e-7);
// }

// TEST_FIXTURE(CartPendulum, ForwardDynamicsADTest){
//   srand((unsigned int) time(0));

//   for(unsigned int trial = 0; trial < 10; trial++) {
// 	VectorNd q = VectorNd::Random(model.q_size);
// 	VectorNd qdot = VectorNd::Random(model.q_size);
// 	VectorNd tau = VectorNd::Random(model.q_size);

// 	unsigned int ndirs = 3 * model.q_size;
// 	MatrixNd x = MatrixNd::Identity(ndirs, ndirs);
// 	MatrixNd q_dirs = x.block(0, 0, model.q_size, ndirs);
// 	MatrixNd qdot_dirs = x.block(model.q_size, 0, model.q_size, ndirs);
// 	MatrixNd tau_dirs = x.block(2*model.q_size, 0, model.q_size, ndirs);

// 	std::vector<SpatialVector> f_ext (model.mBodies.size(),SpatialVector::Zero(model.q_size));
// 	for (int i = 0; i < model.mBodies.size(); ++i) {
// 	  f_ext[i]=SpatialVector::Zero(model.q_size);
// 	}

// 	VectorNd qddot (VectorNd::Zero(model.q_size));
// 	MatrixNd ad_qddot  = MatrixNd::Zero(model.qdot_size, ndirs);
// 	MatrixNd fd_qddot  = MatrixNd::Zero(model.qdot_size, ndirs);


// 	fd_ForwardDynamics(model, q, q_dirs, qdot, qdot_dirs, tau, tau_dirs, qddot, fd_qddot,&f_ext);

// 	ad_ForwardDynamics(model, ad_model, q, q_dirs, qdot, qdot_dirs, tau, tau_dirs, qddot, ad_qddot, &f_ext);

// 	CHECK_ARRAY_CLOSE (fd_qddot.data(), ad_qddot.data(), fd_qddot.cols()*fd_qddot.rows(), 1e-7);
//   }
//   }

// TEST_FIXTURE (CartPendulum, ForwardDynamicsCholesky) {
//   // set nominal values

//   q = VectorNd::Random(model.q_size);
//   qdot = VectorNd::Random(model.q_size);
//   qddot = VectorNd::Random(model.q_size);
//   tau = VectorNd::Random(model.q_size);


//   std::vector<SpatialVector> f_ext (model.mBodies.size(),SpatialVector::Zero(model.q_size));
//   for (int i = 0; i < model.mBodies.size(); ++i) {
// 	f_ext[i]=SpatialVector::Random(model.q_size);
//   }


//   VectorNd qddot_ref (VectorNd::Zero (model.q_size));

//   ForwardDynamicsCholesky(model,q,qdot,tau,qddot,&f_ext);

//   ForwardDynamics(model,q,qdot,tau,qddot_ref,&f_ext);

//   CHECK_ARRAY_CLOSE (qddot, qddot_ref, model.q_size, TEST_PREC);
//   /*
// 	cout << "qddot_ref: " << endl << qddot_ref << endl;
// 	cout << "qddot_test: " << endl << qddot << endl;
// 	cout << "error qddot: " << endl << qddot_ref - qddot << endl;
//   */
// }


// TEST_FIXTURE(CartPendulum, ForwardDynamicsCholeskyADTest){

//   VectorNd q = VectorNd::Random(model.q_size);
//   VectorNd qdot = VectorNd::Random(model.q_size);
//   VectorNd tau = VectorNd::Random(model.q_size);
//   VectorNd qddot = VectorNd::Zero(model.q_size);
//   VectorNd qddot_ref = VectorNd::Zero(model.q_size);

//   unsigned int ndirs = 3 * model.q_size;
//   MatrixNd x = MatrixNd::Identity(ndirs, ndirs);
//   MatrixNd q_dirs = x.block(0, 0, model.q_size, ndirs);
//   MatrixNd qdot_dirs = x.block(model.q_size, 0, model.q_size, ndirs);
//   MatrixNd tau_dirs = x.block(2*model.q_size, 0, model.q_size, ndirs);


//   std::vector<SpatialVector> f_ext (model.mBodies.size(),SpatialVector::Zero(model.q_size));
//   for (int i = 0; i < model.mBodies.size(); ++i) {
// 	f_ext[i]=SpatialVector::Zero(model.q_size);
//   }

//   MatrixNd ad_qddot  = MatrixNd::Zero(model.q_size, ndirs);
//   MatrixNd fd_qddot  = MatrixNd::Zero(model.q_size, ndirs);

//   fd_ForwardDynamics(model,
// 			 q,
// 			 q_dirs,
// 			 qdot,
// 			 qdot_dirs,
// 			 tau,
// 			 tau_dirs,
// 			 qddot_ref,
// 			 fd_qddot,
// 			 &f_ext);

//   ad_ForwardDynamicsCholesky(model,
// 				 ad_model,
// 				 q,
// 				 q_dirs,
// 				 qdot,
// 				 qdot_dirs,
// 				 tau,
// 				 tau_dirs,
// 				 qddot,
// 				 ad_qddot,
// 				 &f_ext);

//   CHECK_ARRAY_CLOSE (qddot_ref.data(), qddot.data(), qddot.size(), 1e-7);
//   CHECK_ARRAY_CLOSE (fd_qddot.data(), ad_qddot.data(), fd_qddot.cols()*fd_qddot.rows(), 1e-7);
// }




// TEST_FIXTURE(CartPendulum, InverseDynamicsADTest_with_external_forces){
// 	for(unsigned int i = 0; i < model.qdot_size; i++){
// 		q[i] = (i+1)*(0.897878435);
// 		qdot[i] = (i+1)*(0.27563682);
// 		qddot[i] = (i+1)*(0.06565644564455);
// 	}

// 	MatrixNd q_dirs 	= MatrixNd::Identity (model.qdot_size, model.qdot_size);
// 	MatrixNd qdot_dirs 	= MatrixNd::Identity (model.qdot_size, model.qdot_size);
// 	MatrixNd qddot_dirs = MatrixNd::Identity (model.qdot_size, model.qdot_size);

// 	std::vector<SpatialVector> f_ext (model.mBodies.size(),SpatialVector::Zero(model.q_size));
// 	for (int i = 0; i < model.mBodies.size(); ++i) {
// 	  f_ext[i]=SpatialVector::Random(model.q_size);
// 	}


// 	VectorNd tau_ref (tau);

// 	MatrixNd ad_tau  = MatrixNd::Zero(model.qdot_size, model.qdot_size);
// 	MatrixNd fd_tau  = MatrixNd::Zero(model.qdot_size, model.qdot_size);

// 	ad_InverseDynamics(model,
// 			ad_model,
// 			q,
// 			q_dirs,
// 			qdot,
// 			qdot_dirs,
// 			qddot,
// 			qddot_dirs,
// 			tau,
// 			ad_tau,
// 			&f_ext);


// 	fd_InverseDynamics(model,
// 			q,
// 			q_dirs,
// 			qdot,
// 			qdot_dirs,
// 			qddot,
// 			qddot_dirs,
// 			tau_ref,
// 			fd_tau,
// 			&f_ext);

// 	CHECK_ARRAY_CLOSE (fd_tau.data(), ad_tau.data(), fd_tau.cols()*fd_tau.rows(), 1e-7);
// 	CHECK_ARRAY_CLOSE (tau_ref.data(), tau.data(), tau_ref.rows(), 1e-7);

// }

// TEST_FIXTURE(CartPendulum, ForwardDynamicsADTest_with_external_forces){
//   srand((unsigned int) time(0));

//   for(unsigned int trial = 0; trial < 10; trial++) {
// 	VectorNd q = VectorNd::Random(model.q_size);
// 	VectorNd qdot = VectorNd::Random(model.q_size);
// 	VectorNd tau = VectorNd::Random(model.q_size);

// 	unsigned int ndirs = 3 * model.q_size;
// 	MatrixNd x = MatrixNd::Identity(ndirs, ndirs);
// 	MatrixNd q_dirs = x.block(0, 0, model.q_size, ndirs);
// 	MatrixNd qdot_dirs = x.block(model.q_size, 0, model.q_size, ndirs);
// 	MatrixNd tau_dirs = x.block(2*model.q_size, 0, model.q_size, ndirs);

// 	std::vector<SpatialVector> f_ext (model.mBodies.size(),SpatialVector::Zero(model.q_size));
// 	for (int i = 0; i < model.mBodies.size(); ++i) {
// 	  f_ext[i]=SpatialVector::Random(model.q_size);
// 	}


// 	VectorNd qddot (VectorNd::Zero(model.q_size));
// 	MatrixNd ad_qddot  = MatrixNd::Zero(model.qdot_size, ndirs);
// 	MatrixNd fd_qddot  = MatrixNd::Zero(model.qdot_size, ndirs);


// 	fd_ForwardDynamics(model, q, q_dirs, qdot, qdot_dirs, tau, tau_dirs, qddot, fd_qddot,&f_ext);

// 	ad_ForwardDynamics(model, ad_model, q, q_dirs, qdot, qdot_dirs, tau, tau_dirs, qddot, ad_qddot, &f_ext);

// 	CHECK_ARRAY_CLOSE (fd_qddot.data(), ad_qddot.data(), fd_qddot.cols()*fd_qddot.rows(), 1e-7);
//   }
//   }



// TEST_FIXTURE(CartPendulum, ForwardDynamicsCholeskyADTest_with_external_forces){

//   VectorNd q = VectorNd::Random(model.q_size);
//   VectorNd qdot = VectorNd::Random(model.q_size);
//   VectorNd tau = VectorNd::Random(model.q_size);
//   VectorNd qddot = VectorNd::Zero(model.q_size);
//   VectorNd qddot_ref = VectorNd::Zero(model.q_size);

//   unsigned int ndirs = 3 * model.q_size;
//   MatrixNd x = MatrixNd::Identity(ndirs, ndirs);
//   MatrixNd q_dirs = x.block(0, 0, model.q_size, ndirs);
//   MatrixNd qdot_dirs = x.block(model.q_size, 0, model.q_size, ndirs);
//   MatrixNd tau_dirs = x.block(2*model.q_size, 0, model.q_size, ndirs);


//   std::vector<SpatialVector> f_ext (model.mBodies.size(),SpatialVector::Zero(model.q_size));
//   for (int i = 0; i < model.mBodies.size(); ++i) {
// 	f_ext[i]=SpatialVector::Random(model.q_size);
//   }

//   MatrixNd ad_qddot  = MatrixNd::Zero(model.q_size, ndirs);
//   MatrixNd fd_qddot  = MatrixNd::Zero(model.q_size, ndirs);

//   fd_ForwardDynamics(model,
// 			 q,
// 			 q_dirs,
// 			 qdot,
// 			 qdot_dirs,
// 			 tau,
// 			 tau_dirs,
// 			 qddot_ref,
// 			 fd_qddot,
// 			 &f_ext);

//   ad_ForwardDynamicsCholesky(model,
// 				 ad_model,
// 				 q,
// 				 q_dirs,
// 				 qdot,
// 				 qdot_dirs,
// 				 tau,
// 				 tau_dirs,
// 				 qddot,
// 				 ad_qddot,
// 				 &f_ext);

//   CHECK_ARRAY_CLOSE (qddot_ref.data(), qddot.data(), qddot.size(), 1e-7);
//   CHECK_ARRAY_CLOSE (fd_qddot.data(), ad_qddot.data(), fd_qddot.cols()*fd_qddot.rows(), 1e-7);
// }
