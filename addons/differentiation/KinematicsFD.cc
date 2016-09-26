/*
 * RBDL - Rigid Body Dynamics Library
 * Copyright (c) 2011-2015 Martin Felis <martin.felis@iwr.uni-heidelberg.de>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#include "rbdl/rbdl_mathutils.h"
#include "rbdl/Logging.h"

#include "rbdl/Model.h"
#include "rbdl/Kinematics.h"

#include "KinematicsAD.h"
#include "KinematicsFD.h"

using std::vector;

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------
namespace FD {
// -----------------------------------------------------------------------------

using namespace RigidBodyDynamics::Math;

RBDL_DLLAPI Vector3d CalcBodyToBaseCoordinatesSingleFunc (
        Model &model,
        const VectorNd &q,
        const MatrixNd &q_dirs,
        unsigned int body_id,
        const Vector3d &point_body_coordinates,
        MatrixNd &out) {
    Vector3d ref = CalcBodyToBaseCoordinatesSingleFunc (model, q, body_id, point_body_coordinates);

    unsigned int ndirs = q_dirs.cols();
    double const h = 1e-8;
    for (unsigned int idir = 0; idir < ndirs; idir++) {
        VectorNd q_dir = q_dirs.block(0, idir, model.q_size, 1);
        Vector3d res_hd = CalcBodyToBaseCoordinatesSingleFunc (model, q + h * q_dir, body_id, point_body_coordinates);
//        Vector3d res_hd_rbdl = CalcBodyToBaseCoordinates (model, q + h * q_dir, body_id, point_body_coordinates);
        out.block<3,1>(0, idir) = (res_hd - ref) / h;
    }

    return ref;
}

RBDL_DLLAPI
Vector3d CalcBodyToBaseCoordinates (
    Model & model,
    ADModel & ad_model,
    VectorNd const & q,
    MatrixNd const & q_dirs,
    unsigned int body_id,
    Vector3d const & point_body_coordinates,
    MatrixNd & fd_body_to_base_coordinates
) {
  Vector3d ret = Vector3d::Zero();

  // update the Kinematics if necessary
  UpdateKinematicsCustom (model, &q, NULL, NULL);

  if (body_id >= model.fixed_body_discriminator) {
    unsigned int fbody_id = body_id - model.fixed_body_discriminator;
    unsigned int parent_id = model.mFixedBodies[fbody_id].mMovableParent;

    Matrix3d fixed_rotation = model.mFixedBodies[fbody_id].mParentTransform.E.transpose();
    Vector3d fixed_position = model.mFixedBodies[fbody_id].mParentTransform.r;

    Matrix3d parent_body_rotation = model.X_base[parent_id].E.transpose();
    Vector3d parent_body_position = model.X_base[parent_id].r;

    ret = parent_body_position
        + parent_body_rotation * (
          fixed_position + fixed_rotation * (point_body_coordinates)
          );
  } else {
    Matrix3d body_rotation = model.X_base[body_id].E.transpose();
    Vector3d body_position = model.X_base[body_id].r;

    ret = body_position + body_rotation * point_body_coordinates;
  }

  unsigned int ndirs =  q_dirs.cols();
  Vector3d temp = Vector3d::Zero();

  assert (fd_body_to_base_coordinates.cols() == ndirs);

  double const h = 1e-8;

  for (unsigned int idir = 0; idir < ndirs; ++idir) {
    VectorNd Q_dir  =  q_dirs.col(idir);

    // update the Kinematics if necessary
    VectorNd Q_temp = q + h*Q_dir;
    UpdateKinematicsCustom (model, &Q_temp, NULL, NULL);

    if (body_id >= model.fixed_body_discriminator) {
      std::cerr << "fixed bodies not supported yet." << std::endl;
      abort();
      unsigned int fbody_id = body_id - model.fixed_body_discriminator;
      unsigned int parent_id = model.mFixedBodies[fbody_id].mMovableParent;

      Matrix3d fixed_rotation = model.mFixedBodies[fbody_id].mParentTransform.E.transpose();
      Vector3d fixed_position = model.mFixedBodies[fbody_id].mParentTransform.r;

      Matrix3d parent_body_rotation = model.X_base[parent_id].E.transpose();
      Vector3d parent_body_position = model.X_base[parent_id].r;

      temp = parent_body_position
          + parent_body_rotation
          * (fixed_position + fixed_rotation * (point_body_coordinates));
    } else {
      Matrix3d body_rotation = model.X_base[body_id].E.transpose();
      Vector3d body_position = model.X_base[body_id].r;

      temp = body_position + body_rotation * point_body_coordinates;
    }

    // derivative evaluation
    fd_body_to_base_coordinates.col(idir) = (temp - ret) / h;
  }
  return ret;
}

RBDL_DLLAPI Vector3d CalcBaseToBodyCoordinates (
    Model & model,
    VectorNd const & q,
    MatrixNd const & q_dirs,
    unsigned body_id,
    Vector3d const & base_point_position,
    MatrixNd const & base_point_position_dirs,
    MatrixNd & fd_base_to_body_coordinates
) {
  unsigned ndirs = q_dirs.cols();
  assert(ndirs == base_point_position_dirs.cols());

  Vector3d b2b = CalcBaseToBodyCoordinates (model, q, body_id,
      base_point_position);
  double h = 1e-8;
  for (unsigned idir = 0; idir < ndirs; idir++) {
    Vector3d b2bh = CalcBaseToBodyCoordinates(model, q + h * q_dirs.col(idir),
        body_id, base_point_position + h * base_point_position_dirs.col(idir));
    fd_base_to_body_coordinates.col(idir) = (b2bh - b2b) / h;
  }
  return b2b;
}

RBDL_DLLAPI Matrix3d CalcBodyWorldOrientation (
        Model & model,
        VectorNd const & q,
        MatrixNd const & q_dirs,
        const unsigned int body_id,
        vector<Matrix3d> & fd_derivative
) {
    unsigned ndirs = q_dirs.cols();
    double h = 1e-8;
    Matrix3d ref = CalcBodyWorldOrientation(model, q, body_id, true);
    for (unsigned idir = 0; idir < ndirs; idir++) {
        VectorNd q_dir  = q_dirs.block(0, idir, model.q_size, 1);
        Matrix3d res_hd = CalcBodyWorldOrientation(model, q + h * q_dir, body_id, true);
        fd_derivative[idir] = (res_hd - ref) / h;
    }
    return ref;
}

RBDL_DLLAPI Vector3d CalcPointVelocity (
    Model & model,
    VectorNd const & q,
    MatrixNd const & q_dirs,
    VectorNd const & qdot,
    MatrixNd const & qdot_dirs,
    unsigned int body_id,
    Vector3d const & point_position,
    MatrixNd & fd_point_velocity
) {
  unsigned ndirs = q_dirs.cols();
  assert(ndirs == qdot_dirs.cols());
  assert(3     == fd_point_velocity.rows());
  assert(ndirs == fd_point_velocity.cols());
  double h = 1e-8;
  Vector3d ref = CalcPointVelocity(model, q, qdot, body_id, point_position,
      true);
  for (unsigned idir = 0; idir < ndirs; idir++) {
    VectorNd qh = q + h * q_dirs.col(idir);
    VectorNd qdh = qdot + h * qdot_dirs.col(idir);
    fd_point_velocity.col(idir) = (CalcPointVelocity(model, qh, qdh, body_id,
        point_position, true) - ref) / h;
  }
  return ref;
}

RBDL_DLLAPI Vector3d CalcPointAcceleration (
		Model & model,
		VectorNd const & q,
		MatrixNd const & q_dirs,
		VectorNd const & qdot,
		MatrixNd const & qdot_dirs,
		VectorNd const & qddot,
		MatrixNd const & qddot_dirs,
		unsigned int body_id,
		Vector3d const & point_position,
		MatrixNd & fd_derivative
) {
	Vector3d ref = CalcPointAcceleration(model, q, qdot, qddot, body_id,
			point_position, true);
	unsigned int const ndirs = q_dirs.cols();
	double const h = 1e-8;
	for (unsigned idir = 0; idir < ndirs; idir++) {
		VectorNd q_dir     = q_dirs.block(0, idir, model.q_size, 1);
		VectorNd qdot_dir  = qdot_dirs.block(0, idir, model.qdot_size, 1);
		VectorNd qddot_dir = qddot_dirs.block(0, idir, model.qdot_size, 1);
		Vector3d res_hd  = CalcPointAcceleration(model, q + h * q_dir,
				qdot + h * qdot_dir, qddot + h * qddot_dir, body_id,
				point_position, true);
		fd_derivative.block<3,1>(0, idir) = (res_hd - ref) / h;
	}
	return ref;
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
    unsigned int ndirs = Q_dirs.cols();
    ad_model.resize_directions(ndirs);

    // temporary evaluation at current point
    CalcPointJacobian (
        model, Q, body_id, point_position, G, update_kinematics
    );

    double h = 1e-8;
    MatrixNd G_temp = MatrixNd::Zero (3, model.dof_count);
    for (unsigned int idir = 0; idir < ndirs; idir++) {
        VectorNd Q_dir = Q_dirs.block(0, idir, model.q_size, 1);

        // temporary evaluation at perturbed point
        CalcPointJacobian (
            model, Q + h * Q_dir, body_id, point_position, G_temp, update_kinematics
        );

        // calculate finite difference
        G_dirs[idir] = (G_temp - G) / h;
    }
}

// -----------------------------------------------------------------------------
} // Namespace FD
// -----------------------------------------------------------------------------



//    void fd_UpdateKinematics (
//            Model &model,
//            const Math::VectorNd &q,
//            const Math::VectorNd &q_dirs,
//            const Math::VectorNd &qdot,
//            const Math::VectorNd &qdot_dirs,
//            const Math::VectorNd &qddot,
//            const Math::VectorNd &qddot_dirs,
//            std::vector<std::vector<SpatialMatrix> > &fd_X_lambda,
//            std::vector<std::vector<SpatialMatrix> > &fd_X_base,
//            std::vector<std::vector<SpatialVector> > &fd_a,
//            std::vector<std::vector<SpatialVector> > &fd_v,
//            std::vector<std::vector<SpatialVector> > &fd_c
//    ) {
//        unsigned int ndirs = q_dirs.cols();
//        double h = 1.0e-8;

//        std::vector<SpatialMatrix> ref_X_lambda (model.mBodies.size(), SpatialMatrix::Zero());
//        std::vector<SpatialMatrix> ref_X_base (model.mBodies.size(), SpatialMatrix::Zero());
//        std::vector<SpatialVector> ref_a (model.mBodies.size(), SpatialVector::Zero());
//        std::vector<SpatialVector> ref_v (model.mBodies.size(), SpatialVector::Zero());
//        std::vector<SpatialVector> ref_c (model.mBodies.size(), SpatialVector::Zero());

//        // evaluate y(t)
//        UpdateKinematics (model, q, qdot, qddot);
//        for (unsigned int i = 0; i < model.mBodies.size(); i++) {
//            ref_X_lambda[i] = model.X_lambda[i].toMatrix();
//            ref_X_base[i] = model.X_base[i].toMatrix();
//            ref_a[i] = model.a[i];
//            ref_v[i] = model.v[i];
//            ref_c[i] = model.c[i];
//        }

//        for (unsigned int j = 0; j < ndirs; j++) {
//            VectorNd q_dir = q_dirs.block(0, j, model.q_size, 1);
//            VectorNd qdot_dir = qdot_dirs.block(0, j, model.q_size, 1);
//            VectorNd qddot_dir = qddot_dirs.block(0, j, model.q_size, 1);

//            // evaluate y(t+h*d)
//            UpdateKinematics (
//                model,
//                q + h * q_dir,
//                qdot + h * qdot_dir,
//                qddot + h * qddot_dir
//            );

//            for (unsigned int i = 0; i < model.mBodies.size(); i++) {
//                fd_X_lambda[i][j] = (model.X_lambda[i].toMatrix() - ref_X_lambda[i]) / h;
//                fd_X_base[i][j] = (model.X_base[i].toMatrix() - ref_X_base[i]) / h;
//                fd_a[i][j] = (model.a[i] - ref_a[i]) / h;
//                fd_v[i][j] = (model.v[i] - ref_v[i]) / h;
//                fd_c[i][j] = (model.c[i] - ref_c[i]) / h;
//            }
//        }
//    };

//    RBDL_DLLAPI
//    Vector3d fd_CalcPointAcceleration (
//            Model &model,
//            const Math::VectorNd &q,
//            const Math::VectorNd &q_dirs,
//            const Math::VectorNd &qdot,
//            const Math::VectorNd &qdot_dirs,
//            const Math::VectorNd &qddot,
//            const Math::VectorNd &qddot_dirs,
//            unsigned int body_id,
//            const Math::Vector3d &point_position,
//            const Math::MatrixNd &fd_derivative,
//            bool update_kinematics = true
//    ) {
//        // evaluate y(t)
//        Vector3d ref = CalcPointAcceleration (model, q, qdot, qddot, body_id, point_position);

//        unsigned int ndirs = q_dirs.cols();
//        double h = 1.0e-8;

//        for (unsigned int j = 0; j < ndirs; j++) {
//            VectorNd q_dir = q_dirs.block(0, j, model.q_size, 1);
//            VectorNd qdot_dir = qdot_dirs.block(0, j, model.q_size, 1);
//            VectorNd qddot_dir = qddot_dirs.block(0, j, model.q_size, 1);

//            // evaluate y(t+h*d)
//            Vector3d res_hd = CalcPointAcceleration (
//                model,
//                q + h * q_dir,
//                qdot + h * qdot_dir,
//                qddot + h * qddot_dir,
//                body_id, point_position
//            );

//            //fd_derivative.block<3,1>(0, j) = (res_hd - ref) / h;
//        }

//        return ref;
//    };

// -----------------------------------------------------------------------------
} // Namespace RigidBodyDynamics
// -----------------------------------------------------------------------------
