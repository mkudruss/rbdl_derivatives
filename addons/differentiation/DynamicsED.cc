/*
 * RBDL - Rigid Body Dynamics Library
 * Copyright (c) 2011-2015 Martin Felis <martin.felis@iwr.uni-heidelberg.de>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#include <iostream>

#include "JointAD.h"
#include "SpatialAlgebraOperatorsAD.h"

#include "DynamicsED.h"
#include "JointED.h"
#include "SpatialAlgebraOperatorsED.h"

using std::cout;
using std::cerr;
using std::endl;
using std::vector;

using namespace RigidBodyDynamics::Math;

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------
namespace ED {
// -----------------------------------------------------------------------------

RBDL_DLLAPI
void InverseDynamics(
  RigidBodyDynamics::Model &model,
  EDModel &ed_model,
  Math::VectorNd const &q,
  Math::MatrixNd const &q_dirs,
  Math::VectorNd const &qdot,
  Math::MatrixNd const &qdot_dirs,
  Math::VectorNd const &qddot,
  Math::MatrixNd const &qddot_dirs,
  Math::VectorNd &tau,
  Math::MatrixNd &ed_tau,
  std::vector<Math::SpatialVector> const *f_ext,
  std::vector<std::vector<Math::SpatialVector> > const *f_ext_dirs
) {

  const unsigned int ndirs = q_dirs.cols();
  assert(ndirs == qdot_dirs.cols());
  assert(ndirs == qddot_dirs.cols());
  ed_model.resize_directions(ndirs);

  // Reset the velocity of the root body
  // nominal evaluation
  model.v[0].setZero ();
  model.a[0].set (0., 0., 0., -model.gravity[0], -model.gravity[1], -model.gravity[2]);
  // derivative evaluation
  ed_model.v_q[0].setZero ();
  ed_model.v_qdot[0].setZero ();
  ed_model.v_qddot[0].setZero ();
  ed_model.a_q[0].setZero ();
  ed_model.a_qdot[0].setZero ();
  ed_model.a_qddot[0].setZero ();

  for (unsigned int i = 1; i < model.mBodies.size(); i++) {
    unsigned int q_index = model.mJoints[i].q_index;
    unsigned int lambda = model.lambda[i];

    jcalc (model, i, q, qdot);

    // nominal evaluation
    SpatialVector temp = model.X_lambda[i].apply(model.v[model.lambda[i]]);
    model.v[i] = temp + model.v_J[i];
    // derivative evaluation
    // d v[i] / d q
    ed_model.v_q[i].leftCols(ndirs)
        = crossm(temp)*model.S[i]*q_dirs.row(model.mJoints[i].q_index)
        + model.X_lambda[i].toMatrix()*ed_model.v_q[model.lambda[i]].leftCols(ndirs);
      // TODO exploit structure from S
//      = crossm(temp)*model.S[i]*qdot_dirs.row(model.mJoints[i].q_index)
//      + model.X_lambda[i].toMatrix()*ed_model.v_q[lambda].leftCols(ndirs);
    // d v[i] / d qdot
    ed_model.v_qdot[i].leftCols(ndirs)
        = model.X_lambda[i].toMatrix()*ed_model.v_qdot[model.lambda[i]].leftCols(ndirs)
                + model.S[i]*qdot_dirs.row(model.mJoints[i].q_index);
//      = model.X_lambda[i].toMatrix()*ed_model.v_qdot[lambda].leftCols(ndirs)
//      // + model.X_lambda[i].apply(model.v[lambda])
//      // d v_J(q, qdot) / d qdot * qdot_dir = S_i(q) * qdot_dir = ed_model
//      + model.S[i]*qdot_dirs.row(model.mJoints[i].q_index);
      // model.c[i] = model.c_J[i] + crossm(model.v[i],model.v_J[i]);
    // d v[i] / d qddot = 0
    ed_model.v_qddot[i].leftCols(ndirs).setZero();

    // nominal evaluation
    model.c[i] = model.c_J[i] + crossm(model.v[i],model.v_J[i]);
    // derivative evaluation
    // d c[i] / d q
    ed_model.c_q[i].leftCols(ndirs)
      = -crossm(model.v_J[i]) * ed_model.v_q[i].leftCols(ndirs);
    // d c[i] / d qdot
    ed_model.c_qdot[i].leftCols(ndirs)
      = crossm(model.v[i])*model.S[i]*qdot_dirs.row(model.mJoints[i].q_index)
        -crossm(model.v_J[i]) * ed_model.v_qdot[i].leftCols(ndirs);
    // d c[i] / d qddot = 0
    ed_model.c_qddot[i].leftCols(ndirs).setZero();

    if(model.mJoints[i].mJointType != JointTypeCustom){
      if (model.mJoints[i].mDoFCount == 1) {
        // nominal evaluation
        temp = model.X_lambda[i].apply(model.a[model.lambda[i]]);
        model.a[i] = temp + model.c[i] + model.S[i] * qddot[q_index];
        // derivative evaluation
        // d a[i] / d q
        ed_model.a_q[i].leftCols(ndirs)
          = crossm(temp)*model.S[i]*q_dirs.row(model.mJoints[i].q_index)
          + model.X_lambda[i].toMatrix()*ed_model.a_q[lambda].leftCols(ndirs)
          + ed_model.c_q[i].leftCols(ndirs);
          // + model.S[i] * qddot[q_index];
        // d a[i] / d qdot
        ed_model.a_qdot[i].leftCols(ndirs)
          = model.X_lambda[i].toMatrix()*ed_model.a_qdot[lambda].leftCols(ndirs)
          + ed_model.c_qdot[i].leftCols(ndirs);
        // d a[i] / d qddot
        ed_model.a_qddot[i].leftCols(ndirs)
          = model.X_lambda[i].toMatrix()* ed_model.a_qddot[lambda].leftCols(ndirs)
          + ed_model.c_qddot[i].leftCols(ndirs)
          + model.S[i] * qddot_dirs.row(model.mJoints[i].q_index).leftCols(ndirs);

      } else if (model.mJoints[i].mDoFCount == 3) {
        cerr << "Multi-dof not supported." << endl;
        abort();
      }
    }else if(model.mJoints[i].mJointType == JointTypeCustom){
      cerr << __FILE__ << " " << __LINE__
           << ": Custom joints not supported." << endl;    abort();
      abort();
    }

    if (!model.mBodies[i].mIsVirtual) {
      // nominal evaluation
      ed_model.h[i] = model.I[i] * model.v[i];
      // derivative evaluation
      // d a[i] / d q
      ed_model.h_q[i].leftCols(ndirs)
        = model.I[i].toMatrix()*ed_model.v_q[i].leftCols(ndirs);
      // d a[i] / d qdot
      ed_model.h_qdot[i].leftCols(ndirs)
        = model.I[i].toMatrix()*ed_model.v_qdot[i].leftCols(ndirs);
      // d a[i] / d qddot
      ed_model.h_qddot[i].leftCols(ndirs)
        = model.I[i].toMatrix()*ed_model.v_qddot[i].leftCols(ndirs);

      // nominal evaluation
      model.f[i] = model.I[i] * model.a[i] + crossf(model.v[i], ed_model.h[i]);
      // derivative evaluation

      // d a[i] / d q
      ed_model.f_q[i].leftCols(ndirs)
        = model.I[i].toMatrix() * ed_model.a_q[i].leftCols(ndirs)
           + crossf_rhs(ed_model.h[i]).transpose() * ed_model.v_q[i].leftCols(ndirs)
           + crossf(model.v[i]) * ed_model.h_q[i].leftCols(ndirs);
      // d a[i] / d qdot
      ed_model.f_qdot[i].leftCols(ndirs)
        = model.I[i].toMatrix() * ed_model.a_qdot[i].leftCols(ndirs)
          + crossf_rhs(ed_model.h[i]).transpose() * ed_model.v_qdot[i].leftCols(ndirs)
          + crossf(model.v[i]) * ed_model.h_qdot[i].leftCols(ndirs);
      // d a[i] / d qddot
      ed_model.f_qddot[i].leftCols(ndirs)
        = model.I[i].toMatrix() * ed_model.a_qddot[i].leftCols(ndirs)
          + crossf_rhs(ed_model.h[i]).transpose() * ed_model.v_qddot[i].leftCols(ndirs)
          + crossf(model.v[i]) * ed_model.h_qddot[i].leftCols(ndirs);

    } else {
      // nominal evaluation
      model.f[i].setZero();
      // derivative evaluation
      // d a[i] / d q
      ed_model.f_q[i].leftCols(ndirs).setZero();
      // d a[i] / d qdot
      ed_model.f_qdot[i].leftCols(ndirs).setZero();
      // d a[i] / d qddot
      ed_model.f_qddot[i].leftCols(ndirs).setZero();
    }
  }

  if (f_ext != NULL) {
    cerr << __FILE__ << " " << __LINE__
         << ": External forces are not allowed." << endl;    abort();
    abort();
  }

  for (unsigned int i = model.mBodies.size() - 1; i > 0; i--) {
    if(model.mJoints[i].mJointType != JointTypeCustom){
      if (model.mJoints[i].mDoFCount == 1) {
        const unsigned int q_index = model.mJoints[i].q_index;
        // nominal evaluation
        tau[q_index] = model.S[i].dot(model.f[i]);
        // derivative evaluation
        // d tau [i] = d tau [i] / d q + d tau [i] / d qdot + d tau [i] / d qddot
        ed_tau.row(q_index)
          = model.S[i].transpose() * ed_model.f_q[i].leftCols(ndirs)
          + model.S[i].transpose() * ed_model.f_qdot[i].leftCols(ndirs)
          + model.S[i].transpose() * ed_model.f_qddot[i].leftCols(ndirs);

      } else if (model.mJoints[i].mDoFCount == 3) {
        cerr << "Multi-dof not supported." << endl;
        abort();
      }
    } else if (model.mJoints[i].mJointType == JointTypeCustom) {
      cerr << __FILE__ << " " << __LINE__
           << ": Custom joints not supported." << endl;    abort();
      abort();
    }

    if (model.lambda[i] != 0) {
      // nominal evaluation
      SpatialVector temp = model.X_lambda[i].applyTranspose(model.f[i]);
      model.f[model.lambda[i]] = model.f[model.lambda[i]] + temp;
      // derivative evaluation
      // d a[i] / d q
      ed_model.f_q[model.lambda[i]].leftCols(ndirs)
        += model.X_lambda[i].toMatrixTranspose()
          * crossf_rhs(model.f[i]).transpose()
          * model.S[i]*q_dirs.row(model.mJoints[i].q_index).leftCols(ndirs)
        + model.X_lambda[i].toMatrixTranspose()*ed_model.f_q[i].leftCols(ndirs);
      // d a[i] / d qdot
      ed_model.f_qdot[model.lambda[i]].leftCols(ndirs)
        += model.X_lambda[i].toMatrixTranspose()*ed_model.f_qdot[i].leftCols(ndirs);
      // d a[i] / d qddot
      ed_model.f_qddot[model.lambda[i]].leftCols(ndirs)
        += model.X_lambda[i].toMatrixTranspose()*ed_model.f_qddot[i].leftCols(ndirs);
    }
  }

  return;
}

RBDL_DLLAPI void NonlinearEffects (
    Model & model,
    EDModel & ed_model,
    const Math::VectorNd & q,
    const Math::MatrixNd & q_dirs,
    const Math::VectorNd & qdot,
    const Math::MatrixNd & qdot_dirs,
    Math::VectorNd & tau,
    Math::MatrixNd & ed_tau
) {
  LOG << "-------- " << __func__ << " --------" << std::endl;
  const unsigned int ndirs = q_dirs.cols();
  assert(ndirs == qdot_dirs.cols());
  // assert(ndirs == qddot_dirs.cols());
  ed_model.resize_directions(ndirs);

  SpatialVector spatial_gravity (0., 0., 0., -model.gravity[0], -model.gravity[1], -model.gravity[2]);

  // Reset the velocity of the root body
  // nominal evaluation
  model.v[0].setZero();
  model.a[0] = spatial_gravity;
  model.f[0].setZero();
  // derivative evaluation
  ed_model.v_q[0].leftCols(ndirs).setZero ();
  ed_model.v_qdot[0].leftCols(ndirs).setZero ();
  ed_model.v_qddot[0].leftCols(ndirs).setZero ();
  ed_model.a_q[0].leftCols(ndirs).setZero ();
  ed_model.a_qdot[0].leftCols(ndirs).setZero ();
  ed_model.a_qddot[0].leftCols(ndirs).setZero ();
  ed_model.f[0].leftCols(ndirs).setZero();

  for (unsigned int i = 1; i < model.mJointUpdateOrder.size(); i++) {
    jcalc (model, model.mJointUpdateOrder[i], q, qdot);
  }

  for (unsigned int i = 1; i < model.mBodies.size(); i++) {
    if (model.lambda[i] == 0) {
      // nominal evaluation
      model.v[i] = model.v_J[i];
      // derivative evaluation
      // d v[i] / d q
      ed_model.v_q[i].leftCols(ndirs).setZero();
      // d v[i] / d qdot
      // TODO just fill correct row instead of multiplying
      ed_model.v_qdot[i].leftCols(ndirs) = model.S[i]*qdot_dirs.row(model.mJoints[i].q_index);

      // nominal evaluation
      model.a[i] = model.X_lambda[i].apply(spatial_gravity);
      // derivative evaluation
      // d a[i] / d q
      ed_model.a_q[i].leftCols(ndirs)
        = crossm(model.a[i])*model.S[i]*q_dirs.row(model.mJoints[i].q_index);
      // d a[i] / d qdot
      ed_model.a_qdot[i].leftCols(ndirs).setZero();
    } else {
      // nominal evaluation
      SpatialVector temp = model.X_lambda[i].apply(model.v[model.lambda[i]]);
      model.v[i] = temp + model.v_J[i];
      // derivative evaluation
      // d v[i] / d q
      ed_model.v_q[i].leftCols(ndirs)
        = crossm(temp)*model.S[i]*q_dirs.row(model.mJoints[i].q_index)
        + model.X_lambda[i].toMatrix()*ed_model.v_q[model.lambda[i]].leftCols(ndirs);
      // d v[i] / d qdot
      ed_model.v_qdot[i].leftCols(ndirs)
        = model.X_lambda[i].toMatrix()*ed_model.v_qdot[model.lambda[i]].leftCols(ndirs)
        + model.S[i]*qdot_dirs.row(model.mJoints[i].q_index);

      // nominal evaluation
      model.c[i] = model.c_J[i] + crossm(model.v[i],model.v_J[i]);
      // derivative evaluation
      // d c[i] / d q
      ed_model.c_q[i].leftCols(ndirs) = -crossm(model.v_J[i])*ed_model.v_q[i].leftCols(ndirs);
      // d c[i] / d qdot
      ed_model.c_qdot[i].leftCols(ndirs)
        = crossm(model.v[i])*model.S[i]*qdot_dirs.row(model.mJoints[i].q_index)
        - crossm(model.v_J[i])*ed_model.v_qdot[i].leftCols(ndirs);

      // nominal evaluation
      temp = model.X_lambda[i].apply(model.a[model.lambda[i]]);
      model.a[i] = temp + model.c[i];
      // derivative evaluation
      // d a[i] / d q
      ed_model.a_q[i].leftCols(ndirs)
        = crossm(temp)*model.S[i]*q_dirs.row(model.mJoints[i].q_index)
        + model.X_lambda[i].toMatrix()*ed_model.a_q[model.lambda[i]].leftCols(ndirs)
        + ed_model.c_q[i].leftCols(ndirs);
      // d a[i] / d qdot
      ed_model.a_qdot[i].leftCols(ndirs)
        = model.X_lambda[i].toMatrix()*ed_model.a_qdot[model.lambda[i]].leftCols(ndirs)
        + ed_model.c_qdot[i].leftCols(ndirs);
    }

    if (!model.mBodies[i].mIsVirtual) {
      // nominal evaluation
      ed_model.h[i] = model.I[i]*model.v[i];
      // derivative evaluation
      // d h[i] / d q
      ed_model.h_q[i].leftCols(ndirs)
        = model.I[i].toMatrix()*ed_model.v_q[i].leftCols(ndirs);
      // d h[i] / d qdot
      ed_model.h_qdot[i].leftCols(ndirs)
        = model.I[i].toMatrix()*ed_model.v_qdot[i].leftCols(ndirs);

      Math::MatrixNd cross_rhs_hi_T = crossf_rhs(ed_model.h[i]).transpose();

      // nominal evaluation
      model.f[i] = model.I[i] * model.a[i]
        + crossf(model.v[i], ed_model.h[i]);
      // derivative evaluation
      // d f[i] / d q

      ed_model.f_q[i].leftCols(ndirs)
        = model.I[i].toMatrix()*ed_model.a_q[i].leftCols(ndirs)
                + cross_rhs_hi_T*ed_model.v_q[i].leftCols(ndirs) //+ b
                // !!! Note by PM !!!
                // Replaced
                //- crossf(ed_model.h[i])*ed_model.v_q[i].leftCols(ndirs)
                // by d because spatial force cross product is not skew-symmetric
                + crossf(model.v[i])*ed_model.h_q[i].leftCols(ndirs);



      // d f[i] / d qdot
      ed_model.f_qdot[i].leftCols(ndirs)
        = model.I[i].toMatrix()*ed_model.a_qdot[i].leftCols(ndirs)
          + crossf_rhs(ed_model.h[i]).transpose()*ed_model.v_qdot[i].leftCols(ndirs) //d
          // !!! Note by PM !!!
          // Replaced
          // - crossf(ed_model.h[i])*ed_model.v_qdot[i].leftCols(ndirs)
          // by d because spatial force cross product is not skew-symmetric
          + crossf(model.v[i])*ed_model.h_qdot[i].leftCols(ndirs);
    } else {
      // nominal evaluation
      model.f[i].setZero();
      // derivative evaluation
      // d f[i] / d q
      ed_model.f_q[i].leftCols(ndirs).setZero();
      // d f[i] / d qdot
      ed_model.f_qdot[i].leftCols(ndirs).setZero();
    }
  }

  for (unsigned int i = model.mBodies.size() - 1; i > 0; i--) {
    if(model.mJoints[i].mJointType != JointTypeCustom){
      if (model.mJoints[i].mDoFCount == 1) {
        const unsigned int q_index = model.mJoints[i].q_index;
        // nominal evaluation
        tau[q_index] = model.S[i].dot(model.f[i]);
      // derivative evaluation
      // d tau [i] = d tau [i] / d q + d tau [i] / d qdot
        ed_tau.row(q_index)
          = model.S[i].transpose() * ed_model.f_q[i].leftCols(ndirs)
          + model.S[i].transpose() * ed_model.f_qdot[i].leftCols(ndirs);

      } else if (model.mJoints[i].mDoFCount == 3) {
        cerr << "Multi-dof not supported." << endl;
        abort();
      }
    } else if(model.mJoints[i].mJointType == JointTypeCustom) {
      cerr << __FILE__ << " " << __LINE__
           << ": Custom joints not supported." << endl;    abort();
      abort();
    }

    if (model.lambda[i] != 0) {
      // nominal evaluation
      SpatialVector temp = model.X_lambda[i].applyTranspose(model.f[i]);
      model.f[model.lambda[i]] += temp;
      // derivative evaluation

      Eigen::MatrixXd temp_cross2(6, ndirs);
      temp_cross2 = model.X_lambda[i].toMatrixTranspose()
          * crossf_rhs(model.f[i]).transpose()
          * model.S[i]*q_dirs.row(model.mJoints[i].q_index).leftCols(ndirs);

//      Eigen::MatrixXd temp_cross(6, ndirs);
//      for (unsigned int jj = 0; jj < ndirs; jj++) {
//        temp_cross.col(jj) = model.X_lambda[i].applyTranspose(crossf(model.S[i]*q_dirs(model.mJoints[i].q_index, jj), model.f[i]));
//      }
//      if (model.lambda[i] == 1 && model.mBodies.size() > 30) {
//        std::cout << "i = " << i << "  lbd = " << 1 << std::endl;
//        std::cout << temp_cross.col(0).transpose() << "\n" << std::endl;
//        std::cout << (model.X_lambda[i].toMatrixTranspose()*ed_model.f_q[i].col(0)).transpose() << "\n" << std::endl;
//        std::cout << (model.X_lambda[i].toMatrixTranspose()*ed_model.f_qdot[i].col(0)).transpose() << "\n" << std::endl;
//      }
//      if (model.mBodies.size() > 30) {
//          for (unsigned int jj = 0; jj < ndirs; jj++) {
//            temp_cross.col(jj) = model.X_lambda[i].applyTranspose(crossf(model.S[i]*q_dirs(model.mJoints[i].q_index, jj), model.f[i]));
//          }
//      }

      // d f[model.lambda[i]] / d q
      ed_model.f_q[model.lambda[i]].leftCols(ndirs)
        += temp_cross2
        //- crossf(temp)*model.S[i]*q_dirs.row(model.mJoints[i].q_index)
        + model.X_lambda[i].toMatrixTranspose()*ed_model.f_q[i].leftCols(ndirs);
      // d f[model.lambda[i]] / d qdot
      ed_model.f_qdot[model.lambda[i]].leftCols(ndirs)
        += model.X_lambda[i].toMatrixTranspose()*ed_model.f_qdot[i].leftCols(ndirs);
    }
  }
}


/*
RBDL_DLLAPI
void InverseDynamics(
  RigidBodyDynamics::Model &model,
  EDModel &ed_model,
  Math::VectorNd const &q,
  Math::MatrixNd const &q_dirs,
  Math::VectorNd const &qdot,
  Math::MatrixNd const &qdot_dirs,
  Math::VectorNd const &qddot,
  Math::MatrixNd const &qddot_dirs,
  Math::VectorNd &tau,
  Math::MatrixNd &tau_dirs,
  std::vector<Math::SpatialVector> const *f_ext,
  std::vector<std::vector<Math::SpatialVector> > const *f_ext_dirs
) {

  // receive number of directions
  const unsigned ndirs = q_dirs.cols();
  assert(ndirs == qdot_dirs.cols());
  assert(ndirs == qddot_dirs.cols());

  // resize container if necessary
  ed_model.resize_directions(ndirs);

  // Reset the velocity of the root body
  // derivative evaluation
  ed_model.v[0].leftCols(ndirs).setZero();
  // nominal evaluation
  model.v[0].setZero();

  // Set the root body acceleration to gravity
  // derivative evaluation
  ed_model.a[0].leftCols(ndirs).setZero();
  // nominal evaluation
  model.a[0].set (0., 0., 0., -model.gravity[0], -model.gravity[1], -model.gravity[2]);

  // RNEA forward propagation
  for (unsigned int i = 1; i < model.mBodies.size(); i++) {
    unsigned int q_index = model.mJoints[i].q_index;
    unsigned int lambda = model.lambda[i];

    // derivative evaluation
    ED::jcalc (model, ed_model, i, q, q_dirs, qdot, qdot_dirs);
    // nominal evaluation
    // NOTE joints are already calculated in ED::jcalc
    // jcalc (model, i, q, qdot);

    // derivative evaluation
    Math::ED::X_apply_v_add_u (
      model.v[i], ed_model.v[i],
      model.X_lambda[i], ed_model.X_lambda[i],
      model.v[lambda], ed_model.v[lambda],
      model.v_J[i], ed_model.v_J[i],
      ndirs
    );
    // nominal evaluation
    // NOTE full action is performed in AD code
    // model.v[i] = model.X_lambda[i].apply(model.v[lambda])
    //   + model.v_J[i];

    // derivative evaluation
    Math::ED::crossm(
      model.c[i], ed_model.c[i],
      model.v[i], ed_model.v[i],
      model.v_J[i], ed_model.v_J[i],
      ndirs
    );
    ed_model.c[i].leftCols(ndirs) += ed_model.c_J[i].leftCols(ndirs);
    // nominal evaluation
    model.c[i] += model.c_J[i]; // + crossm (model.v[i], model.v_J[i]);

    if(model.mJoints[i].mJointType != JointTypeCustom){
      if (model.mJoints[i].mDoFCount == 1) {
        // derivative evaluation
        Math::ED::X_apply_v_add_u (
          model.a[i], ed_model.a[i],
          model.X_lambda[i], ed_model.X_lambda[i],
          model.a[lambda], ed_model.a[lambda],
          model.c[i], ed_model.c[i],
          ndirs
        );
        int vidx = -1;
        for (int idx = 0; idx < 6; idx++) {
          if (model.S[i][idx] == 1.) {
            vidx = idx;
          }
        }
        ed_model.a[i].leftCols(ndirs).row(vidx) += qddot_dirs.row(q_index);
        // nominal evaluation
        model.a[i] += model.S[i] * qddot(q_index);
        // model.a[i] =  model.X_lambda[i].apply(model.a[lambda])
        //   + model.c[i]
        //   + model.S[i] * qddot[q_index];
      } else if (model.mJoints[i].mDoFCount == 3) {
        cerr << __FILE__ << " " << __LINE__ << ":"
             << " Multi-DoF joint not supported." << endl;
        abort();
      }
    }else if(model.mJoints[i].mJointType == JointTypeCustom){
      cerr << __FILE__ << " " << __LINE__ << ":"
           << " Custom joint not supported." << endl;
      abort();
    }

    if (!model.mBodies[i].mIsVirtual) {
      // derivative evaluation
      // NOTE inertia model.I[i] is constant
      // NOTE save ndir times computation of Iv
      // compute model.f = I * model.v
      //      ed_model.f = I * ed_model.v

      SpatialVector Iv;
      Math::ED::SRBI_apply_v (
        Iv, ed_model.Iv,
        model.I[i],
        model.v[i], ed_model.v[i],
        ndirs
      );
      // compute model.f = v x* Iv
      //      ed_model.f = v x* ed_model.Iv
      Math::ED::crossf (
        model.f[i], ed_model.f[i],
        model.v[i], ed_model.v[i],
        Iv, ed_model.Iv,
        ndirs
      );
      // compute model.f += I * model.a
      //      ed_model.f += I * ed_model.a
      Math::ED::add_SRBI_apply_v (
        model.f[i], ed_model.f[i],
        model.I[i],
        model.a[i], ed_model.a[i],
        ndirs
      );
      // nominal evaluation
      // model.f[i] += model.I[i] * model.a[i];// + crossf (model.v[i], Iv);

      for(unsigned int j = 0; j < ndirs; j++) {
        ed_model.f[i].col(j) = model.I[i] * ed_model.a[i].col(j)
          + crossf(ed_model.v[i].col(j),model.I[i] * model.v[i])
          + crossf(model.v[i],model.I[i] * ed_model.v[i].col(j));
      }
      // nominal evaluation
      model.f[i] = model.I[i] * model.a[i]
          + crossf(model.v[i],model.I[i] * model.v[i]);

    } else {
      // derivative evaluation
      ed_model.f[i].leftCols(ndirs).setZero();
      // nominal evaluation
      model.f[i].setZero();
    }
  }

  if (f_ext != NULL) {
    cerr << __FILE__ << " " << __LINE__ << ":"
         << " external forces are not supported." << endl;
    abort();
  }

  for (unsigned int i = model.mBodies.size() - 1; i > 0; i--) {
    if(model.mJoints[i].mJointType != JointTypeCustom){
      if (model.mJoints[i].mDoFCount == 1) {
        // derivative evaluation
        tau_dirs.row(model.mJoints[i].q_index)
          = model.S[i].transpose()*ed_model.f[i].leftCols(ndirs);
            // NOTE ed_model.S[i][idir] = 0 for all joint types
            // NOTE not true for multidof joints
            // + ed_model.S[i][idir].dot(model.f[i]);
        // nominal evaluation
        tau[model.mJoints[i].q_index] = model.S[i].dot(model.f[i]);
      } else if (model.mJoints[i].mDoFCount == 3) {
        cerr << __FILE__ << " " << __LINE__ << ":"
             << " Multi-DoF joint not supported." << endl;
        abort();
      }
    } else if (model.mJoints[i].mJointType == JointTypeCustom) {
      cerr << __FILE__ << " " << __LINE__ << ":"
           << " Custom joint not supported." << endl;
      abort();
    }

    if (model.lambda[i] != 0) {
      // derivative evaluation
      Math::ED::inplace_X_applyTranspose_f (
        model.f[model.lambda[i]], ed_model.f[model.lambda[i]],
        model.X_lambda[i], ed_model.X_lambda[i],
        model.f[i], ed_model.f[i],
        ndirs
      );
      // for (unsigned idir = 0; idir < ndirs; idir++) {
      //   ed_model.f[model.lambda[i]][idir] =
      // }
      // // nominal evaluation
      // model.f[model.lambda[i]] = model.f[model.lambda[i]]
      //   + model.X_lambda[i].applyTranspose(model.f[i]);
    }
  }

  return;
}
*/

RBDL_DLLAPI
void CompositeRigidBodyAlgorithm (
  Model &model,
  EDModel &ed_model,
  Math::VectorNd const & q,
  Math::MatrixNd const & q_dirs,
  Math::MatrixNd & H,
  std::vector<Math::MatrixNd> & H_dirs,
  bool update_kinematics
)
{
  assert (H.rows() == model.dof_count && H.cols() == model.dof_count);

  // get number if directions
  const size_t ndirs = q_dirs.cols();
  ed_model.resize_directions(ndirs);

  for (unsigned int i = 1; i < model.mBodies.size(); i++)
  {
    if (update_kinematics)
    {
      jcalc_X_lambda_S (model, i, q);
    }
    // nominal evaluation
    model.Ic[i] = model.I[i];
    // derivative evaluation
    for (unsigned idir = 0; idir < ndirs; idir++)
    {
      ed_model.Ic[i][idir].setZero();
    }
  }

  for (unsigned int i = model.mBodies.size() - 1; i > 0; i--)
  {
    if (model.lambda[i] != 0)
    {
      SpatialRigidBodyInertia temp
        = model.X_lambda[i].applyTranspose(model.Ic[i]);
      // nominal evaluation
      model.Ic[model.lambda[i]] = model.Ic[model.lambda[i]] + temp;
      // derivative evaluation
      for (unsigned idir = 0; idir < ndirs; idir++) {
      //   // TODO make this fast
        SpatialRigidBodyInertia rbi;
        rbi.createFromMatrix(ed_model.Ic[i][idir]);

        ed_model.Ic[model.lambda[i]][idir]
          = ed_model.Ic[model.lambda[i]][idir]
          + crossf(model.S[i]*q_dirs(model.mJoints[i].q_index, idir))*temp.toMatrix()
          - temp.toMatrix()*crossm(model.S[i]*q_dirs(model.mJoints[i].q_index,idir))
          + model.X_lambda[i].applyTranspose(rbi).toMatrix()
        ;
      }
    }

    unsigned int dof_index_i = model.mJoints[i].q_index;

    if (model.mJoints[i].mDoFCount == 1
        && model.mJoints[i].mJointType != JointTypeCustom)
    {
      // nominal evaluation
      SpatialVector F             = model.Ic[i] * model.S[i];
      H(dof_index_i, dof_index_i) = model.S[i].dot(F);
      // derivative evaluation
      for (unsigned idir = 0; idir < ndirs; idir++)
      {
        ed_model.F[i].col(idir) = ed_model.Ic[i][idir] * model.S[i];
        H_dirs[idir](dof_index_i, dof_index_i)
          = model.S[i].dot(ed_model.F[i].col(idir));
      }

      unsigned int j = i;
      unsigned int dof_index_j = dof_index_i;

      while (model.lambda[j] != 0) {
        // nominal evaluation
        F = model.X_lambda[j].applyTranspose(F);
        // derivative evaluation
        ed_model.F[i].leftCols(ndirs)
          = model.X_lambda[i].toMatrixTranspose()*crossf_rhs(F).transpose()*model.S[j]*q_dirs.row(model.mJoints[j].q_index)
          + model.X_lambda[j].toMatrixTranspose()*ed_model.F[i].leftCols(ndirs);

        j = model.lambda[j];
        dof_index_j = model.mJoints[j].q_index;

        if(model.mJoints[j].mJointType != JointTypeCustom)
        {
          if (model.mJoints[j].mDoFCount == 1)
          {
            // nominal evaluation
            H(dof_index_i,dof_index_j) = F.dot(model.S[j]);
            H(dof_index_j,dof_index_i) = H(dof_index_i,dof_index_j);
            // derivative evaluation
            for (unsigned idir = 0; idir < ndirs; idir++) {
              H_dirs[idir](dof_index_i, dof_index_j)
                = ed_model.F[i].col(idir).dot(model.S[j]);
              H_dirs[idir](dof_index_j, dof_index_i)
                = H_dirs[idir](dof_index_i, dof_index_j);
            }
          }
          else if (model.mJoints[j].mDoFCount == 3
          ) {
            cerr << __FILE__ << " " << __LINE__ << ":"
                 << "Multi-DoF joint not supported." << endl;
            abort();
          }
        }
        else if (model.mJoints[j].mJointType == JointTypeCustom
        ){
          cerr << __FILE__ << " " << __LINE__ << ":"
               << " Custom joint not supported." << endl;
          abort();
        }
      }
    }
    else if (model.mJoints[i].mDoFCount == 3
        && model.mJoints[i].mJointType != JointTypeCustom
    ) {
      cerr << __FILE__ << " " << __LINE__ << ":"
           << "Multi-DoF joint not supported." << endl;
      abort();
    }
    else if (model.mJoints[i].mJointType == JointTypeCustom
    ) {
      cerr << __FILE__ << " " << __LINE__ << ":"
           << " Custom joint not supported." << endl;
      abort();
    }

  }
}

// -----------------------------------------------------------------------------
} // ED
// -----------------------------------------------------------------------------
} // RigidBodyDynamics
// -----------------------------------------------------------------------------
