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
      /*
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
      */

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
      ed_model.f[i].setZero();
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
        model.f[i], ed_model.f[i]
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

// -----------------------------------------------------------------------------
} // ED
// -----------------------------------------------------------------------------
} // RigidBodyDynamics
// -----------------------------------------------------------------------------
