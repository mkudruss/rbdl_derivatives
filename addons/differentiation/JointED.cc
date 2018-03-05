/*
 * RBDL - Rigid Body Dynamics Library
 * Copyright (c) 2011-2015 Martin Felis <martin.felis@iwr.uni-heidelberg.de>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#include <iostream>
#include <limits>
#include <assert.h>

#include "rbdl/Logging.h"

#include "JointED.h"
#include "ModelED.h"
#include "ModelAD.h"
#include "SpatialAlgebraOperatorsED.h"

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------
namespace ED {
// -----------------------------------------------------------------------------

using std::cerr;
using std::endl;
using namespace RigidBodyDynamics::Math;

RBDL_DLLAPI void jcalc (
    RigidBodyDynamics::Model &model,
    RigidBodyDynamics::EDModel &ad_model,
    unsigned int joint_id,
    const VectorNd &q,
    const MatrixNd &q_dirs,
    const VectorNd &qdot,
    const MatrixNd &qdot_dirs
) {
  LOG << "-------- " << __func__ << " --------" << std::endl;

  // exception if we calculate it for the root body
  assert (joint_id > 0);

  // receive number of directions
  const unsigned ndirs = q_dirs.cols();
  assert(ndirs == qdot_dirs.cols());
  // assert(ndirs == qddot_dirs.cols());

  // resize container if necessary
  ad_model.resize_directions(ndirs);

  if (model.mJoints[joint_id].mJointType == JointTypeRevoluteX) {

    // derivative evaluation
    for (unsigned idir = 0; idir < ndirs; idir++) {
      // TODO: Most entries of q_dirs are ZERO !!!
      // TODO implement structure exploitation of zero entries?
      ad_model.X_J[joint_id][idir] = Math::ED::Xrotx (
        q[model.mJoints[joint_id].q_index],
        q_dirs(model.mJoints[joint_id].q_index, idir)
      );
    }
    // nominal evaluation
    model.X_J[joint_id] = Xrotx (q[model.mJoints[joint_id].q_index]);

    // derivative evaluation
    for (unsigned idir = 0; idir < ndirs; idir++) {
      ad_model.v_J[joint_id][idir][0] = qdot_dirs(model.mJoints[joint_id].q_index, idir);
    }
    // nominal evaluation
    model.v_J[joint_id][0] = qdot[model.mJoints[joint_id].q_index];

  } else if (model.mJoints[joint_id].mJointType == JointTypeRevoluteY) {

    // derivative evaluation
    for (unsigned idir = 0; idir < ndirs; idir++) {
      // TODO: Most entries of q_dirs are ZERO !!!
      // TODO implement structure exploitation of zero entries?
      ad_model.X_J[joint_id][idir] = Math::ED::Xroty (
        q[model.mJoints[joint_id].q_index],
        q_dirs(model.mJoints[joint_id].q_index, idir)
      );
    }
    // nominal evaluation
    model.X_J[joint_id] = Xroty (q[model.mJoints[joint_id].q_index]);

    // derivative evaluation
    for (unsigned idir = 0; idir < ndirs; idir++) {
      ad_model.v_J[joint_id][idir][1] = qdot_dirs(model.mJoints[joint_id].q_index, idir);
    }
    // nominal evaluation
    model.v_J[joint_id][1] = qdot[model.mJoints[joint_id].q_index];

  } else if (model.mJoints[joint_id].mJointType == JointTypeRevoluteZ) {

    // derivative evaluation
    Math::ED::Xrotz(
      model.X_J[joint_id], ad_model.X_J[joint_id],
      q[model.mJoints[joint_id].q_index],
      q_dirs.row(model.mJoints[joint_id].q_index),
      ndirs
    );
    // for (unsigned idir = 0; idir < ndirs; idir++) {
    //   // TODO: Most entries of q_dirs are ZERO !!!
    //   // TODO implement structure exploitation of zero entries?
    //   ad_model.X_J[joint_id][idir] = Math::ED::Xrotz (
    //     q[model.mJoints[joint_id].q_index],
    //     q_dirs(model.mJoints[joint_id].q_index, idir)
    //   );
    // }
    // // nominal evaluation
    // model.X_J[joint_id] = Xrotz (q[model.mJoints[joint_id].q_index]);

    // derivative evaluation
    for (unsigned idir = 0; idir < ndirs; idir++) {
      ad_model.v_J[joint_id][idir][2] = qdot_dirs(model.mJoints[joint_id].q_index, idir);
    }
    // nominal evaluation
    model.v_J[joint_id][2] = qdot[model.mJoints[joint_id].q_index];

  } else if (model.mJoints[joint_id].mDoFCount == 1 &&
      model.mJoints[joint_id].mJointType != JointTypeCustom) {
      cerr << __FILE__ << " " << __LINE__ << ":"
           << " MOEP!." << endl;
      abort();

      model.X_J[joint_id] = jcalc_XJ (model, joint_id, q);
      model.v_J[joint_id] =
        model.S[joint_id] * qdot[model.mJoints[joint_id].q_index];
  } else if (model.mJoints[joint_id].mJointType == JointTypeSpherical) {
    cerr << __FILE__ << " " << __LINE__ << ":"
         << " Multi-DoF joint not supported." << endl;
    abort();
  } else if (model.mJoints[joint_id].mJointType == JointTypeEulerZYX) {
    cerr << __FILE__ << " " << __LINE__ << ":"
         << " Multi-DoF joint not supported." << endl;
    abort();
  } else if (model.mJoints[joint_id].mJointType == JointTypeEulerXYZ) {
    cerr << __FILE__ << " " << __LINE__ << ":"
         << " Multi-DoF joint not supported." << endl;
    abort();
  } else if (model.mJoints[joint_id].mJointType == JointTypeEulerYXZ) {
    cerr << __FILE__ << " " << __LINE__ << ":"
         << " Multi-DoF joint not supported." << endl;
    abort();
  } else if(model.mJoints[joint_id].mJointType == JointTypeTranslationXYZ){
    cerr << __FILE__ << " " << __LINE__ << ":"
         << " Multi-DoF joint not supported." << endl;
    abort();
  } else if (model.mJoints[joint_id].mJointType == JointTypeCustom) {
    cerr << __FILE__ << " " << __LINE__ << ":"
         << " Custom joint not supported." << endl;
    abort();
  } else {
    std::cerr << "Error: invalid joint type " << model.mJoints[joint_id].mJointType << " at id " << joint_id << std::endl;
    abort();
  }

  // derivative evaluation
  for (unsigned idir = 0; idir < ndirs; idir++) {
    ad_model.X_lambda[joint_id][idir] = ad_model.X_J[joint_id][idir] * model.X_T[joint_id];
  }
  // nominal evaluation
  model.X_lambda[joint_id] = model.X_J[joint_id] * model.X_T[joint_id];

  return;
}

RBDL_DLLAPI void jcalc_X_lambda_S (
  RigidBodyDynamics::Model &model,
  ADModel &ad_model,
  unsigned int joint_id,
  const Math::VectorNd & q,
  const Math::MatrixNd & q_dirs
) {
  LOG << "-------- " << __func__ << " --------" << std::endl;

  // exception if we calculate it for the root body
  assert (joint_id > 0);

  // receive number of directions
  const unsigned ndirs = q_dirs.cols();
  // assert(ndirs == qdot_dirs.cols());

  // resize container if necessary
  ad_model.resize_directions(ndirs);

  if (model.mJoints[joint_id].mJointType == JointTypeRevoluteX) {
    std::cerr << "Error: invalid joint type " << model.mJoints[joint_id].mJointType << " at id " << joint_id << std::endl;
    abort();
    // derivative evaluation
    for (unsigned idir = 0; idir < ndirs; idir++) {
      // TODO: Most entries of q_dirs are ZERO !!!
      // TODO implement structure exploitation of zero entries?
      // NOTE model.X_T is constant
      ad_model.X_lambda[joint_id][idir] = Math::ED::Xrotx (
          q[model.mJoints[joint_id].q_index],
          q_dirs(model.mJoints[joint_id].q_index, idir)
        ) * model.X_T[joint_id];
    }
    // nominal evaluation
    model.X_lambda[joint_id] =
      Xrotx (q[model.mJoints[joint_id].q_index]) * model.X_T[joint_id];

    // derivative evaluation
    // NOTE for 1-DoF joints ad_model.S = 0
    // nominal evaluation
    model.S[joint_id] = model.mJoints[joint_id].mJointAxes[0];

  } else if (model.mJoints[joint_id].mJointType == JointTypeRevoluteY) {
    std::cerr << "Error: invalid joint type " << model.mJoints[joint_id].mJointType << " at id " << joint_id << std::endl;
    abort();

    // derivative evaluation
    for (unsigned idir = 0; idir < ndirs; idir++) {
      // TODO: Most entries of q_dirs are ZERO !!!
      // TODO implement structure exploitation of zero entries?
      // NOTE model.X_T is constant
      ad_model.X_lambda[joint_id][idir] = Math::ED::Xroty (
          q[model.mJoints[joint_id].q_index],
          q_dirs(model.mJoints[joint_id].q_index, idir)
        ) * model.X_T[joint_id];
    }
    // nominal evaluation
    model.X_lambda[joint_id] =
      Xroty (q[model.mJoints[joint_id].q_index]) * model.X_T[joint_id];

    // derivative evaluation
    // NOTE for 1-DoF joints ad_model.S = 0
    // nominal evaluation
    model.S[joint_id] = model.mJoints[joint_id].mJointAxes[0];

  } else if (model.mJoints[joint_id].mJointType == JointTypeRevoluteZ) {
    // derivative evaluation
    // NOTE evaluations of cos and sin can be saved
    // NOTE evaluations of cos and sin can be saved
    Math::ED::Xrotz(
      model.X_lambda[joint_id], ad_model.X_lambda[joint_id],
      q[model.mJoints[joint_id].q_index],
      q_dirs.row(model.mJoints[joint_id].q_index),
      model.X_T[joint_id],
      ndirs
    );
    // nominal evaluation
    // model.X_lambda[joint_id] =
    //   Xrotz (q[model.mJoints[joint_id].q_index]) * model.X_T[joint_id];

    // derivative evaluation
    // NOTE for 1-DoF joints ad_model.S = 0
    // nominal evaluation
    model.S[joint_id] = model.mJoints[joint_id].mJointAxes[0];

  } else if (model.mJoints[joint_id].mDoFCount == 1
      && model.mJoints[joint_id].mJointType != JointTypeCustom){
    cerr << __FILE__ << " " << __LINE__ << ":"
         << " this joint is not supported." << endl;
    abort();

    // derivative evaluation
    for (unsigned idir = 0; idir < ndirs; idir++) {
      // TODO: Most entries of q_dirs are ZERO !!!
      // TODO implement structure exploitation of zero entries?
      // NOTE model.X_T is constant
      ad_model.X_lambda[joint_id][idir] =
        jcalc_XJ (model, joint_id, q) * model.X_T[joint_id];
    }
    // nominal evaluation
    model.X_lambda[joint_id] =
      jcalc_XJ (model, joint_id, q) * model.X_T[joint_id];

    // Set the joint axis
    // derivative evaluation
    // NOTE for 1-DoF joints ad_model.S = 0
    // nominal evaluation
    model.S[joint_id] = model.mJoints[joint_id].mJointAxes[0];

  } else if (model.mJoints[joint_id].mJointType == JointTypeSpherical) {
    cerr << __FILE__ << " " << __LINE__ << ":"
         << " Multi-DoF joint not supported." << endl;
    abort();
  } else if (model.mJoints[joint_id].mJointType == JointTypeEulerZYX) {
    cerr << __FILE__ << " " << __LINE__ << ":"
         << " Multi-DoF joint not supported." << endl;
    abort();
  } else if (model.mJoints[joint_id].mJointType == JointTypeEulerXYZ) {
    cerr << __FILE__ << " " << __LINE__ << ":"
         << " Multi-DoF joint not supported." << endl;
    abort();
  } else if (model.mJoints[joint_id].mJointType == JointTypeEulerYXZ ) {
    cerr << __FILE__ << " " << __LINE__ << ":"
         << " Multi-DoF joint not supported." << endl;
    abort();
  } else if (model.mJoints[joint_id].mJointType == JointTypeTranslationXYZ) {
    cerr << __FILE__ << " " << __LINE__ << ":"
         << " Multi-DoF joint not supported." << endl;
    abort();
  } else if (model.mJoints[joint_id].mJointType == JointTypeCustom) {
    cerr << __FILE__ << " " << __LINE__ << ":"
         << " Custom joint not supported." << endl;
    abort();
  } else {
    std::cerr << "Error: invalid joint type!" << std::endl;
    abort();
  }

  return;
}


RBDL_DLLAPI
Math::SpatialTransform jcalc_XJ (
  RigidBodyDynamics::Model &model,
  ADModel &ad_model,
  unsigned int joint_id,
  unsigned int idir,
  const Math::VectorNd &q,
  const Math::MatrixNd &q_dirs
) {
  // exception if we calculate it for the root body
  assert (joint_id > 0);

  if (model.mJoints[joint_id].mDoFCount == 1
      && model.mJoints[joint_id].mJointType != JointTypeCustom) {
    if (model.mJoints[joint_id].mJointType == JointTypeRevolute) {
      std::cerr << "In: " << __func__ << "this joint type is not supported." << std::endl;
      abort();
      return Xrot (q[model.mJoints[joint_id].q_index], Vector3d (
            model.mJoints[joint_id].mJointAxes[0][0],
            model.mJoints[joint_id].mJointAxes[0][1],
            model.mJoints[joint_id].mJointAxes[0][2]
            ));
    } else if (model.mJoints[joint_id].mJointType == JointTypePrismatic) {
      return Math::ED::Xtrans (
        // r
        Vector3d (
            model.mJoints[joint_id].mJointAxes[0][3]
            * q[model.mJoints[joint_id].q_index],
            model.mJoints[joint_id].mJointAxes[0][4]
            * q[model.mJoints[joint_id].q_index],
            model.mJoints[joint_id].mJointAxes[0][5]
            * q[model.mJoints[joint_id].q_index]
        ),
        // r_dir
        Vector3d (
            model.mJoints[joint_id].mJointAxes[0][3]
            * q_dirs(model.mJoints[joint_id].q_index, idir),
            model.mJoints[joint_id].mJointAxes[0][4]
            * q_dirs(model.mJoints[joint_id].q_index, idir),
            model.mJoints[joint_id].mJointAxes[0][5]
            * q_dirs(model.mJoints[joint_id].q_index, idir)
        )
      );
    }
  }
  std::cerr << "Error: invalid joint type!" << std::endl;
  abort();
  return SpatialTransform();
}

// -----------------------------------------------------------------------------
} /* ED */
// -----------------------------------------------------------------------------
} /* RigidBodyDynamics */
// -----------------------------------------------------------------------------
