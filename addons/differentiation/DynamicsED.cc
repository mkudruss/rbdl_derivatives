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
  EDModel &ad_model,
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

  LOG << "-------- " << __func__ << " --------" << std::endl;

  // receive number of directions
  const unsigned ndirs = q_dirs.cols();
  assert(ndirs == qdot_dirs.cols());
  assert(ndirs == qddot_dirs.cols());

  // resize container if necessary
  ad_model.resize_directions(ndirs);

  // Reset the velocity of the root body
  // derivative evaluation
  for (unsigned idir = 0; idir < ndirs; idir++) {
    ad_model.v[0][idir].setZero();
  }
  // nominal evaluation
  model.v[0].setZero();

  // Set the root body acceleration to gravity
  // derivative evaluation
  for (unsigned idir = 0; idir < ndirs; idir++) {
    ad_model.a[0][idir].setZero();
  }
  // nominal evaluation
  model.a[0].set (0., 0., 0., -model.gravity[0], -model.gravity[1], -model.gravity[2]);

  // RNEA forward propagation
  for (unsigned int i = 1; i < model.mBodies.size(); i++) {
    unsigned int q_index = model.mJoints[i].q_index;
    unsigned int lambda = model.lambda[i];

    // derivative evaluation
    // TODO check that code!
    ED::jcalc (model, ad_model, i, q, q_dirs, qdot, qdot_dirs);
    // nominal evaluation
    // NOTE joints are already calculated in ED::jcalc
    // jcalc (model, i, q, qdot);

    // derivative evaluation
    // TODO X v + v inplace
    // for (unsigned idir = 0; idir < ndirs; idir++) {
    //   ad_model.v[i][idir] = ad_model.X_lambda[i][idir].apply(model.v[lambda])
    //     + model.X_lambda[i].apply(ad_model.v[lambda][idir])
    //     + ad_model.v_J[i][idir];
    // }
    // nominal evaluation
    model.v[i] = model.X_lambda[i].apply(model.v[lambda]) + model.v_J[i];

    // derivative evaluation
    // for (unsigned idir = 0; idir < ndirs; idir++) {
    //   ad_model.c[i][idir] = ad_model.c_J[i][idir]
    //     + Math::ED::crossm (
    //       model.v[i], ad_model.v[i][idir],
    //       model.v_J[i], ad_model.v_J[i][idir]
    //     );
    // }
    // nominal evaluation
    model.c[i] = model.c_J[i] + crossm (model.v[i],model.v_J[i]);

    if(model.mJoints[i].mJointType != JointTypeCustom){
      if (model.mJoints[i].mDoFCount == 1) {
        // derivative evaluation
        // TODO X v + v inplace
        // for (unsigned idir = 0; idir < ndirs; idir++) {
        //   ad_model.a[i][idir] =  ad_model.X_lambda[i][idir].apply(model.a[lambda])
        //     + model.X_lambda[i].apply(ad_model.a[lambda][idir])
        //     + ad_model.c[i][idir]
        //     // NOTE ad_model.S[i][idir] = 0 for all joint types
        //     // NOTE not true for multidof joints
        //     // + ad_model.S[i][idir] * qddot[q_index];
        //     + model.S[i] * qddot_dirs(q_index, idir);
        // }
        // nominal evaluation
        model.a[i] =  model.X_lambda[i].apply(model.a[lambda])
          + model.c[i]
          + model.S[i] * qddot[q_index];
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
      // NOTE save ndir times computation of Iv
      SpatialVector Iv = model.I[i] * model.v[i];
      // TODO I a + v* I v
      // for (unsigned idir = 0; idir < ndirs; idir++) {
      //   ad_model.f[i][idir]
      //     // NOTE inertia model.I[i] is constant
      //     = model.I[i] * ad_model.a[i][idir]
      //     + Math::ED::crossf (
      //       model.v[i], ad_model.v[i][idir],
      //       Iv, model.I[i] * ad_model.v[i][idir]
      //     );
      // }
      // nominal evaluation
      model.f[i] = model.I[i] * model.a[i]
        + crossf (model.v[i], Iv);
    } else {
      for (unsigned idir = 0; idir < ndirs; idir++) {
        ad_model.f[i][idir].setZero();
      }
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
        for (unsigned idir = 0; idir < ndirs; idir++) {
          tau_dirs(model.mJoints[i].q_index, idir)
            = model.S[i].dot(ad_model.f[i][idir]);
            // NOTE ad_model.S[i][idir] = 0 for all joint types
            // NOTE not true for multidof joints
            // + ad_model.S[i][idir].dot(model.f[i]);
        }
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
      // TODO
      // for (unsigned idir = 0; idir < ndirs; idir++) {
      //   ad_model.f[model.lambda[i]][idir] = ad_model.f[model.lambda[i]][idir]
      //     + ad_model.X_lambda[i][idir].applyTranspose(model.f[i])
      //     + model.X_lambda[i].applyTranspose(ad_model.f[i][idir]);
      // }
      // nominal evaluation
      model.f[model.lambda[i]] = model.f[model.lambda[i]]
        + model.X_lambda[i].applyTranspose(model.f[i]);
    }
  }

  return;
}

/*
RBDL_DLLAPI
void NonlinearEffects (
    Model & model,
    ADModel & ad_model,
    const Math::VectorNd & q,
    const Math::MatrixNd & q_dirs,
    const Math::VectorNd & qdot,
    const Math::MatrixNd & qdot_dirs,
    Math::VectorNd & tau,
    Math::MatrixNd & tau_dirs
) {
  LOG << "-------- " << __func__ << " --------" << std::endl;

  // receive number of directions
  const unsigned ndirs = q_dirs.cols();
  assert(ndirs == qdot_dirs.cols());
  assert(ndirs == qddot_dirs.cols());

  SpatialVector spatial_gravity (0., 0., 0., -model.gravity[0], -model.gravity[1], -model.gravity[2]);

  // Reset the velocity of the root body
  // derivative evaluation
  for (unsigned idir = 0; idir < ndirs; idir++) {
    ad_model.v[0][idir].setZero();
  }
  // nominal evaluation
  model.v[0].setZero();

  // Set the root body acceleration to gravity
  // derivative evaluation
  for (unsigned idir = 0; idir < ndirs; idir++) {
    ad_model.a[0][idir].setZero();
  }
  // nominal evaluation
  model.a[0] = spatial_gravity;

  for (unsigned int i = 1; i < model.mJointUpdateOrder.size(); i++) {
    // derivative evaluation
    // TODO check that code!
    ED::jcalc (model, ad_model, model.mJointUpdateOrder[i], q, q_dirs, qdot, qdot_dirs);
    // nominal evaluation
    // NOTE joints are already calculated in ED::jcalc
    // jcalc (model, model.mJointUpdateOrder[i], q, qdot);
  }

  for (unsigned int i = 1; i < model.mBodies.size(); i++) {
    if (model.lambda[i] == 0) {
      // derivative evaluation
      for (unsigned idir = 0; idir < ndirs; idir++) {
        ad_model.v[i][idir] = ad_model.v_J[i][idir];
      }
      // nominal evaluation
      model.v[i] = model.v_J[i];

      // derivative evaluation
      for (unsigned idir = 0; idir < ndirs; idir++) {
        // NOTE spatial_gravity is constant
        ad_model.a[i][idir] = ad_model.X_lambda[i][idir].apply(spatial_gravity);
      }
      // nominal evaluation
      model.a[i] = model.X_lambda[i].apply(spatial_gravity);

    } else {

      // derivative evaluation
      for (unsigned idir = 0; idir < ndirs; idir++) {
        ad_model.v[i][idir]
          = ad_model.X_lambda[i][idir].apply(model.v[model.lambda[i]])
          + model.X_lambda[i].apply(ad_model.v[model.lambda[i]][idir])
          + ad_model.v_J[i][idir];
      }
      // nominal evaluation
      model.v[i] = model.X_lambda[i].apply(model.v[model.lambda[i]]) + model.v_J[i];

      // derivative evaluation
      for (unsigned idir = 0; idir < ndirs; idir++) {
        ad_model.c[i][idir]
          = ad_model.c_J[i][idir]
          + Math::ED::crossm (
            model.v[i], ad_model.v[i][idir],
            model.v_J[i], ad_model.v_J[i][idir]
          );
      }
      // nominal evaluation
      model.c[i] = model.c_J[i] + crossm(model.v[i],model.v_J[i]);

      // derivative evaluation
      for (unsigned idir = 0; idir < ndirs; idir++) {
        ad_model.a[i][idir]
          = ad_model.X_lambda[i][idir].apply(model.a[model.lambda[i]])
          + model.X_lambda[i].apply(ad_model.a[model.lambda[i]][idir])
          + ad_model.c[i][idir];
      }
      // nominal evaluation
      model.a[i] = model.X_lambda[i].apply(model.a[model.lambda[i]])
        + model.c[i];

    }

    if (!model.mBodies[i].mIsVirtual) {
      // derivative evaluation
      // NOTE save ndir times computation if Iv
      SpatialVector Iv = model.I[i] * model.v[i];
      for (unsigned idir = 0; idir < ndirs; idir++) {
        ad_model.f[i][idir]
          // NOTE inertia model.I[i] is constant
          = model.I[i] * ad_model.a[i][idir]
          + Math::ED::crossf (
            model.v[i], ad_model.v[i][idir],
            Iv, model.I[i] * ad_model.v[i][idir]
          );
      }
      // nominal evaluation
      model.f[i] = model.I[i] * model.a[i] + crossf(model.v[i], Iv);
    } else {
      for (unsigned idir = 0; idir < ndirs; idir++) {
        ad_model.f[i][idir].setZero();
      }
      // nominal evaluation
      model.f[i].setZero();
    }
  }

  for (unsigned int i = model.mBodies.size() - 1; i > 0; i--) {
    if(model.mJoints[i].mJointType != JointTypeCustom){
      if (model.mJoints[i].mDoFCount == 1) {
        // derivative evaluation
        for (unsigned idir = 0; idir < ndirs; idir++) {
          tau_dirs(model.mJoints[i].q_index, idir)
            = model.S[i].dot(ad_model.f[i][idir]);
            // NOTE ad_model.S[i][idir] = 0 for all joint types
            // NOTE not true for multidof joints
            // + ad_model.S[i][idir].dot(model.f[i]);
        }
        // nominal evaluation
        tau[model.mJoints[i].q_index]
          = model.S[i].dot(model.f[i]);
      } else if (model.mJoints[i].mDoFCount == 3) {
        cerr << __FILE__ << " " << __LINE__ << ":"
             << " Multi-DoF joint not supported." << endl;
        abort();
      }
    } else if(model.mJoints[i].mJointType == JointTypeCustom) {
      cerr << __FILE__ << " " << __LINE__ << ":"
           << " Custom joint not supported." << endl;
      abort();
    }

    if (model.lambda[i] != 0) {
      // derivative evaluation
      for (unsigned idir = 0; idir < ndirs; idir++) {
        ad_model.f[model.lambda[i]][idir] = ad_model.f[model.lambda[i]][idir]
          + ad_model.X_lambda[i][idir].applyTranspose(model.f[i])
          + model.X_lambda[i].applyTranspose(ad_model.f[i][idir]);
      }
      // nominal evaluation
      model.f[model.lambda[i]] = model.f[model.lambda[i]]
        + model.X_lambda[i].applyTranspose(model.f[i]);
    }
  }

  return;
}

RBDL_DLLAPI
void CompositeRigidBodyAlgorithm (
  Model &model,
  ADModel &ad_model,
  Math::VectorNd const & q,
  Math::MatrixNd const & q_dirs,
  Math::MatrixNd & H,
  std::vector<Math::MatrixNd> & H_dirs,
  bool update_kinematics
) {
  LOG << "-------- " << __func__ << " --------" << std::endl;

  assert (H.rows() == model.dof_count && H.cols() == model.dof_count);

  // receive number of directions
  const unsigned ndirs = q_dirs.cols();

  // resize container if necessary
  ad_model.resize_directions(ndirs);

  for (unsigned int i = 1; i < model.mBodies.size(); i++) {
    if (update_kinematics) {
      // derivative evaluation
      ED::jcalc_X_lambda_S (model, ad_model, i, q, q_dirs);
      // nominal evaluation
      // NOTE nominal evaluation is part of ED::jcalc_X_lambda_S
      // jcalc_X_lambda_S (model, i, q);
    }

    // derivative evaluation
    for (unsigned idir = 0; idir < ndirs; idir++) {
      // NOTE model.I[i] is constant
      ad_model.Ic[i][idir].setZero();
    }
    // nominal evaluation
    model.Ic[i] = model.I[i];
  }

  for (unsigned int i = model.mBodies.size() - 1; i > 0; i--) {
    if (model.lambda[i] != 0) {
      // derivative evaluation
      // NOTE AD code also takes care of nominal evaluation
      // NOTE applyTranspose computes X_lambda[j]^T * Ic * X_lambda[j]
      //      therefore special product rule is applied, respectively
      Math::ED::plusApplyTranspose(
        model.Ic[model.lambda[i]], ad_model.Ic[model.lambda[i]],
        model.X_lambda[i], ad_model.X_lambda[i],
        model.Ic[i], ad_model.Ic[i],
        ndirs
      );
      // nominal evaluation
      // NOTE applyTranspose computes X_lambda[j]^T * Ic * X_lambda[j]
      // model.Ic[model.lambda[i]] = model.Ic[model.lambda[i]]
      //   + model.X_lambda[i].applyTranspose(model.Ic[i]);
    }

    unsigned int dof_index_i = model.mJoints[i].q_index;

    if (model.mJoints[i].mDoFCount == 1
        && model.mJoints[i].mJointType != JointTypeCustom) {

      // derivative evaluation
      for (unsigned idir = 0; idir < ndirs; idir++) {
        // TODO is this necessary?
        ad_model.F[i][idir] = ad_model.Ic[i][idir] * model.S[i];
      }
      // nominal evaluation
      SpatialVector F = model.Ic[i] * model.S[i];

      // derivative evaluation
      for (unsigned idir = 0; idir < ndirs; idir++) {
        // NOTE model.S[i] is constant
        H_dirs[idir](dof_index_i, dof_index_i) = model.S[i].dot(ad_model.F[i][idir]);
      }
      // nominal evaluation
      H(dof_index_i, dof_index_i) = model.S[i].dot(F);

      unsigned int j = i;
      unsigned int dof_index_j = dof_index_i;

      // NOTE everything in the while loop has significant impact
      while (model.lambda[j] != 0) {
        // derivative evaluation
        for (unsigned idir = 0; idir < ndirs; idir++) {
          // NOTE applyTranspose requires special product rule
          // ad_model.F[i][idir]
          //   = ad_model.X_lambda[j][idir].applyTranspose(F)
          //   + model.X_lambda[j].applyTranspose(ad_model.F[i][idir]);
          Math::ED::applyTranspose(
            ad_model.F[i][idir],
            model.X_lambda[j], ad_model.X_lambda[j][idir],
            F,                 ad_model.F[i][idir]
          );
        }
        // nominal evaluation
        F = model.X_lambda[j].applyTranspose(F);

        j = model.lambda[j];
        dof_index_j = model.mJoints[j].q_index;

        if(model.mJoints[j].mJointType != JointTypeCustom) {
          if (model.mJoints[j].mDoFCount == 1) {

            // derivative evaluation
            for (unsigned idir = 0; idir < ndirs; idir++) {
              H_dirs[idir](dof_index_i,dof_index_j)
                = ad_model.F[i][idir].dot(model.S[j]);
              H_dirs[idir](dof_index_j,dof_index_i) = H_dirs[idir](dof_index_i,dof_index_j);
            }
            // nominal evaluation
            H(dof_index_i,dof_index_j) = F.dot(model.S[j]);
            H(dof_index_j,dof_index_i) = H(dof_index_i,dof_index_j);

          } else if (model.mJoints[j].mDoFCount == 3) {
            cerr << __FILE__ << " " << __LINE__ << ":"
                 << " Multi-DoF joint not supported." << endl;
            abort();
          }
        } else if (model.mJoints[j].mJointType == JointTypeCustom){
          cerr << __FILE__ << " " << __LINE__ << ":"
               << " Custom joint not supported." << endl;
          abort();
        }
      }
    } else if (model.mJoints[i].mDoFCount == 3
        && model.mJoints[i].mJointType != JointTypeCustom) {
      cerr << __FILE__ << " " << __LINE__ << ":"
           << " Multi-DoF joint not supported." << endl;
      abort();
    } else if (model.mJoints[i].mJointType == JointTypeCustom) {
      cerr << __FILE__ << " " << __LINE__ << ":"
           << " Custom joint not supported." << endl;
      abort();
    }
  }

  return;
}
*/

// -----------------------------------------------------------------------------
} // ED
// -----------------------------------------------------------------------------
} // RigidBodyDynamics
// -----------------------------------------------------------------------------
