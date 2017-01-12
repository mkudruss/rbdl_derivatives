#include "rbdl_utilsAD.h"
#include "rbdl_mathutilsAD.h"
#include "KinematicsAD.h"

using std::fill_n;
using std::vector;

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------
namespace Utils {
// -----------------------------------------------------------------------------
namespace AD {
// -----------------------------------------------------------------------------

using namespace Math;

RBDL_DLLAPI void CalcCenterOfMass (
    Model & model,
    ADModel & ad_model,
    VectorNd const & q,
    MatrixNd const & q_dirs,
    VectorNd const & qdot,
    MatrixNd const & qdot_dirs,
    double & mass,
    Vector3d & com,
    MatrixNd & ad_com,
    Vector3d * com_velocity,
    MatrixNd * ad_com_velocity,
    Vector3d * angular_momentum,
    MatrixNd * ad_angular_momentum,
    bool update_kinematics) {
  int ndirs = q_dirs.cols();
  assert(ndirs == qdot_dirs.cols());
  assert(ndirs == ad_com.cols());

  if (com_velocity) {
    assert(q_dirs.cols() == ad_com_velocity->cols());
  }
  if (angular_momentum) {
    assert(q_dirs.cols() == ad_angular_momentum->cols());
  }

  if (update_kinematics) {
    RigidBodyDynamics::AD::UpdateKinematicsCustom (model, ad_model, &q,
                                                   &q_dirs, &qdot, &qdot_dirs, NULL, NULL);
  }

  for (size_t i = 1; i < model.mBodies.size(); i++) {
    // derivative evaluation
    for(int idir = 0; idir < ndirs; idir++) {
      ad_model.Ic[i][idir].setZero();
    }
    // nominal evaluation
    model.Ic[i] = model.I[i];

    // derivative evaluation
    for(int idir = 0; idir < ndirs; idir++) {
      ad_model.hc[i][idir] = model.Ic[i].toMatrix() * ad_model.v[i][idir];
      // + ad_model.Ic[i][idir] * model.v[i];
      // summand is skipped because it is always zero
    }
    // nominal evaluation
    model.hc[i] = model.Ic[i].toMatrix() * model.v[i];
  }

  SpatialRigidBodyInertia Itot (0., Vector3dZero, Matrix3dZero);
  /// TODO: Consider replacing ad_Itot by vector<SpatialMatrix>
  vector<SpatialRigidBodyInertia> ad_Itot(ndirs, Itot);

  SpatialVector htot (SpatialVector::Zero(6));
  vector<SpatialVector> ad_htot(ndirs, htot);

  for (size_t i = model.mBodies.size() - 1; i > 0; i--) {
    unsigned int lambda = model.lambda[i];
    if (lambda != 0) {
      // derivative evaluation
      for (int idir = 0; idir < ndirs; idir++) {
        ad_model.Ic[lambda][idir] = ad_model.Ic[lambda][idir]
            + model.X_lambda[i].toMatrixTranspose() * ad_model.Ic[i][idir] * model.X_lambda[i].toMatrix()
            + model.X_lambda[i].toMatrixTranspose() * model.Ic[i].toMatrix() * ad_model.X_lambda[i][idir]
            + ad_model.X_lambda[i][idir].transpose() * model.Ic[i].toMatrix() * model.X_lambda[i].toMatrix();
      }
      // nominal evaluation
      model.Ic[lambda] = model.Ic[lambda] + model.X_lambda[i].applyTranspose (model.Ic[i]);

      // derivative evaluation
      for (int idir = 0; idir < ndirs; idir++) {
        ad_model.hc[lambda][idir] = ad_model.hc[lambda][idir]
            + model.X_lambda[i].toMatrixTranspose() * ad_model.hc[i][idir]
            + ad_model.X_lambda[i][idir].transpose() * model.hc[i];
      }
      // nominal evaluation
      model.hc[lambda] = model.hc[lambda] + model.X_lambda[i].applyTranspose (model.hc[i]);
    } else {
      // derivative evaluation
      for (int idir = 0; idir < ndirs; idir++) {
        SpatialMatrix m =
            model.X_lambda[i].toMatrixTranspose() * ad_model.Ic[i][idir] * model.X_lambda[i].toMatrix()
            + model.X_lambda[i].toMatrixTranspose() * model.Ic[i].toMatrix() * ad_model.X_lambda[i][idir]
            + ad_model.X_lambda[i][idir].transpose() * model.Ic[i].toMatrix() * model.X_lambda[i].toMatrix();

        SpatialRigidBodyInertia m2i;
        m2i.createFromMatrix(m);
        ad_Itot[idir] = ad_Itot[idir] + m2i;
      }
      // nominal evaluation
      Itot = Itot + model.X_lambda[i].applyTranspose (model.Ic[i]);

      // derivative evaluation
      for (int idir = 0; idir < ndirs; idir++) {
        ad_htot[idir] +=
            model.X_lambda[i].toMatrixTranspose() * ad_model.hc[i][idir]
            + ad_model.X_lambda[i][idir].transpose() * model.hc[i];
      }
      // nominal evaluation
      htot = htot + model.X_lambda[i].applyTranspose (model.hc[i]);
    }
  }

  std::cout << "AD HTOT 1 = " << std::endl;
  for (int i = 0; i < ndirs; i++) {
    std::cout << ad_htot[i].transpose() << std::endl;
  }

  mass = Itot.m;
  com = Itot.h / mass;

  for (int idir = 0; idir < ndirs; idir++) {
    ad_com.block<3,1>(0, idir) = ad_Itot[idir].h / mass;
  }

  LOG << "mass = " << mass << " com = " << com.transpose() << " htot = " << htot.transpose() << std::endl;

  if (com_velocity) {
    // derivative evaluation
    for (int idir = 0; idir < ndirs; idir++) {
      ad_com_velocity->block<3, 1>(0, idir) = Vector3d(
            ad_htot[idir][3] / mass,
          ad_htot[idir][4] / mass,
          ad_htot[idir][5] / mass);
    }
    // nominal evaluation
    *com_velocity = Vector3d (
          htot[3] / mass,
        htot[4] / mass,
        htot[5] / mass);
  }

  if (angular_momentum) {
    // derivative evaluation
    for (int idir = 0; idir < ndirs; idir++) {
      ad_htot[idir] = Xtrans(com).toMatrixAdjoint() * ad_htot[idir]
          + Math::AD::Xtrans(com, ad_com.block<3,1>(0, idir)).adjoint() * htot;
    }
    // nominal evaluation
    htot = Xtrans (com).applyAdjoint (htot);

    // derivative evaluation
    for (int idir = 0; idir < ndirs; idir++) {
      ad_angular_momentum->block<3,1>(0, idir) = ad_htot[idir].block<3, 1>(0, 0);
      // Vector3d(ad_htot[idir][0], ad_htot[idir][1], ad_htot[idir][2]);
    }
    // nominal evaluation
    angular_momentum->set (htot[0], htot[1], htot[2]);
  }

//  std::cout << "AD HTOT 2 = " << std::endl;
//  for (int i = 0; i < ndirs; i++) {
//    std::cout << ad_htot[i].transpose() << std::endl;
//  }
}

RBDL_DLLAPI double CalcPotentialEnergy (
    Model & model,
    ADModel & ad_model,
    VectorNd const & q,
    MatrixNd const & q_dirs,
    MatrixNd & ad_pote,
    bool update_kinematics)
{
  int ndirs = q_dirs.cols();
  assert(ad_pote.cols() == ndirs);
  assert(ad_pote.rows() == 1);

  double mass = 0;
  Vector3d com;
  MatrixNd ad_com(3, ndirs);
  CalcCenterOfMass (model, ad_model, q, q_dirs,
                    VectorNd::Zero (model.qdot_size), /// TODO: Cache in model? (model.zeroVecNq
                    MatrixNd::Zero (model.qdot_size, q_dirs.cols()), /// TODO: Cache in model? (if q_dirs.size() == nq : use model.zeroMatNqNq)
                    mass, com, ad_com, NULL, NULL, NULL, NULL, update_kinematics);

  Vector3d g = -Vector3d(model.gravity[0], model.gravity[1], model.gravity[2]);

  LOG << "pot_energy: " << " mass = " << mass << " com = " << com.transpose() << std::endl;

  // derivative value
  ad_pote = mass * g.transpose() * ad_com;
  // nominal value
  return mass * com.dot(g);
}

RBDL_DLLAPI double CalcKineticEnergy (
    Model & model,
    ADModel & ad_model,
    VectorNd const & q,
    MatrixNd const & q_dirs,
    VectorNd const & qdot,
    MatrixNd const & qdot_dirs,
    MatrixNd & ad_kine,
    bool update_kinematics) {
  int ndirs = q_dirs.cols();
  assert(ndirs == qdot_dirs.cols());
  assert(ndirs == ad_kine.cols());
  assert(1 == ad_kine.rows());

  if (update_kinematics) {
    RigidBodyDynamics::AD::UpdateKinematicsCustom(model, ad_model,
                                                  &q, &q_dirs, &qdot, &qdot_dirs, 0, 0);
  }

  ad_kine.setZero();
  double kine = 0.;
  for (size_t i = 1; i < model.mBodies.size(); i++) {
    SpatialVector Iv = model.I[i] * model.v[i];
    // derivative value
    for (int idir = 0; idir < ndirs; idir++) {
      ad_kine.block<1,1>(0, idir) += .5 * (
            ad_model.v[i][idir].transpose() * Iv
            + Iv.transpose() * ad_model.v[i][idir]);
    }
    // nominal value
    kine += 0.5 * model.v[i].transpose() * Iv;
  }
  return kine;
}

// -----------------------------------------------------------------------------
} // namespace RigidBodyDynamics
// -----------------------------------------------------------------------------
} // namespace Utils
// -----------------------------------------------------------------------------
} // namespace AD
// -----------------------------------------------------------------------------


