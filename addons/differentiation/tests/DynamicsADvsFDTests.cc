#include <UnitTest++.h>

#include <iostream>

#include "rbdl/Model.h"
#include "rbdl/Dynamics.h"
#include "rbdl/rbdl_mathutils.h"

#include "Fixtures.h"
#include "ModelAD.h"
#include "DynamicsAD.h"
#include "DynamicsFD.h"
#include "rbdl_mathutilsAD.h"

#include "ModelCheckADvsFD.h"

using namespace std;
using namespace RigidBodyDynamics;
using namespace RigidBodyDynamics::Math;

const double TEST_PREC = 1.0e-8;

// -----------------------------------------------------------------------------

template <typename T>
void ForwardDynamicsADTestTemplate(
    T & obj,
    unsigned int trial_count,
    double array_close_prec) {
  Model ad_model = obj.model;
  Model fd_model = obj.model;
  ADModel ad_d_model = obj.ad_model;
  ADModel fd_d_model = obj.ad_model;
  srand(666);

  for(unsigned int trial = 0; trial < trial_count; trial++) {
    VectorNd q = VectorNd::Random(obj.model.q_size);
    VectorNd qdot = VectorNd::Random(obj.model.q_size);
    VectorNd tau = VectorNd::Random(obj.model.q_size);

    unsigned int ndirs = 3 * obj.model.q_size;
    MatrixNd x = MatrixNd::Identity(ndirs, ndirs);
    MatrixNd q_dirs = x.block(0, 0, obj.model.q_size, ndirs);
    MatrixNd qdot_dirs = x.block(obj.model.q_size, 0, obj.model.q_size, ndirs);
    MatrixNd tau_dirs = x.block(2 * obj.model.q_size, 0, obj.model.q_size, ndirs);

    vector<SpatialVector> f_ext (
          obj.model.mBodies.size(),
          SpatialVector::Zero()
          );

    vector<vector<SpatialVector>> f_ext_dirs (
          obj.model.mBodies.size(),
          vector<SpatialVector>(
            ndirs,
            SpatialVector::Zero()));

    if (trial != 0) {
      for (unsigned i = 0; i < f_ext.size(); i++) {
        f_ext[i].setRandom();
      }
      for (unsigned i = 0; i < f_ext.size(); i++) {
        for (unsigned idir = 0; idir < f_ext_dirs[i].size(); idir++) {
          f_ext_dirs[i][idir].setRandom();
        }
      }
    }

    VectorNd ad_qddot (VectorNd::Zero(obj.model.q_size));
    VectorNd fd_qddot (VectorNd::Zero(obj.model.q_size));
    MatrixNd ad_dqddot  = MatrixNd::Zero(obj.model.qdot_size, ndirs);
    MatrixNd fd_dqddot  = MatrixNd::Zero(obj.model.qdot_size, ndirs);

    AD::ForwardDynamics(
          ad_model,
          ad_d_model,
          q, q_dirs,
          qdot, qdot_dirs,
          tau, tau_dirs,
          ad_qddot, ad_dqddot,
          &f_ext, &f_ext_dirs);

    FD::ForwardDynamics(
          fd_model,
          &fd_d_model,
          q, q_dirs,
          qdot, qdot_dirs,
          tau, tau_dirs,
          fd_qddot, fd_dqddot,
          &f_ext, &f_ext_dirs);

    checkModelsADvsFD(
          ndirs,
          ad_model, ad_d_model,
          fd_model, fd_d_model);

    CHECK_ARRAY_CLOSE(ad_qddot.data(), fd_qddot.data(), obj.model.q_size,
                      array_close_prec);

    CHECK_ARRAY_CLOSE(fd_dqddot.data(), ad_dqddot.data(),
                      fd_dqddot.cols() * fd_dqddot.rows(), array_close_prec);
  }
}

TEST_FIXTURE(CartPendulum, CartPendulumForwardDynamicsADTest){
  ForwardDynamicsADTestTemplate(*this, 10, 1e-5);
}

TEST_FIXTURE(Arm2DofX, Arm2DofXForwardDynamicsADTest){
  ForwardDynamicsADTestTemplate(*this, 10, 1e-5);
}

TEST_FIXTURE(Arm2DofZ, Arm2DofZForwardDynamicsADTest){
  ForwardDynamicsADTestTemplate(*this, 10, 1e-5);
}

TEST_FIXTURE(Arm3DofXZYp, Arm3DofXZYpForwardDynamicsADTest){
  ForwardDynamicsADTestTemplate(*this, 10, 1e-5);
}

TEST_FIXTURE(Arm3DofXZZp, Arm3DofXZZpForwardDynamicsADTest){
  ForwardDynamicsADTestTemplate(*this, 10, 1e-5);
}

// -----------------------------------------------------------------------------

template<typename T>
void InverseDynamicsADTestTemplate(
    T & obj,
    unsigned trial_count,
    double array_close_prec) {
  Model ad_model      = obj.model;
  Model fd_model      = obj.model;
  ADModel ad_d_model = obj.ad_model;
  ADModel fd_d_model = obj.ad_model;
  VectorNd & q     = obj.q;
  VectorNd & qdot  = obj.qdot;
  VectorNd & qddot = obj.qddot;
  VectorNd & tau   = obj.tau;
  unsigned ndirs = ad_model.qdot_size;

  for (unsigned trial = 0; trial < trial_count; trial++) {
    q.setRandom();
    qdot.setRandom();
    qddot.setRandom();

    MatrixNd q_dirs 	= MatrixNd::Random (ad_model.qdot_size, ad_model.qdot_size);
    MatrixNd qdot_dirs 	= MatrixNd::Random (ad_model.qdot_size, ad_model.qdot_size);
    MatrixNd qddot_dirs = MatrixNd::Random (ad_model.qdot_size, ad_model.qdot_size);

    std::vector<SpatialVector> f_ext (
          ad_model.mBodies.size(),
          SpatialVector::Zero());

    vector<vector<SpatialVector>> f_ext_dirs (
          obj.model.mBodies.size(),
          vector<SpatialVector>(
            ndirs,
            SpatialVector::Zero()));

    for (unsigned i = 0; i < f_ext.size(); i++) {
      f_ext[i].setRandom();
    }
    for (unsigned i = 0; i < f_ext.size(); i++) {
      for (unsigned idir = 0; idir < f_ext_dirs[i].size(); idir++) {
        f_ext_dirs[i][idir].setRandom();
      }
    }

    VectorNd tau_ref (tau);

    MatrixNd ad_tau  = MatrixNd::Random(ad_model.qdot_size, ad_model.qdot_size);
    MatrixNd fd_tau  = MatrixNd::Random(ad_model.qdot_size, ad_model.qdot_size);

    AD::InverseDynamics(
          ad_model, ad_d_model,
          q, q_dirs,
          qdot, qdot_dirs,
          qddot, qddot_dirs,
          tau, ad_tau,
          &f_ext, &f_ext_dirs);

    FD::InverseDynamics(
          fd_model, &fd_d_model,
          q, q_dirs,
          qdot, qdot_dirs,
          qddot, qddot_dirs,
          tau_ref, fd_tau,
          &f_ext, &f_ext_dirs);

    checkModelsADvsFD(ndirs, ad_model, ad_d_model, fd_model, fd_d_model);

    CHECK_ARRAY_CLOSE (tau_ref.data(), tau.data(), tau_ref.rows(), array_close_prec);
    CHECK_ARRAY_CLOSE (fd_tau.data(), ad_tau.data(), fd_tau.cols()*fd_tau.rows(), array_close_prec);
  }
}

TEST_FIXTURE(CartPendulum, CartPendulumInverseDynamicsADTest){
  InverseDynamicsADTestTemplate<CartPendulum>(*this, 10, 1e-5);
}

TEST_FIXTURE(Arm2DofX, Arm2DofXInverseDynamicsADTest) {
  InverseDynamicsADTestTemplate(*this, 10, 1e-5);
}

TEST_FIXTURE(Arm2DofZ, Arm2DofZInverseDynamicsADTest) {
  InverseDynamicsADTestTemplate(*this, 10, 1e-5);
}

TEST_FIXTURE(Arm3DofXZYp, Arm3DofXZYpInverseDynamicsADTest) {
  InverseDynamicsADTestTemplate(*this, 10, 1e-5);
}

TEST_FIXTURE(Arm3DofXZZp, Arm3DofXZYpInverseDynamicsADTest) {
  InverseDynamicsADTestTemplate(*this, 10, 1e-5);
}

// -----------------------------------------------------------------------------

template<typename T>
void NonlinearEffectsADTestTemplate(
    T & obj,
    unsigned trial_count,
    double array_close_prec) {
  Model ad_model      = obj.model;
  Model fd_model      = obj.model;
  ADModel ad_d_model = obj.ad_model;
  ADModel fd_d_model = obj.ad_model;
  VectorNd q(obj.model.dof_count);
  VectorNd qdot(obj.model.dof_count);
  VectorNd tau(obj.model.dof_count);

  unsigned ndirs = ad_model.qdot_size;

  for (unsigned trial = 0; trial < trial_count; trial++) {
    q.setRandom();
    qdot.setRandom();

    MatrixNd q_dirs     = MatrixNd::Random (ad_model.qdot_size, ad_model.qdot_size);
    MatrixNd qdot_dirs  = MatrixNd::Random (ad_model.qdot_size, ad_model.qdot_size);

    VectorNd tau_nom (tau);
    VectorNd ad_d_tau_nom (tau);
    VectorNd fd_d_tau_nom (tau);
    MatrixNd ad_tau_der = MatrixNd::Zero(ad_model.qdot_size, ad_model.qdot_size);
    MatrixNd fd_tau_der = MatrixNd::Zero(ad_model.qdot_size, ad_model.qdot_size);

    NonlinearEffects(ad_model, q, qdot, tau_nom);

    AD::NonlinearEffects(ad_model, ad_d_model, q, q_dirs, qdot, qdot_dirs,
                         ad_d_tau_nom, ad_tau_der);

    FD::NonlinearEffects(fd_model, &fd_d_model, q, q_dirs, qdot, qdot_dirs,
                         fd_d_tau_nom, fd_tau_der);

    checkModelsADvsFD(ndirs, ad_model, ad_d_model, fd_model, fd_d_model);

    CHECK_ARRAY_CLOSE (tau_nom.data(), ad_d_tau_nom.data(), tau_nom.rows(), array_close_prec);
    CHECK_ARRAY_CLOSE (tau_nom.data(), fd_d_tau_nom.data(), tau_nom.rows(), array_close_prec);

    CHECK_ARRAY_CLOSE (fd_tau_der.data(), ad_tau_der.data(),
                       fd_tau_der.cols() * fd_tau_der.rows(), array_close_prec);
  }
}

TEST_FIXTURE( CartPendulum, CartPendulumNonlinearEffectsADTest) {
  NonlinearEffectsADTestTemplate(*this, 10, 1e-6);
}

TEST_FIXTURE( Arm2DofX, Arm2DofXNonlinearEffectsADTest) {
  NonlinearEffectsADTestTemplate(*this, 10, 1e-6);
}

TEST_FIXTURE( Arm2DofZ, Arm2DofZNonlinearEffectsADTest) {
  NonlinearEffectsADTestTemplate(*this, 10, 1e-6);
}

TEST_FIXTURE( Arm3DofXZYp, Arm3DofXZYpNonlinearEffectsADTest) {
  NonlinearEffectsADTestTemplate(*this, 10, 1e-5);
}

TEST_FIXTURE( Arm3DofXZZp, Arm3DofXZZpNonlinearEffectsADTest) {
  NonlinearEffectsADTestTemplate(*this, 10, 1e-5);
}

TEST_FIXTURE (FixedBase6DoF, FixedBase6DoFNonlinearEffectsADTest) {
  // add contacts and bind them to constraint set
  constraint_set.AddContactConstraint (contact_body_id, Vector3d (1., 0., 0.), contact_normal);
  constraint_set.AddContactConstraint (contact_body_id, Vector3d (0., 1., 0.), contact_normal);
  constraint_set.Bind (model);
  ad_constraint_set = ADConstraintSet(constraint_set, model.dof_count);
  NonlinearEffectsADTestTemplate(*this, 10, 1e-5);
}

// -----------------------------------------------------------------------------

template<typename T>
void CompositeRigidBodyAlgorithmADTestTemplate(
    T & obj,
    unsigned trial_count,
    double array_close_prec) {
  Model model = obj.model;
  Model fd_model = obj.model;
  Model ad_model = obj.model;
  ADModel ad_d_model(model);
  ADModel fd_d_model = ad_d_model;
  for (unsigned trial = 0; trial < trial_count; trial++) {
    unsigned ndirs = model.qdot_size;
    VectorNd q = VectorNd::Random(model.qdot_size);
    MatrixNd q_dirs = MatrixNd::Random (model.qdot_size, ndirs);
    MatrixNd H = MatrixNd::Zero(q.size(),q.size());
    MatrixNd ad_H = MatrixNd::Random(model.q_size, ndirs);
    MatrixNd fd_H = MatrixNd::Random(model.q_size, ndirs);
    vector<MatrixNd> fd_H_dirs(ndirs);
    vector<MatrixNd> ad_H_dirs(ndirs, MatrixNd::Zero(model.q_size,model.q_size));

    CompositeRigidBodyAlgorithm(model, q, H, true);
    FD::CompositeRigidBodyAlgorithm(fd_model, &fd_d_model, q, q_dirs, fd_H, fd_H_dirs);
    AD::CompositeRigidBodyAlgorithm(ad_model, ad_d_model, q, q_dirs, ad_H, ad_H_dirs);

    checkModelsADvsFD(ndirs,
                      ad_model, ad_d_model,
                      fd_model, fd_d_model);

    CHECK_ARRAY_CLOSE (H.data(), ad_H.data(),
                       model.q_size * model.q_size,
                       array_close_prec);
    CHECK_ARRAY_CLOSE (H.data(), fd_H.data(),
                       model.q_size * model.q_size,
                       array_close_prec);
    for (size_t idir = 0; idir < model.qdot_size; idir++) {
      CHECK_ARRAY_CLOSE (fd_H_dirs[idir].data(), ad_H_dirs[idir].data(),
                         model.q_size * model.q_size,
                         array_close_prec);
    }
  }
}

TEST_FIXTURE( CartPendulum, CartPendulumCompositeRigidBodyAlgorithmADTest) {
  CompositeRigidBodyAlgorithmADTestTemplate(*this, 10, 1e-5);
}

TEST_FIXTURE( Arm2DofX, Arm2DofXCompositeRigidBodyAlgorithmADTest) {
  CompositeRigidBodyAlgorithmADTestTemplate(*this, 10, 1e-5);
}

TEST_FIXTURE( Arm2DofZ, Arm2DofZCompositeRigidBodyAlgorithmADTest) {
  CompositeRigidBodyAlgorithmADTestTemplate(*this, 10, 1e-5);
}

TEST_FIXTURE( Arm3DofXZYp, Arm3DofXZYpCompositeRigidBodyAlgorithmADTest) {
  CompositeRigidBodyAlgorithmADTestTemplate(*this, 10, 1e-5);
}

TEST_FIXTURE( Arm3DofXZZp, Arm3DofXZZpCompositeRigidBodyAlgorithmADTest) {
  CompositeRigidBodyAlgorithmADTestTemplate(*this, 10, 1e-5);
}

// -----------------------------------------------------------------------------
