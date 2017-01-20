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
void ForwardDynamicsADTestTemplate(T & obj, unsigned int numTrials,
                                   double CHECK_ARRAY_PREC = 1e-7) {
  Model ad_model = obj.model;
  Model fd_model = obj.model;
  ADModel ad_d_model = obj.ad_model;
  ADModel fd_d_model = obj.ad_model;
  srand(666);

  for(unsigned int trial = 0; trial < numTrials; trial++) {
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
          fd_d_model,
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
                      CHECK_ARRAY_PREC);

    CHECK_ARRAY_CLOSE(fd_dqddot.data(), ad_dqddot.data(),
                      fd_dqddot.cols() * fd_dqddot.rows(), CHECK_ARRAY_PREC);
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
void InverseDynamicsADTestTemplate(T & obj, double CHECK_ARRAY_PREC = 1e-7) {
  Model & model      = obj.model;
  ADModel & ad_model = obj.ad_model;
  VectorNd & q     = obj.q;
  VectorNd & qdot  = obj.qdot;
  VectorNd & qddot = obj.qddot;
  VectorNd & tau   = obj.tau;

  for(unsigned int i = 0; i < model.qdot_size; i++) {
    q[i] = (i+1)*(0.897878435);
    qdot[i] = (i+1)*(0.27563682);
    qddot[i] = (i+1)*(0.06565644564455);
  }

  MatrixNd q_dirs 	= MatrixNd::Identity (model.qdot_size, model.qdot_size);
  MatrixNd qdot_dirs 	= MatrixNd::Identity (model.qdot_size, model.qdot_size);
  MatrixNd qddot_dirs = MatrixNd::Identity (model.qdot_size, model.qdot_size);

  std::vector<SpatialVector> f_ext (
        model.mBodies.size(),
        SpatialVector::Zero());

  VectorNd tau_ref (tau);

  MatrixNd ad_tau  = MatrixNd::Zero(model.qdot_size, model.qdot_size);
  MatrixNd fd_tau  = MatrixNd::Zero(model.qdot_size, model.qdot_size);

  RigidBodyDynamics::AD::InverseDynamics(
        model, ad_model,
        q, q_dirs,
        qdot, qdot_dirs,
        qddot, qddot_dirs,
        tau, ad_tau,
        &f_ext
        );

  RigidBodyDynamics::FD::InverseDynamics(
        model,
        q, q_dirs,
        qdot, qdot_dirs,
        qddot, qddot_dirs,
        tau_ref, fd_tau,
        &f_ext
        );

  // NOTE: 1e-8 is not enough to make test pass
  CHECK_ARRAY_CLOSE (tau_ref.data(), tau.data(), tau_ref.rows(), CHECK_ARRAY_PREC);
  CHECK_ARRAY_CLOSE (fd_tau.data(), ad_tau.data(), fd_tau.cols()*fd_tau.rows(), CHECK_ARRAY_PREC);
}

TEST_FIXTURE(CartPendulum, CartPendulumInverseDynamicsADTest){
  InverseDynamicsADTestTemplate<CartPendulum>(*this);
}

//TEST_FIXTURE(Arm2DofX, Arm2DofXInverseDynamicsADTest) {
//  InverseDynamicsADTestTemplate(*this, 1e-6);
//}

//TEST_FIXTURE(Arm2DofZ, Arm2DofZInverseDynamicsADTest) {
//  InverseDynamicsADTestTemplate(*this, 1e-6);
//}

//TEST_FIXTURE(Arm3DofXZYp, Arm3DofXZYpInverseDynamicsADTest) {
//  InverseDynamicsADTestTemplate(*this, 1e-5);
//}

//TEST_FIXTURE(Arm3DofXZZp, Arm3DofXZYpInverseDynamicsADTest) {
//    InverseDynamicsADTestTemplate(*this, 1e-6);
//}

// -----------------------------------------------------------------------------

template<typename T>
void NonlinearEffectsADTestTemplate(T & obj, double CHECK_ARRAY_PREC = 1e-7) {
  Model & model      = obj.model;
  ADModel & ad_model = obj.ad_model;
  VectorNd & q     = obj.q;
  VectorNd & qdot  = obj.qdot;
  VectorNd & tau   = obj.tau;

  for(unsigned int i = 0; i < model.qdot_size; i++) {
      q[i]     = (i+1)*(0.897878435);
      qdot[i]  = (i+1)*(0.27563682);
  }

  MatrixNd q_dirs     = MatrixNd::Identity (model.qdot_size, model.qdot_size);
  MatrixNd qdot_dirs  = MatrixNd::Identity (model.qdot_size, model.qdot_size);

  VectorNd tau_nom (tau);
  VectorNd ad_tau_nom (tau);
  VectorNd fd_tau_nom (tau);
  MatrixNd ad_tau_der = MatrixNd::Zero(model.qdot_size, model.qdot_size);
  MatrixNd fd_tau_der = MatrixNd::Zero(model.qdot_size, model.qdot_size);

  NonlinearEffects(model, q, qdot, tau_nom);

  RigidBodyDynamics::AD::NonlinearEffects(model, ad_model, q, q_dirs,
      qdot, qdot_dirs, ad_tau_nom, ad_tau_der);

  RigidBodyDynamics::FD::NonlinearEffects(model, q, q_dirs, qdot, qdot_dirs,
      fd_tau_nom, fd_tau_der);

  CHECK_ARRAY_CLOSE (tau_nom.data(), ad_tau_nom.data(), tau_nom.rows(), CHECK_ARRAY_PREC);
  CHECK_ARRAY_CLOSE (tau_nom.data(), fd_tau_nom.data(), tau_nom.rows(), CHECK_ARRAY_PREC);
  CHECK_ARRAY_CLOSE (fd_tau_der.data(), ad_tau_der.data(),
                     fd_tau_der.cols() * fd_tau_der.rows(), CHECK_ARRAY_PREC);
}

//TEST_FIXTURE( CartPendulum, CartPendulumNonlinearEffectsADTest) {
//    NonlinearEffectsADTestTemplate(*this);
//}

//TEST_FIXTURE( Arm2DofX, Arm2DofXNonlinearEffectsADTest) {
//    NonlinearEffectsADTestTemplate(*this, 1e-6);
//}

//TEST_FIXTURE( Arm2DofZ, Arm2DofZNonlinearEffectsADTest) {
//    NonlinearEffectsADTestTemplate(*this, 1e-6);
//}

//TEST_FIXTURE( Arm3DofXZYp, Arm3DofXZYpNonlinearEffectsADTest) {
//    NonlinearEffectsADTestTemplate(*this, 1e-5);
//}

//TEST_FIXTURE( Arm3DofXZZp, Arm3DofXZZpNonlinearEffectsADTest) {
//    NonlinearEffectsADTestTemplate(*this, 1e-5);
//}

// -----------------------------------------------------------------------------

template<typename T>
void CompositeRigidBodyAlgorithmADTestTemplate(T & obj) {
    Model & model = obj.model;

    MatrixNd q_dirs = MatrixNd::Identity (model.qdot_size, model.qdot_size);

    MatrixNd H = MatrixNd::Zero(model.q_size,model.q_size);
    std::vector<MatrixNd> fd_out;
    std::vector<MatrixNd> ad_out(q_dirs.cols(),MatrixNd::Zero(model.q_size,model.q_size));
    ADModel ad_model(model);

    VectorNd q = VectorNd::Random(model.qdot_size);
    RigidBodyDynamics::FD::CompositeRigidBodyAlgorithm(model, q, q_dirs, fd_out);
    RigidBodyDynamics::AD::CompositeRigidBodyAlgorithm(model, ad_model, q, q_dirs, H, ad_out);

    MatrixNd inertia_test = MatrixNd::Zero(q.size(),q.size());
    CompositeRigidBodyAlgorithm(model, q, inertia_test, true);
    for (size_t nIdx = 0; nIdx < fd_out.size(); nIdx++) {
        // cout << "Dir " << nIdx << endl;
        // cout << "fd_CRBA: " << endl << fd_out[nIdx] << std::endl << std::endl;
        // cout << "ad_CRBA: " << endl << ad_out[nIdx] << std::endl << std::endl;
        // cout << "ad_CRBA error (fd,ad)" << endl << (fd_out[nIdx] - ad_out[nIdx]) << std::endl;
        CHECK_ARRAY_CLOSE (inertia_test.data(), H.data(), model.q_size*model.q_size, TEST_PREC);
        CHECK_ARRAY_CLOSE (fd_out[nIdx].data(), ad_out[nIdx].data(), model.q_size*model.q_size, 1e-7);
    }
}

//TEST_FIXTURE( CartPendulum, CartPendulumCompositeRigidBodyAlgorithmADTest) {
//    CompositeRigidBodyAlgorithmADTestTemplate(*this);
//}

//TEST_FIXTURE( Arm2DofX, Arm2DofXCompositeRigidBodyAlgorithmADTest) {
//    CompositeRigidBodyAlgorithmADTestTemplate(*this);
//}

//TEST_FIXTURE( Arm2DofZ, Arm2DofZCompositeRigidBodyAlgorithmADTest) {
//    CompositeRigidBodyAlgorithmADTestTemplate(*this);
//}

//TEST_FIXTURE( Arm3DofXZYp, Arm3DofXZYpCompositeRigidBodyAlgorithmADTest) {
//    CompositeRigidBodyAlgorithmADTestTemplate(*this);
//}

//TEST_FIXTURE( Arm3DofXZZp, Arm3DofXZZpCompositeRigidBodyAlgorithmADTest) {
//    CompositeRigidBodyAlgorithmADTestTemplate(*this);
//}

// -----------------------------------------------------------------------------
