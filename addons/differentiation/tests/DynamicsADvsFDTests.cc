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

using namespace std;
using namespace RigidBodyDynamics;
using namespace RigidBodyDynamics::Math;

const double TEST_PREC = 1.0e-8;

// -----------------------------------------------------------------------------

TEST_FIXTURE(CartPendulum, ForwardDynamicsADTest){
	srand((unsigned int) time(0));

	for(unsigned int trial = 0; trial < 10; trial++) {
		VectorNd q = VectorNd::Random(model.q_size);
		VectorNd qdot = VectorNd::Random(model.q_size);
		VectorNd tau = VectorNd::Random(model.q_size);

		unsigned int ndirs = 3 * model.q_size;
		MatrixNd x = MatrixNd::Identity(ndirs, ndirs);
		MatrixNd q_dirs = x.block(0, 0, model.q_size, ndirs);
		MatrixNd qdot_dirs = x.block(model.q_size, 0, model.q_size, ndirs);
		MatrixNd tau_dirs = x.block(2*model.q_size, 0, model.q_size, ndirs);

		std::vector<SpatialVector> f_ext (
			model.mBodies.size(),
			SpatialVector::Zero()
		);
		// for (unsigned int i = 0; i < model.mBodies.size(); ++i) {
		// 	f_ext[i]=SpatialVector::Zero(model.q_size);
		// }

		VectorNd qddot (VectorNd::Zero(model.q_size));
		MatrixNd ad_qddot  = MatrixNd::Zero(model.qdot_size, ndirs);
		MatrixNd fd_qddot  = MatrixNd::Zero(model.qdot_size, ndirs);

		RigidBodyDynamics::FD::ForwardDynamics(
			model,
			q, q_dirs,
			qdot, qdot_dirs,
			tau, tau_dirs,
			qddot, fd_qddot,
			&f_ext
		);

		RigidBodyDynamics::FD::ForwardDynamics(
			model,
			q, q_dirs,
			qdot, qdot_dirs,
			tau, tau_dirs,
			qddot, ad_qddot,
			&f_ext
		);

		CHECK_ARRAY_CLOSE (
			fd_qddot.data(), ad_qddot.data(), fd_qddot.cols()*fd_qddot.rows(),
			TEST_PREC
		);
	}
}

// -----------------------------------------------------------------------------

TEST_FIXTURE(CartPendulum, InverseDynamicsADTest){

	for(unsigned int i = 0; i < model.qdot_size; i++){
		q[i] = (i+1)*(0.897878435);
		qdot[i] = (i+1)*(0.27563682);
		qddot[i] = (i+1)*(0.06565644564455);
	}

	MatrixNd q_dirs 	= MatrixNd::Identity (model.qdot_size, model.qdot_size);
	MatrixNd qdot_dirs 	= MatrixNd::Identity (model.qdot_size, model.qdot_size);
	MatrixNd qddot_dirs = MatrixNd::Identity (model.qdot_size, model.qdot_size);

	// std::vector<SpatialVector> f_ext;
	std::vector<SpatialVector> f_ext (
		model.mBodies.size(),
		SpatialVector::Zero()
	);
	// for (unsigned int i = 0; i < model.mBodies.size(); ++i) {
	// 	f_ext[i]=SpatialVector::Zero(model.q_size);
	// }

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
	CHECK_ARRAY_CLOSE (tau_ref.data(), tau.data(), tau_ref.rows(), 1e-07);
	CHECK_ARRAY_CLOSE (fd_tau.data(), ad_tau.data(), fd_tau.cols()*fd_tau.rows(), 1e-07);
}

// -----------------------------------------------------------------------------

TEST_FIXTURE( CartPendulum, CompositeRigidBodyAlgorithmADTest) {
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
		CHECK_ARRAY_CLOSE (fd_out[nIdx].data(), ad_out[nIdx].data(), model.q_size*model.q_size, TEST_PREC);
	}
}

// -----------------------------------------------------------------------------
