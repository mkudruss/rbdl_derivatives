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
