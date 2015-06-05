#include <UnitTest++.h>

#include <iostream>

#include "rbdl/rbdl_mathutils.h"
#include "rbdl/Logging.h"

#include "rbdl/Model.h"
#include "rbdl/Kinematics.h"

using namespace std;
using namespace RigidBodyDynamics;
using namespace RigidBodyDynamics::Math;

const double TEST_PREC = 1.0e-12;

struct CartPendulum {
	CartPendulum () {
		using namespace RigidBodyDynamics;
		using namespace RigidBodyDynamics::Math;

		ClearLogOutput();
		model = new Model;

		double pendulum_length = 0.4;

		body_cart = Body (1., Vector3d (0., 0., 0.), Vector3d (0.1, 0.2, 0.5));
		body_pendulum = Body (0.4, Vector3d (0., 0., 0.5 * pendulum_length), Vector3d (0.1, 0.2, 0.5));

		joint_cart = Joint (SpatialVector (0., 0., 0., 1., 0., 0.));
		joint_pendulum = Joint (SpatialVector (0., 1., 0., 0., 0., 0.));

		id_cart = model->AddBody (0, Xtrans(Vector3d (0., 0., 0.)), joint_cart, body_cart, "cart");
		id_pendulum = model->AddBody (id_cart, Xtrans(Vector3d (0., 0., 0.)), joint_pendulum, body_pendulum, "pendulum");

		q = VectorNd::Constant ((size_t) model->dof_count, 0.);
		qdot = VectorNd::Constant ((size_t) model->dof_count, 0.);
		qddot = VectorNd::Constant ((size_t) model->dof_count, 0.);
		tau = VectorNd::Constant ((size_t) model->dof_count, 0.);

		ClearLogOutput();
	}
	~CartPendulum () {
		delete model;
	}
	
	RigidBodyDynamics::Model *model;

	unsigned int id_cart, id_pendulum;
	RigidBodyDynamics::Body body_cart, body_pendulum;
	RigidBodyDynamics::Joint joint_cart, joint_pendulum;

	RigidBodyDynamics::Math::VectorNd q;
	RigidBodyDynamics::Math::VectorNd qdot;
	RigidBodyDynamics::Math::VectorNd qddot;
	RigidBodyDynamics::Math::VectorNd tau;
};

TEST_FIXTURE ( CartPendulum, CartPendulumSimple ) {
//	CHECK_ARRAY_CLOSE (v_fixed_body.data(), v_body.data(), 6, TEST_PREC);
}
