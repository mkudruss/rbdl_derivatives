#include "rbdl/rbdl.h"

#include "ModelAD.h"

struct CartPendulum {
    CartPendulum () {
        using namespace RigidBodyDynamics;
        using namespace RigidBodyDynamics::Math;

        ClearLogOutput();

        double cart_w = 0.5;
        double cart_d = 0.2;
        double cart_h = 0.2;
        double cart_m = 100.0;
        Vector3d cart_com (0., 0., 0.);

        double pend_l = 0.5;
        double pend_r = 0.1;
        double pend_m = 1.0;
        Vector3d pend_com (0., 0., pend_l);

        body_cart = Body (cart_m, cart_com, Matrix3d (
                    1. / 12. * cart_m * (cart_h*cart_h + cart_d * cart_d), 0., 0.,
                    0., 1. / 12. * cart_m * (cart_w*cart_w + cart_h * cart_h), 0.,
                    0., 0., 1. / 12. * cart_m * (cart_w*cart_w + cart_d * cart_d)
                    ));
        body_pendulum = Body (pend_m, pend_com, Matrix3d(
                    2. / 5. * pend_m * pend_r, 0., 0.,
                    0., 2. / 5. * pend_m * pend_r, 0.,
                    0., 0., 2. / 5. * pend_m * pend_r
                ));

        joint_cart = Joint (SpatialVector (0., 0., 0., 1., 0., 0.));
        joint_pendulum = Joint (SpatialVector (0., 1., 0., 0., 0., 0.));

        id_cart = model.AddBody (0, Xtrans(Vector3d (0., 0., 0.)), joint_cart, body_cart, "cart");
        id_pendulum = model.AddBody (id_cart, Xtrans(Vector3d (0., 0., 0.)), joint_pendulum, body_pendulum, "pendulum");

        q = VectorNd::Constant ((size_t) model.dof_count, 0.);
        qdot = VectorNd::Constant ((size_t) model.dof_count, 0.);
        qddot = VectorNd::Constant ((size_t) model.dof_count, 0.);
        tau = VectorNd::Constant ((size_t) model.dof_count, 0.);

        ad_model = ADModel(model);

		body_point = RigidBodyDynamics::Math::Vector3d (0., 0., pend_l);

        ClearLogOutput();
    }

    RigidBodyDynamics::Model model;
    ADModel ad_model;

    unsigned int id_cart, id_pendulum;
    RigidBodyDynamics::Body body_cart, body_pendulum;
    RigidBodyDynamics::Joint joint_cart, joint_pendulum;

    RigidBodyDynamics::Math::VectorNd q;
    RigidBodyDynamics::Math::VectorNd qdot;
    RigidBodyDynamics::Math::VectorNd qddot;
    RigidBodyDynamics::Math::VectorNd tau;

	RigidBodyDynamics::Math::Vector3d body_point;
};

struct Arm2DofX {
	Arm2DofX() {
		using namespace RigidBodyDynamics;
		using namespace RigidBodyDynamics::Math;

		ClearLogOutput();

		double arm_l = 0.5;
		double arm_r = 0.1;
		double arm_m = 3.5;
		Matrix3d arm_I = .5 * arm_m * Matrix3d(
					.5 * (arm_r * arm_r + arm_l * arm_l / 3.), 0., 0.,
					arm_r * arm_r, 0., 0.,
					.5 * (arm_r * arm_r + arm_l * arm_l / 3.), 0., 0.);
		Vector3d arm_com (0., .5 * arm_l, 0.);

		body_distal    = Body(arm_m, arm_com, arm_I);
		body_proximal  = Body(arm_m, arm_com, arm_I);

		joint_distal   = Joint(SpatialVector(1., 0., 0., 0., 0., 0.));
		joint_proximal = Joint(SpatialVector(1., 0., 0., 0., 0., 0.));

		id_distal      = model.AddBody(0, Xtrans(Vector3d(0., 0., 0.)),
									   joint_distal, body_distal, "distal");

		id_proximal    = model.AddBody(id_distal, Xtrans(Vector3d(0., 0., 0.)),
									   joint_proximal, body_proximal,
									   "proximal");

		q     = VectorNd::Constant((size_t) model.dof_count, 0.);
		qdot  = VectorNd::Constant((size_t) model.dof_count, 0.);
		qddot = VectorNd::Constant((size_t) model.dof_count, 0.);
		tau   = VectorNd::Constant((size_t) model.dof_count, 0.);

		ad_model = ADModel(model);

		body_point = RigidBodyDynamics::Math::Vector3d (0., arm_l, 0.);

		ClearLogOutput();
	}

	RigidBodyDynamics::Model model;
	ADModel ad_model;

	unsigned int id_distal, id_proximal;
	RigidBodyDynamics::Body  body_distal, body_proximal;
	RigidBodyDynamics::Joint joint_distal, joint_proximal;

	RigidBodyDynamics::Math::VectorNd q;
	RigidBodyDynamics::Math::VectorNd qdot;
	RigidBodyDynamics::Math::VectorNd qddot;
	RigidBodyDynamics::Math::VectorNd tau;

	RigidBodyDynamics::Math::Vector3d body_point;
};

struct Arm2DofZ {
	Arm2DofZ() {
		using namespace RigidBodyDynamics;
		using namespace RigidBodyDynamics::Math;

		ClearLogOutput();

		double arm_l = 0.5;
		double arm_r = 0.1;
		double arm_m = 3.5;
		Matrix3d arm_I = .5 * arm_m * Matrix3d(
					.5 * (arm_r * arm_r + arm_l * arm_l / 3.), 0., 0.,
					arm_r * arm_r, 0., 0.,
					.5 * (arm_r * arm_r + arm_l * arm_l / 3.), 0., 0.);
		Vector3d arm_com (0., .5 * arm_l, 0.);

		body_distal    = Body(arm_m, arm_com, arm_I);
		body_proximal  = Body(arm_m, arm_com, arm_I);

		joint_distal   = Joint(SpatialVector(0., 0., 1., 0., 0., 0.));
		joint_proximal = Joint(SpatialVector(0., 0., 1., 0., 0., 0.));

		id_distal      = model.AddBody(0, Xtrans(Vector3d(0., 0., 0.)),
									   joint_distal, body_distal, "distal");

		id_proximal    = model.AddBody(id_distal, Xtrans(Vector3d(0., 0., 0.)),
									   joint_proximal, body_proximal,
									   "proximal");

		q     = VectorNd::Constant((size_t) model.dof_count, 0.);
		qdot  = VectorNd::Constant((size_t) model.dof_count, 0.);
		qddot = VectorNd::Constant((size_t) model.dof_count, 0.);
		tau   = VectorNd::Constant((size_t) model.dof_count, 0.);

		ad_model = ADModel(model);

		body_point = RigidBodyDynamics::Math::Vector3d (0., arm_l, 0.);

		ClearLogOutput();
	}

	RigidBodyDynamics::Model model;
	ADModel ad_model;

	unsigned int id_distal, id_proximal;
	RigidBodyDynamics::Body  body_distal, body_proximal;
	RigidBodyDynamics::Joint joint_distal, joint_proximal;

	RigidBodyDynamics::Math::VectorNd q;
	RigidBodyDynamics::Math::VectorNd qdot;
	RigidBodyDynamics::Math::VectorNd qddot;
	RigidBodyDynamics::Math::VectorNd tau;

	RigidBodyDynamics::Math::Vector3d body_point;
};
