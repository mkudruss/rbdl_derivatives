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
