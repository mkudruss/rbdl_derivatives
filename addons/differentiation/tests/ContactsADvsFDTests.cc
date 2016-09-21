#include <UnitTest++.h>

#include <iostream>

#include "rbdl/Model.h"
#include "rbdl/Kinematics.h"
#include "rbdl/rbdl_mathutils.h"

#include "ContactsAD.h"
#include "ContactsFD.h"

#include "Fixtures.h"

using namespace std;
using namespace RigidBodyDynamics;
using namespace RigidBodyDynamics::Math;

double const TEST_PREC = 1.0e-8;

// -----------------------------------------------------------------------------

template <typename T>
void CalcContactJacobianTemplate(T & obj) {
    // rename for convenience
    Model& model = *(obj.model);
    ADModel& ad_model = *(obj.ad_model);

    ConstraintSet& cs = obj.constraint_set;
    ADConstraintSet& ad_cs = obj.ad_constraint_set;

    // set up input quantities
    int const nq = model.dof_count;
    int const ndirs = nq;

    bool update_kinematics = true;
    VectorNd q = VectorNd::Zero(nq);
    MatrixNd q_dirs = MatrixNd::Identity(nq, nq);

    // set up no output quantities
    MatrixNd G = MatrixNd::Zero (3, nq);
    MatrixNd G_ad = MatrixNd::Zero (3, nq);
    MatrixNd G_fd = MatrixNd::Zero (3, nq);

    // set up derivative output quantities
    vector<MatrixNd> derivative_ad (ndirs, G_ad);
    vector<MatrixNd> derivative_fd (ndirs, G_fd);

    // call nominal version
    RigidBodyDynamics::CalcContactJacobian(model, q, cs, G, update_kinematics);

    // call AD version
    RigidBodyDynamics::AD::CalcContactJacobian(
        model, ad_model,
        q, q_dirs,
        cs, ad_cs,
        G, derivative_ad,
        update_kinematics
    );

    // call FD version
    RigidBodyDynamics::FD::CalcContactJacobian(
        model, ad_model,
        q, q_dirs,
        cs, ad_cs,
        G, derivative_fd,
        update_kinematics
    );

    // compare nominal results
    CHECK_ARRAY_CLOSE(G_ad.data(), G.data(), G.size(), TEST_PREC);
    CHECK_ARRAY_CLOSE(G_fd.data(), G.data(), G.size(), TEST_PREC);
    CHECK_ARRAY_CLOSE(G_ad.data(), G_fd.data(), G.size(), TEST_PREC);

    for (unsigned idir = 0; idir < ndirs; idir++) {
        CHECK_ARRAY_CLOSE(
            derivative_ad[idir].data(),
            derivative_fd[idir].data(),
            derivative_ad[idir].size(),
            TEST_PREC
        );
    }
}

TEST_FIXTURE (FixedBase6DoF, FixedBase6DoFCalcContactJacobian) {

    // add contacts and bind them to constraint set
    constraint_set.AddConstraint (contact_body_id, Vector3d (1., 0., 0.), contact_normal);
    constraint_set.AddConstraint (contact_body_id, Vector3d (0., 1., 0.), contact_normal);
    constraint_set.Bind (*model);

/*
    VectorNd QDDot_lagrangian = VectorNd::Constant (model->mBodies.size() - 1, 0.);

    ClearLogOutput();

    ForwardDynamicsContactsDirect (*model, Q, QDot, Tau, constraint_set_lagrangian, QDDot_lagrangian);
    ForwardDynamicsContactsKokkevis (*model, Q, QDot, Tau, constraint_set, QDDot_contacts);

    point_accel_lagrangian = CalcPointAcceleration (*model, Q, QDot, QDDot_lagrangian, contact_body_id, contact_point, true);
    point_accel_contacts = CalcPointAcceleration (*model, Q, QDot, QDDot_contacts, contact_body_id, contact_point, true);

    // check whether FDContactsLagrangian and FDContacts match
    CHECK_ARRAY_CLOSE (
            constraint_set_lagrangian.force.data(),
            constraint_set.force.data(),
            constraint_set.size(), TEST_PREC
            );

    // check whether the point accelerations match
    CHECK_ARRAY_CLOSE (point_accel_lagrangian.data(), point_accel_contacts.data(), 3, TEST_PREC);

    // check whether the generalized accelerations match
    CHECK_ARRAY_CLOSE (QDDot_lagrangian.data(), QDDot_contacts.data(), QDDot_lagrangian.size(), TEST_PREC);
    */
    CalcContactJacobianTemplate(*this);
}


/*
TEST_FIXTURE ( CartPendulum, CartPendulumCalcBodyWorldOrientation) {
    CalcBodyWorldOrientationTemplate(*this);
}

TEST_FIXTURE ( Arm2DofX, Arm2DofXCalcBodyWorldOrientation) {
    CalcBodyWorldOrientationTemplate(*this);
}

TEST_FIXTURE ( Arm2DofZ, Arm2DofZCalcBodyWorldOrientation) {
    CalcBodyWorldOrientationTemplate(*this);
}

TEST_FIXTURE ( Arm3DofXZYp, Arm3DofXZYpCalcBodyWorldOrientation) {
    CalcBodyWorldOrientationTemplate(*this);
}

TEST_FIXTURE ( Arm3DofXZZp, Arm3DofXZZpCalcBodyWorldOrientation) {
    CalcBodyWorldOrientationTemplate(*this);
}
*/

// -----------------------------------------------------------------------------
