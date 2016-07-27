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

/*
template <typename T>
void CalcPointJacobianTemplate(T & obj) {
    // rename for convenience
    Model& model = obj.model;
    ADModel& ad_model = obj.ad_model;

    ConstraintSet& cs = obj.cs;
    ADConstraintSet& ad_cs = obj.ad_cs;

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
    CalcPointJacobian (model, q, cs, G, update_kinematics);

    // int nTrials = 0;
    // do {
        // call AD version
        AD::CalcPointJacobian (
            model, ad_model,
            q, q_dirs,
            cs, ad_cs,
            G_ad, ad_derivative,
            update_kinematics
        );

        // call FD version
        FD::CalcPointJacobian (
            model, ad_model,
            q, q_dirs,
            cs, ad_cs,
            G_fd, fd_derivative,
            update_kinematics
        );

        // compare nominal results
        CHECK_ARRAY_CLOSE(G_ad.data(), G.data(), G.size(), TEST_PREC);
        CHECK_ARRAY_CLOSE(G_fd.data(), G.data(), G.size(), TEST_PREC);
        CHECK_ARRAY_CLOSE(G_ad.data(), G_fd.data(), G.size(), TEST_PREC);

        for (unsigned idir = 0; idir < ndirs; idir++) {
            // CHECK_ARRAY_CLOSE(ad_E.data(), fd_E.data(), 9, 1e-6);
            // for (int idir = 0; idir < ndirs; idir++) {
            //     CHECK_ARRAY_CLOSE(fd_derivative[idir].data(), ad_derivative[idir].data(), 9, 1e-6);
            // }
        }
        // q = VectorNd::Random(nq);
    // } while(nTrials++ < 10);
}

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
