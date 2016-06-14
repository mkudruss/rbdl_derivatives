#include <UnitTest++.h>

#include <iostream>

#include "rbdl/Model.h"
#include "rbdl/Kinematics.h"
#include "rbdl/rbdl_mathutils.h"

#include "KinematicsAD.h"
#include "KinematicsFD.h"

#include "Fixtures.h"

using namespace std;
using namespace RigidBodyDynamics;
using namespace RigidBodyDynamics::Math;

const double TEST_PREC = 1.0e-8;

template <typename T>
void CalcBodyWorldOrientationTemplate(T & obj) {
    Model & model = obj.model;
    ADModel & ad_model = obj.ad_model;
    int const nq = model.dof_count;

    VectorNd q = VectorNd::Zero(nq);
    MatrixNd q_dirs = MatrixNd::Identity(nq, nq);

    int ndirs = nq;
    vector<Matrix3d> ad_derivative(ndirs);
    vector<Matrix3d> fd_derivative(ndirs);

    int nTrials = 0;
    do {

        for (unsigned i = 1; i < model.mBodies.size(); i++) {
            unsigned id = model.mBodyNameMap[model.GetBodyName(i)];
            Matrix3d ad_E = RigidBodyDynamics::AD::CalcBodyWorldOrientation(model, ad_model, q, q_dirs, id, ad_derivative);
            Matrix3d fd_E = RigidBodyDynamics::FD::CalcBodyWorldOrientation(model, q, q_dirs, id, fd_derivative);
            CHECK_ARRAY_CLOSE(ad_E.data(), fd_E.data(), 9, 1e-6);
            for (int idir = 0; idir < ndirs; idir++) {
                CHECK_ARRAY_CLOSE(fd_derivative[idir].data(), ad_derivative[idir].data(), 9, 1e-6);
            }
        }
        q = VectorNd::Random(nq);
    } while(nTrials++ < 10);
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

