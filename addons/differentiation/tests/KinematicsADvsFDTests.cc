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

// -----------------------------------------------------------------------------

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

// -----------------------------------------------------------------------------

template <typename T>
void CalcPointAccelerationTemplate(T & obj) {
	Model   & model    = obj.model;
	ADModel & ad_model = obj.ad_model;
	int const nq       = model.dof_count;

	VectorNd q    = VectorNd::Zero(nq);
	VectorNd qd   = VectorNd::Zero(nq);
	VectorNd qdd  = VectorNd::Zero(nq);
	MatrixNd q_dirs   = MatrixNd::Zero(nq, 3 * nq);
	MatrixNd qd_dirs  = MatrixNd::Zero(nq, 3 * nq);
	MatrixNd qdd_dirs = MatrixNd::Zero(nq, 3 * nq);

	q_dirs.block(0, 0, nq, nq)        = MatrixNd::Identity(nq, nq);
	qd_dirs.block(0, nq, nq, nq)      = MatrixNd::Identity(nq, nq);
	qdd_dirs.block(0, 2 * nq, nq, nq) = MatrixNd::Identity(nq, nq);

	int ndirs = 3 * nq;
	MatrixNd ad_derivative(3, ndirs);
	MatrixNd fd_derivative(3, ndirs);

	Vector3d reference_point = Vector3dZero;
	int nTrials = 0;
	do {
		for (unsigned i = 1; i < model.mBodies.size(); i++) {
			unsigned id = model.mBodyNameMap[model.GetBodyName(i)];
			Vector3d ad_a = RigidBodyDynamics::AD::CalcPointAcceleration(model,
					ad_model, q, q_dirs, qd, qd_dirs, qdd, qdd_dirs, id,
					reference_point, ad_derivative);
			Vector3d fd_a = RigidBodyDynamics::FD::CalcPointAcceleration(model,
					q, q_dirs, qd, qd_dirs, qdd, qdd_dirs, id,
					reference_point, fd_derivative);
			CHECK_ARRAY_CLOSE(ad_a.data(), fd_a.data(), 3, 1e-12);
			CHECK_ARRAY_CLOSE(fd_derivative.data(), ad_derivative.data(),
							  3 * ndirs, 1e-7);
		}
		q   = VectorNd::Random(nq);
		qd  = VectorNd::Random(nq);
		qdd = VectorNd::Random(nq);
	} while(nTrials++ < 10);
}


TEST_FIXTURE ( CartPendulum, CartPendulumCalcPointAcceleration) {
	CalcPointAccelerationTemplate(*this);
}

TEST_FIXTURE ( Arm2DofX, Arm2DofXCalcPointAcceleration) {
	CalcPointAccelerationTemplate(*this);
}

TEST_FIXTURE ( Arm2DofZ, Arm2DofZCalcPointAcceleration) {
	CalcPointAccelerationTemplate(*this);
}

TEST_FIXTURE ( Arm3DofXZYp, Arm3DofXZYpCalcPointAcceleration) {
	CalcPointAccelerationTemplate(*this);
}

TEST_FIXTURE ( Arm3DofXZZp, Arm3DofXZZpCalcPointAcceleration) {
	CalcPointAccelerationTemplate(*this);
}

// -----------------------------------------------------------------------------


