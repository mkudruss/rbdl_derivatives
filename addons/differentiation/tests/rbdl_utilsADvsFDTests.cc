#include <UnitTest++.h>

#include <iomanip>
#include <iostream>

#include "Fixtures.h"
#include "rbdl_utilsFD.h"
#include "rbdl_utilsAD.h"

using namespace RigidBodyDynamics;
using namespace RigidBodyDynamics::Math;
using namespace RigidBodyDynamics::Utils;

const double TEST_PREC = 1.0e-12;

TEST_FIXTURE ( CartPendulum, CartPendulumCalcCenterOfMass) {
    q[0] = 0.3;
    q[1] = -0.2;
    qdot[0] = .1;
    qdot[1] = -.15;

    int nrows = model.dof_count;
    int ndirs = 2 * model.dof_count;

    MatrixNd q_dirs = MatrixNd::Zero(nrows, ndirs);
    MatrixNd qdot_dirs = MatrixNd::Zero(nrows, ndirs);

    q_dirs.block(0, 0, nrows, model.dof_count)
            = MatrixNd::Identity(nrows, model.dof_count);
    q_dirs.block(0, model.dof_count, nrows, model.dof_count)
            = MatrixNd::Zero(nrows, model.dof_count);
    qdot_dirs.block(0, 0, nrows, model.dof_count)
            = MatrixNd::Zero(nrows, model.dof_count);
    qdot_dirs.block(0, model.dof_count, nrows, model.dof_count)
            = MatrixNd::Identity(nrows, model.dof_count);

    double   ad_mass;
    Vector3d ad_com;
    MatrixNd ad_d_com(3, q_dirs.cols());
    Vector3d ad_comVelocity(0., 0., 0.);
    MatrixNd ad_d_comVelocity(3, ndirs);
    Vector3d ad_angMomentum(0., 0., 0.);
    MatrixNd ad_d_angMomentum(3, ndirs);
    Utils::AD::CalcCenterOfMass(model, ad_model, q, q_dirs, qdot, qdot_dirs, ad_mass,
                        ad_com, ad_d_com,
                        &ad_comVelocity, &ad_d_comVelocity,
                        &ad_angMomentum, &ad_d_angMomentum, true);

    double   fd_mass;
    Vector3d fd_com;
    MatrixNd fd_d_com(3, ndirs);
    Vector3d fd_comVelocity(0., 0., 0.);
    MatrixNd fd_d_comVelocity(3, ndirs);
    Vector3d fd_angMomentum(0., 0., 0.);
    MatrixNd fd_d_angMomentum(3, ndirs);
    FD::CalcCenterOfMass(model, q, q_dirs, qdot, qdot_dirs, fd_mass, fd_com,
                         fd_d_com, &fd_comVelocity, &fd_d_comVelocity,
                         &fd_angMomentum, &fd_d_angMomentum);

    CHECK_EQUAL(ad_mass, fd_mass);
    CHECK_ARRAY_CLOSE(ad_com.data(), fd_com.data(), 3, TEST_PREC);
    CHECK_ARRAY_CLOSE(ad_d_com.data(), fd_d_com.data(), 3 * ndirs, 1e-8);
    CHECK_ARRAY_CLOSE(ad_comVelocity.data(), fd_comVelocity.data(), 3, 1e-8);
    CHECK_ARRAY_CLOSE(ad_d_comVelocity.data(), fd_d_comVelocity.data(), 3 * ndirs, 1e-8);
    CHECK_ARRAY_CLOSE(ad_angMomentum.data(), fd_angMomentum.data(), 3, 1e-8);
    CHECK_ARRAY_CLOSE(ad_d_angMomentum.data(), fd_d_angMomentum.data(), 3 * ndirs, 1e-8);
}

TEST_FIXTURE ( CartPendulum, CartPendulumCalcPotentialEnergy) {
    q[0] = 0.3;
    q[1] = -0.2;
    qdot[0] = .1;
    qdot[1] = -.15;

    int nrows = model.dof_count;
    int ndirs = 2 * model.dof_count;

    MatrixNd q_dirs = MatrixNd::Zero(nrows, ndirs);
    MatrixNd qdot_dirs = MatrixNd::Zero(nrows, ndirs);

    q_dirs.block(0, 0, nrows, model.dof_count)
            = MatrixNd::Identity(nrows, model.dof_count);
    q_dirs.block(0, model.dof_count, nrows, model.dof_count)
            = MatrixNd::Zero(nrows, model.dof_count);
    qdot_dirs.block(0, 0, nrows, model.dof_count)
            = MatrixNd::Zero(nrows, model.dof_count);
    qdot_dirs.block(0, model.dof_count, nrows, model.dof_count)
            = MatrixNd::Identity(nrows, model.dof_count);

    MatrixNd fd_pote = MatrixNd::Zero(1, ndirs);
    FD::CalcPotentialEnergy(model, q, q_dirs, fd_pote);

    MatrixNd ad_pote = MatrixNd::Zero(1, ndirs);
    Utils::AD::CalcPotentialEnergy(model, ad_model, q, q_dirs, ad_pote, true);

    CHECK_ARRAY_CLOSE(fd_pote.data(), ad_pote.data(), 1 * ndirs, 1e-8);
}

TEST_FIXTURE ( CartPendulum, CartPendulumCalcKineticEnergy) {
    q[0] = 0.3;
    q[1] = -0.2;
    qdot[0] = .1;
    qdot[1] = -.15;

    int nrows = model.dof_count;
    int ndirs = 2 * model.dof_count;

    MatrixNd q_dirs = MatrixNd::Zero(nrows, ndirs);
    MatrixNd qdot_dirs = MatrixNd::Zero(nrows, ndirs);

    q_dirs.block(0, 0, nrows, model.dof_count)
            = MatrixNd::Identity(nrows, model.dof_count);
    q_dirs.block(0, model.dof_count, nrows, model.dof_count)
            = MatrixNd::Zero(nrows, model.dof_count);
    qdot_dirs.block(0, 0, nrows, model.dof_count)
            = MatrixNd::Zero(nrows, model.dof_count);
    qdot_dirs.block(0, model.dof_count, nrows, model.dof_count)
            = MatrixNd::Identity(nrows, model.dof_count);

    MatrixNd fd_kine = MatrixNd::Zero(1, ndirs);
    Utils::FD::CalcKineticEnergy(model, q, q_dirs, qdot, qdot_dirs, fd_kine);

    MatrixNd ad_kine = MatrixNd::Zero(1, ndirs);
    Utils::AD::CalcKineticEnergy(model, ad_model, q, q_dirs, qdot, qdot_dirs, ad_kine, true);

    CHECK_ARRAY_CLOSE(fd_kine.data(), ad_kine.data(), 1 * ndirs, 1e-6);
}

TEST_FIXTURE( Arm2Dof, Arm2DofCalcCenterOfMass)
{
    q[0]    =  M_PI / 4.;
    q[1]    = -M_PI / 3.;
    qdot[0] = -.2;
    qdot[1] =  .11;

    int nrows = model.dof_count;
    int ndirs = 2 * model.dof_count;

    MatrixNd q_dirs = MatrixNd::Zero(nrows, ndirs);
    MatrixNd qdot_dirs = MatrixNd::Zero(nrows, ndirs);

    q_dirs.block(0, 0, nrows, model.dof_count)
            = MatrixNd::Identity(nrows, model.dof_count);
    q_dirs.block(0, model.dof_count, nrows, model.dof_count)
            = MatrixNd::Zero(nrows, model.dof_count);
    qdot_dirs.block(0, 0, nrows, model.dof_count)
            = MatrixNd::Zero(nrows, model.dof_count);
    qdot_dirs.block(0, model.dof_count, nrows, model.dof_count)
            = MatrixNd::Identity(nrows, model.dof_count);

    double   ad_mass;
    Vector3d ad_com;
    MatrixNd ad_d_com(3, q_dirs.cols());
    Vector3d ad_comVelocity(0., 0., 0.);
    MatrixNd ad_d_comVelocity(3, ndirs);
    Vector3d ad_angMomentum(0., 0., 0.);
    MatrixNd ad_d_angMomentum(3, ndirs);
    Utils::AD::CalcCenterOfMass(model, ad_model, q, q_dirs, qdot, qdot_dirs, ad_mass,
                        ad_com, ad_d_com,
                        &ad_comVelocity, &ad_d_comVelocity,
                        &ad_angMomentum, &ad_d_angMomentum, true);

    double   fd_mass;
    Vector3d fd_com;
    MatrixNd fd_d_com(3, ndirs);
    Vector3d fd_comVelocity(0., 0., 0.);
    MatrixNd fd_d_comVelocity(3, ndirs);
    Vector3d fd_angMomentum(0., 0., 0.);
    MatrixNd fd_d_angMomentum(3, ndirs);
    FD::CalcCenterOfMass(model, q, q_dirs, qdot, qdot_dirs, fd_mass, fd_com,
                         fd_d_com, &fd_comVelocity, &fd_d_comVelocity,
                         &fd_angMomentum, &fd_d_angMomentum);

    CHECK_EQUAL(ad_mass, fd_mass);
    CHECK_ARRAY_CLOSE(ad_com.data(), fd_com.data(), 3, TEST_PREC);
    CHECK_ARRAY_CLOSE(ad_d_com.data(), fd_d_com.data(), 3 * ndirs, 1e-8);
    CHECK_ARRAY_CLOSE(ad_comVelocity.data(), fd_comVelocity.data(), 3, 1e-8);
    CHECK_ARRAY_CLOSE(ad_d_comVelocity.data(), fd_d_comVelocity.data(), 3 * ndirs, 1e-8);
    CHECK_ARRAY_CLOSE(ad_angMomentum.data(), fd_angMomentum.data(), 3, 1e-8);
    CHECK_ARRAY_CLOSE(ad_d_angMomentum.data(), fd_d_angMomentum.data(), 3 * ndirs, 1e-8);
}
