#include <iostream>
#include <unittest++/UnitTest++.h>


#include "rbdl/Model.h"
#include "rbdl/Dynamics.h"
#include "rbdl/rbdl_mathutils.h"

#include "Fixtures.h"
#include "Human36Fixture.h"
#include "ModelAD.h"
#include "ModelED.h"
#include "DynamicsAD.h"
#include "DynamicsED.h"
#include "DynamicsFD.h"
#include "DynamicsFDC.h"
#include "rbdl_mathutilsAD.h"

#include "ModelCheckADvsFD.h"

using namespace std;
using namespace RigidBodyDynamics;
using namespace RigidBodyDynamics::Math;

static const double TEST_PREC = 1.0e-8;
static const double EPS = std::sqrt(numeric_limits<double>::epsilon());

// -----------------------------------------------------------------------------
/*! \brief We want to compute the derivative of d X(q(t)) / dt
*/
template <typename T>
void analytic_derivative_of_spatial_transform(
  T & obj,
  std::function<SpatialTransform (const double &)> get_transform,
  const SpatialVector & S,
  const double TEST_PRECISION=1e-8
) {
  const VectorNd Q = VectorNd::Random(1);
  const VectorNd Qdot = VectorNd::Random(1);

  // nominal evaluation
  SpatialMatrix X = get_transform(Q(0)).toMatrix();

  // derivative evaluation
  SpatialMatrix X_fd = get_transform(Q(0) + EPS * Qdot(0)).toMatrix();
  X_fd = (X_fd - X) / EPS;

  // analytic derivative
  SpatialMatrix X_ad = -crossm(S*Qdot(0))*X;

  // print results
  // cout << "X_fd: " << endl << X_fd.transpose() << endl;
  // cout << "X_ad: " << endl << X_ad.transpose() << endl;
  // cout << "error(max): " << endl << (X_ad - X_fd).cwiseAbs().transpose()
  //      << endl << " (" << (X_ad - X_fd).cwiseAbs().maxCoeff() << ")" << endl;
  // cout << endl;

  // test result
  CHECK_ARRAY_CLOSE(X_fd.data(), X_ad.data(), X.size(), TEST_PRECISION);
}

TEST(Xrotx_analytic_derivative_of_spatial_transform){
  srand (421337);
  SpatialVector S = SpatialVector::Zero();
  S(0) = 1.0;
  analytic_derivative_of_spatial_transform(*this, Xrotx, S, 1e-8);
}

TEST(Xroty_analytic_derivative_of_spatial_transform){
  srand (421337);
  SpatialVector S = SpatialVector::Zero();
  S(1) = 1.0;
  analytic_derivative_of_spatial_transform(*this, Xroty, S, 1e-8);
}

TEST(Xrotz_analytic_derivative_of_spatial_transform){
  srand (421337);
  SpatialVector S = SpatialVector::Zero();
  S(2) = 1.0;
  analytic_derivative_of_spatial_transform(*this, Xrotz, S, 1e-8);
}

TEST(Xtransx_analytic_derivative_of_spatial_transform){
  srand (421337);
  SpatialVector S = SpatialVector::Zero();
  S(3) = 1.0;
  std::function<SpatialTransform (const double &)> Xtransx
    = [](const double &q) { return Xtrans(Vector3d(q, 0., 0.)); };
  analytic_derivative_of_spatial_transform(*this, Xtransx, S, 1e-8);
}

TEST(Xtransy_analytic_derivative_of_spatial_transform){
  srand (421337);
  SpatialVector S = SpatialVector::Zero();
  S(4) = 1.0;
  std::function<SpatialTransform (const double &)> Xtransy
    = [](const double &q) { return Xtrans(Vector3d(0., q, 0.)); };
  analytic_derivative_of_spatial_transform(*this, Xtransy, S, 1e-8);
}

TEST(Xtransz_analytic_derivative_of_spatial_transform){
  srand (421337);
  SpatialVector S = SpatialVector::Zero();
  S(5) = 1.0;
  std::function<SpatialTransform (const double &)> Xtransz
    = [](const double &q) { return Xtrans(Vector3d(0., 0., q)); };
  analytic_derivative_of_spatial_transform(*this, Xtransz, S, 1e-8);
}

// -----------------------------------------------------------------------------
/*! \brief We want to compute the derivative of d X(q(t)) / dt * v
*/
template <typename T>
void analytic_derivative_of_spatial_transform_times_v(
  T & obj,
  std::function<SpatialTransform (const double &)> get_transform,
  const SpatialVector & S,
  const double TEST_PRECISION=1e-8
) {
  const VectorNd Q = VectorNd::Random(1);
  const VectorNd Qdot = VectorNd::Random(1);

  // nominal evaluation
  SpatialTransform X = get_transform(Q(0));
  SpatialVector v = SpatialVector::Random();
  SpatialVector res = X.apply(v);

  // derivative evaluation
  SpatialTransform X_fd = get_transform(Q(0) + EPS * Qdot(0));
  SpatialVector res_fd = (X_fd.apply(v) - res) / EPS;

  // analytic derivative
  SpatialVector res_ad = -crossm(S*Qdot(0), res);

  // print results
  // cout << "res_fd:      " << res_fd.transpose() << endl;
  // cout << "res_ad:      " << res_ad.transpose() << endl;
  // cout << "error(max):  " << (res_ad - res_fd).cwiseAbs().transpose()
  //      << " (" << (res_ad - res_fd).cwiseAbs().maxCoeff() << ")" << endl;
  // cout << endl;

  // test result
  CHECK_ARRAY_CLOSE(res_fd.data(), res_ad.data(), res.size(), TEST_PRECISION);
}

TEST(Xrotx_analytic_derivative_of_spatial_transform_times_v){
  srand (421337);
  SpatialVector S = SpatialVector::Zero();
  S(0) = 1.0;
  analytic_derivative_of_spatial_transform_times_v(*this, Xrotx, S, 1e-8);
}

TEST(Xroty_analytic_derivative_of_spatial_transform_times_v){
  srand (421337);
  SpatialVector S = SpatialVector::Zero();
  S(1) = 1.0;
  analytic_derivative_of_spatial_transform_times_v(*this, Xroty, S, 1e-7);
}

TEST(Xrotz_analytic_derivative_of_spatial_transform_times_v){
  srand (421337);
  SpatialVector S = SpatialVector::Zero();
  S(2) = 1.0;
  analytic_derivative_of_spatial_transform_times_v(*this, Xrotz, S, 1e-8);
}

TEST(Xtransx_analytic_derivative_of_spatial_transform_times_v){
  srand (421337);
  SpatialVector S = SpatialVector::Zero();
  S(3) = 1.0;
  std::function<SpatialTransform (const double &)> Xtransx
    = [](const double &q) { return Xtrans(Vector3d(q, 0., 0.)); };
  analytic_derivative_of_spatial_transform_times_v(*this, Xtransx, S, 1e-8);
}

TEST(Xtransy_analytic_derivative_of_spatial_transform_times_v){
  srand (421337);
  SpatialVector S = SpatialVector::Zero();
  S(4) = 1.0;
  std::function<SpatialTransform (const double &)> Xtransy
    = [](const double &q) { return Xtrans(Vector3d(0., q, 0.)); };
  analytic_derivative_of_spatial_transform_times_v(*this, Xtransy, S, 1e-8);
}

TEST(Xtransz_analytic_derivative_of_spatial_transform_times_v){
  srand (421337);
  SpatialVector S = SpatialVector::Zero();
  S(5) = 1.0;
  std::function<SpatialTransform (const double &)> Xtransz
    = [](const double &q) { return Xtrans(Vector3d(0., 0., q)); };
  analytic_derivative_of_spatial_transform_times_v(*this, Xtransz, S, 1e-8);
}

// -----------------------------------------------------------------------------
/*! \brief We want to compute the derivative of d X(q(t))^* / dt * f
*/
template <typename T>
void analytic_derivative_of_spatial_transform_transpose_times_f(
  T & obj,
  std::function<SpatialTransform (const double &)> get_transform,
  const SpatialVector & S,
  const double TEST_PRECISION=1e-8
) {
  const VectorNd Q = VectorNd::Random(1);
  const VectorNd Qdot = VectorNd::Random(1);

  // nominal evaluation
  SpatialTransform X = get_transform(Q(0));
  SpatialVector f = SpatialVector::Random();
  SpatialVector res = X.applyTranspose(f);

  // derivative evaluation
  SpatialTransform X_fd = get_transform(Q(0) + EPS * Qdot(0));
  SpatialVector res_fd = (X_fd.applyTranspose(f) - res) / EPS;

  // analytic derivative
  SpatialVector res_ad = crossf(S*Qdot(0), res);

  // print results
  // cout << "res_fd:      " << res_fd.transpose() << endl;
  // cout << "res_ad:      " << res_ad.transpose() << endl;
  // cout << "error(max):  " << (res_ad - res_fd).cwiseAbs().transpose()
  //      << " (" << (res_ad - res_fd).cwiseAbs().maxCoeff() << ")" << endl;
  // cout << endl;

  // test result
  CHECK_ARRAY_CLOSE(res_fd.data(), res_ad.data(), res.size(), TEST_PRECISION);
}

TEST(Xrotx_analytic_derivative_of_spatial_transform_transpose_times_f){
  srand (421337);
  SpatialVector S = SpatialVector::Zero();
  S(0) = 1.0;
  analytic_derivative_of_spatial_transform_transpose_times_f(*this, Xrotx, S, 1e-8);
}

TEST(Xroty_analytic_derivative_of_spatial_transform_transpose_times_f){
  srand (421337);
  SpatialVector S = SpatialVector::Zero();
  S(1) = 1.0;
  analytic_derivative_of_spatial_transform_transpose_times_f(*this, Xroty, S, 1e-7);
}

TEST(Xrotz_analytic_derivative_of_spatial_transform_transpose_times_f){
  srand (421337);
  SpatialVector S = SpatialVector::Zero();
  S(2) = 1.0;
  analytic_derivative_of_spatial_transform_transpose_times_f(*this, Xrotz, S, 1e-8);
}

TEST(Xtransx_analytic_derivative_of_spatial_transform_transpose_times_f){
  srand (421337);
  SpatialVector S = SpatialVector::Zero();
  S(3) = 1.0;
  std::function<SpatialTransform (const double &)> Xtransx
    = [](const double &q) { return Xtrans(Vector3d(q, 0., 0.)); };
  analytic_derivative_of_spatial_transform_transpose_times_f(*this, Xtransx, S, 1e-8);
}

TEST(Xtransy_analytic_derivative_of_spatial_transform_transpose_times_f){
  srand (421337);
  SpatialVector S = SpatialVector::Zero();
  S(4) = 1.0;
  std::function<SpatialTransform (const double &)> Xtransy
    = [](const double &q) { return Xtrans(Vector3d(0., q, 0.)); };
  analytic_derivative_of_spatial_transform_transpose_times_f(*this, Xtransy, S, 1e-8);
}

TEST(Xtransz_analytic_derivative_of_spatial_transform_transpose_times_f){
  srand (421337);
  SpatialVector S = SpatialVector::Zero();
  S(5) = 1.0;
  std::function<SpatialTransform (const double &)> Xtransz
    = [](const double &q) { return Xtrans(Vector3d(0., 0., q)); };
  analytic_derivative_of_spatial_transform_transpose_times_f(*this, Xtransz, S, 1e-8);
}


// -----------------------------------------------------------------------------
/*! \brief We want to compute the derivative of d ^B I(q(t)) / dt
*/

template <typename T>
void analytic_derivative_of_spatial_inertia(
  T & obj,
  std::function<SpatialTransform (const double &)> get_transform,
  const SpatialVector & S,
  const double TEST_PRECISION=1e-8
) {
  const VectorNd Q = VectorNd::Random(1);
  const VectorNd Qdot = VectorNd::Random(1);

  VectorNd v = VectorNd::Random(1);
  const double m = v(0);
  const Vector3d h = Vector3d::Random();
  const Matrix3d Ic = Matrix3d::Random();
  const SpatialRigidBodyInertia IC = SpatialRigidBodyInertia(m, h, Ic);

  // nominal evaluation
  SpatialTransform X = get_transform(Q(0));
  SpatialMatrix I = X.applyTranspose(IC).toMatrix();

  // derivative evaluation
  SpatialTransform X_fd = get_transform(Q(0) + EPS * Qdot(0));
  SpatialMatrix I_fd = (X_fd.applyTranspose(IC).toMatrix() - I) / EPS;

  // analytic derivative
  SpatialMatrix I_ad = crossf(S*Qdot(0))*I - I*crossm(S*Qdot(0));

  // print results
  // cout << "I_fd: " << endl << I_fd.transpose() << endl;
  // cout << "I_ad: " << endl << I_ad.transpose() << endl;
  // cout << "error(max): " << endl << (I_ad - I_fd).cwiseAbs().transpose()
  //      << endl << " (" << (I_ad - I_fd).cwiseAbs().maxCoeff() << ")" << endl;
  // cout << endl;

  // test Iult
  CHECK_ARRAY_CLOSE(I_fd.data(), I_ad.data(), I_fd.size(), TEST_PRECISION);
}

TEST(Xrotx_analytic_derivative_of_spatial_inertia){
  srand (421337);
  SpatialVector S = SpatialVector::Zero();
  S(0) = 1.0;
  analytic_derivative_of_spatial_inertia(*this, Xrotx, S, 1e-7);
}

TEST(Xroty_analytic_derivative_of_spatial_inertia){
  srand (421337);
  SpatialVector S = SpatialVector::Zero();
  S(1) = 1.0;
  analytic_derivative_of_spatial_inertia(*this, Xroty, S, 1e-7);
}

TEST(Xrotz_analytic_derivative_of_spatial_inertia){
  srand (421337);
  SpatialVector S = SpatialVector::Zero();
  S(2) = 1.0;
  analytic_derivative_of_spatial_inertia(*this, Xrotz, S, 1e-7);
}

TEST(Xtransx_analytic_derivative_of_spatial_inertia){
  srand (421337);
  SpatialVector S = SpatialVector::Zero();
  S(3) = 1.0;
  std::function<SpatialTransform (const double &)> Xtransx
    = [](const double &q) { return Xtrans(Vector3d(q, 0., 0.)); };
  analytic_derivative_of_spatial_inertia(*this, Xtransx, S, 1e-8);
}

TEST(Xtransy_analytic_derivative_of_spatial_inertia){
  srand (421337);
  SpatialVector S = SpatialVector::Zero();
  S(4) = 1.0;
  std::function<SpatialTransform (const double &)> Xtransy
    = [](const double &q) { return Xtrans(Vector3d(0., q, 0.)); };
  analytic_derivative_of_spatial_inertia(*this, Xtransy, S, 1e-8);
}

TEST(Xtransz_analytic_derivative_of_spatial_inertia){
  srand (421337);
  SpatialVector S = SpatialVector::Zero();
  S(5) = 1.0;
  std::function<SpatialTransform (const double &)> Xtransz
    = [](const double &q) { return Xtrans(Vector3d(0., 0., q)); };
  analytic_derivative_of_spatial_transform_times_v(*this, Xtransz, S, 1e-8);
}

// -----------------------------------------------------------------------------
/*! \brief We want to compute the derivative of d ^B I(q(t)) / dt
*/

// SpatialRigidBodyInertia derivative_of_X_applyTranspose_rbi(
//   const SpatialRigidBodyInertia& I,
//   const SpatialVector& S,
//   const SpatialDirection& Qdirs
// ) {
//   return SpatialRigidBodyInertia();
// }

Matrix3d cross(Vector3d r) {
  return Matrix3d (
        0., -r[2],  r[1],
      r[2],    0., -r[0],
     -r[1],  r[0],    0.
  );
};


template <typename T>
void efficient_analytic_derivative_of_spatial_inertia(
  T & obj,
  std::function<SpatialTransform (const double &)> get_transform,
  const SpatialVector & S,
  const double TEST_PRECISION=1e-8
) {
  const VectorNd Q = VectorNd::Random(1);
  const VectorNd Qdot = VectorNd::Random(1);

  VectorNd v = VectorNd::Random(1);
  const double m = v(0);
  const Vector3d h = Vector3d::Random();
  const Matrix3d Ic = Matrix3d::Random();
  const SpatialRigidBodyInertia IC = SpatialRigidBodyInertia(m, h, Ic);

  // nominal evaluation
  SpatialTransform X = get_transform(Q(0));
  SpatialRigidBodyInertia I = X.applyTranspose(IC);

  // derivative evaluation
  SpatialTransform X_fd = get_transform(Q(0) + EPS * Qdot(0));
  SpatialRigidBodyInertia I_fd = (X_fd.applyTranspose(IC) - I) * (1.0 / EPS);

  // analytic derivative
  // SpatialMatrix I_ad = crossf(S*Qdot(0))*I - I*crossm(S*Qdot(0));
  SpatialMatrix I_an = crossf(S*Qdot(0))*I.toMatrix() - I.toMatrix()*crossm(S*Qdot(0));

  // separate vector into 3D components
  Vector3d w = S.head(3)*Qdot(0);
  Vector3d v0 = S.tail(3)*Qdot(0);

  SpatialRigidBodyInertia I_ad
    = SpatialRigidBodyInertia (
        0,
        w.cross(I.h) + I.m * v0,
        Matrix3d(
          -I.Iyx*w[2] + I.Izx*w[1] - I.Iyx*w[2] + I.Izx*w[1] + 2.*(I.h[1]*v0[1] + I.h[2]*v0[2]),
           I.Ixx*w[2] - I.Izx*w[0] - I.Iyy*w[2] + I.Izy*w[1] -     I.h[0]*v0[1] - I.h[1]*v0[0] ,
          -I.Ixx*w[1] + I.Iyx*w[0] - I.Izy*w[2] + I.Izz*w[1] -     I.h[0]*v0[2] - I.h[2]*v0[0] ,

           I.Ixx*w[2] - I.Iyy*w[2] + I.Izy*w[1] - I.Izx*w[0] -     I.h[0]*v0[1] - I.h[1]*v0[0] ,
           I.Iyx*w[2] + I.Iyx*w[2] - I.Izy*w[0] - I.Izy*w[0] + 2.*(I.h[0]*v0[0] + I.h[2]*v0[2]),
           I.Izx*w[2] - I.Iyx*w[1] + I.Iyy*w[0] - I.Izz*w[0] -     I.h[1]*v0[2] - I.h[2]*v0[1] ,

          -I.Ixx*w[1] + I.Iyx*w[0] - I.Izy*w[2] + I.Izz*w[1] -     I.h[0]*v0[2] - I.h[2]*v0[0] ,
          -I.Iyx*w[1] + I.Iyy*w[0] + I.Izx*w[2] - I.Izz*w[0] -     I.h[1]*v0[2] - I.h[2]*v0[1] ,
          -I.Izx*w[1] + I.Izy*w[0] - I.Izx*w[1] + I.Izy*w[0] + 2.*(I.h[0]*v0[0] + I.h[1]*v0[1])
        )
    );


  const SpatialMatrix I_err = (I_ad.toMatrix() - I_an).cwiseAbs();
  const double error = I_err.maxCoeff();

  if (error > 1e-7) {
    cout << "S:" << S.transpose() << endl;
    cout << "w:" << w.transpose() << endl;
    cout << "v0:" << v0.transpose() << endl;

    cout << "I: " << endl << I.toMatrix() << endl;
    cout << "I_fd: " << endl << I_fd.toMatrix() << endl;
    cout << "I_an: " << endl << I_an << endl << endl;
    cout << "I_ad: " << endl << I_ad.toMatrix() << endl << endl;
    cout << "error(max): " << " (" << error << ")" << endl
      << endl << I_err << endl;
  }

  // test result
  CHECK_ARRAY_CLOSE(
    I_fd.toMatrix().data(), I_ad.toMatrix().data(),
    I_fd.toMatrix().size(), TEST_PRECISION
  );
}

TEST(Xrotx_efficient_analytic_derivative_of_spatial_inertia){
  srand (421337);
  SpatialVector S = SpatialVector::Zero();
  S(0) = 1.0;
  efficient_analytic_derivative_of_spatial_inertia(*this, Xrotx, S, 1e-7);
}

TEST(Xroty_efficient_analytic_derivative_of_spatial_inertia){
  srand (421337);
  SpatialVector S = SpatialVector::Zero();
  S(1) = 1.0;
  efficient_analytic_derivative_of_spatial_inertia(*this, Xroty, S, 1e-7);
}

TEST(Xrotz_efficient_analytic_derivative_of_spatial_inertia){
  srand (421337);
  SpatialVector S = SpatialVector::Zero();
  S(2) = 1.0;
  efficient_analytic_derivative_of_spatial_inertia(*this, Xrotz, S, 1e-7);
}

TEST(Xtransx_efficient_analytic_derivative_of_spatial_inertia){
  srand (421337);
  SpatialVector S = SpatialVector::Zero();
  S(3) = 1.0;
  std::function<SpatialTransform (const double &)> Xtransx
    = [](const double &q) { return Xtrans(Vector3d(q, 0., 0.)); };
  efficient_analytic_derivative_of_spatial_inertia(*this, Xtransx, S, 1e-8);
}

TEST(Xtransy_efficient_analytic_derivative_of_spatial_inertia){
  srand (421337);
  SpatialVector S = SpatialVector::Zero();
  S(4) = 1.0;
  std::function<SpatialTransform (const double &)> Xtransy
    = [](const double &q) { return Xtrans(Vector3d(0., q, 0.)); };
  efficient_analytic_derivative_of_spatial_inertia(*this, Xtransy, S, 1e-8);
}

TEST(Xtransz_efficient_analytic_derivative_of_spatial_inertia){
  srand (421337);
  SpatialVector S = SpatialVector::Zero();
  S(5) = 1.0;
  std::function<SpatialTransform (const double &)> Xtransz
    = [](const double &q) { return Xtrans(Vector3d(0., 0., q)); };
  efficient_analytic_derivative_of_spatial_inertia(*this, Xtransz, S, 1e-7);
}

// -----------------------------------------------------------------------------

