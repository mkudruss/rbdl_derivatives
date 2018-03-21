#ifndef RBDL_SPATIALALGEBRAOPERATORS_ED_H
#define RBDL_SPATIALALGEBRAOPERATORS_ED_H

#include <rbdl/rbdl.h>
#include <rbdl/SpatialAlgebraOperators.h>

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------
namespace Math {
// -----------------------------------------------------------------------------
namespace ED {
// -----------------------------------------------------------------------------

using namespace RigidBodyDynamics::Math;

typedef Eigen::Matrix<double, 6, Eigen::Dynamic> SpatialDirection;

/*
inline void X_apply_v_plus_w (
  SpatialVector &res,
  SpatialDirection &res_dir,
  const SpatialTransform &X, const std::vector<SpatialTransform> &X_dir,
  const SpatialVector &v, const SpatialDirection &v_dir,
  const SpatialVector &w
) {
  Vector3d v_rxw (
      v[3] - X.r[1]*v[2] + X.r[2]*v[1],
      v[4] - X.r[2]*v[0] + X.r[0]*v[2],
      v[5] - X.r[0]*v[1] + X.r[1]*v[0]
  );

  const unsigned int ndirs = res_dir.cols();
  for (unsigned idir = 0; idir < ndirs; idir++) {
    Vector3d v_rxw_dir (
        v_dir(3, idir) - X.r[1]*v_dir(2, idir) + X.r[2]*v_dir(1, idir),
        v_dir(4, idir) - X.r[2]*v_dir(0, idir) + X.r[0]*v_dir(2, idir),
        v_dir(5, idir) - X.r[0]*v_dir(1, idir) + X.r[1]*v_dir(0, idir)
    );

      res_dir(0, idir) = X_dir[idir].E(0,0) * v[0]     + X_dir[idir].E(0,1) * v[1]     + X_dir[idir].E(0,2) * v[2]
             + X.E(0,0)     * v_dir(0, idir) + X.E(0,1)     * v_dir(1, idir) + X.E(0,2)     * v_dir(2, idir);

      res_dir(1, idir) = X_dir[idir].E(1,0) * v[0]     + X_dir[idir].E(1,1) * v[1]     + X_dir[idir].E(1,2) * v[2]
             + X.E(1,0)     * v_dir(0, idir) + X.E(1,1)     * v_dir(1, idir) + X.E(1,2)     * v_dir(2, idir);

      res_dir(2, idir) = X_dir[idir].E(2,0) * v[0]     + X_dir[idir].E(2,1) * v[1]     + X_dir[idir].E(2,2) * v[2]
             + X.E(2,0)     * v_dir(0, idir) + X.E(2,1)     * v_dir(1, idir) + X.E(2,2)     * v_dir(2, idir);

      res_dir(3, idir) = X_dir[idir].E(0,0) * v_rxw[0]     + X_dir[idir].E(0,1) * v_rxw[1]     + X_dir[idir].E(0,2) * v_rxw[2]
             + X.E(0,0)     * v_rxw_dir[0] + X.E(0,1)     * v_rxw_dir[1] + X.E(0,2)     * v_rxw_dir[2];

      res_dir(4, idir) = X_dir[idir].E(1,0) * v_rxw[0]     + X_dir[idir].E(1,1) * v_rxw[1]     + X_dir[idir].E(1,2) * v_rxw[2]
             + X.E(1,0)     * v_rxw_dir[0] + X.E(1,1)     * v_rxw_dir[1] + X.E(1,2)     * v_rxw_dir[2];

      res_dir(5, idir) = X_dir[idir].E(2,0) * v_rxw[0]     + X_dir[idir].E(2,1) * v_rxw[1]     + X_dir[idir].E(2,2) * v_rxw[2]
             + X.E(2,0)     * v_rxw_dir[0] + X.E(2,1)     * v_rxw_dir[1] + X.E(2,2)     * v_rxw_dir[2];
    }

    // nominal evaluation
    res[0] = X.E(0,0) * v[0] + X.E(0,1) * v[1] + X.E(0,2) * v[2];
    res[1] = X.E(1,0) * v[0] + X.E(1,1) * v[1] + X.E(1,2) * v[2];
    res[2] = X.E(2,0) * v[0] + X.E(2,1) * v[1] + X.E(2,2) * v[2];
    res[3] = X.E(0,0) * v_rxw[0] + X.E(0,1) * v_rxw[1] + X.E(0,2) * v_rxw[2];
    res[4] = X.E(1,0) * v_rxw[0] + X.E(1,1) * v_rxw[1] + X.E(1,2) * v_rxw[2];
    res[5] = X.E(2,0) * v_rxw[0] + X.E(2,1) * v_rxw[1] + X.E(2,2) * v_rxw[2];
}
*/

/** Same as res = X_d * v + X * v_d.
 *
 * \returns (E * w, - E * rxw + E * v)
 */
inline void X_apply_v (
  SpatialVector &res,
  SpatialDirection &res_dir,
  const SpatialTransform &X, const std::vector<SpatialTransform> &X_dir,
  const SpatialVector &v, const SpatialDirection &v_dir
) {

  Vector3d v_rxw (
      v[3] - X.r[1] * v[2] + X.r[2] * v[1],
      v[4] - X.r[2] * v[0] + X.r[0] * v[2],
      v[5] - X.r[0] * v[1] + X.r[1] * v[0]
      );

  const unsigned int ndirs = res_dir.cols();
  for (unsigned idir = 0; idir < ndirs; idir++) {
    Vector3d w_dir = v_dir.col(idir).segment<3>(0);
    res_dir.col(idir).segment<3>(0) =
        X.E * w_dir
        + X_dir[idir].E * v.segment<3>(0);
    res_dir.col(idir).segment<3>(3) =
        X.E * (v_dir.col(idir).segment<3>(3)
                - X.r.cross(w_dir)
                - X_dir[idir].r.cross(v.segment<3>(0)))
        + X_dir[idir].E * v_rxw;
  }

  res = SpatialVector(
      X.E(0,0) * v[0] + X.E(0,1) * v[1] + X.E(0,2) * v[2],
      X.E(1,0) * v[0] + X.E(1,1) * v[1] + X.E(1,2) * v[2],
      X.E(2,0) * v[0] + X.E(2,1) * v[1] + X.E(2,2) * v[2],
      X.E(0,0) * v_rxw[0] + X.E(0,1) * v_rxw[1] + X.E(0,2) * v_rxw[2],
      X.E(1,0) * v_rxw[0] + X.E(1,1) * v_rxw[1] + X.E(1,2) * v_rxw[2],
      X.E(2,0) * v_rxw[0] + X.E(2,1) * v_rxw[1] + X.E(2,2) * v_rxw[2]
      );
}

/** Same as f += X^T * f.
 */
inline void inplace_X_applyT_f(
  SpatialVector &res,
  SpatialDirection &res_dir,
  const SpatialTransform &X, const std::vector<SpatialTransform> &X_dir,
  const SpatialVector &f_sp, const SpatialDirection &f_sp_dir
) {
  Vector3d E_T_f (
    X.E(0,0) * f_sp[3] + X.E(1,0) * f_sp[4] + X.E(2,0) * f_sp[5],
    X.E(0,1) * f_sp[3] + X.E(1,1) * f_sp[4] + X.E(2,1) * f_sp[5],
    X.E(0,2) * f_sp[3] + X.E(1,2) * f_sp[4] + X.E(2,2) * f_sp[5]
  );

  const unsigned int ndirs = res_dir.cols();
  for (unsigned idir = 0; idir < ndirs; idir++) {
    Vector3d E_T_f_dir (
      X_dir[idir].E(0,0) * f_sp[3] + X_dir[idir].E(1,0) * f_sp[4] + X_dir[idir].E(2,0) * f_sp[5] + X.E(0,0) * f_sp_dir(3, idir) + X.E(1,0) * f_sp_dir(4, idir) + X.E(2,0) * f_sp_dir(5, idir),
      X_dir[idir].E(0,1) * f_sp[3] + X_dir[idir].E(1,1) * f_sp[4] + X_dir[idir].E(2,1) * f_sp[5] + X.E(0,1) * f_sp_dir(3, idir) + X.E(1,1) * f_sp_dir(4, idir) + X.E(2,1) * f_sp_dir(5, idir),
      X_dir[idir].E(0,2) * f_sp[3] + X_dir[idir].E(1,2) * f_sp[4] + X_dir[idir].E(2,2) * f_sp[5] + X.E(0,2) * f_sp_dir(3, idir) + X.E(1,2) * f_sp_dir(4, idir) + X.E(2,2) * f_sp_dir(5, idir)
    );

    res_dir(0, idir) += X_dir[idir].E(0,0) * f_sp[0]     + X_dir[idir].E(1,0) * f_sp[1]     + X_dir[idir].E(2,0) * f_sp[2]     - X_dir[idir].r[2] * E_T_f[1]     + X_dir[idir].r[1] * E_T_f[2]
           + X.E(0,0)     * f_sp_dir(0, idir) + X.E(1,0)     * f_sp_dir(1, idir) + X.E(2,0)     * f_sp_dir(2, idir) - X.r[2]     * E_T_f_dir[1] + X.r[1]     * E_T_f_dir[2];
    res_dir(1, idir) += X_dir[idir].E(0,1) * f_sp[0]     + X_dir[idir].E(1,1) * f_sp[1]     + X_dir[idir].E(2,1) * f_sp[2]     + X_dir[idir].r[2] * E_T_f[0]     - X_dir[idir].r[0] * E_T_f[2]
           + X.E(0,1)     * f_sp_dir(0, idir) + X.E(1,1)     * f_sp_dir(1, idir) + X.E(2,1)     * f_sp_dir(2, idir) + X.r[2]     * E_T_f_dir[0] - X.r[0]     * E_T_f_dir[2];
    res_dir(2, idir) += X_dir[idir].E(0,2) * f_sp[0]     + X_dir[idir].E(1,2) * f_sp[1]     + X_dir[idir].E(2,2) * f_sp[2]     - X_dir[idir].r[1] * E_T_f[0]     + X_dir[idir].r[0] * E_T_f[1]
           + X.E(0,2)     * f_sp_dir(0, idir) + X.E(1,2)     * f_sp_dir(1, idir) + X.E(2,2)     * f_sp_dir(2, idir) - X.r[1]     * E_T_f_dir[0] + X.r[0]     * E_T_f_dir[1];
    res_dir(3, idir) += E_T_f_dir[0];
    res_dir(4, idir) += E_T_f_dir[1];
    res_dir(5, idir) += E_T_f_dir[2];
  }

  res[0] += X.E(0,0) * f_sp[0] + X.E(1,0) * f_sp[1] + X.E(2,0) * f_sp[2] - X.r[2] * E_T_f[1] + X.r[1] * E_T_f[2];
  res[1] += X.E(0,1) * f_sp[0] + X.E(1,1) * f_sp[1] + X.E(2,1) * f_sp[2] + X.r[2] * E_T_f[0] - X.r[0] * E_T_f[2];
  res[2] += X.E(0,2) * f_sp[0] + X.E(1,2) * f_sp[1] + X.E(2,2) * f_sp[2] - X.r[1] * E_T_f[0] + X.r[0] * E_T_f[1];
  res[3] += E_T_f [0];
  res[4] += E_T_f [1];
  res[5] += E_T_f [2];

  return;
}

/** Same as X^T * f.
 */
inline void X_applyT_f(
  SpatialVector &res,
  SpatialDirection &res_dir,
  const SpatialTransform &X, const std::vector<SpatialTransform> &X_dir,
  const SpatialVector &f_sp, const SpatialDirection &f_sp_dir
) {
  Vector3d E_T_f (
    X.E(0,0) * f_sp[3] + X.E(1,0) * f_sp[4] + X.E(2,0) * f_sp[5],
    X.E(0,1) * f_sp[3] + X.E(1,1) * f_sp[4] + X.E(2,1) * f_sp[5],
    X.E(0,2) * f_sp[3] + X.E(1,2) * f_sp[4] + X.E(2,2) * f_sp[5]
  );

  const unsigned int ndirs = res_dir.cols();
  for (unsigned idir = 0; idir < ndirs; idir++) {
    Vector3d E_T_f_dir (
      X_dir[idir].E(0,0) * f_sp[3] + X_dir[idir].E(1,0) * f_sp[4] + X_dir[idir].E(2,0) * f_sp[5] + X.E(0,0) * f_sp_dir(3, idir) + X.E(1,0) * f_sp_dir(4, idir) + X.E(2,0) * f_sp_dir(5, idir),
      X_dir[idir].E(0,1) * f_sp[3] + X_dir[idir].E(1,1) * f_sp[4] + X_dir[idir].E(2,1) * f_sp[5] + X.E(0,1) * f_sp_dir(3, idir) + X.E(1,1) * f_sp_dir(4, idir) + X.E(2,1) * f_sp_dir(5, idir),
      X_dir[idir].E(0,2) * f_sp[3] + X_dir[idir].E(1,2) * f_sp[4] + X_dir[idir].E(2,2) * f_sp[5] + X.E(0,2) * f_sp_dir(3, idir) + X.E(1,2) * f_sp_dir(4, idir) + X.E(2,2) * f_sp_dir(5, idir)
    );

    res_dir(0, idir) = X_dir[idir].E(0,0) * f_sp[0]     + X_dir[idir].E(1,0) * f_sp[1]     + X_dir[idir].E(2,0) * f_sp[2]     - X_dir[idir].r[2] * E_T_f[1]     + X_dir[idir].r[1] * E_T_f[2]
           + X.E(0,0)     * f_sp_dir(0, idir) + X.E(1,0)     * f_sp_dir(1, idir) + X.E(2,0)     * f_sp_dir(2, idir) - X.r[2]     * E_T_f_dir[1] + X.r[1]     * E_T_f_dir[2];
    res_dir(1, idir) = X_dir[idir].E(0,1) * f_sp[0]     + X_dir[idir].E(1,1) * f_sp[1]     + X_dir[idir].E(2,1) * f_sp[2]     + X_dir[idir].r[2] * E_T_f[0]     - X_dir[idir].r[0] * E_T_f[2]
           + X.E(0,1)     * f_sp_dir(0, idir) + X.E(1,1)     * f_sp_dir(1, idir) + X.E(2,1)     * f_sp_dir(2, idir) + X.r[2]     * E_T_f_dir[0] - X.r[0]     * E_T_f_dir[2];
    res_dir(2, idir) = X_dir[idir].E(0,2) * f_sp[0]     + X_dir[idir].E(1,2) * f_sp[1]     + X_dir[idir].E(2,2) * f_sp[2]     - X_dir[idir].r[1] * E_T_f[0]     + X_dir[idir].r[0] * E_T_f[1]
           + X.E(0,2)     * f_sp_dir(0, idir) + X.E(1,2)     * f_sp_dir(1, idir) + X.E(2,2)     * f_sp_dir(2, idir) - X.r[1]     * E_T_f_dir[0] + X.r[0]     * E_T_f_dir[1];
    res_dir(3, idir) = E_T_f_dir[0];
    res_dir(4, idir) = E_T_f_dir[1];
    res_dir(5, idir) = E_T_f_dir[2];
  }

  res[0] = X.E(0,0) * f_sp[0] + X.E(1,0) * f_sp[1] + X.E(2,0) * f_sp[2] - X.r[2] * E_T_f[1] + X.r[1] * E_T_f[2];
  res[1] = X.E(0,1) * f_sp[0] + X.E(1,1) * f_sp[1] + X.E(2,1) * f_sp[2] + X.r[2] * E_T_f[0] - X.r[0] * E_T_f[2];
  res[2] = X.E(0,2) * f_sp[0] + X.E(1,2) * f_sp[1] + X.E(2,2) * f_sp[2] - X.r[1] * E_T_f[0] + X.r[0] * E_T_f[1];
  res[3] = E_T_f [0];
  res[4] = E_T_f [1];
  res[5] = E_T_f [2];

  return;
}


/** Same as X^* I X^{-1}
*/
// void apply (
//   SpatialRigidBodyInertia & res,
//   const SpatialTransform &X,
//   const SpatialTransform &X_dir,
//   const SpatialRigidBodyInertia &rbi,
//   const SpatialRigidBodyInertia &rbi_dir
// ) {
//   std::cerr << __FILE__ << " " << __LINE__ << ":"
//        << " not yet implemented." << std::endl;
//   abort();
//   // return SpatialRigidBodyInertia (
//   //     rbi.m,
//   //     E * (rbi.h - rbi.m * r),
//   //     E *
//   //     (
//   //      Matrix3d (
//   //        rbi.Ixx, rbi.Iyx, rbi.Izx,
//   //        rbi.Iyx, rbi.Iyy, rbi.Izy,
//   //        rbi.Izx, rbi.Izy, rbi.Izz
//   //        )
//   //      + VectorCrossMatrix (r) * VectorCrossMatrix (rbi.h)
//   //      + (VectorCrossMatrix(rbi.h - rbi.m * r) * VectorCrossMatrix (r))
//   //     )
//   //     * E.transpose()
//   //     );
// }

/** Same as X^T I X
*/

inline void applyTranspose (
  SpatialRigidBodyInertia & res,
  std::vector<SpatialRigidBodyInertia> &res_dirs,
  const SpatialTransform &X,
  const std::vector<SpatialTransform> &X_dirs,
  const SpatialRigidBodyInertia &rbi,
  const std::vector<SpatialRigidBodyInertia> &rbi_dirs,
  const unsigned int &ndirs
) {

  // intermediate value
  Matrix3d I = rbi.getInertia();

  // nominal evaluation
  Vector3d E_T_mr = X.E.transpose() * rbi.h + rbi.m * X.r;

  // nominal evaluation
  Matrix3d Z
    = X.E.transpose() * I * X.E
    - VectorCrossMatrix(X.r) * VectorCrossMatrix (X.E.transpose() * rbi.h)
    - VectorCrossMatrix (E_T_mr) * VectorCrossMatrix (X.r);

  // derivative evaluation
  for (unsigned idir = 0; idir < ndirs; idir++) {
    // intermediate value
    Vector3d E_T_mr_dirs
      = X_dirs[idir].E.transpose() * rbi.h
      + X.E.transpose() * rbi_dirs[idir].h
      + rbi_dirs[idir].m * X.r
      + rbi.m * X_dirs[idir].r;

    // intermediate value
    Matrix3d Z_dirs
      = X_dirs[idir].E.transpose() * I * X.E
      + X.E.transpose() * rbi_dirs[idir].getInertia() * X.E
      + X.E.transpose() * I * X_dirs[idir].E

      - VectorCrossMatrix(X_dirs[idir].r) * VectorCrossMatrix (X.E.transpose() * rbi.h)
      - VectorCrossMatrix(X.r) * (
        VectorCrossMatrix (X_dirs[idir].E.transpose() * rbi.h)
        + VectorCrossMatrix (X.E.transpose() * rbi_dirs[idir].h)
      )

      - VectorCrossMatrix (E_T_mr_dirs) * VectorCrossMatrix (X.r)
      - VectorCrossMatrix (E_T_mr) * VectorCrossMatrix (X_dirs[idir].r)
    ;
    res_dirs[idir]
      = SpatialRigidBodyInertia (rbi_dirs[idir].m, E_T_mr_dirs, Z_dirs);
  }

  // nominal evaluation
  res = SpatialRigidBodyInertia (rbi.m, E_T_mr, Z);

  return;
}

inline void plusApplyTranspose (
  SpatialRigidBodyInertia & res,
  std::vector<SpatialRigidBodyInertia> &res_dirs,
  const SpatialTransform &X,
  const std::vector<SpatialTransform> &X_dirs,
  const SpatialRigidBodyInertia &rbi,
  const std::vector<SpatialRigidBodyInertia> &rbi_dirs,
  const unsigned int &ndirs
) {

  // intermediate value
  Matrix3d I = rbi.getInertia();

  // nominal evaluation
  Vector3d E_T_mr = X.E.transpose() * rbi.h + rbi.m * X.r;

  // nominal evaluation
  Matrix3d Z
    = X.E.transpose() * I * X.E
    - VectorCrossMatrix(X.r) * VectorCrossMatrix (X.E.transpose() * rbi.h)
    - VectorCrossMatrix (E_T_mr) * VectorCrossMatrix (X.r);

  // derivative evaluation
  for (unsigned idir = 0; idir < ndirs; idir++) {
    // intermediate value
    Vector3d E_T_mr_dirs
      = X_dirs[idir].E.transpose() * rbi.h
      + X.E.transpose() * rbi_dirs[idir].h
      + rbi_dirs[idir].m * X.r
      + rbi.m * X_dirs[idir].r;

    // intermediate value
    Matrix3d Z_dirs
      = X_dirs[idir].E.transpose() * I * X.E
      + X.E.transpose() * rbi_dirs[idir].getInertia() * X.E
      + X.E.transpose() * I * X_dirs[idir].E

      - VectorCrossMatrix(X_dirs[idir].r) * VectorCrossMatrix (X.E.transpose() * rbi.h)
      - VectorCrossMatrix(X.r) * (
        VectorCrossMatrix (X_dirs[idir].E.transpose() * rbi.h)
        + VectorCrossMatrix (X.E.transpose() * rbi_dirs[idir].h)
      )

      - VectorCrossMatrix (E_T_mr_dirs) * VectorCrossMatrix (X.r)
      - VectorCrossMatrix (E_T_mr) * VectorCrossMatrix (X_dirs[idir].r)
    ;
    res_dirs[idir] = res_dirs[idir]
      + SpatialRigidBodyInertia (rbi_dirs[idir].m, E_T_mr_dirs, Z_dirs);
  }

  // nominal evaluation
  res = res + SpatialRigidBodyInertia (rbi.m, E_T_mr, Z);

  return;
}


inline void Xrotx (
  SpatialTransform &res,
  std::vector<SpatialTransform> &res_dirs,
  const double &rot,
  const VectorNd &rot_dirs,
  const SpatialTransform &X,
  const unsigned int &ndirs
) {
  double s, c, s_t, c_t;
  // NOTE sin and cos are slow to compute?
  s = sin (rot);
  c = cos (rot);

  res.E = Matrix3d (
    1., 0., 0.,
    0., c, s,
    0., -s, c
  );
  res.r = X.r; // + X.E.transpose() * Vector3d (0., 0., 0.);

  for (unsigned idir = 0; idir < ndirs; idir++) {
    s_t = s*rot_dirs[idir];
    c_t = c*rot_dirs[idir];
    res_dirs[idir].E = Matrix3d (
        0.,   0.,   0.,
        0., -s_t,  c_t,
        0., -c_t, -s_t
    ) * X.E;
    res_dirs[idir].r = X.E * Vector3d (0., 0., 0.);
  }

  return;
}

inline void Xrotx (
  SpatialTransform &res,
  std::vector<SpatialTransform> &res_dirs,
  const double &rot,
  const VectorNd &rot_dirs,
  const unsigned int &ndirs
) {
  double s, c, s_t, c_t;
  // NOTE sin and cos are slow to compute?
  s = sin (rot);
  c = cos (rot);

  res.E = Matrix3d (
    1., 0., 0.,
    0., c, s,
    0., -s, c
  );
  res.r.setZero();

  for (unsigned idir = 0; idir < ndirs; idir++) {
    s_t = s*rot_dirs[idir];
    c_t = c*rot_dirs[idir];
    res_dirs[idir].E = Matrix3d (
        0.,   0.,   0.,
        0., -s_t,  c_t,
        0., -c_t, -s_t
    );
    res_dirs[idir].r.setZero();
  }

  return;
}

inline void Xroty (
  SpatialTransform &res,
  std::vector<SpatialTransform> &res_dirs,
  const double &rot,
  const VectorNd &rot_dirs,
  const unsigned int &ndirs
) {
  double s, c, s_t, c_t;
  // NOTE sin and cos are slow to compute?
  s = sin (rot);
  c = cos (rot);

  res.E = Matrix3d (
    c, 0., -s,
    0., 1., 0.,
    s, 0., c
  );
  res.r.setZero();

  for (unsigned idir = 0; idir < ndirs; idir++) {
    s_t = s*rot_dirs[idir];
    c_t = c*rot_dirs[idir];
    res_dirs[idir].E = Matrix3d (
     -s_t, 0., -c_t,
       0., 0.,   0.,
      c_t, 0., -s_t
    );
    res_dirs[idir].r.setZero();
  }

  return;
}

inline void Xroty (
  SpatialTransform &res,
  std::vector<SpatialTransform> &res_dirs,
  const double &rot,
  const VectorNd &rot_dirs,
  const SpatialTransform &X,
  const unsigned int &ndirs
) {
  double s, c, s_t, c_t;
  // NOTE sin and cos are slow to compute?
  s = sin (rot);
  c = cos (rot);

  res.E = Matrix3d (
    c, 0., -s,
    0., 1., 0.,
    s, 0., c
  );
  res.r = X.r; // + X.E.transpose() * Vector3d (0., 0., 0.);

  for (unsigned idir = 0; idir < ndirs; idir++) {
    s_t = s*rot_dirs[idir];
    c_t = c*rot_dirs[idir];
    res_dirs[idir].E = Matrix3d (
     -s_t, 0., -c_t,
       0., 0.,   0.,
      c_t, 0., -s_t
    ) * X.E;
    res_dirs[idir].r = X.E * Vector3d (0., 0., 0.);
  }

  return;
}


inline void Xrotz (
  SpatialTransform &res,
  std::vector<SpatialTransform> &res_dirs,
  const double &rot,
  const VectorNd &rot_dirs,
  const unsigned int &ndirs
) {
  double s, c, s_t, c_t;
  // NOTE sin and cos are slow to compute?
  s = sin (rot);
  c = cos (rot);

  res.E = Matrix3d (
    c, s, 0.,
    -s, c, 0.,
    0., 0., 1.
  );
  res.r.setZero();

  for (unsigned idir = 0; idir < ndirs; idir++) {
    s_t = s*rot_dirs[idir];
    c_t = c*rot_dirs[idir];
    res_dirs[idir].E = Matrix3d (
      -s_t,  c_t, 0.,
      -c_t, -s_t, 0.,
        0.,   0., 0.
    );
    res_dirs[idir].r.setZero();
  }

  return;
}

inline void Xrotz (
  SpatialTransform &res,
  std::vector<SpatialTransform> &res_dirs,
  const double &rot,
  const VectorNd &rot_dirs,
  const SpatialTransform &X,
  const unsigned int &ndirs
) {
  double s, c, s_t, c_t;
  // NOTE sin and cos are slow to compute?
  s = sin (rot);
  c = cos (rot);

  res.E = Matrix3d (
    c, s, 0.,
    -s, c, 0.,
    0., 0., 1.
  );
  res.r = X.r; // + X.E.transpose() * Vector3d (0., 0., 0.);

  for (unsigned idir = 0; idir < ndirs; idir++) {
    s_t = s*rot_dirs[idir];
    c_t = c*rot_dirs[idir];
    res_dirs[idir].E = Matrix3d (
      -s_t,  c_t, 0.,
      -c_t, -s_t, 0.,
        0.,   0., 0.
    ) * X.E;
    res_dirs[idir].r = X.E * Vector3d (0., 0., 0.);
  }

  return;
}

inline SpatialVector crossm (
    const SpatialVector &v1, const SpatialVector &v1_dirs,
    const SpatialVector &v2, const SpatialVector &v2_dirs
    ) {
  return SpatialVector (
      // nominal -v1[2] * v2[1] + v1[1] * v2[2],
      - v1_dirs[2] * v2[1] - v1[2] * v2_dirs[1]
      + v1_dirs[1] * v2[2] + v1[1] * v2_dirs[2],
      // nominal v1[2] * v2[0] - v1[0] * v2[2],
      + v1_dirs[2] * v2[0] + v1[2] * v2_dirs[0]
      - v1_dirs[0] * v2[2] - v1[0] * v2_dirs[2],
      // nominal -v1[1] * v2[0] + v1[0] * v2[1],
      - v1_dirs[1] * v2[0] - v1[1] * v2_dirs[0]
      + v1_dirs[0] * v2[1] + v1[0] * v2_dirs[1],
      // nominal -v1[5] * v2[1] + v1[4] * v2[2] - v1[2] * v2[4] + v1[1] * v2[5],
      - v1_dirs[5] * v2[1] - v1[5] * v2_dirs[1]
      + v1_dirs[4] * v2[2] + v1[4] * v2_dirs[2]
      - v1_dirs[2] * v2[4] - v1[2] * v2_dirs[4]
      + v1_dirs[1] * v2[5] + v1[1] * v2_dirs[5],
      // nominal v1[5] * v2[0] - v1[3] * v2[2] + v1[2] * v2[3] - v1[0] * v2[5],
      + v1_dirs[5] * v2[0] + v1[5] * v2_dirs[0]
      - v1_dirs[3] * v2[2] - v1[3] * v2_dirs[2]
      + v1_dirs[2] * v2[3] + v1[2] * v2_dirs[3]
      - v1_dirs[0] * v2[5] - v1[0] * v2_dirs[5],
      // nominal -v1[4] * v2[0] + v1[3] * v2[1] - v1[1] * v2[3] + v1[0] * v2[4]
      - v1_dirs[4] * v2[0] - v1[4] * v2_dirs[0]
      + v1_dirs[3] * v2[1] + v1[3] * v2_dirs[1]
      - v1_dirs[1] * v2[3] - v1[1] * v2_dirs[3]
      + v1_dirs[0] * v2[4] + v1[0] * v2_dirs[4]
      );
}

inline SpatialVector crossf (
  const SpatialVector &v1, const SpatialVector &v1_dirs,
  const SpatialVector &v2, const SpatialVector &v2_dirs
) {
  return SpatialVector (
      // nominal -v1[2] * v2[1] + v1[1] * v2[2] - v1[5] * v2[4] + v1[4] * v2[5],
      - v1_dirs[2] * v2[1] - v1[2] * v2_dirs[1]
      + v1_dirs[1] * v2[2] + v1[1] * v2_dirs[2]
      - v1_dirs[5] * v2[4] - v1[5] * v2_dirs[4]
      + v1_dirs[4] * v2[5] + v1[4] * v2_dirs[5],
      // nominal v1[2] * v2[0] - v1[0] * v2[2] + v1[5] * v2[3] - v1[3] * v2[5],
      + v1_dirs[2] * v2[0] + v1[2] * v2_dirs[0]
      - v1_dirs[0] * v2[2] - v1[0] * v2_dirs[2]
      + v1_dirs[5] * v2[3] + v1[5] * v2_dirs[3]
      - v1_dirs[3] * v2[5] - v1[3] * v2_dirs[5],
      // nominal -v1[1] * v2[0] + v1[0] * v2[1] - v1[4] * v2[3] + v1[3] * v2[4],
      - v1_dirs[1] * v2[0] - v1[1] * v2_dirs[0]
      + v1_dirs[0] * v2[1] + v1[0] * v2_dirs[1]
      - v1_dirs[4] * v2[3] - v1[4] * v2_dirs[3]
      + v1_dirs[3] * v2[4] + v1[3] * v2_dirs[4],
      // nominal - v1[2] * v2[4] + v1[1] * v2[5],
      - v1_dirs[2] * v2[4] - v1[2] * v2_dirs[4]
      + v1_dirs[1] * v2[5] + v1[1] * v2_dirs[5],
      // nominal + v1[2] * v2[3] - v1[0] * v2[5],
      + v1_dirs[2] * v2[3] + v1[2] * v2_dirs[3]
      - v1_dirs[0] * v2[5] - v1[0] * v2_dirs[5],
      // nominal - v1[1] * v2[3] + v1[0] * v2[4]
      - v1_dirs[1] * v2[3] - v1[1] * v2_dirs[3]
      + v1_dirs[0] * v2[4] + v1[0] * v2_dirs[4]
      );
}


/*
inline SpatialMatrix crossm (const SpatialVector &v) {
  return SpatialMatrix (
      0,  -v[2],  v[1],         0,          0,         0,
      v[2],          0, -v[0],         0,          0,         0,
      -v[1],   v[0],         0,         0,          0,         0,
      0,  -v[5],  v[4],         0,  -v[2],  v[1],
      v[5],          0, -v[3],  v[2],          0, -v[0],
      -v[4],   v[3],         0, -v[1],   v[0],         0
      );
}

inline SpatialVector crossm (const SpatialVector &v1, const SpatialVector &v2) {
  return SpatialVector (
      -v1[2] * v2[1] + v1[1] * v2[2],
      v1[2] * v2[0] - v1[0] * v2[2],
      -v1[1] * v2[0] + v1[0] * v2[1],
      -v1[5] * v2[1] + v1[4] * v2[2] - v1[2] * v2[4] + v1[1] * v2[5],
      v1[5] * v2[0] - v1[3] * v2[2] + v1[2] * v2[3] - v1[0] * v2[5],
      -v1[4] * v2[0] + v1[3] * v2[1] - v1[1] * v2[3] + v1[0] * v2[4]
      );
}

inline SpatialMatrix crossf (const SpatialVector &v) {
  return SpatialMatrix (
      0,  -v[2],  v[1],         0,  -v[5],  v[4],
      v[2],          0, -v[0],  v[5],          0, -v[3],
      -v[1],   v[0],         0, -v[4],   v[3],         0,
      0,          0,         0,         0,  -v[2],  v[1],
      0,          0,         0,  v[2],          0, -v[0],
      0,          0,         0, -v[1],   v[0],         0
      );
}
*/

// -----------------------------------------------------------------------------
} // namespace ED
// -----------------------------------------------------------------------------
} // namespace Math
// -----------------------------------------------------------------------------
} // namespace RigidBodyDynamics
// -----------------------------------------------------------------------------

#endif
