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

/**
 *
 */

// enum SpatialTransformType {
//   SpatialTransformType_XEr = 0,
//   SpatialTransformType_X0r,
//   SpatialTransformType_XE0,
//   SpatialTransformType_X00
// };

struct RBDL_DLLAPI SpatialTransformDot {
  SpatialTransformDot() :
    E (Matrix3d::Identity(3,3)),
    Erx (Matrix3d::Identity(3,3))
  {}
  SpatialTransformDot (const Matrix3d &E_, const Matrix3d &Erx_) :
    E (E_),
    Erx (Erx_)
  {}

  void setZero() {
    E.setZero();
    Erx.setZero();
  }

  static SpatialTransformDot Zero() {
    SpatialTransformDot res;
    res.setZero();
    return res;
  }

  Matrix3d E;
  Matrix3d Erx;
};

  /** Same as X * v.
   *
   * \returns (E * w, - E * rxw + E * v)
   */
  // SpatialVector apply (const SpatialVector &v_sp) const {
  //   Vector3d v_rxw (
  //       v_sp[3] - r[1]*v_sp[2] + r[2]*v_sp[1],
  //       v_sp[4] - r[2]*v_sp[0] + r[0]*v_sp[2],
  //       v_sp[5] - r[0]*v_sp[1] + r[1]*v_sp[0]
  //       );
  //   return SpatialVector (
  //       E(0,0) * v_sp[0] + E(0,1) * v_sp[1] + E(0,2) * v_sp[2],
  //       E(1,0) * v_sp[0] + E(1,1) * v_sp[1] + E(1,2) * v_sp[2],
  //       E(2,0) * v_sp[0] + E(2,1) * v_sp[1] + E(2,2) * v_sp[2],
  //       E(0,0) * v_rxw[0] + E(0,1) * v_rxw[1] + E(0,2) * v_rxw[2],
  //       E(1,0) * v_rxw[0] + E(1,1) * v_rxw[1] + E(1,2) * v_rxw[2],
  //       E(2,0) * v_rxw[0] + E(2,1) * v_rxw[1] + E(2,2) * v_rxw[2]
  //       );
  // }

  /** Same as X^T * f.
   *
   * \returns (E^T * n + rx * E^T * f, E^T * f)
   */
  // SpatialVector applyTranspose (const SpatialVector &f_sp) const {
  //   Vector3d E_T_f (
  //       E(0,0) * f_sp[3] + E(1,0) * f_sp[4] + E(2,0) * f_sp[5],
  //       E(0,1) * f_sp[3] + E(1,1) * f_sp[4] + E(2,1) * f_sp[5],
  //       E(0,2) * f_sp[3] + E(1,2) * f_sp[4] + E(2,2) * f_sp[5]
  //       );

  //   return SpatialVector (
  //       E(0,0) * f_sp[0] + E(1,0) * f_sp[1] + E(2,0) * f_sp[2] - r[2] * E_T_f[1] + r[1] * E_T_f[2],
  //       E(0,1) * f_sp[0] + E(1,1) * f_sp[1] + E(2,1) * f_sp[2] + r[2] * E_T_f[0] - r[0] * E_T_f[2],
  //       E(0,2) * f_sp[0] + E(1,2) * f_sp[1] + E(2,2) * f_sp[2] - r[1] * E_T_f[0] + r[0] * E_T_f[1],
  //       E_T_f [0],
  //       E_T_f [1],
  //       E_T_f [2]
  //       );
  // }

  /** Same as X^* I X^{-1}
  */
  // SpatialRigidBodyInertia apply (const SpatialRigidBodyInertia &rbi) {
  //   return SpatialRigidBodyInertia (
  //       rbi.m,
  //       E * (rbi.h - rbi.m * r),
  //       E *
  //       (
  //        Matrix3d (
  //          rbi.Ixx, rbi.Iyx, rbi.Izx,
  //          rbi.Iyx, rbi.Iyy, rbi.Izy,
  //          rbi.Izx, rbi.Izy, rbi.Izz
  //          )
  //        + VectorCrossMatrix (r) * VectorCrossMatrix (rbi.h)
  //        + (VectorCrossMatrix(rbi.h - rbi.m * r) * VectorCrossMatrix (r))
  //       )
  //       * E.transpose()
  //       );
  // }

  /** Same as X^T I X
  */
  // SpatialRigidBodyInertia applyTranspose (const SpatialRigidBodyInertia &rbi) const {
  //   Vector3d E_T_mr = E.transpose() * rbi.h + rbi.m * r;
  //   return SpatialRigidBodyInertia (
  //       rbi.m,
  //       E_T_mr,
  //       E.transpose() *
  //       Matrix3d (
  //         rbi.Ixx, rbi.Iyx, rbi.Izx,
  //         rbi.Iyx, rbi.Iyy, rbi.Izy,
  //         rbi.Izx, rbi.Izy, rbi.Izz
  //         ) * E
  //       - VectorCrossMatrix(r) * VectorCrossMatrix (E.transpose() * rbi.h)
  //       - VectorCrossMatrix (E_T_mr) * VectorCrossMatrix (r));
  // }

  // SpatialVector applyAdjoint (const SpatialVector &f_sp) {
  //   Vector3d En_rxf = E * (Vector3d (f_sp[0], f_sp[1], f_sp[2]) - r.cross(Vector3d (f_sp[3], f_sp[4], f_sp[5])));
  //   //    Vector3d En_rxf = E * (Vector3d (f_sp[0], f_sp[1], f_sp[2]) - r.cross(Eigen::Map<Vector3d> (&(f_sp[3]))));

  //   return SpatialVector (
  //       En_rxf[0],
  //       En_rxf[1],
  //       En_rxf[2],
  //       E(0,0) * f_sp[3] + E(0,1) * f_sp[4] + E(0,2) * f_sp[5],
  //       E(1,0) * f_sp[3] + E(1,1) * f_sp[4] + E(1,2) * f_sp[5],
  //       E(2,0) * f_sp[3] + E(2,1) * f_sp[4] + E(2,2) * f_sp[5]
  //       );
  // }

  // SpatialMatrix toMatrix () const {
  //   Matrix3d _Erx =
  //     E * Matrix3d (
  //         0., -r[2], r[1],
  //         r[2], 0., -r[0],
  //         -r[1], r[0], 0.
  //         );
  //   SpatialMatrix result;
  //   result.block<3,3>(0,0) = E;
  //   result.block<3,3>(0,3) = Matrix3d::Zero(3,3);
  //   result.block<3,3>(3,0) = -_Erx;
  //   result.block<3,3>(3,3) = E;

  //   return result;
  // }

  // SpatialMatrix toMatrixAdjoint () const {
  //   Matrix3d _Erx =
  //     E * Matrix3d (
  //         0., -r[2], r[1],
  //         r[2], 0., -r[0],
  //         -r[1], r[0], 0.
  //         );
  //   SpatialMatrix result;
  //   result.block<3,3>(0,0) = E;
  //   result.block<3,3>(0,3) = -_Erx;
  //   result.block<3,3>(3,0) = Matrix3d::Zero(3,3);
  //   result.block<3,3>(3,3) = E;

  //   return result;
  // }

  // SpatialMatrix toMatrixTranspose () const {
  //   Matrix3d _Erx =
  //     E * Matrix3d (
  //         0., -r[2], r[1],
  //         r[2], 0., -r[0],
  //         -r[1], r[0], 0.
  //         );
  //   SpatialMatrix result;
  //   result.block<3,3>(0,0) = E.transpose();
  //   result.block<3,3>(0,3) = -_Erx.transpose();
  //   result.block<3,3>(3,0) = Matrix3d::Zero(3,3);
  //   result.block<3,3>(3,3) = E.transpose();

  //   return result;
  // }

  // SpatialTransform inverse() const {
  //   return SpatialTransform (
  //       E.transpose(),
  //       - E * r
  //       );
  // }

  // SpatialTransform operator* (const SpatialTransform &XT) const {
  //   return SpatialTransform (E * XT.E, XT.r + XT.E.transpose() * r);
  // }

  // void operator*= (const SpatialTransform &XT) {
  //   r = XT.r + XT.E.transpose() * r;
  //   E *= XT.E;
  // }

  // void setZero() {
  //   E.setZero();
  //   r.setZero();
  // }

  // static SpatialTransform Zero() {
  //   SpatialTransform res;
  //   res.setZero();
  //   return res;
  // }

  // Matrix3d E;
  // Vector3d r;
// };

/** Same as res = X_d * v + X * v_d.
 *
 * \returns (E * w, - E * rxw + E * v)
 */
inline void apply (
  SpatialVector & res,
  const SpatialTransform &X, const SpatialTransform &X_dir,
  const SpatialVector &v, const SpatialVector &v_dir
) {
    Vector3d v_rxw (
        v[3] - X_dir.r[1]*v[2] + X_dir.r[2]*v[1],
        v[4] - X_dir.r[2]*v[0] + X_dir.r[0]*v[2],
        v[5] - X_dir.r[0]*v[1] + X_dir.r[1]*v[0]
    );
    Vector3d v_rxw_dir (
        v_dir[3] - X.r[1]*v_dir[2] + X.r[2]*v_dir[1],
        v_dir[4] - X.r[2]*v_dir[0] + X.r[0]*v_dir[2],
        v_dir[5] - X.r[0]*v_dir[1] + X.r[1]*v_dir[0]
    );

    // Assign results
    res[0] = X_dir.E(0,0) * v[0]     + X_dir.E(0,1) * v[1]     + X_dir.E(0,2) * v[2]
           + X.E(0,0)     * v_dir[0] + X.E(0,1)     * v_dir[1] + X.E(0,2)     * v_dir[2];

    res[1] = X_dir.E(1,0) * v[0]     + X_dir.E(1,1) * v[1]     + X_dir.E(1,2) * v[2]
           + X.E(1,0)     * v_dir[0] + X.E(1,1)     * v_dir[1] + X.E(1,2)     * v_dir[2];

    res[2] = X_dir.E(2,0) * v[0]     + X_dir.E(2,1) * v[1]     + X_dir.E(2,2) * v[2]
           + X.E(2,0)     * v_dir[0] + X.E(2,1)     * v_dir[1] + X.E(2,2)     * v_dir[2];

    res[3] = X_dir.E(0,0) * v_rxw[0]     + X_dir.E(0,1) * v_rxw[1]     + X_dir.E(0,2) * v_rxw[2]
           + X.E(0,0)     * v_rxw_dir[0] + X.E(0,1)     * v_rxw_dir[1] + X.E(0,2)     * v_rxw_dir[2];

    res[4] = X_dir.E(1,0) * v_rxw[0]     + X_dir.E(1,1) * v_rxw[1]     + X_dir.E(1,2) * v_rxw[2]
           + X.E(1,0)     * v_rxw_dir[0] + X.E(1,1)     * v_rxw_dir[1] + X.E(1,2)     * v_rxw_dir[2];

    res[5] = X_dir.E(2,0) * v_rxw[0]     + X_dir.E(2,1) * v_rxw[1]     + X_dir.E(2,2) * v_rxw[2]
           + X.E(2,0)     * v_rxw_dir[0] + X.E(2,1)     * v_rxw_dir[1] + X.E(2,2)     * v_rxw_dir[2];
}

/** Same as X^T * f.
 *
 * \returns (E^T * n + rx * E^T * f, E^T * f)
 */
inline void applyTranspose (
  SpatialVector & res,
  const SpatialTransform &X, const SpatialTransform &X_dir,
  const SpatialVector &f_sp, const SpatialVector &f_sp_dir
) {
  Vector3d E_T_f (
    X.E(0,0) * f_sp[3] + X.E(1,0) * f_sp[4] + X.E(2,0) * f_sp[5],
    X.E(0,1) * f_sp[3] + X.E(1,1) * f_sp[4] + X.E(2,1) * f_sp[5],
    X.E(0,2) * f_sp[3] + X.E(1,2) * f_sp[4] + X.E(2,2) * f_sp[5]
  );

  Vector3d E_T_f_dir (
    X_dir.E(0,0) * f_sp[3] + X_dir.E(1,0) * f_sp[4] + X_dir.E(2,0) * f_sp[5] + X.E(0,0) * f_sp_dir[3] + X.E(1,0) * f_sp_dir[4] + X.E(2,0) * f_sp_dir[5],
    X_dir.E(0,1) * f_sp[3] + X_dir.E(1,1) * f_sp[4] + X_dir.E(2,1) * f_sp[5] + X.E(0,1) * f_sp_dir[3] + X.E(1,1) * f_sp_dir[4] + X.E(2,1) * f_sp_dir[5],
    X_dir.E(0,2) * f_sp[3] + X_dir.E(1,2) * f_sp[4] + X_dir.E(2,2) * f_sp[5] + X.E(0,2) * f_sp_dir[3] + X.E(1,2) * f_sp_dir[4] + X.E(2,2) * f_sp_dir[5]
  );

  res[0] = X_dir.E(0,0) * f_sp[0]     + X_dir.E(1,0) * f_sp[1]     + X_dir.E(2,0) * f_sp[2]     - X_dir.r[2] * E_T_f[1]     + X_dir.r[1] * E_T_f[2]
         + X.E(0,0)     * f_sp_dir[0] + X.E(1,0)     * f_sp_dir[1] + X.E(2,0)     * f_sp_dir[2] - X.r[2]     * E_T_f_dir[1] + X.r[1]     * E_T_f_dir[2];
  res[1] = X_dir.E(0,1) * f_sp[0]     + X_dir.E(1,1) * f_sp[1]     + X_dir.E(2,1) * f_sp[2]     + X_dir.r[2] * E_T_f[0]     - X_dir.r[0] * E_T_f[2]
         + X.E(0,1)     * f_sp_dir[0] + X.E(1,1)     * f_sp_dir[1] + X.E(2,1)     * f_sp_dir[2] + X.r[2]     * E_T_f_dir[0] - X.r[0]     * E_T_f_dir[2];
  res[2] = X_dir.E(0,2) * f_sp[0]     + X_dir.E(1,2) * f_sp[1]     + X_dir.E(2,2) * f_sp[2]     - X_dir.r[1] * E_T_f[0]     + X_dir.r[0] * E_T_f[1]
         + X.E(0,2)     * f_sp_dir[0] + X.E(1,2)     * f_sp_dir[1] + X.E(2,2)     * f_sp_dir[2] - X.r[1]     * E_T_f_dir[0] + X.r[0]     * E_T_f_dir[1];
  res[3] = E_T_f_dir [0];
  res[4] = E_T_f_dir [1];
  res[5] = E_T_f_dir [2];
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


inline SpatialTransform Xrotx (
    const double &xrot,
    const double &xrot_dirs
) {
  double s, c;
  s = sin (xrot)*xrot_dirs;
  c = cos (xrot)*xrot_dirs;
  return SpatialTransform (
      Matrix3d (
        0., 0., 0.,
        0., -s,  c,
        0., -c, -s
        ),
      Vector3d (0., 0., 0.)
      );
}

inline SpatialTransform Xroty (
  const double &yrot,
  const double &yrot_dirs
) {
  double s, c;
  s = sin (yrot)*yrot_dirs;
  c = cos (yrot)*yrot_dirs;
  return SpatialTransform (
      Matrix3d (
        -s, 0., -c,
        0., 0., 0.,
         c, 0., -s
        ),
      Vector3d (0., 0., 0.)
      );
}

inline void Xrotz (
  SpatialTransform &res,
  std::vector<SpatialTransformDot> &res_dirs,
  const double &zrot,
  const VectorNd &zrot_dirs,
  const unsigned int &ndirs
) {
  double s, c, s_t, c_t;
  // NOTE sin and cos are slow to compute?
  s = sin (zrot);
  c = cos (zrot);

  res.E = Matrix3d (
    c, s, 0.,
    -s, c, 0.,
    0., 0., 1.
  );
  res.r = Vector3d (0., 0., 0.);

  for (unsigned idir = 0; idir < ndirs; idir++) {
    s_t = s*zrot_dirs[idir];
    c_t = c*zrot_dirs[idir];
    res_dirs[idir].E = Matrix3d (
      -s_t,  c_t, 0.,
      -c_t, -s_t, 0.,
        0.,   0., 0.
    );
    res_dirs[idir].Erx.setZero();
  }

  return;
}

inline void Xrotz (
  SpatialTransform &res,
  std::vector<SpatialTransform> &res_dirs,
  const double &zrot,
  const VectorNd &zrot_dirs,
  const SpatialTransform &X,
  const unsigned int &ndirs
) {
  double s, c, s_t, c_t;
  // NOTE sin and cos are slow to compute?
  s = sin (zrot);
  c = cos (zrot);

  res.E = Matrix3d (
    c, s, 0.,
    -s, c, 0.,
    0., 0., 1.
  );
  res.r = X.r; // + X.E.transpose() * Vector3d (0., 0., 0.);

  for (unsigned idir = 0; idir < ndirs; idir++) {
    s_t = s*zrot_dirs[idir];
    c_t = c*zrot_dirs[idir];
    res_dirs[idir].E = Matrix3d (
      -s_t,  c_t, 0.,
      -c_t, -s_t, 0.,
        0.,   0., 0.
    ) * X.E;
    res_dirs[idir].r = X.E * Vector3d (0., 0., 0.);
  }

  return;
}

inline SpatialTransform Xrotz (
  const double &zrot,
  const double &zrot_dirs
) {
  double s, c;
  s = sin (zrot)*zrot_dirs;
  c = cos (zrot)*zrot_dirs;
  return SpatialTransform (
      Matrix3d (
        -s,  c, 0.,
        -c, -s, 0.,
        0., 0., 0.
        ),
      Vector3d (0., 0., 0.)
      );
}

inline SpatialTransform Xtrans (
  const Vector3d &r,
  const Vector3d &r_dirs
) {
  return SpatialTransform (
      Matrix3d::Zero(3,3),
      r_dirs
      );
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
