#ifndef RBDL_SPATIALALGEBRAOPERATORS_AD_H
#define RBDL_SPATIALALGEBRAOPERATORS_AD_H

#include <rbdl/rbdl.h>
#include <rbdl/SpatialAlgebraOperators.h>

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------
namespace AD {
// -----------------------------------------------------------------------------

using namespace RigidBodyDynamics::Math;

inline void mulSTST(
    SpatialTransform const & lhs,
    std::vector<SpatialTransform> const & lhs_dirs,
    SpatialTransform const & rhs,
    std::vector<SpatialTransform> const & rhs_dirs,
    SpatialTransform & res,
    std::vector<SpatialTransform> & res_dirs) {
  unsigned const ndirs = lhs_dirs.size();

  assert(ndirs == rhs_dirs.size());
  assert(ndirs == res_dirs.size());

  // derivative code
  for (unsigned idir = 0; idir < ndirs; idir++) {
    res_dirs[idir].E = lhs.E * rhs_dirs[idir].E + lhs_dirs[idir].E * rhs.E;
    res_dirs[idir].r =
        rhs_dirs[idir].r
        + rhs.E.transpose() * lhs_dirs[idir].r
        + rhs_dirs[idir].E.transpose() * lhs.r;
  }

  // nominal code
  res.E = lhs.E * rhs.E;
  res.r = rhs.r + rhs.E.transpose() * lhs.r;
}

inline void mulSTST(
    SpatialTransform const & lhs,
    std::vector<SpatialTransform> const & lhs_dirs,
    SpatialTransform const & rhs,
    SpatialTransform & res,
    std::vector<SpatialTransform> & res_dirs) {
  unsigned const ndirs = lhs_dirs.size();

  assert(ndirs == res_dirs.size());

  // derivative code
  for (unsigned idir = 0; idir < ndirs; idir++) {
    res_dirs[idir].E = lhs_dirs[idir].E * rhs.E;
    res_dirs[idir].r = rhs.E.transpose() * lhs_dirs[idir].r;
  }

  // nominal code
  res.E = lhs.E * rhs.E;
  res.r = rhs.r + rhs.E.transpose() * lhs.r;
}


inline void applySTSV(
    unsigned ndirs,
    SpatialTransform const & st,
    std::vector<SpatialTransform> const & st_dirs,
    SpatialVector const & sv,
    SpatialVector & res,
    std::vector<SpatialVector> & res_dirs) {
  assert(ndirs <= st_dirs.size());
  assert(ndirs <= res_dirs.size());

  Vector3d v_rxw (
      sv[3] - st.r[1] * sv[2] + st.r[2] * sv[1],
      sv[4] - st.r[2] * sv[0] + st.r[0] * sv[2],
      sv[5] - st.r[0] * sv[1] + st.r[1] * sv[0]
      );

  for (unsigned idir = 0; idir < ndirs; idir++) {
    res_dirs[idir].segment<3>(0) =
        st_dirs[idir].E * sv.segment<3>(0);
    res_dirs[idir].segment<3>(3) =
        st.E * (- st_dirs[idir].r.cross(sv.segment<3>(0)))
        + st_dirs[idir].E * v_rxw;
  }

  res = SpatialVector(
      st.E(0,0) * sv[0] + st.E(0,1) * sv[1] + st.E(0,2) * sv[2],
      st.E(1,0) * sv[0] + st.E(1,1) * sv[1] + st.E(1,2) * sv[2],
      st.E(2,0) * sv[0] + st.E(2,1) * sv[1] + st.E(2,2) * sv[2],
      st.E(0,0) * v_rxw[0] + st.E(0,1) * v_rxw[1] + st.E(0,2) * v_rxw[2],
      st.E(1,0) * v_rxw[0] + st.E(1,1) * v_rxw[1] + st.E(1,2) * v_rxw[2],
      st.E(2,0) * v_rxw[0] + st.E(2,1) * v_rxw[1] + st.E(2,2) * v_rxw[2]
      );
}

inline void applySTSV(
    unsigned ndirs,
    SpatialTransform const & st,
    std::vector<SpatialTransform> const & st_dirs,
    SpatialVector const & sv,
    std::vector<SpatialVector> const & sv_dirs,
    SpatialVector & res,
    std::vector<SpatialVector> & res_dirs) {
  assert(ndirs <= st_dirs.size());
  assert(ndirs <= sv_dirs.size());
  assert(ndirs <= res_dirs.size());

  Vector3d v_rxw (
      sv[3] - st.r[1] * sv[2] + st.r[2] * sv[1],
      sv[4] - st.r[2] * sv[0] + st.r[0] * sv[2],
      sv[5] - st.r[0] * sv[1] + st.r[1] * sv[0]
      );

  for (unsigned idir = 0; idir < ndirs; idir++) {
    res_dirs[idir].segment<3>(0) =
        st.E * sv_dirs[idir].segment<3>(0)
        + st_dirs[idir].E * sv.segment<3>(0);
    res_dirs[idir].segment<3>(3) =
        st.E * (sv_dirs[idir].segment<3>(3)
                - st.r.cross(sv_dirs[idir].segment<3>(0))
                - st_dirs[idir].r.cross(sv.segment<3>(0)))
        + st_dirs[idir].E * v_rxw;
  }

  res = SpatialVector(
      st.E(0,0) * sv[0] + st.E(0,1) * sv[1] + st.E(0,2) * sv[2],
      st.E(1,0) * sv[0] + st.E(1,1) * sv[1] + st.E(1,2) * sv[2],
      st.E(2,0) * sv[0] + st.E(2,1) * sv[1] + st.E(2,2) * sv[2],
      st.E(0,0) * v_rxw[0] + st.E(0,1) * v_rxw[1] + st.E(0,2) * v_rxw[2],
      st.E(1,0) * v_rxw[0] + st.E(1,1) * v_rxw[1] + st.E(1,2) * v_rxw[2],
      st.E(2,0) * v_rxw[0] + st.E(2,1) * v_rxw[1] + st.E(2,2) * v_rxw[2]
      );
}

inline void applyTransposeAD (
    unsigned ndirs,
    SpatialTransform const & st,
    std::vector<SpatialTransform> const & st_dirs,
    SpatialRigidBodyInertia const & srbi,
    std::vector<SpatialRigidBodyInertia> const & srbi_dirs,
    SpatialRigidBodyInertia & res,
    std::vector<SpatialRigidBodyInertia> & res_dirs) {

  Math::MatrixNd E_T_mr_dirs(3, ndirs);
  Math::Vector3d E_T_mr;

  // derivative code
  for (unsigned i = 0; i < ndirs; i++) {
    E_T_mr_dirs.col(i) =
        st_dirs[i].E.transpose() * srbi.h + st.E.transpose() * srbi_dirs[i].h
        + srbi.m * st_dirs[i].r + srbi_dirs[i].m * st.r;
  }
  // nominal code
  E_T_mr = st.E.transpose() * srbi.h + srbi.m * st.r;

  Math::Matrix3d Ixyz = srbi.getInertia();
  Math::Matrix3d res_Ixyz_1 = st.E.transpose() * Ixyz * st.E;

  //  // fd derivative code
  //  for (unsigned i = 0; i < ndirs; i++) {
  //    Matrix3d res_Ixyz_1_i_ad =
  //        st.E.transpose() * Ixyz * st_dirs[i].E
  //        + st.E.transpose() * srbi_dirs[i].getInertia() * st.E
  //        + st_dirs[i].E.transpose() * Ixyz * st.E;
  //    double h = 1e-8;
  //    Matrix3d res_Ixyz_1_i_h_ad =
  //        (st.E + h * st_dirs[i].E).transpose()
  //        * (Ixyz  + h * srbi_dirs[i].getInertia())
  //        * (st.E + h * st_dirs[i].E);
  //    Matrix3d res_Ixyz_1_i_fd = //(res_Ixyz_1_i_h_ad - st.E.transpose() * Ixyz) / h;
  //        (res_Ixyz_1_i_h_ad - res_Ixyz_1) / h;
  //    cout << res_Ixyz_1_i_ad << endl << endl << res_Ixyz_1_i_fd << endl << endl;
  //    cout << (res_Ixyz_1_i_ad - res_Ixyz_1_i_fd).norm() << endl;
  //    CHECK_ARRAY_CLOSE(res_Ixyz_1_i_ad.data(), res_Ixyz_1_i_fd.data(), 9, 1e-6);
  //  }

  // derivative code
  for (unsigned i = 0; i < ndirs; i++) {
    res_dirs[i] = Math::SpatialRigidBodyInertia (
          srbi_dirs[i].m,
          E_T_mr_dirs.col(i),
          st.E.transpose() * Ixyz * st_dirs[i].E
          + st.E.transpose() * srbi_dirs[i].getInertia() * st.E
          + st_dirs[i].E.transpose() * Ixyz * st.E

          - Math::VectorCrossMatrix(st.r) * Math::VectorCrossMatrix (st.E.transpose() * srbi_dirs[i].h)
          - Math::VectorCrossMatrix(st.r) * Math::VectorCrossMatrix (st_dirs[i].E.transpose() * srbi.h)
          - Math::VectorCrossMatrix(st_dirs[i].r) * Math::VectorCrossMatrix (st.E.transpose() * srbi.h)

          - Math::VectorCrossMatrix (E_T_mr) * Math::VectorCrossMatrix (st_dirs[i].r)
          - Math::VectorCrossMatrix (E_T_mr_dirs.col(i)) * Math::VectorCrossMatrix (st.r)
          );
  }

  // nominal code
  res = Math::SpatialRigidBodyInertia (
        srbi.m,
        E_T_mr,
        res_Ixyz_1
        - Math::VectorCrossMatrix(st.r) * Math::VectorCrossMatrix (st.E.transpose() * srbi.h)
        - Math::VectorCrossMatrix (E_T_mr) * Math::VectorCrossMatrix (st.r));
}

// -----------------------------------------------------------------------------
} // namespace AD
// -----------------------------------------------------------------------------
} // namespace RigidBodyDynamics
// -----------------------------------------------------------------------------


#endif
