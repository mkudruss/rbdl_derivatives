#include "ContactsAD.h"
#include "DynamicsAD.h"

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------

ADConstraintSet::ADConstraintSet(ConstraintSet CS, int dof_count) {
  ndirs = 4 * dof_count;

  G.resize(ndirs, CS.G);
}

// -----------------------------------------------------------------------------
namespace AD {
// -----------------------------------------------------------------------------

using namespace RigidBodyDynamics::Math;
using namespace std;

RBDL_DLLAPI
void CalcContactJacobian(
        Model &model,
        ADModel &ad_model,
        const Math::VectorNd &Q,
        const Math::MatrixNd &Q_dirs,
        const ConstraintSet &CS,
        ADConstraintSet &ad_CS,
        Math::MatrixNd &G,
        std::vector<Math::MatrixNd> &G_dirs,
        bool update_kinematics
) {
    unsigned int ndirs = Q_dirs.cols();
    ad_model.resize_directions(ndirs);

    if (update_kinematics) {
        // derivative evaluation
        // UpdateKinematicsCustom (
        //     model, ad_model,
        //     &Q, &Q_dirs,
        //     NULL, NULL,
        //     NULL, NULL
        // );
        // nominal evaluation
        // NOTE kinematics are already updated in the AD version
        // UpdateKinematicsCustom (model, &Q, NULL, NULL);
    }

    unsigned int i,j;

    // variables to check whether we need to recompute G
    unsigned int prev_body_id = 0;
    Math::Vector3d prev_body_point = Math::Vector3d::Zero();
    // derivative evaluation
    std::vector<Math::MatrixNd> Gi_dirs (
        ndirs, Math::MatrixNd (3, model.dof_count)
    );
    // nominal evaluation
    Math::MatrixNd Gi (3, model.dof_count);

    for (i = 0; i < CS.size(); i++) {
        // only compute the matrix Gi if actually needed
        if (prev_body_id != CS.body[i] || prev_body_point != CS.point[i]) {
            Gi.setZero();
            // derivative evaluation
            AD::CalcPointJacobian(
                model, ad_model,
                Q, Q_dirs,
                CS.body[i], CS.point[i],
                Gi, Gi_dirs,
                false
            );
            // nominal evaluation
            // NOTE nominal Jacobian is already computed by the AD version
            // CalcPointJacobian (model, Q, CS.body[i], CS.point[i], Gi, false);

            prev_body_id = CS.body[i];
            prev_body_point = CS.point[i];
        }

        for (j = 0; j < model.dof_count; j++) {
            // derivative evaluation
            for (unsigned int idirs = 0; idirs < ndirs; idirs++) {
                Math::Vector3d gaxis_dirs (
                    Gi_dirs[idirs](0,j),
                    Gi_dirs[idirs](1,j),
                    Gi_dirs[idirs](2,j)
                );
                G_dirs[idirs](i,j) = gaxis_dirs.transpose() * CS.normal[i];
            }
            // nominal evaluation
            Math::Vector3d gaxis (Gi(0,j), Gi(1,j), Gi(2,j));
            // nominal evaluation
            G(i,j) = gaxis.transpose() * CS.normal[i];
        }
    }
}

RBDL_DLLAPI
void SolveContactSystemDirect (
    const MatrixNd &H,
    const vector<MatrixNd> & H_dirs,
    const MatrixNd &G,
    const vector<MatrixNd> & G_dirs,
    const VectorNd & c,
    const MatrixNd & c_dirs,
    const VectorNd & gamma,
    const MatrixNd & gamma_dirs,
    MatrixNd &A,
    vector<MatrixNd> & A_dirs,
    VectorNd & b,
    MatrixNd & b_dirs,
    VectorNd & x,
    MatrixNd & x_ad,
    LinearSolver & linear_solver
) {
  int ndirs = H_dirs.size();
  assert(ndirs == G_dirs.size());
  assert(ndirs == c_dirs.cols());
  assert(ndirs == gamma_dirs.cols());
  assert(ndirs == A_dirs.size());
  assert(ndirs == b_dirs.cols());
  assert(ndirs == x_ad.cols());

  // derivative construction
  for (int i = 0; i < ndirs; i++) {
    A_dirs[i].block(0, 0, c.rows(), c.rows())            = H_dirs[i];
    A_dirs[i].block(0, c.rows(), c.rows(), gamma.rows()) =
        G_dirs[i].transpose();
    A_dirs[i].block(c.rows(), 0, gamma.rows(), c.rows()) = G_dirs[i];
    b_dirs.block(0, i, c.rows(), 1)                      = c_dirs.col(i);
    b_dirs.block(c.rows(), i, gamma.rows(), 1)           = gamma_dirs.col(i);
  }
  // nominal construction
  //    Build the system: Copy H
  A.block(0, 0, c.rows(), c.rows()) = H;
  //    Copy G and G^T
  A.block(0, c.rows(), c.rows(), gamma.rows()) = G.transpose();
  A.block(c.rows(), 0, gamma.rows(), c.rows()) = G;
  //    Build the system: Copy -C + \tau
  b.block(0, 0, c.rows(), 1) = c;
  b.block(c.rows(), 0, gamma.rows(), 1) = gamma;

  LOG << "A = " << std::endl << A << std::endl;
  LOG << "b = " << std::endl << b << std::endl;

  switch (linear_solver) {
    case (LinearSolverPartialPivLU) :
#ifdef RBDL_USE_SIMPLE_MATH
      // SimpleMath does not have a LU solver so just use its QR solver
      x = A.householderQr().solve(b);
      /// TODO: implement for simple math
#else
      {
        Eigen::PartialPivLU<MatrixNd::PlainObject> A_LU = A.partialPivLu();
        // nominal evaluation
        x = A_LU.solve(b);
        // derivative evaluation
        for (int i = 0; i < ndirs; i++) {
          x_ad.col(i) = A_LU.solve(b_dirs.col(i) - A_dirs[i] * x);
        }
      }
#endif
      break;
    case (LinearSolverColPivHouseholderQR) :
      {
        Eigen::ColPivHouseholderQR<MatrixNd::PlainObject> A_CPQR =
            A.colPivHouseholderQr();
        // nominal evaluation
        x = A_CPQR.solve(b);
        // derivative evaluation
        for (int i = 0; i < ndirs; i++) {
          x_ad.col(i) = A_CPQR.solve(b_dirs.col(i) - A_dirs[i] * x);
        }
      }
      break;
    case (LinearSolverHouseholderQR) :
      {
        Eigen::HouseholderQR<MatrixNd::PlainObject> A_QR = A.householderQr();
        // nominal evaluation
        x = A_QR.solve(b);
        // derivative evaluation;
        for (int i = 0; i < ndirs; i++) {
          x_ad.col(i) = A_QR.solve(b_dirs.col(i) - A_dirs[i] * x);
        }
      }
      break;
    default:
      LOG << "Error: Invalid linear solver: " << linear_solver << std::endl;
      assert (0);
      break;
  }
  LOG << "x = " << std::endl << x << std::endl;
}

RBDL_DLLAPI
void ComputeContactImpulsesDirect (
    Model & model,
    ADModel & ad_model,
    const VectorNd & q,
    const MatrixNd & q_dirs,
    const VectorNd & qDotMinus,
    const MatrixNd & qDotMinus_dirs,
    ConstraintSet & CS,
    ADConstraintSet & ad_CS,
    VectorNd & qDotPlus,
    MatrixNd & ad_qDotPlus
) {
  int ndirs = q_dirs.size();
  assert(ndirs + qDotMinus_dirs.size());

  // Compute H
  UpdateKinematicsCustom(model, ad_model, &q, &q_dirs, 0, 0, 0, 0);

  vector<MatrixNd> H_ad(ndirs, CS.H);
  CompositeRigidBodyAlgorithm(model, ad_model, q, q_dirs, CS.H, H_ad, false);
  /// TODO: Check if false makes sense as last argument here.

  // Compute G
  CalcContactJacobian (model, ad_model, q, q_dirs, CS, ad_CS, CS.G, ad_CS.G, false);
  /// TODO: Check if false makes sense as last argument here.

  VectorNd c   = CS.H * qDotMinus;
  MatrixNd c_ad(CS.H.rows(), ndirs);
  for (int i = 0; i < ndirs; i++) {
    c_ad.block(0, i, CS.H.rows(), 1) = CS.H * qDotMinus_dirs.col(i) + H_ad[i] * qDotMinus;
  }

  MatrixNd         gamma_dirs = MatrixNd::Zero(CS.v_plus.rows(), ndirs);
  vector<MatrixNd> A_dirs(ndirs, MatrixNd::Zero(CS.A.rows(), CS.A.cols()));
  MatrixNd         b_dirs = MatrixNd::Zero(CS.b.rows(), ndirs);
  MatrixNd         x_ad   = MatrixNd::Zero(CS.x.rows(), ndirs);
  SolveContactSystemDirect (CS.H, H_ad, CS.G, ad_CS.G, c, c_ad, CS.v_plus,
                            gamma_dirs, CS.A, A_dirs, CS.b, b_dirs,
                            CS.x, x_ad, CS.linear_solver);
  /// TODO: make gamma_dirs, A_dirs, b_dirs, x_ad members of ADConstraintSet

  // derivative evaluation
  ad_qDotPlus = x_ad.block(0, 0, model.dof_count, ndirs);
  // nominal evaluation
  //    Copy back QDotPlus
  qDotPlus = CS.x.segment(0, model.dof_count);

  // Copy back constraint impulses
  for (unsigned int i = 0; i < CS.size(); i++) {
    CS.impulse[i] = CS.x[model.dof_count + i];
  }
  /// TODO: compute derivative of impulse
}

// -----------------------------------------------------------------------------
} // namespace AD
// -----------------------------------------------------------------------------
} // namespace RigidBodyDynamics
// -----------------------------------------------------------------------------
