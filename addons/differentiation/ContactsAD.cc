#include "ContactsAD.h"

#include <rbdl/rbdl_math.h>
#include <rbdl/rbdl_mathutils.h>

//#include <Eigen/LU>

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
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
        const ADConstraintSet &ad_CS,
        Math::MatrixNd &G,
        Math::MatrixNd &G_dirs,
        bool update_kinematics = true
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
    Math::MatrixNd Gi (3, model.dof_count);

    for (i = 0; i < CS.size(); i++) {
        // only compute the matrix Gi if actually needed
        if (prev_body_id != CS.body[i] || prev_body_point != CS.point[i]) {
            Gi.setZero();
            // nominal evaluation
            CalcPointJacobian (model, Q, CS.body[i], CS.point[i], Gi, false);
            prev_body_id = CS.body[i];
            prev_body_point = CS.point[i];
        }

        for (j = 0; j < model.dof_count; j++) {
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
    const VectorNd &c,
    const vector<MatrixNd> & c_dirs,
    const VectorNd &gamma,
    const vector<MatrixNd> & gamma_dirs,
    VectorNd &qddot,
    VectorNd &lambda,
    MatrixNd &A,
    vector<MatrixNd> & A_dirs,
    VectorNd &b,
    vector<MatrixNd> & b_dirs,
    VectorNd &x,
    vector<VectorNd> & ad_x,
    LinearSolver &linear_solver
) {
  int ndirs = H_dirs.size();
  assert(ndirs == G_dirs.size());
  assert(ndirs == c_dirs.size());
  assert(ndirs == gamma_dirs.size());
  assert(ndirs == A_dirs.size());
  assert(ndirs == b_dirs.size());

  // derivative construction
  for (int i = 0; i < ndirs; i++) {
    A_dirs[i].block(0, 0, c.rows(), c.rows())            = H_dirs[i];
    A_dirs[i].block(0, c.rows(), c.rows(), gamma.rows()) = G_dirs[i].transpose();
    A_dirs[i].block(c.rows(), 0, gamma.rows(), c.rows()) = G_dirs[i];
    b_dirs[i].block(0, 0, c.rows(), 1)                   = c_dirs[i];
    b_dirs[i].block(c.rows(), 0, gamma.rows(), 1)        = gamma_dirs[i];
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
#else
      {
        Eigen::PartialPivLU<MatrixNd::PlainObject> A_LU = A.partialPivLu();
        // nominal evaluation
        x = A_LU.solve(b);
        // derivative evaluation
        for (int i = 0; i < ndirs; i++) {
          ad_x[i] = A_LU.solve(b_dirs[i] - A_dirs[i] * x);
        }
      }
#endif
      break;
    case (LinearSolverColPivHouseholderQR) :
      {
        Eigen::ColPivHouseholderQR<MatrixNd::PlainObject> A_CPQR = A.colPivHouseholderQr();
        // nominal evaluation
        x = A_CPQR.solve(b);
        // derivative evaluation
        for (int i = 0; i < ndirs; i++) {
          ad_x[i] = A_CPQR.solve(b_dirs[i] - A_dirs[i] * x);
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
          ad_x[i] = A_QR.solve(b_dirs[i] - A_dirs[i] * x);
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


// -----------------------------------------------------------------------------
} // namespace AD
// -----------------------------------------------------------------------------
} // namespace RigidBodyDynamics
// -----------------------------------------------------------------------------
