#include "ConstraintsAD.h"
#include "DynamicsAD.h"
#include "SpatialAlgebraOperatorsAD.h"

#include <iomanip>

using namespace std;

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------

using namespace RigidBodyDynamics::Math;

ADConstraintSet::ADConstraintSet(const ConstraintSet & CS, int dof_count) {
  ndirs = 4 * dof_count;

  G.resize(ndirs,  MatrixNd::Zero(CS.G.rows(), CS.G.cols()));
  Gi.resize(ndirs, MatrixNd::Zero(3, dof_count));
  A.resize(ndirs, MatrixNd::Zero(CS.A.rows(), CS.A.cols()));
  H.resize(ndirs, MatrixNd::Zero(CS.H.rows(), CS.H.cols()));
  b             = MatrixNd::Zero(CS.b.rows(), ndirs);
  v_plus        = MatrixNd::Zero(CS.v_plus.rows(), ndirs);
  x             = MatrixNd::Zero(CS.x.rows(), ndirs);
  impulse       = MatrixNd::Zero(CS.impulse.rows(), ndirs);
  QDDot_0       = MatrixNd::Zero(CS.QDDot_0.rows(), ndirs);
  C             = MatrixNd::Zero(CS.C.rows(), ndirs);
  gamma         = MatrixNd::Zero(CS.gamma.rows(), ndirs);
  force         = MatrixNd::Zero(CS.force.rows(), ndirs);
}

// -----------------------------------------------------------------------------
namespace AD {
// -----------------------------------------------------------------------------

RBDL_DLLAPI
void CalcContactJacobian(
    Model &model,
    ADModel &ad_model,
    const Math::VectorNd &q,
    const Math::MatrixNd &q_dirs,
    ConstraintSet & CS,
    ADConstraintSet & ad_CS,
    Math::MatrixNd &G,
    vector<Math::MatrixNd> &G_dirs,
    bool update_kinematics) {
  unsigned int ndirs = q_dirs.cols();
  assert(ndirs <= G_dirs.size());

  if (update_kinematics) {
    // derivative evaluation
    UpdateKinematicsCustom (
          model, ad_model,
          &q, &q_dirs,
          NULL, NULL,
          NULL, NULL
          );
    // nominal evaluation
    // NOTE kinematics are already updated in the AD version
    // UpdateKinematicsCustom (model, &Q, NULL, NULL);
  }

  // variables to check whether we need to recompute G.
  ConstraintSet::ConstraintType prev_constraint_type
    = ConstraintSet::ConstraintTypeLast;
  unsigned int prev_body_id_1 = 0;
  // unsigned int prev_body_id_2 = 0;
  SpatialTransform prev_body_X_1;
  SpatialTransform prev_body_X_2;

  for (unsigned int i = 0; i < CS.contactConstraintIndices.size(); i++) {
    const unsigned int c = CS.contactConstraintIndices[i];

    // only compute the matrix Gi if actually needed
    if (prev_constraint_type != CS.constraintType[c]
        || prev_body_id_1 != CS.body[c]
        || prev_body_X_1.r != CS.point[c]) {

      // Compute the jacobian for the point.
      for (unsigned idir = 0; idir < ndirs; idir++) {
        ad_CS.Gi[idir].setZero();
      }
      CS.Gi.setZero();
      CalcPointJacobian (model, ad_model, q, q_dirs,
                         CS.body[c], CS.point[c],
                         CS.Gi, ad_CS.Gi,
                         false);

      prev_constraint_type = ConstraintSet::ContactConstraint;

      // Update variables for optimization check.
      prev_body_id_1 = CS.body[c];
      prev_body_X_1 = Xtrans(CS.point[c]);
    }

    for(unsigned int j = 0; j < model.dof_count; j++) {
      MatrixNd gaxis_dirs(3, ndirs);
      for (unsigned idir = 0; idir < ndirs; idir++) {
        gaxis_dirs.col(idir) = ad_CS.Gi[idir].col(j);
      }
      Vector3d gaxis (CS.Gi(0,j), CS.Gi(1,j), CS.Gi(2,j));
      for (unsigned idir = 0; idir < ndirs; idir++) {
        G_dirs[idir](c,j) = gaxis_dirs.col(idir).transpose() * CS.normal[c];
      }
      G(c,j) = gaxis.transpose() * CS.normal[c];
    }
  }

//  unsigned int i,j;

//  // variables to check whether we need to recompute G
//  unsigned int prev_body_id = 0;
//  Vector3d prev_body_point = Vector3d::Zero();
//  // derivative evaluation
//  vector<MatrixNd> Gi_dirs (ndirs, MatrixNd::Zero (3, model.dof_count));
//  // nominal evaluation
//  MatrixNd Gi (3, model.dof_count);

//  for (i = 0; i < CS.size(); i++) {
//    // only compute the matrix Gi if actually needed
//    if (prev_body_id != CS.body[i] || prev_body_point != CS.point[i]) {
//      Gi.setZero();
//      // derivative evaluation
//      AD::CalcPointJacobian(
//            model, ad_model,
//            q, q_dirs,
//            CS.body[i], CS.point[i],
//            Gi, Gi_dirs,
//            false
//            );
//      // nominal evaluation
//      // NOTE nominal Jacobian is already computed by the AD version
//      // CalcPointJacobian (model, Q, CS.body[i], CS.point[i], Gi, false);

//      prev_body_id = CS.body[i];
//      prev_body_point = CS.point[i];
//    }

//    for (j = 0; j < model.dof_count; j++) {
//      // derivative evaluation
//      for (unsigned idir = 0; idir < ndirs; idir++) {
//        Vector3d gaxis_dirs (
//              Gi_dirs[idir](0,j),
//              Gi_dirs[idir](1,j),
//              Gi_dirs[idir](2,j));
//            G_dirs[idir](i,j) = gaxis_dirs.transpose() * CS.normal[i];
//      }
//      // nominal evaluation
//      Math::Vector3d gaxis (Gi(0,j), Gi(1,j), Gi(2,j));
//      // nominal evaluation
//      G(i,j) = gaxis.transpose() * CS.normal[i];
//    }
//  }
}

RBDL_DLLAPI
void CalcContactSystemVariables (
    Model &model,
    ADModel &ad_model,
    const VectorNd  &q,
    const MatrixNd  &q_dirs,
    const VectorNd  &qdot,
    const MatrixNd  &qdot_dirs,
    const VectorNd  &tau,
    const MatrixNd  &tau_dirs,
    ConstraintSet   &CS,
    ADConstraintSet &ad_CS) {
  int ndirs = q_dirs.cols();
  assert(ndirs == qdot_dirs.cols());
  assert(ndirs == tau_dirs.cols());

  // Compute C
  NonlinearEffects (model, ad_model, q, q_dirs, qdot, qdot_dirs, CS.C, ad_CS.C);
  assert (CS.H.cols() == model.dof_count && CS.H.rows() == model.dof_count);

  // Compute H
  CompositeRigidBodyAlgorithm (model, ad_model, q, q_dirs, CS.H, ad_CS.H, false);
  /// TODO: Check if false makes sense as last argument here.

  // Compute G
  // We have to update model.X_base as they are not automatically computed
  // by NonlinearEffects()
  for (unsigned int i = 1; i < model.mBodies.size(); i++) {
    unsigned int lambda = model.lambda[i];
    // derivative evaluation
//    for (int idir = 0; idir < ndirs; idir++) {
//      ad_model.X_base[i][idir] =
//          ad_model.X_lambda[i][idir] * model.X_base[lambda] //.toMatrix()
//          + model.X_lambda[i] .toMatrix() * ad_model.X_base[lambda][idir];
//    }

    mulSTST (ndirs,
             model.X_lambda[i], ad_model.X_lambda[i],
             model.X_base[lambda], ad_model.X_base[lambda],
             model.X_base[i], ad_model.X_base[i]);

    // nominal evaluation
    model.X_base[i] = model.X_lambda[i] * model.X_base[model.lambda[i]];
  }
  CalcContactJacobian(model, ad_model, q, q_dirs, CS, ad_CS, CS.G, ad_CS.G, false);
  /// TODO: Check if false makes sense as last argument here.

  // Compute gamma
  unsigned int prev_body_id = 0;
  Vector3d prev_body_point = Vector3d::Zero();
  Vector3d gamma_i = Vector3d::Zero();
  MatrixNd ad_gamma_i(gamma_i.rows(), ndirs);
  // derivative code
  ad_CS.QDDot_0.setZero();
  // nominal code
  CS.QDDot_0.setZero();
  UpdateKinematicsCustom (model, ad_model, 0, 0, 0, 0, &CS.QDDot_0, &ad_CS.QDDot_0);
  for (unsigned int i = 0; i < CS.size(); i++) {
    // only compute point accelerations when necessary
    if (prev_body_id != CS.body[i] || prev_body_point != CS.point[i]) {
      // derivative and nominal  code
      gamma_i = CalcPointAcceleration(model, ad_model, q, q_dirs,
          qdot, qdot_dirs, CS.QDDot_0, ad_CS.QDDot_0, CS.body[i], CS.point[i],
          ad_gamma_i, false);
      /// TODO: Check if false makes sense as last argument here.
      // nominal code
      prev_body_id = CS.body[i];
      prev_body_point = CS.point[i];
    }

    // we also substract ContactData[i].acceleration such that the contact
    // point will have the desired acceleration
    // derivative code
    for (int idir = 0; idir < ndirs; idir++) {
      ad_CS.gamma(i, idir) = - CS.normal[i].dot(ad_gamma_i.col(idir));
    }
    // nominal code
    CS.gamma[i] = CS.acceleration[i] - CS.normal[i].dot(gamma_i);
  }
  /// TODO: Implement derivative of for loop, esp. CalcPointAcceleration
  ///
  ///
}

RBDL_DLLAPI
void ForwardDynamicsContactsDirect (
    Model   & model,
    ADModel & ad_model,
    const VectorNd & q,
    const MatrixNd & q_dirs,
    const VectorNd & qdot,
    const MatrixNd & qdot_dirs,
    const VectorNd & tau,
    const MatrixNd & tau_dirs,
    ConstraintSet & CS,
    ADConstraintSet & ad_CS,
    VectorNd & qddot,
    MatrixNd & ad_qddot
    ) {
  LOG << "-------- " << __func__ << " --------" << std::endl;
  int ndirs = q_dirs.cols();
  assert(ndirs == qdot_dirs.cols());
  assert(ndirs == tau_dirs.cols());

  // derivative and nominal code
  CalcContactSystemVariables (model, ad_model, q, q_dirs, qdot, qdot_dirs,
                              tau, tau_dirs, CS, ad_CS);

  // derivative and nominal code
  VectorNd c    = tau - CS.C;
  MatrixNd ad_c(c.rows(), ndirs);
  for (int idir = 0; idir < ndirs; idir++) {
    ad_c.col(idir) = tau_dirs.col(idir) - ad_CS.C.col(idir);
  }
  SolveContactSystemDirect (CS.H, ad_CS.H, CS.G, ad_CS.G, c, ad_c,
                            CS.gamma, ad_CS.gamma, CS.A, ad_CS.A, CS.b, ad_CS.b,
                            CS.x, ad_CS.x, CS.linear_solver, ndirs);

  // Copy back QDDot
  // derivative code
  ad_qddot = ad_CS.x.block(0, 0, model.dof_count, ndirs);
  // nominal code
  qddot = CS.x.segment(0, model.dof_count);

  // Copy back contact forces
  // derivative code
  ad_CS.force = -ad_CS.x.block(model.dof_count, 0, CS.size(), ndirs);
  // nominal code
  CS.force    = -CS.x.segment(model.dof_count, CS.size());
}

RBDL_DLLAPI
void ComputeContactImpulsesDirect (
    Model & model,
    ADModel & ad_model,
    const VectorNd & q,
    const MatrixNd & q_dirs,
    const VectorNd & qdot_minus,
    const MatrixNd & qdot_minus_dirs,
    ConstraintSet & CS,
    ADConstraintSet & ad_CS,
    VectorNd & qdot_plus,
    MatrixNd & ad_qdot_plus
) {
  int ndirs = q_dirs.cols();
  assert(ndirs == qdot_minus_dirs.cols());
  assert(ndirs == ad_qdot_plus.cols());

  // Compute H
  UpdateKinematicsCustom(model, ad_model, &q, &q_dirs, 0, 0, 0, 0);

//  vector<MatrixNd> H_ad(ndirs, CS.H);
  CompositeRigidBodyAlgorithm(model, ad_model, q, q_dirs, CS.H, ad_CS.H, false);
  /// TODO: Check if false makes sense as last argument here.

  // Compute G
  CalcContactJacobian (model, ad_model, q, q_dirs, CS, ad_CS, CS.G, ad_CS.G, false);
  /// TODO: Check if false makes sense as last argument here.

  VectorNd c = CS.H * qdot_minus;
  MatrixNd c_ad(CS.H.rows(), ndirs);
  for (int i = 0; i < ndirs; i++) {
    c_ad.block(0, i, CS.H.rows(), 1) = CS.H * qdot_minus_dirs.col(i) + ad_CS.H[i] * qdot_minus;
  }

  SolveContactSystemDirect (CS.H, ad_CS.H, CS.G, ad_CS.G, c, c_ad, CS.v_plus,
                            ad_CS.v_plus, CS.A, ad_CS.A, CS.b, ad_CS.b,
                            CS.x, ad_CS.x, CS.linear_solver, ndirs);

  // derivative evaluation
  ad_qdot_plus = ad_CS.x.block(0, 0, model.dof_count, ndirs);
  // nominal evaluation
  //    Copy back QDotPlus
  qdot_plus = CS.x.segment(0, model.dof_count);

  // derivative evaluation
  ad_CS.impulse = ad_CS.x.block(model.dof_count, 0, CS.size(), ndirs);
  // nominal evaluation
  //    Copy back constraint impulses
  CS.impulse = CS.x.segment(model.dof_count, CS.size());
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
    LinearSolver & linear_solver,
    unsigned ndirs) {
  assert(ndirs <= H_dirs.size());
  assert(ndirs <= G_dirs.size());
  assert(ndirs <= c_dirs.cols());
  assert(ndirs <= gamma_dirs.cols());
  assert(ndirs <= A_dirs.size());
  assert(ndirs <= b_dirs.cols());
  assert(ndirs <= x_ad.cols());

  // derivative construction
  for (unsigned idir = 0; idir < ndirs; idir++) {
    A_dirs[idir].block(0, 0, c.rows(), c.rows())            = H_dirs[idir];
    A_dirs[idir].block(0, c.rows(), c.rows(), gamma.rows()) =
        G_dirs[idir].transpose();
    A_dirs[idir].block(c.rows(), 0, gamma.rows(), c.rows()) = G_dirs[idir];
    b_dirs.block(0, idir, c.rows(), 1)                      = c_dirs.col(idir);
    b_dirs.block(c.rows(), idir, gamma.rows(), 1)           = gamma_dirs.col(idir);
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
        // nominal code
        x = A_LU.solve(b);
        // derivative code
        for (unsigned idir = 0; idir < ndirs; idir++) {
          x_ad.col(idir) = A_LU.solve(b_dirs.col(idir) - A_dirs[idir] * x);
        }
      }
#endif
      break;
    case (LinearSolverColPivHouseholderQR) :
      {
        Eigen::ColPivHouseholderQR<MatrixNd::PlainObject> A_CPQR =
            A.colPivHouseholderQr();
        // nominal code
        x = A_CPQR.solve(b);
        // derivative code
        for (unsigned idir = 0; idir < ndirs; idir++) {
          x_ad.col(idir) = A_CPQR.solve(b_dirs.col(idir) - A_dirs[idir] * x);
        }
      }
      break;
    case (LinearSolverHouseholderQR) :
      {
        Eigen::HouseholderQR<MatrixNd::PlainObject> A_QR = A.householderQr();
        // nominal evaluation
        x = A_QR.solve(b);
        // derivative evaluation
        for (unsigned idir = 0; idir < ndirs; idir++) {
          x_ad.col(idir) = A_QR.solve(b_dirs.col(idir) - A_dirs[idir] * x);
        }
      }
      break;
//    case (LinearSolverLLT) :
//      {
//        Eigen::LLT<MatrixNd::PlainObject> A_LLT = A.llt();
//        // nominal code
//        x = A_LLT.solve(b);
//        // derivative code
//        for (int i = 0; i < ndirs; i++) {
//          x_ad.col(i) = A_LLT.solve(b_dirs.col(i) - A_dirs[i] * x);
//        }
//      }
//      break;
/// TODO: Add support for this solver also.
    default:
      LOG << "Error: Invalid linear solver: " << linear_solver << std::endl;
      cerr << __FILE__ << " " << __LINE__
           << ": Error: Invalid linear solver: " << linear_solver << std::endl;
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
