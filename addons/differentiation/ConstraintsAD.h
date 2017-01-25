#ifndef RBDL_CONTACTS_AD_H
#define RBDL_CONTACTS_AD_H

#include <rbdl/rbdl_math.h>
#include <rbdl/rbdl_mathutils.h>

#include "rbdl/Constraints.h"

#include "ModelAD.h"
#include "KinematicsAD.h"

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------

struct ADConstraintSet {
  int ndirs;

  std::vector<Math::MatrixNd> G;
  std::vector<Math::MatrixNd> Gi;
  std::vector<Math::MatrixNd> A;
  std::vector<Math::MatrixNd> H;
  Math::MatrixNd              b;
  Math::MatrixNd              v_plus;
  Math::MatrixNd              x;
  Math::MatrixNd              impulse;
  Math::MatrixNd              QDDot_0;
  Math::MatrixNd              C;
  Math::MatrixNd              gamma;
  Math::MatrixNd              force;
  Math::MatrixNd              err;
  Math::MatrixNd              errd;

  ADConstraintSet() {}
  ADConstraintSet(const ConstraintSet &CS, int dof_count);
};

// -----------------------------------------------------------------------------
namespace AD {
// -----------------------------------------------------------------------------

RBDL_DLLAPI void CalcConstraintsPositionError (
    Model &model,
    ADModel &ad_model,
    const Math::VectorNd &q,
    const Math::MatrixNd &q_dirs,
    ConstraintSet &cs,
    ADConstraintSet &ad_cs,
    Math::VectorNd &err,
    Math::MatrixNd &ad_err,
    bool update_kinematics
    );

RBDL_DLLAPI
void CalcContactJacobian(
		Model   & model,
		ADModel & ad_model,
		const Math::VectorNd & q,
		const Math::MatrixNd & q_dirs,
    ConstraintSet  & CS,
		ADConstraintSet & ad_CS,
		Math::MatrixNd  & G,
		std::vector<Math::MatrixNd> & G_dirs,
		bool update_kinematics = true
		);

RBDL_DLLAPI void CalcConstrainedSystemVariables (
    Model &model,
    ADModel &ad_model,
    const Math::VectorNd  &q,
    const Math::MatrixNd  &q_dirs,
    const Math::VectorNd  &qdot,
    const Math::MatrixNd  &qdot_dirs,
    const Math::VectorNd  &tau,
    const Math::MatrixNd  &tau_dirs,
    ConstraintSet   &CS,
    ADConstraintSet &ad_CS);

RBDL_DLLAPI
void ForwardDynamicsContactsDirect (
    Model   & model,
    ADModel & ad_model,
    const Math::VectorNd & q,
    const Math::MatrixNd & q_dirs,
    const Math::VectorNd & qdot,
    const Math::MatrixNd & qdot_dirs,
    const Math::VectorNd & tau,
    const Math::MatrixNd & tau_dirs,
    ConstraintSet   & CS,
    ADConstraintSet & ad_CS,
    Math::VectorNd  & qddot,
    Math::MatrixNd  & ad_qddot
    );

/*
RBDL_DLLAPI
void ForwardDynamicsContactsRangeSpaceSparse (
				Model &model,
				const Math::VectorNd &Q,
				const Math::VectorNd &QDot,
				const Math::VectorNd &Tau,
				ConstraintSet &CS,
				Math::VectorNd &QDDot
				);

RBDL_DLLAPI
void ForwardDynamicsContactsNullSpace (
				Model &model,
				const Math::VectorNd &Q,
				const Math::VectorNd &QDot,
				const Math::VectorNd &Tau,
				ConstraintSet &CS,
				Math::VectorNd &QDDot
				);

RBDL_DLLAPI
void ForwardDynamicsContactsKokkevis (
				Model &model,
				const Math::VectorNd &Q,
				const Math::VectorNd &QDot,
				const Math::VectorNd &Tau,
				ConstraintSet &CS,
				Math::VectorNd &QDDot
				);

*/

RBDL_DLLAPI
void ComputeContactImpulsesDirect(
		Model & model,
		ADModel & ad_model,
		const Math::VectorNd & q,
		const Math::MatrixNd & q_dirs,
		const Math::VectorNd & qdot_minus,
		const Math::MatrixNd & qdot_minus_dirs,
		ConstraintSet   & CS,
		ADConstraintSet & ad_CS,
		Math::VectorNd  & qdot_plus,
		Math::MatrixNd  & ad_qdot_plus
		);

/*
RBDL_DLLAPI
void ComputeContactImpulsesRangeSpaceSparse (
				Model &model,
				const Math::VectorNd &Q,
				const Math::VectorNd &QDotMinus,
				ConstraintSet &CS,
				Math::VectorNd &QDotPlus
				);

RBDL_DLLAPI
void ComputeContactImpulsesNullSpace (
				Model &model,
				const Math::VectorNd &Q,
				const Math::VectorNd &QDotMinus,
				ConstraintSet &CS,
				Math::VectorNd &QDotPlus
				);
*/

RBDL_DLLAPI
void SolveContactSystemDirect (
    const Math::MatrixNd &H,
    const std::vector<Math::MatrixNd> & H_dirs,
    const Math::MatrixNd &G,
    const std::vector<Math::MatrixNd> & G_dirs,
    const Math::VectorNd & c,
    const Math::MatrixNd & c_dirs,
    const Math::VectorNd & gamma,
    const Math::MatrixNd & gamma_dirs,
    Math::MatrixNd & A,
    std::vector<Math::MatrixNd> & A_dirs,
    Math::VectorNd & b,
    Math::MatrixNd & b_dirs,
    Math::VectorNd & x,
    Math::MatrixNd & x_ad,
    Math::LinearSolver & linear_solver,
    unsigned ndirs
    );

/*
RBDL_DLLAPI
void SolveContactSystemRangeSpaceSparse (
				Model &model,
				Math::MatrixNd &H,
				const Math::MatrixNd &G,
				const Math::VectorNd &c,
				const Math::VectorNd &gamma,
				Math::VectorNd &qddot,
				Math::VectorNd &lambda,
				Math::MatrixNd &K,
				Math::VectorNd &a,
				Math::LinearSolver linear_solver
				);

RBDL_DLLAPI
void SolveContactSystemNullSpace (
				Math::MatrixNd &H,
				const Math::MatrixNd &G,
				const Math::VectorNd &c,
				const Math::VectorNd &gamma,
				Math::VectorNd &qddot,
				Math::VectorNd &lambda,
				Math::MatrixNd &Y,
				Math::MatrixNd &Z,
				Math::VectorNd &qddot_y,
				Math::VectorNd &qddot_z,
				Math::LinearSolver &linear_solver
				);
*/

// -----------------------------------------------------------------------------
} // namespace AD
// -----------------------------------------------------------------------------
} // namespace RigidBodyDynamics
// -----------------------------------------------------------------------------

/* RBDL_CONTACTS_AD_H */
#endif
