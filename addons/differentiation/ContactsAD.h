#ifndef RBDL_CONTACTS_AD_H
#define RBDL_CONTACTS_AD_H

#include <rbdl/rbdl_math.h>
#include <rbdl/rbdl_mathutils.h>

#include "rbdl/Contacts.h"

#include "ModelAD.h"
#include "KinematicsAD.h"

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------

struct ADConstraintSet {
    ADConstraintSet() {};
};

// -----------------------------------------------------------------------------
namespace AD {
// -----------------------------------------------------------------------------

RBDL_DLLAPI
void CalcContactJacobian(
        Model &model,
        const Math::VectorNd &Q,
        const ConstraintSet &CS,
        Math::MatrixNd &G,
        bool update_kinematics = true
        );

/*
RBDL_DLLAPI
void CalcContactSystemVariables (
        Model &model,
        const Math::VectorNd &Q,
        const Math::VectorNd &QDot,
        const Math::VectorNd &Tau,
        ConstraintSet &CS
        );

RBDL_DLLAPI
void ForwardDynamicsContactsDirect (
        Model &model,
        const Math::VectorNd &Q,
        const Math::VectorNd &QDot,
        const Math::VectorNd &Tau,
        ConstraintSet &CS,
        Math::VectorNd &QDDot
        );

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

RBDL_DLLAPI
void ComputeContactImpulsesDirect (
        Model &model,
        const Math::VectorNd &Q,
        const Math::VectorNd &QDotMinus,
        ConstraintSet &CS,
        Math::VectorNd &QDotPlus
        );

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

RBDL_DLLAPI
void SolveContactSystemDirect (
        Math::MatrixNd &H,
        const Math::MatrixNd &G,
        const Math::VectorNd &c,
        const Math::VectorNd &gamma,
        Math::VectorNd &qddot,
        Math::VectorNd &lambda,
        Math::MatrixNd &A,
        Math::VectorNd &b,
        Math::VectorNd &x,
        Math::LinearSolver &linear_solver
        );

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
