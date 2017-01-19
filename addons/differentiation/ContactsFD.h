#ifndef RBDL_CONTACTS_FD_H
#define RBDL_CONTACTS_FD_H

#include <rbdl/rbdl_math.h>
#include <rbdl/rbdl_mathutils.h>

#include "rbdl/Constraints.h"

#include "KinematicsAD.h"
#include "ContactsAD.h"

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------
namespace FD {
// -----------------------------------------------------------------------------

RBDL_DLLAPI
void CalcConstraintsJacobian (
    Model &model,
    ADModel &ad_model,
    const Math::VectorNd &Q,
    const Math::MatrixNd &Q_dirs,
    ConstraintSet &CS,
    ADConstraintSet &ad_CS,
    Math::MatrixNd &G,
    std::vector<Math::MatrixNd> &G_dirs
    );

RBDL_DLLAPI
void ComputeConstraintImpulsesDirect (
    Model & model,
    const Math::VectorNd & q,
    const Math::MatrixNd & q_dirs,
    const Math::VectorNd & qdot_minus,
    const Math::MatrixNd & qdot_minus_dirs,
    ConstraintSet   & CS,
    Math::VectorNd  & qdot_plus,
    Math::MatrixNd  & fd_qdot_plus
    );


RBDL_DLLAPI
void ComputeConstraintImpulsesDirect (
    Model & model,
    const Math::VectorNd & q,
    const Math::MatrixNd & q_dirs,
    const Math::VectorNd & qdot_minus,
    const Math::MatrixNd & qdot_minus_dirs,
    ConstraintSet & CS,
    Math::MatrixNd & fd_b,
    std::vector<Math::MatrixNd> & fd_A,
    Math::VectorNd & qdot_plus,
    Math::MatrixNd & fd_qdot_plus
    );

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
    Math::MatrixNd & x_fd,
    Math::LinearSolver & linear_solver,
    int ndirs
    );


// -----------------------------------------------------------------------------
} // namespace FD
// -----------------------------------------------------------------------------
} // namespace RigidBodyDynamics
// -----------------------------------------------------------------------------

/* RBDL_CONTACTS_FD_H */
#endif