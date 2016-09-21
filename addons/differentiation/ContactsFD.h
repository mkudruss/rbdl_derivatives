#ifndef RBDL_CONTACTS_FD_H
#define RBDL_CONTACTS_FD_H

#include <rbdl/rbdl_math.h>
#include <rbdl/rbdl_mathutils.h>

#include "rbdl/Contacts.h"

#include "KinematicsAD.h"
#include "ContactsAD.h"

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------
namespace FD {
// -----------------------------------------------------------------------------

RBDL_DLLAPI
void CalcContactJacobian(
    Model &model,
    const Math::VectorNd &Q,
    const Math::VectorNd &Q_dirs,
    const ConstraintSet &CS,
    ADConstraintSet &ad_CS,
    Math::MatrixNd &G,
    std::vector<Math::MatrixNd> &G_dirs,
    bool update_kinematics = true
    );

RBDL_DLLAPI
void ComputeContactImpulsesDirect(Model & model,
    const Math::VectorNd & q,
    const Math::MatrixNd & q_dirs,
    const Math::VectorNd & qdot_minus,
    const Math::MatrixNd & qdot_minus_dirs,
    ConstraintSet   & CS,
    Math::VectorNd  & qdot_plus,
    Math::MatrixNd  & ad_qdot_plus
    );

// -----------------------------------------------------------------------------
} // namespace FD
// -----------------------------------------------------------------------------
} // namespace RigidBodyDynamics
// -----------------------------------------------------------------------------

/* RBDL_CONTACTS_FD_H */
#endif
