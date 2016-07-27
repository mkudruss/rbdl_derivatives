#ifndef RBDL_CONTACTS_FD_H
#define RBDL_CONTACTS_FD_H

#include <rbdl/rbdl_math.h>
#include <rbdl/rbdl_mathutils.h>

#include "rbdl/Contacts.h"

#include "ModelAD.h"
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
        ADModel &ad_model,
        const Math::VectorNd &Q,
        const Math::VectorNd &Q_dirs,
        const ConstraintSet &CS,
        const ADConstraintSet &ad_CS,
        Math::MatrixNd &G,
        Math::MatrixNd &G_dirs,
        bool update_kinematics = true
        );

// -----------------------------------------------------------------------------
} // namespace FD
// -----------------------------------------------------------------------------
} // namespace RigidBodyDynamics
// -----------------------------------------------------------------------------

/* RBDL_CONTACTS_FD_H */
#endif
