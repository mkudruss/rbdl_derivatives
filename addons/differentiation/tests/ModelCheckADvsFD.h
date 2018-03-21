#ifndef MODELCHECKADVSFD_h
#define MODELCHECKADVSFD_h

#include <rbdl/rbdl.h>
#include "ModelAD.h"
#include "ModelED.h"
#include "ConstraintsAD.h"

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------

void checkModelsADvsFD (
    unsigned ndirs,
    Model const & ad_model,
    ADModel const & ad_d_model,
    Model const & fd_model,
    ADModel const & fd_d_model);

void checkModelsADvsED (
    unsigned ndirs,
    Model const & ad_model,
    ADModel const & ad_d_model,
    Model const & ed_model,
    EDModel const & ed_d_model);

void checkConstraintSetsADvsFD (
    unsigned ndirs,
    ConstraintSet const & ad_cs,
    ADConstraintSet const & ad_d_cs,
    ConstraintSet const & fd_cs,
    ADConstraintSet const & fd_d_cs);

// -----------------------------------------------------------------------------
} // namespace RigidBodyDynamics
// -----------------------------------------------------------------------------

#endif // MODELCHECKADVSFD_h
