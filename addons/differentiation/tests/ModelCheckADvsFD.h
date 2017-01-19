#ifndef MODELCHECKADVSFD_h
#define MODELCHECKADVSFD_h

#include <rbdl/rbdl.h>
#include "ModelAD.h"

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------

void checkModelsADvsFD (
    unsigned ndirs,
    Model const & ad_model,
    ADModel const & ad_d_model,
    Model const & fd_model,
    ADModel const & fd_d_model);

// -----------------------------------------------------------------------------
} // namespace RigidBodyDynamics
// -----------------------------------------------------------------------------

#endif // MODELCHECKADVSFD_h
