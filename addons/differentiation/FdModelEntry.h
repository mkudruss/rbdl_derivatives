#ifndef FDMODELENTRY_h
#define FDMODELENTRY_h

#include "rbdl/rbdl.h"
#include "ModelAD.h"

namespace RigidBodyDynamics {

void computeFDEntry(
    Model const & model,
    Model const & modelh,
    double h,
    int idir,
    ADModel & fd_model);

}

#endif // FDMODELENTRY_h
