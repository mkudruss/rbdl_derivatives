#include "ContactsFD.h"

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------
namespace FD {
// -----------------------------------------------------------------------------

/*
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
) {
    unsigned int ndirs = q_dirs.cols();
    ad_model.resize_directions(ndirs);

    if (update_kinematics) {
        // derivative evaluation
        UpdateKinematicsCustom (model, ad_model, &Q, &Q_dirs, NULL, NULL);
        // nominal evaluation
        // NOTE kinematics are already updated in the AD version
        // UpdateKinematicsCustom (model, &Q, NULL, NULL);
    }

    unsigned int i,j;

    // variables to check whether we need to recompute G
    unsigned int prev_body_id = 0;
    Vector3d prev_body_point = Vector3d::Zero();
    MatrixNd Gi (3, model.dof_count);

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
            Vector3d gaxis (Gi(0,j), Gi(1,j), Gi(2,j));
            // nominal evaluation
            G(i,j) = gaxis.transpose() * CS.normal[i];
        }
    }
}

*/

// -----------------------------------------------------------------------------
} // namespace FD
// -----------------------------------------------------------------------------
} // namespace RigidBodyDynamics
// -----------------------------------------------------------------------------
