#include "ContactsAD.h"

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------
namespace AD {
// -----------------------------------------------------------------------------

RBDL_DLLAPI
void CalcContactJacobian(
        Model &model,
        ADModel &ad_model,
        const Math::VectorNd &Q,
        const Math::VectorNd &Q_dirs,
        const ConstraintSet &CS,
        ADConstraintSet &ad_CS,
        Math::MatrixNd &G,
        std::vector<Math::MatrixNd> &G_dirs,
        bool update_kinematics
) {
    unsigned int ndirs = Q_dirs.cols();
    ad_model.resize_directions(ndirs);

    if (update_kinematics) {
        // derivative evaluation
        // UpdateKinematicsCustom (
        //     model, ad_model,
        //     &Q, &Q_dirs,
        //     NULL, NULL,
        //     NULL, NULL
        // );
        // nominal evaluation
        // NOTE kinematics are already updated in the AD version
        // UpdateKinematicsCustom (model, &Q, NULL, NULL);
    }

    unsigned int i,j;

    // variables to check whether we need to recompute G
    unsigned int prev_body_id = 0;
    Math::Vector3d prev_body_point = Math::Vector3d::Zero();
    // derivative evaluation
    std::vector<Math::MatrixNd> Gi_dirs (
        ndirs, Math::MatrixNd (3, model.dof_count)
    );
    // nominal evaluation
    Math::MatrixNd Gi (3, model.dof_count);

    for (i = 0; i < CS.size(); i++) {
        // only compute the matrix Gi if actually needed
        if (prev_body_id != CS.body[i] || prev_body_point != CS.point[i]) {
            Gi.setZero();
            // derivative evaluation
            AD::CalcPointJacobian(
                model, ad_model,
                Q, Q_dirs,
                CS.body[i], CS.point[i],
                Gi, Gi_dirs,
                false
            );
            // nominal evaluation
            // NOTE nominal Jacobian is already computed by the AD version
            // CalcPointJacobian (model, Q, CS.body[i], CS.point[i], Gi, false);

            prev_body_id = CS.body[i];
            prev_body_point = CS.point[i];
        }

        for (j = 0; j < model.dof_count; j++) {
            // derivative evaluation
            for (unsigned int idirs = 0; idirs < ndirs; idirs++) {
                Math::Vector3d gaxis_dirs (
                    Gi_dirs[idirs](0,j),
                    Gi_dirs[idirs](1,j),
                    Gi_dirs[idirs](2,j)
                );
                G_dirs[idirs](i,j) = gaxis_dirs.transpose() * CS.normal[i];
            }
            // nominal evaluation
            Math::Vector3d gaxis (Gi(0,j), Gi(1,j), Gi(2,j));
            G(i,j) = gaxis.transpose() * CS.normal[i];
        }
    }
}

// -----------------------------------------------------------------------------
} // namespace AD
// -----------------------------------------------------------------------------
} // namespace RigidBodyDynamics
// -----------------------------------------------------------------------------
