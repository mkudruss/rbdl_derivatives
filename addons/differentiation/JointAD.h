#ifndef RBDL_JOINT_AD_h
#define RBDL_JOINT_AD_h

#include <rbdl/Model.h>

#include "ModelAD.h"


namespace RigidBodyDynamics {
	RBDL_DLLAPI
		void ad_jcalc (
			Model &model,
			ADModel &ad_model,
			unsigned int joint_id,
			const VectorNd &q,
			const MatrixNd &q_dirs,
			const VectorNd &qdot,
			const MatrixNd &qdot_dirs
		);

	RBDL_DLLAPI
	void ad_jcalc_X_lambda_S (
		Model &model,
		ADModel &ad_model,
		unsigned int joint_id,
		const VectorNd &q,
		const VectorNd &q_dirs
		);
}

#endif
