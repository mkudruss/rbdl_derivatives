/*
 * RBDL - Rigid Body Dynamics Library
 * Copyright (c) 2011-2015 Martin Felis <martin.felis@iwr.uni-heidelberg.de>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */


#include "DynamicsAD.h"

using namespace RigidBodyDynamics::Math;

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------
namespace AD {
// -----------------------------------------------------------------------------

RBDL_DLLAPI
void ForwardDynamics (
	Model& model,
	ADModel& ad_model,
	const VectorNd& q,
	const MatrixNd& q_dirs,
	const VectorNd& qdot,
	const MatrixNd& qdot_dirs,
	const VectorNd& tau,
	const MatrixNd& tau_dirs,
	VectorNd& qddot,
	MatrixNd& ad_qddot,
	std::vector<SpatialVector>* f_ext
) {
	SpatialVector spatial_gravity (0., 0., 0., model.gravity[0], model.gravity[1], model.gravity[2]);

	unsigned int ndirs = q_dirs.cols();
	ad_model.resize_directions(ndirs);

	unsigned int i = 0;

	// Reset the velocity of the root body
	model.v[0].setZero();

	for (i = 1; i < model.mBodies.size(); i++) {
		unsigned int lambda = model.lambda[i];

		ad_jcalc(model, ad_model, i, q, q_dirs, qdot, qdot_dirs);

		if (lambda != 0) {
			for(unsigned int j = 0; j < ndirs; j++) {
				ad_model.X_base[i][j] = ad_model.X_lambda[i][j] * model.X_base[lambda].toMatrix()
					+ model.X_lambda[i].toMatrix() * ad_model.X_base[i][j];
			}
			model.X_base[i] = model.X_lambda[i] * model.X_base[lambda];
		} else {
			for(unsigned int j = 0; j < ndirs; j++) {
				ad_model.X_base[i][j] = ad_model.X_lambda[i][j];
			}
			model.X_base[i] = model.X_lambda[i];
		}

		for(unsigned int j = 0; j < ndirs; j++) {
			ad_model.v[i][j] = ad_model.X_lambda[i][j] * model.v[lambda]
				+ model.X_lambda[i].apply(ad_model.v[lambda][j])
				+ ad_model.v_J[i][j];
		}
		model.v[i] = model.X_lambda[i].apply(model.v[lambda]) + model.v_J[i];

		for(unsigned int j = 0; j < ndirs; j++) {
			ad_model.c[i][j] = ad_model.c_J[i][j]
				+ crossm(ad_model.v[i][j],model.v_J[i])
				+ crossm(model.v[i],ad_model.v_J[i][j]);
		}
		model.c[i] = model.c_J[i] + crossm(model.v[i],model.v_J[i]);

		for(unsigned int j = 0; j < ndirs; j++) {
				ad_model.IA[i][j].setZero();
		}
		model.I[i].setSpatialMatrix(model.IA[i]);

		for(unsigned int j = 0; j < ndirs; j++) {
			ad_model.pA[i][j] = crossf(ad_model.v[i][j],model.I[i] * model.v[i])
				+ crossf(model.v[i],model.I[i] * ad_model.v[i][j]);;
		}
		model.pA[i] = crossf(model.v[i],model.I[i] * model.v[i]);

		if (f_ext != NULL && (*f_ext)[i] != SpatialVectorZero) {
			for(unsigned int j = 0; j < ndirs; j++) {
				SpatialMatrix ad_X_base_force(ad_model.X_base[i][j]);
				ad_X_base_force.block<3,3>(3,0) = Matrix3d::Zero();
				ad_X_base_force.block<3,3>(0,3) = ad_model.X_base[i][j].block<3,3>(3,0);
				ad_model.pA[i][j] -= ad_X_base_force * (*f_ext)[i];
			}
			model.pA[i] -= model.X_base[i].toMatrixAdjoint() * (*f_ext)[i];
		}
	}

	for (i = model.mBodies.size() - 1; i > 0; i--) {
		unsigned int q_index = model.mJoints[i].q_index;

		if (model.mJoints[i].mDoFCount == 3) {
			cerr << "Multi-dof not supported." << endl;
			abort();
//          model.multdof3_U[i] = model.IA[i] * model.multdof3_S[i];
// #ifdef EIGEN_CORE_H
//          model.multdof3_Dinv[i] = (model.multdof3_S[i].transpose() * model.multdof3_U[i]).inverse().eval();
// #else
//          model.multdof3_Dinv[i] = (model.multdof3_S[i].transpose() * model.multdof3_U[i]).inverse();
// #endif
//          Vector3d tau_temp (Tau[q_index], Tau[q_index + 1], Tau[q_index + 2]);
//
//          model.multdof3_u[i] = tau_temp - model.multdof3_S[i].transpose() * model.pA[i];
//
// //           LOG << "multdof3_u[" << i << "] = " << model.multdof3_u[i].transpose() << std::endl;
//          unsigned int lambda = model.lambda[i];
//          if (lambda != 0) {
//              SpatialMatrix Ia = model.IA[i] - model.multdof3_U[i] * model.multdof3_Dinv[i] * model.multdof3_U[i].transpose();
//              SpatialVector pa = model.pA[i] + Ia * model.c[i] + model.multdof3_U[i] * model.multdof3_Dinv[i] * model.multdof3_u[i];
// #ifdef EIGEN_CORE_H
//              model.IA[lambda].noalias() += model.X_lambda[i].toMatrixTranspose() * Ia * model.X_lambda[i].toMatrix();
//              model.pA[lambda].noalias() += model.X_lambda[i].applyTranspose(pa);
// #else
//              model.IA[lambda] += model.X_lambda[i].toMatrixTranspose() * Ia * model.X_lambda[i].tadoMatrix();
//              model.pA[lambda] += model.X_lambda[i].applyTranspose(pa);
// #endif
//              LOG << "pA[" << lambda << "] = " << model.pA[lambda].transpose() << std::endl;
//          }
		}
		else {
			for(unsigned int j = 0; j < ndirs; j++) {
				ad_model.U[i][j] = model.IA[i] * ad_model.S[i][j];
				//  + ad_model.IA[i][j] * model.S[i];
			}
			model.U[i] = model.IA[i] * model.S[i];

			for(unsigned int j = 0; j < ndirs; j++) {
				ad_model.d(i,j) = ad_model.S[i][j].dot(model.U[i]) + model.S[i].dot(ad_model.U[i][j]);
			}
			model.d[i] = model.S[i].dot(model.U[i]);

			for(unsigned int j = 0; j < ndirs; j++) {
				ad_model.u(i,j) = tau_dirs(q_index,j)
					- ad_model.S[i][j].dot(model.pA[i])
					- model.S[i].dot(ad_model.pA[i][j]);
			}
			model.u[i] = tau[q_index] - model.S[i].dot(model.pA[i]);

			unsigned int lambda = model.lambda[i];
			if (lambda != 0) {
				vector<SpatialMatrix> ad_Ia(ndirs, SpatialMatrix::Zero());
				for(unsigned int j = 0; j < ndirs; j++) {
					ad_Ia[j] = ad_model.IA[i][j]
						- ad_model.U[i][j] * (model.U[i] / model.d[i]).transpose()
						- model.U[i] * (ad_model.U[i][j] / model.d[i]).transpose()
						- model.U[i] * (model.U[i] * ad_model.d(i,j) / model.d[i]).transpose();
				}
				SpatialMatrix Ia = model.IA[i] - model.U[i] * (model.U[i] / model.d[i]).transpose();

				vector<SpatialVector> ad_pa(ndirs, SpatialVector::Zero());
				for(unsigned int j = 0; j < ndirs; j++) {
					ad_pa[j] = ad_model.pA[i][j]
						+ ad_Ia[j] * model.c[i]
						+ Ia * ad_model.c[i][j]
						+ ad_model.U[i][j] * model.u[i] / model.d[i]
						+ model.U[i] * ad_model.u(i,j) / model.d[i]
						+ model.U[i] * model.u[i] * ad_model.d(i,j) / model.d[i];
				}
				SpatialVector pa = model.pA[i] + Ia * model.c[i] + model.U[i] * model.u[i] / model.d[i];
#ifdef EIGEN_CORE_H
				for(unsigned int j = 0; j < ndirs; j++) {
					ad_model.IA[lambda][j].noalias() +=
						ad_model.X_lambda[i][j].transpose() * Ia * model.X_lambda[i].toMatrix()
						+ model.X_lambda[i].toMatrixTranspose() * ad_Ia[j] * model.X_lambda[i].toMatrix()
						+ model.X_lambda[i].toMatrixTranspose() * Ia * ad_model.X_lambda[i][j];
				}
				model.IA[lambda].noalias() += model.X_lambda[i].toMatrixTranspose() * Ia * model.X_lambda[i].toMatrix();

				for(unsigned int j = 0; j < ndirs; j++) {
					ad_model.pA[lambda][j].noalias() += ad_model.X_lambda[i][j].transpose() * pa
						+ model.X_lambda[i].applyTranspose(ad_pa[j]);
				}
				model.pA[lambda].noalias() += model.X_lambda[i].applyTranspose(pa);
#else
				cerr << "Simple math not yet supported." << endl;
				abort();
				//model.IA[lambda] += model.X_lambda[i].toMatrixTranspose() * Ia * model.X_lambda[i].toMatrix();
				//model.pA[lambda] += model.X_lambda[i].applyTranspose(pa);
#endif
			}
		}
	}

	model.a[0] = spatial_gravity * -1.;

	for (i = 1; i < model.mBodies.size(); i++) {
		unsigned int q_index = model.mJoints[i].q_index;
		unsigned int lambda = model.lambda[i];
		SpatialTransform X_lambda = model.X_lambda[i];

		for(unsigned int j = 0; j < ndirs; j++) {
			ad_model.a[i][j] = ad_model.X_lambda[i][j] * model.a[lambda]
				+ X_lambda.apply(ad_model.a[lambda][j])
				+ ad_model.c[i][j];
		}
		model.a[i] = X_lambda.apply(model.a[lambda]) + model.c[i];

		if (model.mJoints[i].mDoFCount == 3) {
			cerr << "Multi-dof joints not supported." << endl;
			abort();
//          Vector3d qdd_temp = model.multdof3_Dinv[i] * (model.multdof3_u[i] - model.multdof3_U[i].transpose() * model.a[i]);
//          qddot[q_index] = qdd_temp[0];
//          qddot[q_index + 1] = qdd_temp[1];
//          qddot[q_index + 2] = qdd_temp[2];
//          model.a[i] = model.a[i] + model.multdof3_S[i] * qdd_temp;
		}
		else {
			for(unsigned int j = 0; j < ndirs; j++) {
				ad_qddot(q_index,j) = (ad_model.d(i,j)/model.d[i]) * (model.u[i] - model.U[i].dot(model.a[i]))
					+ (1./model.d[i]) * (ad_model.u(i,j) - ad_model.U[i][j].dot(model.a[i]) - model.U[i].dot(ad_model.a[i][j]));
			}
			qddot[q_index] = (1./model.d[i]) * (model.u[i] - model.U[i].dot(model.a[i]));

			for(unsigned int j = 0; j < ndirs; j++) {
				ad_model.a[i][j] = ad_model.a[i][j]
					+ ad_model.S[i][j] * qddot[q_index]
					+ model.S[i] * ad_qddot(q_index,j);
			}
			model.a[i] = model.a[i] + model.S[i] * qddot[q_index];
		}
	}
}

RBDL_DLLAPI
void InverseDynamics(
	Model& model,
	ADModel &ad_model,
	const Math::VectorNd& q,
	const Math::MatrixNd& q_dirs,
	const Math::VectorNd& qdot,
	const Math::MatrixNd& qdot_dirs,
	const Math::VectorNd& qddot,
	const Math::MatrixNd& qddot_dirs,
	Math::VectorNd& tau,
	Math::MatrixNd& ad_tau,
	std::vector<Math::SpatialVector> *f_ext
){
	model.v[0].setZero();
	model.a[0].set (0., 0., 0.,
			-model.gravity[0],
			-model.gravity[1],
			-model.gravity[2]);

	unsigned int ndirs = q_dirs.cols();

	for (unsigned int i = 1; i < model.mBodies.size(); i++) {
		unsigned int lambda = model.lambda[i];
		unsigned int q_index = model.mJoints[i].q_index;

		// derivative evaluation
		ad_jcalc (model, ad_model, i, q, q_dirs, qdot, qdot_dirs);
		// nominal evaluation
		// NOTE joints are already calculated in ad_jcalc
		// jcalc (model, ad_model, i, q, q_dirs, qdot, qdot_dirs);

		// Done in ad_jcalc
		// FIXME this is right?
		// for(unsigned int j = 0; j < ndirs; j++) {
		//  ad_model.X_lambda[i][j] = ad_model.X_J[i][j] * model.X_T[i].toMatrix();
		// }

		if (lambda != 0) {
			// derivative evaluation
			for(unsigned int j = 0; j < ndirs; j++) {
				ad_model.X_base[i][j] = ad_model.X_lambda[i][j] * model.X_base[lambda].toMatrix()
					+ model.X_lambda[i].toMatrix() * ad_model.X_base[i][j];
			}
			// nominal evaluation
			model.X_base[i] = model.X_lambda[i] * model.X_base[lambda];
		} else {
			// derivative evaluation
			for(unsigned int j = 0; j < ndirs; j++) {
				ad_model.X_base[i][j] = ad_model.X_lambda[i][j];
			}
			// nominal evaluation
			model.X_base[i] = model.X_lambda[i];
		}

		// derivative evaluation
		for(unsigned int j = 0; j < ndirs; j++) {
			ad_model.v[i][j] = ad_model.X_lambda[i][j] * model.v[lambda]
				+ model.X_lambda[i].apply (ad_model.v[lambda][j])
				+ ad_model.v_J[i][j];
		}
		// nominal evaluation
		model.v[i] = model.X_lambda[i].apply(model.v[lambda]) + model.v_J[i];

		// derivative evaluation
		for(unsigned int j = 0; j < ndirs; j++) {
			ad_model.c[i][j] = ad_model.c_J[i][j]
				+ crossm(ad_model.v[i][j],model.v_J[i])
				+ crossm(model.v[i],ad_model.v_J[i][j]);
		}
		// nominal evaluation
		model.c[i] = model.c_J[i] + crossm(model.v[i],model.v_J[i]);

		if (model.mJoints[i].mDoFCount == 3) {
			cerr << "Multi-dof not supported." << endl;
			abort();
			// nominal evaluation
			model.a[i] = model.X_lambda[i].apply(model.a[lambda])
				+ model.c[i]
				+ model.multdof3_S[i] * Vector3d (qddot[q_index], qddot[q_index + 1], qddot[q_index + 2]);
		} else {
			// derivative evaluation
			for(unsigned int j = 0; j < ndirs; j++) {
				ad_model.a[i][j] = ad_model.X_lambda[i][j] * model.a[lambda] + model.X_lambda[i].apply(ad_model.a[lambda][j])
					+ ad_model.c[i][j]
					+ ad_model.S[i][j] * qddot[q_index] + model.S[i] * qddot_dirs(q_index,j);
			}
			// nominal evaluation
			model.a[i] = model.X_lambda[i].apply(model.a[lambda]) + model.c[i] + model.S[i] * qddot[q_index];
		}

		if (!model.mBodies[i].mIsVirtual) {
			// derivative evaluation
			for(unsigned int j = 0; j < ndirs; j++) {
				ad_model.f[i][j] = model.I[i] * ad_model.a[i][j]
					+ crossf(ad_model.v[i][j],model.I[i] * model.v[i])
					+ crossf(model.v[i],model.I[i] * ad_model.v[i][j]);
			}
			// nominal evaluation
			model.f[i] = model.I[i] * model.a[i] + crossf(model.v[i],model.I[i] * model.v[i]);
		} else {
			// derivative evaluation
			for(unsigned int j = 0; j < ndirs; j++) {
				ad_model.f[i][j].setZero();
			}
			// nominal evaluation
			model.f[i].setZero();
		}

		if (f_ext != NULL && (*f_ext)[i] != SpatialVectorZero) {
			// derivative evaluation
			for(unsigned int j = 0; j < ndirs; j++) {
				SpatialMatrix ad_X_base_force(ad_model.X_base[i][j]);
				ad_X_base_force.block<3,3>(3,0) = Matrix3d::Zero();
				ad_X_base_force.block<3,3>(0,3) = ad_model.X_base[i][j].block<3,3>(3,0);
				ad_model.f[i][j] -= ad_X_base_force * (*f_ext)[i];
			}
			// nominal evaluation
			model.f[i] -= model.X_base[i].toMatrixAdjoint() * (*f_ext)[i];
		}
	}

	for (unsigned int i = model.mBodies.size() - 1; i > 0; i--) {
		if (model.mJoints[i].mDoFCount == 3) {
			cerr << "Multi-dof not supported." << endl;
			abort();
			// nominal evaluation
			tau.block<3,1>(model.mJoints[i].q_index, 0) = model.multdof3_S[i].transpose() * model.f[i];
		} else {
			// derivative evaluation
			for(unsigned int j = 0; j < ndirs; j++) {
				ad_tau(model.mJoints[i].q_index,j) = ad_model.S[i][j].dot(model.f[i])
					+ model.S[i].dot(ad_model.f[i][j]);
			}
			// nominal evaluation
			tau[model.mJoints[i].q_index] = model.S[i].dot(model.f[i]);
		}

		if (model.lambda[i] != 0) {
			// derivative evaluation
			for(unsigned int j = 0; j < ndirs; j++) {
				ad_model.f[model.lambda[i]][j] = ad_model.f[model.lambda[i]][j]
					+ ad_model.X_lambda[i][j].transpose() * model.f[i]
					+ model.X_lambda[i].applyTranspose(ad_model.f[i][j]);
			}
			// nominal evaluation
			model.f[model.lambda[i]] = model.f[model.lambda[i]]
				+ model.X_lambda[i].applyTranspose(model.f[i]);
		}
	}
}

RBDL_DLLAPI
void CompositeRigidBodyAlgorithm (
	Model &model,
	ADModel &ad_model,
	const VectorNd &q,
	const MatrixNd &q_dirs,
	MatrixNd &H,
	std::vector<MatrixNd> &H_ad,
	bool update_kinematics
) {
	LOG << "-------- " << __func__ << " --------" << std::endl;

	assert (H.rows() == model.dof_count && H.cols() == model.dof_count);

	size_t ndirs = q_dirs.cols();
	ad_model.resize_directions(ndirs);

	assert (H.rows() == model.dof_count && H.cols() == model.dof_count);

	for (unsigned int i = 1; i < model.mBodies.size(); i++) {
		if (update_kinematics) {
			// derivative evaluation
			ad_jcalc_X_lambda_S (model, ad_model, i, q, q_dirs);
			// nominal evaluation
			// NOTE nominal evaluation is part of ad_jcalc_X_lambda_S
			// jcalc_X_lambda_S (model, i, q);
		}
		// derivative evaluation
		for (size_t idir = 0; idir < ndirs; idir++) {
			// NOTE Ic is a constant transformation wrt. q
			ad_model.Ic[i][idir] = SpatialMatrix::Zero();
		}
		// nominal evaluation
		model.Ic[i] = model.I[i];
	}

	for (unsigned int i = model.mBodies.size() - 1; i > 0; i--) {
		if (model.lambda[i] != 0) {
			// derivative evaluation
			for(size_t idir = 0; idir < ndirs; idir++) {
				ad_model.Ic[model.lambda[i]][idir] = ad_model.Ic[model.lambda[i]][idir]
					+ model.X_lambda[i].toMatrixTranspose() * ad_model.Ic[i][idir]*model.X_lambda[i].toMatrix()
					+ ad_model.X_lambda[i][idir].transpose() * model.Ic[i].toMatrix()*model.X_lambda[i].toMatrix()
					+ model.X_lambda[i].toMatrixTranspose() * model.Ic[i].toMatrix()*ad_model.X_lambda[i][idir];
			}
			// nominal evaluation
			model.Ic[model.lambda[i]] = model.Ic[model.lambda[i]] + model.X_lambda[i].applyTranspose(model.Ic[i]);
		}

		unsigned int dof_index_i = model.mJoints[i].q_index;

		assert( model.mJoints[i].mDoFCount == 1);
		// if (model.mJoints[i].mDoFCount == 3) {
		//  Matrix63 F_63 = model.Ic[i].toMatrix() * model.multdof3_S[i];
		//  H.block<3,3>(dof_index_i, dof_index_i) = model.multdof3_S[i].toMatrixTranspose() * F_63;

		//  unsigned int j = i;
		//  unsigned int dof_index_j = dof_index_i;

		//  while (model.lambda[j] != 0) {
		//      F_63 = model.X_lambda[j].toMatrixTranspose() * (F_63);
		//      j = model.lambda[j];
		//      dof_index_j = model.mJoints[j].q_index;

		//      if (model.mJoints[j].mDoFCount == 3) {
		//          Matrix3d H_temp2 = F_63.transpose() * (model.multdof3_S[j]);

		//          H.block<3,3>(dof_index_i,dof_index_j) = H_temp2;
		//          H.block<3,3>(dof_index_j,dof_index_i) = H_temp2.transpose();
		//      } else {
		//          Vector3d H_temp2 = F_63.transpose() * (model.S[j]);

		//          H.block<3,1>(dof_index_i,dof_index_j) = H_temp2;
		//          H.block<1,3>(dof_index_j,dof_index_i) = H_temp2.transpose();
		//      }
		//  }
		// } else {

		// derivative evaluation
		for (size_t j = 0; j < ndirs; j++) {
			ad_model.F[i][j] = ad_model.Ic[i][j]*model.S[i]
				+ model.Ic[i].toMatrix()*ad_model.S[i][j];
		}
		// nominal evaluation
		SpatialVector F = model.Ic[i] * model.S[i];

		// derivative evaluation
		for (size_t j = 0; j < ndirs; j++) {
			H_ad[j](dof_index_i, dof_index_i) = ad_model.S[i][j].dot(F)
				+ model.S[i].dot(ad_model.F[i][j]);
		}
		// nominal evaluation
		H(dof_index_i, dof_index_i) = model.S[i].dot(F);

		unsigned int j = i;
		unsigned int dof_index_j = dof_index_i;

		while (model.lambda[j] != 0) {
			// derivative evaluation
			for (size_t ndir = 0; ndir < ndirs; ndir++) {
				ad_model.F[i][ndir] = ad_model.X_lambda[j][ndir].transpose()*F
					+ model.X_lambda[j].toMatrixTranspose()*ad_model.F[i][ndir];
			}
			// nominal evaluation
			F = model.X_lambda[j].applyTranspose(F);
			j = model.lambda[j];
			dof_index_j = model.mJoints[j].q_index;

			assert(model.mJoints[j].mDoFCount == 1);
			// if (model.mJoints[j].mDoFCount == 3) {
			//   Vector3d H_temp2 = (F.transpose() * model.multdof3_S[j]).transpose();
			//   H.block<1,3>(dof_index_i,dof_index_j) = H_temp2.transpose();
			//   H.block<3,1>(dof_index_j,dof_index_i) = H_temp2;
			// } else {
				// derivative evaluation
				for (size_t ndir = 0; ndir < ndirs; ndir++) {
					H_ad[ndir](dof_index_i,dof_index_j) =
						ad_model.F[i][ndir].dot(model.S[j])
						+ F.dot(ad_model.S[j][ndir]);
					H_ad[ndir](dof_index_j,dof_index_i) = H_ad[ndir](dof_index_i,dof_index_j);
				}
				// nominal evaluation
				H(dof_index_i,dof_index_j) = F.dot(model.S[j]);
				H(dof_index_j,dof_index_i) = H(dof_index_i,dof_index_j);
			// }
		}
	}
}

// -----------------------------------------------------------------------------
} /* AD */
// -----------------------------------------------------------------------------
} /* RigidBodyDynamics */
// -----------------------------------------------------------------------------
