/*
 * RBDL - Rigid Body Dynamics Library
 * Copyright (c) 2011-2015 Martin Felis <martin.felis@iwr.uni-heidelberg.de>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#include "DynamicsFD.h"
#include "FdModelEntry.h"

using namespace RigidBodyDynamics::Math;

using std::vector;

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------
namespace FD {
// -----------------------------------------------------------------------------

RBDL_DLLAPI
void ForwardDynamics(
  Model & model,
  ADModel & fd_model,
	const VectorNd& q,
	const MatrixNd& q_dirs,
	const VectorNd& qdot,
	const MatrixNd& qdot_dirs,
	const VectorNd& tau,
	const MatrixNd& tau_dirs,
	VectorNd& qddot,
	MatrixNd& fd_qddot,
    vector<SpatialVector>* f_ext
) {
	assert(q_dirs.cols() == qdot_dirs.cols()
		&& q_dirs.cols() == fd_qddot.cols()
		&& tau_dirs.cols() == q_dirs.cols()
		&& "q_dirs, qdot_dirs, tau_dirs, fd_qddot have different dimensions");

	double h = sqrt(1e-16);
	unsigned int ndirs = q_dirs.cols();

	ForwardDynamics(model, q, qdot, tau, qddot,f_ext);

	VectorNd hd_qddot(qddot);
	VectorNd q_dir(q);
	VectorNd qdot_dir(qdot);
	VectorNd tau_dir(tau);

  for(unsigned idir = 0; idir < ndirs; ++idir) {
    q_dir = q_dirs.col(idir);
    qdot_dir = qdot_dirs.col(idir);
    tau_dir = tau_dirs.col(idir);
    Model modelh = model;
    ForwardDynamics(model, q + h*q_dir, qdot + h*qdot_dir, tau + h*tau_dir, hd_qddot,f_ext);
    fd_qddot.col(idir) = (hd_qddot - qddot) / h;
    computeFDEntry(model, modelh, h, idir, fd_model);
  }
}

RBDL_DLLAPI
void InverseDynamics(
	Model& model,
	const Math::VectorNd& q,
	const Math::MatrixNd& q_dirs,
	const Math::VectorNd& qdot,
	const Math::MatrixNd& qdot_dirs,
	const Math::VectorNd& qddot,
	const Math::MatrixNd& qddot_dirs,
	Math::VectorNd& tau,
	Math::MatrixNd& fd_tau,
	std::vector<Math::SpatialVector> *f_ext
) {

	assert(q_dirs.cols() == qdot_dirs.cols() &&
		q_dirs.cols() == qddot_dirs.cols() &&
		   "q_dirs, qdot_dirs, qddot_dirs have different dimensions");

	assert( fd_tau.cols() == q_dirs.cols() &&
			fd_tau.rows() == tau.rows() &&
			"fd_tau and tau have different dimensions");

	double h = sqrt(1e-16);

	unsigned int ndirs   = q_dirs.cols();

	InverseDynamics(model, q, qdot, qddot, tau, f_ext);

	VectorNd hd_tau(tau);
	VectorNd q_dir(q);
	VectorNd qdot_dir(qdot);
	VectorNd qddot_dir(qddot);

	for (unsigned int i = 0; i < ndirs; i++ ){
			q_dir = q_dirs.block(0,i, model.q_size, 1);
			qdot_dir = qdot_dirs.block(0,i, model.q_size, 1);
			qddot_dir = qddot_dirs.block(0,i, model.q_size, 1);
			InverseDynamics(model,         q + h*q_dir,
										qdot + h*qdot_dir,
									   qddot + h*qddot_dir,
												 hd_tau,  f_ext);

		fd_tau.block(0,i,hd_tau.rows(),1) = (hd_tau - tau) / h;
	}
}

RBDL_DLLAPI
void NonlinearEffects (
    Model & model,
    const Math::VectorNd & q,
    const Math::MatrixNd & q_dirs,
    const Math::VectorNd & qDot,
    const Math::MatrixNd & qDot_dirs,
    Math::VectorNd & tau,
    Math::MatrixNd & fd_tau
) {
  int ndirs = q_dirs.cols();
  assert(ndirs == qDot_dirs.cols());
  assert(ndirs == fd_tau.cols());

  NonlinearEffects (model, q, qDot, tau);
  double h = 1.0e-8;
  for (int idir = 0; idir < ndirs; idir++) {
    VectorNd q_dir  = q_dirs.col(idir);
    VectorNd qd_dir = qDot_dirs.col(idir);
    VectorNd hd_tau(tau.rows());
    NonlinearEffects (model, q + h * q_dir, qDot + h * qd_dir, hd_tau);
    fd_tau.col(idir) = (hd_tau - tau) / h;
  }
}

RBDL_DLLAPI
void CompositeRigidBodyAlgorithm (
	Model &model,
	const VectorNd &q,
	const MatrixNd &q_dirs,
	std::vector<MatrixNd> &fd_out
) {
	MatrixNd inertia_ref = MatrixNd::Zero(q.size(),q.size());
	MatrixNd inertia_fd = MatrixNd::Zero(q.size(),q.size());

	CompositeRigidBodyAlgorithm(model, q, inertia_ref, true);
	unsigned int ndirs = q_dirs.cols();
	double h = 1.0e-8;
	for (unsigned int j = 0; j < ndirs; j++) {
		VectorNd q_dir = q_dirs.block(0,j, model.qdot_size, 1);
		CompositeRigidBodyAlgorithm(model, q+h*q_dir, inertia_fd, true);
		fd_out.push_back( (inertia_fd - inertia_ref) / h);
	}
}

// -----------------------------------------------------------------------------
} /* FD */
// -----------------------------------------------------------------------------
} /* RigidBodyDynamics */
// -----------------------------------------------------------------------------
