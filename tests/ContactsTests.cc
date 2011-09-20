#include <UnitTest++.h>

#include <iostream>

#include "mathutils.h"
#include "Logging.h"

#include "Model.h"
#include "Contacts.h"
#include "Dynamics.h"
#include "Kinematics.h"

using namespace std;
using namespace SpatialAlgebra;
using namespace SpatialAlgebra::Operators;
using namespace RigidBodyDynamics;
using namespace RigidBodyDynamics::Experimental;

const double TEST_PREC = 1.0e-11;

struct ContactsFixture {
	ContactsFixture () {
		ClearLogOutput();
		model = new Model;
		model->Init();

		model->gravity = Vector3d  (0., -9.81, 0.);

		/* 
		 *
		 *          X Contact point (ref child)
		 *          |
		 *    Base  |
		 *   / body |
		 *  O-------*
		 *           \
		 *             Child body
		 */

		// base body (3 DoF)
		base_rot_z = Body (
				0.,
				Vector3d (0., 0., 0.),
				Vector3d (0., 0., 0.)
				);
		joint_base_rot_z = Joint (
				JointTypeRevolute,
				Vector3d (0., 0., 1.)
				);
		base_rot_z_id = model->AddBody (0, Xtrans (Vector3d (0., 0., 0.)), joint_base_rot_z, base_rot_z);

		base_rot_y = Body (
				0.,
				Vector3d (0., 0., 0.),
				Vector3d (0., 0., 0.)
				);
		joint_base_rot_y = Joint (
				JointTypeRevolute,
				Vector3d (0., 1., 0.)
				);
		base_rot_y_id = model->AddBody (base_rot_z_id, Xtrans (Vector3d (0., 0., 0.)), joint_base_rot_y, base_rot_y);

		base_rot_x = Body (
				1.,
				Vector3d (0.5, 0., 0.),
				Vector3d (1., 1., 1.)
				);
		joint_base_rot_x = Joint (
				JointTypeRevolute,
				Vector3d (1., 0., 0.)
				);
		base_rot_x_id = model->AddBody (base_rot_y_id, Xtrans (Vector3d (0., 0., 0.)), joint_base_rot_x, base_rot_x);

		// child body (3 DoF)
		child_rot_z = Body (
				0.,
				Vector3d (0., 0., 0.),
				Vector3d (0., 0., 0.)
				);
		joint_child_rot_z = Joint (
				JointTypeRevolute,
				Vector3d (0., 0., 1.)
				);
		child_rot_z_id = model->AddBody (base_rot_x_id, Xtrans (Vector3d (1., 0., 0.)), joint_child_rot_z, child_rot_z);

		child_rot_y = Body (
				0.,
				Vector3d (0., 0., 0.),
				Vector3d (0., 0., 0.)
				);
		joint_child_rot_y = Joint (
				JointTypeRevolute,
				Vector3d (0., 1., 0.)
				);
		child_rot_y_id = model->AddBody (child_rot_z_id, Xtrans (Vector3d (0., 0., 0.)), joint_child_rot_y, child_rot_y);

		child_rot_x = Body (
				1.,
				Vector3d (0., 0.5, 0.),
				Vector3d (1., 1., 1.)
				);
		joint_child_rot_x = Joint (
				JointTypeRevolute,
				Vector3d (1., 0., 0.)
				);
		child_rot_x_id = model->AddBody (child_rot_y_id, Xtrans (Vector3d (0., 0., 0.)), joint_child_rot_x, child_rot_x);

		Q = VectorNd::Constant (model->mBodies.size() - 1, 0.);
		QDot = VectorNd::Constant (model->mBodies.size() - 1, 0.);
		QDDot = VectorNd::Constant (model->mBodies.size() - 1, 0.);
		Tau = VectorNd::Constant (model->mBodies.size() - 1, 0.);

		contact_body_id = child_rot_x_id;
		contact_point = Vector3d  (1., 0., 0.);
		contact_normal = Vector3d  (0., 1., 0.);

		ClearLogOutput();
	}
	
	~ContactsFixture () {
		delete model;
	}
	Model *model;

	unsigned int base_rot_z_id, base_rot_y_id, base_rot_x_id,
		child_rot_z_id, child_rot_y_id, child_rot_x_id,
		base_body_id;

	Body base_rot_z, base_rot_y, base_rot_x,
		child_rot_z, child_rot_y, child_rot_x;

	Joint joint_base_rot_z, joint_base_rot_y, joint_base_rot_x,
		joint_child_rot_z, joint_child_rot_y, joint_child_rot_x;

	VectorNd Q;
	VectorNd QDot;
	VectorNd QDDot;
	VectorNd Tau;

	unsigned int contact_body_id;
	Vector3d contact_point;
	Vector3d contact_normal;
	std::vector<ContactInfo> contact_data;
};

TEST ( TestForwardDynamicsContactsLagrangianSimple ) {
	Model model;
	model.Init();
	model.gravity = Vector3d  (0., -9.81, 0.);
	Body base_body (1., Vector3d (0., 0., 0.), Vector3d (1., 1., 1.));
	unsigned int base_body_id = model.SetFloatingBaseBody(base_body);

	VectorNd Q = VectorNd::Constant ((size_t) model.dof_count, 0.);
	VectorNd QDot = VectorNd::Constant ((size_t) model.dof_count, 0.);
	VectorNd QDDot = VectorNd::Constant  ((size_t) model.dof_count, 0.);
	VectorNd Tau = VectorNd::Constant ((size_t) model.dof_count, 0.);

	Q[1] = 1.;
	QDot[0] = 1.;
	QDot[3] = -1.;

	unsigned int contact_body_id = base_body_id;
	Vector3d contact_point ( 0., -1., 0.);

	ContactInfo ground_x (contact_body_id, contact_point, Vector3d (1., 0., 0.));
	ContactInfo ground_y (contact_body_id, contact_point, Vector3d (0., 1., 0.));
	ContactInfo ground_z (contact_body_id, contact_point, Vector3d (0., 0., 1.));

	std::vector<ContactInfo> contact_data;

	contact_data.push_back (ground_x);
	contact_data.push_back (ground_y);
	contact_data.push_back (ground_z);

	ClearLogOutput();

	ForwardDynamicsContactsLagrangian (model, Q, QDot, Tau, contact_data, QDDot);

	Vector3d point_acceleration = CalcPointAcceleration (model, Q, QDot, QDDot, contact_body_id, contact_point);

	CHECK_ARRAY_CLOSE (
			Vector3d (0., 0., 0.).data(),
			point_acceleration.data(),
			3,
			TEST_PREC
			);

	// cout << "LagrangianSimple Logoutput Start" << endl;
	// cout << LogOutput.str() << endl;
	// cout << "LagrangianSimple Logoutput End" << endl;

	/*
	unsigned int i;
	for (i = 0; i < contact_data.size(); i++) {
		cout << "cf[" << i << "] = " << contact_data[i].force << endl;
	}

	cout << QDDot << endl;
	*/
}

TEST ( TestForwardDynamicsContactsLagrangianMoving ) {
	Model model;
	model.Init();
	model.gravity = Vector3d  (0., -9.81, 0.);
	Body base_body (1., Vector3d (0., 0., 0.), Vector3d (1., 1., 1.));
	unsigned int base_body_id = model.SetFloatingBaseBody(base_body);

	VectorNd Q = VectorNd::Constant ((size_t) model.dof_count, 0.);
	VectorNd QDot = VectorNd::Constant ((size_t) model.dof_count, 0.);
	VectorNd QDDot = VectorNd::Constant  ((size_t) model.dof_count, 0.);
	VectorNd Tau = VectorNd::Constant ((size_t) model.dof_count, 0.);

	Q[0] = 0.1;
	Q[1] = 0.2;
	Q[2] = 0.3;
	Q[3] = 0.4;
	Q[4] = 0.5;
	Q[5] = 0.6;
	QDot[0] = 1.1;
	QDot[1] = 1.2;
	QDot[2] = 1.3;
	QDot[3] = -1.4;
	QDot[4] = -1.5;
	QDot[5] = -1.6;

	unsigned int contact_body_id = base_body_id;
	Vector3d contact_point ( 0., -1., 0.);

	ContactInfo ground_x (contact_body_id, contact_point, Vector3d (1., 0., 0.));
	ContactInfo ground_y (contact_body_id, contact_point, Vector3d (0., 1., 0.));
	ContactInfo ground_z (contact_body_id, contact_point, Vector3d (0., 0., 1.));

	std::vector<ContactInfo> contact_data;

	contact_data.push_back (ground_x);
	contact_data.push_back (ground_y);
	contact_data.push_back (ground_z);

	ClearLogOutput();

	ForwardDynamicsContactsLagrangian (model, Q, QDot, Tau, contact_data, QDDot);

	Vector3d point_acceleration = CalcPointAcceleration (model, Q, QDot, QDDot, contact_body_id, contact_point);

	CHECK_ARRAY_CLOSE (
			Vector3d (0., 0., 0.).data(),
			point_acceleration.data(),
			3,
			TEST_PREC
			);

	// cout << "LagrangianSimple Logoutput Start" << endl;
	// cout << LogOutput.str() << endl;
	// cout << "LagrangianSimple Logoutput End" << endl;

	/*
	unsigned int i;
	for (i = 0; i < contact_data.size(); i++) {
		cout << "cf[" << i << "] = " << contact_data[i].force << endl;
	}

	cout << QDDot << endl;
	*/
}

TEST_FIXTURE (ContactsFixture, ForwardDynamicsContactsSingleContact) {
	contact_normal.set (0., 1., 0.);
	contact_data.push_back (ContactInfo(contact_body_id, contact_point, contact_normal, 0.));

	Vector3d point_accel_lagrangian;
	double contact_force_lagrangian;
	
	ForwardDynamicsContactsLagrangian (*model, Q, QDot, Tau, contact_data, QDDot);
	point_accel_lagrangian = CalcPointAcceleration (*model, Q, QDot, QDDot, contact_body_id, contact_point);
	contact_force_lagrangian = contact_data[0].force;
		
//	ClearLogOutput();
	ForwardDynamicsContacts (*model, Q, QDot, Tau, contact_data, QDDot);
//	cout << LogOutput.str() << endl;

	Vector3d point_accel_recursive;
	double contact_force_recursive;
	point_accel_recursive = CalcPointAcceleration (*model, Q, QDot, QDDot, contact_body_id, contact_point, true);
	contact_force_recursive = contact_data[0].force;

	CHECK_CLOSE (contact_force_lagrangian, contact_force_recursive, TEST_PREC);
	CHECK_CLOSE (0., contact_normal.dot(point_accel_recursive), TEST_PREC);
	CHECK_ARRAY_CLOSE (point_accel_lagrangian.data(), point_accel_recursive.data(), 3, TEST_PREC);
}

TEST_FIXTURE (ContactsFixture, ForwardDynamicsContactsSingleContactRotated) {
	contact_normal.set (0., 1., 0.);
	contact_data.push_back (ContactInfo(contact_body_id, contact_point, contact_normal, 0.));

	Q[0] = 0.6;
	Q[3] =   M_PI * 0.6;

	Vector3d point_accel_lagrangian;
	double contact_force_lagrangian;
	
	ForwardDynamicsContactsLagrangian (*model, Q, QDot, Tau, contact_data, QDDot);
	point_accel_lagrangian = CalcPointAcceleration (*model, Q, QDot, QDDot, contact_body_id, contact_point);
	contact_force_lagrangian = contact_data[0].force;
		
	ClearLogOutput();
	ForwardDynamicsContacts (*model, Q, QDot, Tau, contact_data, QDDot);
//	cout << LogOutput.str() << endl;

	Vector3d point_accel_recursive;
	double contact_force_recursive;
	point_accel_recursive = CalcPointAcceleration (*model, Q, QDot, QDDot, contact_body_id, contact_point);
	contact_force_recursive = contact_data[0].force;

	CHECK_CLOSE (0., contact_normal.dot(point_accel_recursive), TEST_PREC);
	CHECK_ARRAY_CLOSE (point_accel_lagrangian.data(), point_accel_recursive.data(), 3, TEST_PREC);
	CHECK_CLOSE (contact_force_lagrangian, contact_force_recursive, TEST_PREC);
}

TEST_FIXTURE (ContactsFixture, ForwardDynamicsContactsMultipleContact) {
	contact_data.push_back (ContactInfo(contact_body_id, contact_point, Vector3d (1., 0., 0.), 0.));
	contact_data.push_back (ContactInfo(contact_body_id, contact_point, Vector3d (0., 1., 0.), 0.));

	std::vector<ContactInfo> contact_data_lagrangian (contact_data);
	

	// we rotate the joints so that we have full mobility at the contact
	// point:
	//
	//  O       X (contact point)
	//   \     /
	//    \   /
	//     \ /
	//      *      
	//

	Q[0] = M_PI * 0.25;
	Q[3] = M_PI * 0.5;
	
	ClearLogOutput();
	ForwardDynamicsContacts (*model, Q, QDot, Tau, contact_data, QDDot);
//	cout << LogOutput.str() << endl;

	Vector3d point_accel_c = CalcPointAcceleration (*model, Q, QDot, QDDot, contact_body_id, contact_point);
//	cout << "point_accel_c = " << point_accel_c.transpose() << endl;

	ForwardDynamicsContactsLagrangian (*model, Q, QDot, Tau, contact_data_lagrangian, QDDot);
//	cout << "Lagrangian contact force " << contact_data_lagrangian[0].force << ", " << contact_data_lagrangian[1].force << endl;

	CHECK_CLOSE (contact_data_lagrangian[0].force, contact_data[0].force, TEST_PREC);
	CHECK_CLOSE (contact_data_lagrangian[1].force, contact_data[1].force, TEST_PREC);

	CHECK_CLOSE (0., point_accel_c[0], TEST_PREC);
	CHECK_CLOSE (0., point_accel_c[1], TEST_PREC);
}

TEST_FIXTURE (ContactsFixture, ForwardDynamicsContactsMultipleContactRotating) {
	contact_data.push_back (ContactInfo(contact_body_id, contact_point, Vector3d (1., 0., 0.), 0.));
	contact_data.push_back (ContactInfo(contact_body_id, contact_point, Vector3d (0., 1., 0.), 0.));

	std::vector<ContactInfo> contact_data_lagrangian (contact_data);
	

	// we rotate the joints so that we have full mobility at the contact
	// point:
	//
	//  O       X (contact point)
	//   \     /
	//    \   /
	//     \ /
	//      *      
	//

	Q[0] = M_PI * 0.25;
	Q[3] = M_PI * 0.5;
	QDot[0] = 2.;
	
	ClearLogOutput();
	ForwardDynamicsContacts (*model, Q, QDot, Tau, contact_data, QDDot);
//	cout << LogOutput.str() << endl;

	Vector3d point_accel_c = CalcPointAcceleration (*model, Q, QDot, QDDot, contact_body_id, contact_point);
//	cout << "point_accel_c = " << point_accel_c.transpose() << endl;

	ForwardDynamicsContactsLagrangian (*model, Q, QDot, Tau, contact_data_lagrangian, QDDot);
//	cout << "Lagrangian contact force " << contact_data_lagrangian[0].force << ", " << contact_data_lagrangian[1].force << endl;

	CHECK_CLOSE (contact_data_lagrangian[0].force, contact_data[0].force, TEST_PREC);
	CHECK_CLOSE (contact_data_lagrangian[1].force, contact_data[1].force, TEST_PREC);

	CHECK_CLOSE (0., point_accel_c[0], TEST_PREC);
	CHECK_CLOSE (0., point_accel_c[1], TEST_PREC);
}
