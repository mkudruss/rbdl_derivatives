/** \file Mainpage.h 
 * \mainpage Rigid Body Dynamics Library
 *
 * This is the documentation of a yet to be named rigid body simulation
 * code. So far the code supports forward and inverse dynamics by using the
 * Articulated Body Algorathm and the Newton-Euler algorithm, respectively.
 * Additionally it also containes the Composite Rigid Body Algorithm that
 * computes the joint space inertia matrix.
 *
 * The code is written by <a
 * href="mailto:martin.felis@iwr.uni-heidelberg.de">Martin Felis
 * <martin.felis@iwr.uni-heidelberg.de></a> and heavily inspired by the
 * pseudo code of the book "Rigid Body Dynamics Algorithms" of <a
 * href="http://users.cecs.anu.edu.au/~roy/">Roy Featherstone</a>.
 * 
 * The library uses the configurable math library (which can be found here:
 * <a href="http://www.cmldev.net">http://www.cmldev.net</a>).
 *
 * Documentation of the functions can be found at the documentation page of
 * the namespace RigidBodyDynamics.
 *
 * \section Example An Example
 * 
 * Here is a simple example how one can create a meaningless model and
 * compute the forward dynamics for it:
 * 
 * \code
 *	#include <iostream>
 *
 *	#include <rbdl.h>
 *
 *	using namespace RigidBodyDynamics;
 *
 *	int main (int argc, char* argv[]) {
 *	Model* model = NULL;
 *
 *	unsigned int body_a_id, body_b_id, body_c_id;
 *	Body body_a, body_b, body_c;
 *	Joint joint_a, joint_b, joint_c;
 *
 *	model = new Model();
 *	model->Init();
 *
 *	model->gravity = Vector3d (0., -9.81, 0.);
 *
 *	body_a = Body (1., Vector3d (0.5, 0., 0.0), Vector3d (1., 1., 1.));
 *		joint_a = Joint(
 *		JointTypeRevolute,
 *		Vector3d (0., 0., 1.)
 *	);
 *	
 *	body_a_id = model->AddBody(0, Xtrans(Vector3d(0., 0., 0.)), joint_a, body_a);
 *	
 *	body_b = Body (1., Vector3d (0., 0.5, 0.), Vector3d (1., 1., 1.));
 *		joint_b = Joint (
 *		JointTypeRevolute,
 *		Vector3d (0., 0., 1.)
 *	);
 *	
 *	body_b_id = model->AddBody(body_a_id, Xtrans(Vector3d(1., 0., 0.)), joint_b, body_b);
 *	
 *	body_c = Body (0., Vector3d (0.5, 0., 0.), Vector3d (1., 1., 1.));
 *		joint_c = Joint (
 *		JointTypeRevolute,
 *		Vector3d (0., 0., 1.)
 *	);
 *	
 *	body_c_id = model->AddBody(body_b_id, Xtrans(Vector3d(0., 1., 0.)), joint_c, body_c);
 *
 *	VectorNd Q = VectorNd::Zero (model->dof_count);
 *	VectorNd QDot = VectorNd::Zero (model->dof_count);
 *	VectorNd Tau = VectorNd::Zero (model->dof_count);
 *	VectorNd QDDot = VectorNd::Zero (model->dof_count);
 *
 * 	ForwardDynamics (*model, Q, QDot, Tau, QDDot);
 *
 *	std::cout << QDDot.transpose() << std::endl;
 *
 *	delete model;
 *
 * 	return 0;
 *}
 * \endcode
 *
 * If the library itself is already created, one can compile this example
 * with:
 * \code
 * 	g++ example.cc -I<path to src folder> -lrbdl -L<path to librbdl.a> -o example
 * \endcode
 *
 * Additionally there is a CMakeLists.txt, that can be used to automatically
 * create the makefiles by using <a href="http://www.cmake.org">CMake</a>.
 * It uses the script FindRBDL.cmake which can be used to find the library
 * and include directory of the headers.
 *
 * The FindRBDL.cmake script can use the environment variables RBDL_PATH,
 * RBDL_INCLUDE_PATH, and RBDL_LIBRARY_PATH to find the required files.
 *
 * \section ModelConstruction Construction of Models
 *
 * The construction of models makes use of carefully designed constructors
 * of the classes Body and Joint to ease the process of creating bodies.
 * Adding bodies to the model is done by specifying the parent body by its
 * id, the transformation from the parent origin to the joint origin, the
 * joint specification as an object, and the body itself. These parameters
 * are then fed to the function RigidBodyDynamics::Model::AddBody().
 *
 * \todo [high] check impulse computation
 * \todo [low] use cml for the SpatialAlgebra quantities
 * \todo [low] incorporate GiNaC (http://www.ginac.de) to generate code
 * \todo [low] serialization of the model?
 *
 * \subsection ToDo_done Done:
 * <ul>
 *   <li>[med] add specification for the visualization to the model</li>
 *   <li>[med] get rid of the std::vector<> values in Model</li>
 * </ul>
 */
