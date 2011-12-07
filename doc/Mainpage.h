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
 * The library comes with version 3 of the the
 * <a href="http://eigen.tuxfamily.org/">Eigen</a> math library. More
 * information about it can be found here:
 * <a href="http://eigen.tuxfamily.org/">http://eigen.tuxfamily.org/</a>.
 *
 * Documentation of the functions can be found at the documentation page of
 * the namespace RigidBodyDynamics.
 *
 * \section ModelConstruction Construction of Models
 *
 * The construction of \link RigidBodyDynamics::Model Models \endlink makes
 * use of carefully designed constructors of the classes \link
 * RigidBodyDynamics::Body Body \endlink and \link RigidBodyDynamics::Joint
 * Joint \endlink to ease the process of creating bodies.  Adding bodies to
 * the model is done by specifying the parent body by its id, the
 * transformation from the parent origin to the joint origin, the joint
 * specification as an object, and the body itself. These parameters are
 * then fed to the function RigidBodyDynamics::Model::AddBody().
 *
 * A simple example can be found \ref SimpleExample "here".
 * 
 * \section ModuleOverview API reference separated by functional modules
 * 
 * \li \ref kinematics_group
 * \li \ref dynamics_group
 * \li \ref contacts_group
 *
 * \section Licensing Licensing
 *
 * The library is published under the very permissive zlib free software
 * license which should allow you to use the software wherever you need.
 * Here is the full license text:
 * \verbatim
RBDL - Rigid Body Dynamics Library
Copyright (c) 2011 Martin Felis <martin.felis@iwr.uni-heidelberg.de>

This software is provided 'as-is', without any express or implied
warranty. In no event will the authors be held liable for any damages
arising from the use of this software.

Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it
freely, subject to the following restrictions:

   1. The origin of this software must not be misrepresented; you must not
   claim that you wrote the original software. If you use this software
   in a product, an acknowledgment in the product documentation would be
   appreciated but is not required.

   2. Altered source versions must be plainly marked as such, and must not be
   misrepresented as being the original software.

   3. This notice may not be removed or altered from any source
   distribution.
\endverbatim
 *
 * \todo [low] incorporate GiNaC (http://www.ginac.de) to generate code (some work done in 'ginac' branch)
 * \todo [low] serialization of the model?
 *
 * \subsection ToDo_done Done:
 * <ul>
 *   <li>[high] check impulse computation 2011-06-07</li>
 *   <li>[med] add specification for the visualization to the model</li>
 *   <li>[med] get replace std::vector<> vectors with proper math vectors in Model</li>
 * </ul>
 */
