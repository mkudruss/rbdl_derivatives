#ifndef _BODY_H
#define _BODY_H

#include <mathwrapper.h>
#include <assert.h>
#include <iostream>
#include "Logging.h"

namespace RigidBodyDynamics {

/** \brief Describes all properties of a single body 
 *
 * A Body contains information about mass, the location of its center of
 * mass, and the ineria tensor in the center of mass. This class is
 * designed to use the given information and transform it such that it can
 * directly be used by the spatial algebra.
 */
struct Body {
	Body() :
		mMass (1.),
		mCenterOfMass (0., 0., 0.),
		mSpatialInertia (
				0., 0., 0., 0., 0., 0.,	
				0., 0., 0., 0., 0., 0.,	
				0., 0., 0., 0., 0., 0.,	
				0., 0., 0., 0., 0., 0.,	
				0., 0., 0., 0., 0., 0.,	
				0., 0., 0., 0., 0., 0.
				)
	{ };
	Body(const Body &body) :
		mMass (body.mMass),
		mCenterOfMass (body.mCenterOfMass),
		mSpatialInertia (body.mSpatialInertia)
	{};
	Body& operator= (const Body &body) {
		if (this != &body) {
			mSpatialInertia = body.mSpatialInertia;
			mCenterOfMass = body.mCenterOfMass;
			mMass = body.mMass;
		}

		return *this;
	}
	/** \brief Constructs a body from mass, center of mass and radii of gyration 
	 *
	 * This constructor eases the construction of a new body as all the required
	 * parameters can simply be specified as parameters to the constructor.
	 * These are then used to generate the spatial inertia matrix.
	 *
	 * \param mass the mass of the body
	 * \param com  the position of the center of mass in the bodies coordinates
	 * \param gyration_radii the radii of gyration at the center of mass of the body
	 */
	Body(const double &mass,
			const Vector3d &com,
			const Vector3d &gyration_radii) :
		mMass (mass),
		mCenterOfMass(com) {
			Matrix3d com_cross (
					0., -com[2],  com[1],
					com[2],      0., -com[0],
					-com[1],  com[0],      0.
					);
			Matrix3d parallel_axis;
			parallel_axis = mass * com_cross * cml::transpose(com_cross);

			LOG << "parrallel axis = " << parallel_axis << std::endl;

			Vector3d gr (gyration_radii);
			Matrix3d pa (parallel_axis);
			Matrix3d mcc = mass * com_cross;
			Matrix3d mccT = transpose(mcc);

			mSpatialInertia.set (
					gr[0] + pa(0, 0), pa(0, 1), pa(0, 2), mcc(0, 0), mcc(0, 1), mcc(0, 2),
					pa(1, 0), gr[1] + pa(1, 1), pa(1, 2), mcc(1, 0), mcc(1, 1), mcc(1, 2),
					pa(2, 0), pa(2, 1), gr[2] + pa(2, 2), mcc(2, 0), mcc(2, 1), mcc(2, 2),
					mccT(0, 0), mccT(0, 1), mccT(0, 2), mass, 0., 0.,
					mccT(1, 0), mccT(1, 1), mccT(1, 2), 0., mass, 0.,
					mccT(2, 0), mccT(2, 1), mccT(2, 2), 0., 0., mass
					);
		}

	/** \brief Constructs a body from mass, center of mass, and a 3x3 inertia matrix 
	 *
	 * This constructor eases the construction of a new body as all the required
	 * parameters can simply be specified as parameters to the constructor.
	 * These are then used to generate the spatial inertia matrix.
	 *
	 * \param mass the mass of the body
	 * \param com  the position of the center of mass in the bodies coordinates
	 * \param inertia_C the inertia at the center of mass
	 */
	Body(const double &mass,
			const Vector3d &com,
			const Matrix3d &inertia_C) :
		mMass (mass),
		mCenterOfMass(com) {
			Matrix3d com_cross (
					0., -com[2],  com[1],
					com[2],      0., -com[0],
					-com[1],  com[0],      0.
					);
			Matrix3d parallel_axis;
			parallel_axis = mass * com_cross * cml::transpose(com_cross);

			LOG << "parrallel axis = " << parallel_axis << std::endl;

			Matrix3d pa (parallel_axis);
			Matrix3d mcc = mass * com_cross;
			Matrix3d mccT = transpose(mcc);

			mSpatialInertia.set (
					inertia_C(0,0) + pa(0, 0), inertia_C(0,1) + pa(0, 1), inertia_C(0,2) + pa(0, 2), mcc(0, 0), mcc(0, 1), mcc(0, 2),
					inertia_C(1,0) + pa(1, 0), inertia_C(1,1) + pa(1, 1), inertia_C(1,2) + pa(1, 2), mcc(1, 0), mcc(1, 1), mcc(1, 2),
					inertia_C(2,0) + pa(2, 0), inertia_C(2,1) + pa(2, 1), inertia_C(2,2) + pa(2, 2), mcc(2, 0), mcc(2, 1), mcc(2, 2),
					mccT(0, 0), mccT(0, 1), mccT(0, 2), mass, 0., 0.,
					mccT(1, 0), mccT(1, 1), mccT(1, 2), 0., mass, 0.,
					mccT(2, 0), mccT(2, 1), mccT(2, 2), 0., 0., mass
					);
		}

	~Body() {};

	/// \brief The mass of the body
	double mMass;
	/// \brief The position of the center of mass in body coordinates
	Vector3d mCenterOfMass;
	/// \brief The spatial inertia that contains both mass and inertia information
	SpatialAlgebra::SpatialMatrix mSpatialInertia;
};

}

#endif /* _BODY_H */
