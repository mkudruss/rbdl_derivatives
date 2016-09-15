/*
 * RBDL - Rigid Body Dynamics Library
 * Copyright (c) 2011-2015 Martin Felis <martin.felis@iwr.uni-heidelberg.de>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#ifndef RBDL_MATH_AD_H
#define RBDL_MATH_AD_H

#include <cmath>
#include <limits>

#include <iostream>
#include <assert.h>

#include <rbdl/Model.h>
#include <rbdl/rbdl_mathutils.h>
#include <rbdl/SpatialAlgebraOperators.h>

#include "rbdl/Logging.h"

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------
namespace Math {
// -----------------------------------------------------------------------------
namespace AD {
// -----------------------------------------------------------------------------

RBDL_DLLAPI
inline Matrix3d E_from_Matrix(const SpatialMatrix X) {
    Matrix3d E = Matrix3d::Zero();
    E = X.block<3, 3>(0,0);
    return E;
}

RBDL_DLLAPI
inline Vector3d r_from_Matrix(const SpatialMatrix X) {
    Matrix3d E = E_from_Matrix(X);
    Matrix3d Erx = Matrix3d::Zero();
    Erx = X.block<3,3>(3,0);
    Matrix3d rx = E.transpose() * Erx;
    Vector3d r = Vector3d::Zero();
    r(0) = -rx(2,1);
    r(1) =  rx(2,0);
    r(2) = -rx(1,0);
    return r;
}

RBDL_DLLAPI
inline Vector3d r_from_Matrix(const SpatialMatrix X, const SpatialMatrix X_dirs) {
    Matrix3d E_dirs = E_from_Matrix(X_dirs);
    Matrix3d E = E_from_Matrix(X);

    Matrix3d Erx_dirs = Matrix3d::Zero();
    Matrix3d Erx = Matrix3d::Zero();
    Erx_dirs = X_dirs.block<3,3>(3,0);
    Erx = X.block<3,3>(3,0);

    Matrix3d rx_dirs = E_dirs.transpose() * Erx + E.transpose() * Erx_dirs;
    // NOTE: currently not used
    // Matrix3d rx = E.transpose() * Erx;

    Vector3d r_dirs = Vector3d::Zero();
    r_dirs(0) = -rx_dirs(2,1);
    r_dirs(1) =  rx_dirs(2,0);
    r_dirs(2) = -rx_dirs(1,0);

    return r_dirs;
}

RBDL_DLLAPI
inline SpatialMatrix Xrotx (double const & xrot, double const & xrot_dirs) {
	SpatialMatrix result (SpatialMatrix::Zero(6,6));

	double s, c;
	s = sin (xrot) * xrot_dirs;
	c = cos (xrot) * xrot_dirs;

	Matrix3d E(0., 0.,  0.,
			   0., -s,  c,
			   0., -c, -s);

	result.block<3,3>(0,0) = E;
	result.block<3,3>(3,3) = E;

	return result;
}

RBDL_DLLAPI
inline SpatialMatrix Xroty (const double &yrot, const double &yrot_dirs) {
	SpatialMatrix result (SpatialMatrix::Zero(6,6));

	double s, c;
	s = sin (yrot) * yrot_dirs;
	c = cos (yrot) * yrot_dirs;

	Matrix3d E(
				-s, 0., -c,
				0., 0., 0.,
				c, 0., -s
				);

	result.block<3,3>(0,0) = E;
	result.block<3,3>(3,3) = E;

	return result;
}

RBDL_DLLAPI
inline SpatialMatrix Xrotz (const double &zrot, const double &zrot_dirs) {
	SpatialMatrix result (SpatialMatrix::Zero(6,6));

	double s, c;

	s = sin (zrot) * zrot_dirs;
	c = cos (zrot) * zrot_dirs;
	Matrix3d E(-s,  c, 0.,
			   -c, -s, 0.,
			   0., 0., 0.);

	result.block<3,3>(0,0) = E;
	result.block<3,3>(3,3) = E;

	return result;
}

RBDL_DLLAPI
inline SpatialMatrix Xtrans (const Vector3d &trans, const Vector3d &trans_dirs) {
    SpatialMatrix result (SpatialMatrix::Zero(6,6));

    result(3,1) =  trans_dirs[2];
    result(3,2) = -trans_dirs[1];

    result(4,0) = -trans_dirs[2];
    result(4,2) =  trans_dirs[0];

    result(5,0) =  trans_dirs[1];
    result(5,1) = -trans_dirs[0];

    return result;
}

RBDL_DLLAPI
inline Matrix3d cross3d (const Vector3d &v) {
    return Matrix3d (
                0, -v[2],  v[1],
             v[2],     0, -v[0],
            -v[1],  v[0],     0
            );
}

// NOTE: ad_crossm is equal to crossm applied to the directions
RBDL_DLLAPI
inline SpatialMatrix crossm (const SpatialVector &v) {
    return SpatialMatrix (
                0, -v[2],  v[1],     0,     0,     0,
             v[2],     0, -v[0],     0,     0,     0,
            -v[1],  v[0],     0,     0,     0,     0,
                0, -v[5],  v[4],     0, -v[2],  v[1],
             v[5],     0, -v[3],  v[2],     0, -v[0],
            -v[4],  v[3],     0, -v[1],  v[0],     0
            );
}

RBDL_DLLAPI
inline SpatialVector crossm (
    const SpatialVector &v1, const SpatialVector &v1_dirs,
    const SpatialVector &v2, const SpatialVector &v2_dirs
) {
    return SpatialVector (
            // nominal -v1[2] * v2[1] + v1[1] * v2[2],
            -v1_dirs[2] * v2[1] - v1[2] * v2_dirs[1] + v1_dirs[1] * v2[2] + v1[1] * v2_dirs[2],
            // nominal v1[2] * v2[0] - v1[0] * v2[2],
            v1_dirs[2] * v2[0] + v1[2] * v2_dirs[0] - v1_dirs[0] * v2[2] - v1[0] * v2_dirs[2],
            // nominal -v1[1] * v2[0] + v1[0] * v2[1],
            -v1_dirs[1] * v2[0] - v1[1] * v2_dirs[0] + v1_dirs[0] * v2[1] + v1[0] * v2_dirs[1],
            // nominal -v1[5] * v2[1] + v1[4] * v2[2] - v1[2] * v2[4] + v1[1] * v2[5],
            -v1_dirs[5] * v2[1] - v1[5] * v2_dirs[1] + v1_dirs[4] * v2[2] + v1[4] * v2_dirs[2]
            - v1_dirs[2] * v2[4] - v1[2] * v2_dirs[4] + v1_dirs[1] * v2[5] + v1[1] * v2_dirs[5],
            // nominal v1[5] * v2[0] - v1[3] * v2[2] + v1[2] * v2[3] - v1[0] * v2[5],
            v1_dirs[5] * v2[0] + v1[5] * v2_dirs[0] - v1_dirs[3] * v2[2] - v1[3] * v2_dirs[2]
            + v1_dirs[2] * v2[3] + v1[2] * v2_dirs[3] - v1_dirs[0] * v2[5] - v1[0] * v2_dirs[5],
            // nominal -v1[4] * v2[0] + v1[3] * v2[1] - v1[1] * v2[3] + v1[0] * v2[4]
            -v1_dirs[4] * v2[0] - v1[4] * v2_dirs[0] + v1_dirs[3] * v2[1] + v1[3] * v2_dirs[1]
            - v1_dirs[1] * v2[3] - v1[1] * v2_dirs[3] + v1_dirs[0] * v2[4] + v1[0] * v2_dirs[4]
            );
}

// -----------------------------------------------------------------------------
} /* AD */
// -----------------------------------------------------------------------------
} /* Math */
// -----------------------------------------------------------------------------
} /* RigidBodyDynamics */
// -----------------------------------------------------------------------------

/* RBDL_MATH_AD_H */
#endif
