#ifndef _CMLWRAPPER_H
#define _CMLWRAPPER_H

#include "cml/cml_config.h"
#include "cml/cml.h"

typedef cml::vector<double, cml::dynamic<> > cmlVector;
typedef cml::matrix<double, cml::dynamic<> > cmlMatrix;

typedef cml::vector<double, cml::fixed<3> > Vector3d;
typedef cml::matrix<double, cml::fixed<3,3> > Matrix3d;

typedef cml::vector<double, cml::dynamic<> > VectorNd;
typedef cml::matrix<double, cml::dynamic<> > MatrixNd;

typedef cml::quaternion<double, cml::fixed<>, cml::vector_first, cml::positive_cross> Quaternion;

#include "SpatialAlgebra.h"

/*
typedef SpatialAlgebra::SpatialVector SpatialVector;
typedef SpatialAlgebra::SpatialMatrix SpatialMatrix;
*/

//typedef cml::vector<double, cml::fixed<6> > SpatialVector;
//typedef cml::matrix<double, cml::fixed<6,6> > SpatialMatrix;

#endif /* _CMLWRAPPER_H */
