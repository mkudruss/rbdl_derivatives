#ifndef _MATHWRAPPER_H
#define _MATHWRAPPER_H

#define USE_EIGEN

#ifndef USE_EIGEN

	#include "cml/cml_config.h"
	#include "cml/cml.h"

	typedef cml::vector<double, cml::fixed<3> > Vector3d;
	typedef cml::matrix<double, cml::fixed<3,3> > Matrix3d;

	typedef cml::vector<double, cml::dynamic<> > VectorNd;
	typedef cml::matrix<double, cml::dynamic<> > MatrixNd;

	#include "SpatialAlgebra.h"

#else

	#define EIGEN_DEFAULT_TO_ROW_MAJOR
	#define EIGEN_MATRIX_PLUGIN "MatrixAddons.h"

	#include "Eigen/Dense"
	#include "Eigen/StdVector"

	typedef Eigen::Matrix< double, 3, 1> Vector3d;
	typedef Eigen::Matrix< double, 3, 3> Matrix3d;

	typedef Eigen::VectorXd VectorNd;
	typedef Eigen::MatrixXd MatrixNd;

	namespace SpatialAlgebra {
		typedef Eigen::Matrix< double, 6, 1> SpatialVector;
		typedef Eigen::Matrix< double, 6, 6> SpatialMatrix;
	}

	inline EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(SpatialAlgebra::SpatialVector)
	inline EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(SpatialAlgebra::SpatialMatrix)
#endif

#include "SpatialAlgebraOperators.h"

#endif /* _MATHWRAPPER_H */
