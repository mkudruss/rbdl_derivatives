#include <UnitTest++.h>

#include <iostream>

#include "rbdl/rbdl_mathutils.h"

#include "rbdl_mathutilsAD.h"
#include "rbdl_mathutilsFD.h"

#include "rbdl/Logging.h"

using namespace std;
using namespace RigidBodyDynamics;
using namespace RigidBodyDynamics::Math;

const double TEST_PREC = 1.0e-12;

// -----------------------------------------------------------------------------

TEST (crossm_v1v2_ADvsAnalyticTest) {
    unsigned int ndirs = 20;

    SpatialVector v1 = SpatialVector::Random();
    SpatialVector v2 = SpatialVector::Random();
    MatrixNd v1_dirs = MatrixNd::Random(6, ndirs);
    MatrixNd v2_dirs = MatrixNd::Random(6, ndirs);

    std::vector<SpatialVector> ad_res (ndirs, SpatialVector::Zero());
    std::vector<SpatialVector> an_res (ndirs, SpatialVector::Zero());

    for (unsigned int idirs = 0; idirs < ndirs; ++idirs) {
        SpatialVector v1_dir = v1_dirs.block(0, idirs, 6, 1);
        SpatialVector v2_dir = v2_dirs.block(0, idirs, 6, 1);

        ad_res[idirs] = AD::crossm (v1, v1_dir, v2, v2_dir);
        an_res[idirs] = crossm (v1_dir, v2) + crossm (v1, v2_dir);
        CHECK_ARRAY_CLOSE (
            an_res[idirs].data(), ad_res[idirs].data(), 6, TEST_PREC
        );
    }
}

TEST (crossm_v1v2_ADvsFDTest) {
    double TEST_PREC = 1e-07;
    unsigned int ndirs = 20;

    SpatialVector v1 = SpatialVector::Random();
    SpatialVector v2 = SpatialVector::Random();
    MatrixNd v1_dirs = MatrixNd::Random(6, ndirs);
    MatrixNd v2_dirs = MatrixNd::Random(6, ndirs);

    std::vector<SpatialVector> ad_res (ndirs, SpatialVector::Zero());
    std::vector<SpatialVector> fd_res (ndirs, SpatialVector::Zero());

    for (unsigned int idirs = 0; idirs < ndirs; ++idirs) {
        SpatialVector v1_dir = v1_dirs.block(0, idirs, 6, 1);
        SpatialVector v2_dir = v2_dirs.block(0, idirs, 6, 1);

        ad_res[idirs] = AD::crossm (v1, v1_dir, v2, v2_dir);
        fd_res[idirs] = FD::crossm (v1, v1_dir, v2, v2_dir);
        CHECK_ARRAY_CLOSE (
            fd_res[idirs].data(), ad_res[idirs].data(), 6, TEST_PREC
        );
    }
}

TEST (crossm_v_ADvsAnalyticTest) {
    double TEST_PREC = 1e-08;
    unsigned int ndirs = 20;

    SpatialVector v = SpatialVector::Random();
    MatrixNd v_dirs = MatrixNd::Random(6, ndirs);

    std::vector<SpatialMatrix> ad_res (ndirs, SpatialMatrix::Zero());
    std::vector<SpatialMatrix> an_res (ndirs, SpatialMatrix::Zero());

    for (unsigned int idirs = 0; idirs < ndirs; ++idirs) {
        SpatialVector v_dir = v_dirs.block(0, idirs, 6, 1);

        ad_res[idirs] = AD::crossm (v_dir);
        an_res[idirs] = crossm (v_dir);
        CHECK_ARRAY_CLOSE (an_res[idirs].data(), ad_res[idirs].data(), 6, TEST_PREC);
    }

}

TEST (crossm_v_ADvsFDTest) {
    double TEST_PREC = 1e-08;
    unsigned int ndirs = 20;

    SpatialVector v = SpatialVector::Random();
    MatrixNd v_dirs = MatrixNd::Random(6, ndirs);

    std::vector<SpatialMatrix> ad_res (ndirs, SpatialMatrix::Zero());
    std::vector<SpatialMatrix> fd_res (ndirs, SpatialMatrix::Zero());

    for (unsigned int idirs = 0; idirs < ndirs; ++idirs) {
        SpatialVector v_dir = v_dirs.block(0, idirs, 6, 1);

        ad_res[idirs] = AD::crossm (v_dir);
        fd_res[idirs] = FD::crossm (v, v_dir);
        CHECK_ARRAY_CLOSE (fd_res[idirs].data(), ad_res[idirs].data(), 6, TEST_PREC);
    }
}

// -----------------------------------------------------------------------------
