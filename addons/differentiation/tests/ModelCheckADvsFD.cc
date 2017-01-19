#include <UnitTest++.h>

#include "ModelCheckADvsFD.h"

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------

void checkModelsADvsFD(
    unsigned ndirs,
    Model const & ad_model,
    ADModel const & ad_d_model,
    Model const & fd_model,
    ADModel const & fd_d_model) {

  CHECK_EQUAL(ad_model.a.size(), fd_model.a.size());
  for (unsigned i = 0; i < ad_model.a.size(); i++) {
    CHECK_ARRAY_CLOSE(ad_model.a[i].data(), fd_model.a[i].data(), 6, 1e-7);
  }

  CHECK_EQUAL(ad_model.v.size(), fd_model.v.size());
  for (unsigned i = 0; i < ad_model.v.size(); i++) {
    CHECK_ARRAY_EQUAL(ad_model.v[i].data(), fd_model.v[i].data(), 6);
  }

  CHECK_EQUAL(ad_model.v_J.size(), fd_model.v_J.size());
  for (unsigned i = 0; i < ad_model.v_J.size(); i++) {
    CHECK_ARRAY_EQUAL(ad_model.v_J[i].data(), fd_model.v_J[i].data(), 6);
  }

  CHECK_EQUAL(ad_model.c.size(), fd_model.c.size());
  for (unsigned i = 0; i < ad_model.c.size(); i++) {
    CHECK_ARRAY_EQUAL(ad_model.c[i].data(), fd_model.c[i].data(), 6);
  }

  CHECK_EQUAL(ad_model.c_J.size(), fd_model.c_J.size());
  for (unsigned i = 0; i < ad_model.c_J.size(); i++) {
    CHECK_ARRAY_EQUAL(ad_model.c_J[i].data(), fd_model.c_J[i].data(), 6);
  }

  CHECK_EQUAL(ad_model.X_lambda.size(), fd_model.X_lambda.size());
  for (unsigned i = 0; i < ad_model.X_lambda.size(); i++) {
    CHECK_ARRAY_EQUAL(ad_model.X_lambda[i].E.data(), fd_model.X_lambda[i].E.data(), 9);
    CHECK_ARRAY_EQUAL(ad_model.X_lambda[i].r.data(), fd_model.X_lambda[i].r.data(), 3);
  }

  CHECK_EQUAL(ad_model.X_base.size(), fd_model.X_base.size());
  for (unsigned i = 0; i < ad_model.X_base.size(); i++) {
    CHECK_ARRAY_EQUAL(ad_model.X_base[i].E.data(), fd_model.X_base[i].E.data(), 9);
    CHECK_ARRAY_EQUAL(ad_model.X_base[i].r.data(), fd_model.X_base[i].r.data(), 3);
  }

  CHECK_EQUAL(ad_model.X_J.size(), fd_model.X_J.size());
  for (unsigned i = 0; i < ad_model.X_J.size(); i++) {
    CHECK_ARRAY_EQUAL(ad_model.X_J[i].E.data(), fd_model.X_J[i].E.data(), 9);
    CHECK_ARRAY_EQUAL(ad_model.X_J[i].r.data(), fd_model.X_J[i].r.data(), 3);
  }

  for (unsigned idir = 0; idir < ndirs; idir++) {
    for (unsigned i = 0; i < ad_model.v_J.size(); i++) {
      CHECK_ARRAY_CLOSE(ad_d_model.v_J[i][idir].data(), fd_d_model.v_J[i][idir].data(), 6, 1e-7 );
    }
  }

  for (unsigned idir = 0; idir < ndirs; idir++) {
    for (unsigned i = 0; i < ad_model.c_J.size(); i++) {
      CHECK_ARRAY_CLOSE(ad_d_model.c_J[i][idir].data(), fd_d_model.c_J[i][idir].data(), 6, 1e-7 );
    }
  }

  for (unsigned idir = 0; idir < ndirs; idir++) {
    for (unsigned i = 0; i < ad_model.X_J.size(); i++) {
      CHECK_ARRAY_CLOSE(ad_d_model.X_J[i][idir].E.data(), fd_d_model.X_J[i][idir].E.data(), 9, 1e-7);
      CHECK_ARRAY_CLOSE(ad_d_model.X_J[i][idir].r.data(), fd_d_model.X_J[i][idir].r.data(), 3, 1e-7);
    }
  }

  for (unsigned idir = 0; idir < ndirs; idir++) {
    for (unsigned i = 0; i < ad_model.X_lambda.size(); i++) {
      CHECK_ARRAY_CLOSE(ad_d_model.X_lambda[i][idir].E.data(), fd_d_model.X_lambda[i][idir].E.data(), 9, 1e-7);
      CHECK_ARRAY_CLOSE(ad_d_model.X_lambda[i][idir].r.data(), fd_d_model.X_lambda[i][idir].r.data(), 3, 1e-7);
    }
  }

  for (unsigned idir = 0; idir < ndirs; idir++) {
    for (unsigned i = 0; i < ad_model.X_base.size(); i++) {
      CHECK_ARRAY_CLOSE(ad_d_model.X_base[i][idir].E.data(), fd_d_model.X_base[i][idir].E.data(), 9, 1e-7);
      CHECK_ARRAY_CLOSE(ad_d_model.X_base[i][idir].r.data(), fd_d_model.X_base[i][idir].r.data(), 3, 1e-7);
    }
  }

  for (unsigned idir = 0; idir < ndirs; idir++) {
    for (unsigned i = 0; i < ad_model.v.size(); i++) {
      CHECK_ARRAY_CLOSE(ad_d_model.v[i][idir].data(), fd_d_model.v[i][idir].data(), 6, 1e-7 );
    }
  }

  for (unsigned idir = 0; idir < ndirs; idir++) {
    for (unsigned i = 0; i < ad_model.c.size(); i++) {
      CHECK_ARRAY_CLOSE(ad_d_model.c[i][idir].data(), fd_d_model.c[i][idir].data(), 6, 1e-7 );
    }
  }

  for (unsigned idir = 0; idir < ndirs; idir++) {
    for (unsigned i = 0; i < ad_model.a.size(); i++) {
      CHECK_ARRAY_CLOSE(ad_d_model.a[i][idir].data(), fd_d_model.a[i][idir].data(), 6, 1e-7 );
    }
  }

  for (unsigned idir = 0; idir < ndirs; idir++) {
    for (unsigned i = 0; i < ad_model.U.size(); i++) {
      CHECK_ARRAY_CLOSE(ad_d_model.U[i][idir].data(), fd_d_model.U[i][idir].data(), 6, 1e-7 );
    }
  }
}

// -----------------------------------------------------------------------------
} //namespace RigidBodyDynamics
// -----------------------------------------------------------------------------

