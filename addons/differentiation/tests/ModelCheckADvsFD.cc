#include <UnitTest++.h>

#include "ModelCheckADvsFD.h"

using namespace std;

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------

void checkModelsADvsFD(
    unsigned ndirs,
    Model const & ad_model,
    ADModel const & ad_d_model,
    Model const & fd_model,
    ADModel const & fd_d_model) {
  // nominal check
  CHECK_EQUAL(ad_model.v_J.size(), fd_model.v_J.size());
  for (unsigned i = 0; i < ad_model.v_J.size(); i++) {
    CHECK_ARRAY_EQUAL(ad_model.v_J[i].data(), fd_model.v_J[i].data(), 6);
  }
  // derivative check
  for (unsigned idir = 0; idir < ndirs; idir++) {
    for (unsigned i = 0; i < ad_model.v_J.size(); i++) {
      CHECK_ARRAY_CLOSE(ad_d_model.v_J[i][idir].data(), fd_d_model.v_J[i][idir].data(), 6, 1e-7 );
    }
  }

  // nominal check
  CHECK_EQUAL(ad_model.c_J.size(), fd_model.c_J.size());
  for (unsigned i = 0; i < ad_model.c_J.size(); i++) {
    CHECK_ARRAY_EQUAL(ad_model.c_J[i].data(), fd_model.c_J[i].data(), 6);
  }
  // derivative check
  for (unsigned idir = 0; idir < ndirs; idir++) {
    for (unsigned i = 0; i < ad_model.c_J.size(); i++) {
      CHECK_ARRAY_CLOSE(ad_d_model.c_J[i][idir].data(), fd_d_model.c_J[i][idir].data(), 6, 1e-7 );
    }
  }

  // nominal check
  CHECK_EQUAL(ad_model.X_J.size(), fd_model.X_J.size());
  for (unsigned i = 0; i < ad_model.X_J.size(); i++) {
    CHECK_ARRAY_CLOSE(ad_model.X_J[i].E.data(), fd_model.X_J[i].E.data(), 9, 1e-7);
    CHECK_ARRAY_CLOSE(ad_model.X_J[i].r.data(), fd_model.X_J[i].r.data(), 3, 1e-7);
  }
  // derivative check
  for (unsigned idir = 0; idir < ndirs; idir++) {
    for (unsigned i = 0; i < ad_model.X_J.size(); i++) {
      CHECK_ARRAY_CLOSE(ad_d_model.X_J[i][idir].E.data(), fd_d_model.X_J[i][idir].E.data(), 9, 1e-7);
      CHECK_ARRAY_CLOSE(ad_d_model.X_J[i][idir].r.data(), fd_d_model.X_J[i][idir].r.data(), 3, 1e-7);
    }
  }

  // nominal check
  CHECK_EQUAL(ad_model.Ic.size(), fd_model.Ic.size());
  for (unsigned i = 0; i < ad_model.Ic.size(); i++) {
    CHECK_CLOSE(ad_model.Ic[i].m, fd_model.Ic[i].m, 1e-7);
    CHECK_ARRAY_CLOSE(ad_model.Ic[i].h.data(), fd_model.Ic[i].h.data(), 3, 1e-7);
    CHECK_CLOSE(ad_model.Ic[i].Ixx, fd_model.Ic[i].Ixx, 1e-7);
    CHECK_CLOSE(ad_model.Ic[i].Iyx, fd_model.Ic[i].Iyx, 1e-7);
    CHECK_CLOSE(ad_model.Ic[i].Iyy, fd_model.Ic[i].Iyy, 1e-7);
    CHECK_CLOSE(ad_model.Ic[i].Izx, fd_model.Ic[i].Izx, 1e-7);
    CHECK_CLOSE(ad_model.Ic[i].Izy, fd_model.Ic[i].Izy, 1e-7);
    CHECK_CLOSE(ad_model.Ic[i].Izz, fd_model.Ic[i].Izz, 1e-7);
  }
  // derivative check
  for (unsigned idir = 0; idir < ndirs; idir++) {
    for (unsigned i = 0; i < ad_model.Ic.size(); i++) {
      CHECK_CLOSE(ad_d_model.Ic[i][idir].m, fd_d_model.Ic[i][idir].m, 1e-6);
      CHECK_ARRAY_CLOSE(ad_d_model.Ic[i][idir].h.data(), fd_d_model.Ic[i][idir].h.data(), 3, 1e-6);
      CHECK_CLOSE(ad_d_model.Ic[i][idir].Ixx, fd_d_model.Ic[i][idir].Ixx, 1e-6);
      CHECK_CLOSE(ad_d_model.Ic[i][idir].Iyx, fd_d_model.Ic[i][idir].Iyx, 1e-6);
      CHECK_CLOSE(ad_d_model.Ic[i][idir].Iyy, fd_d_model.Ic[i][idir].Iyy, 1e-6);
      CHECK_CLOSE(ad_d_model.Ic[i][idir].Izx, fd_d_model.Ic[i][idir].Izx, 1e-6);
      CHECK_CLOSE(ad_d_model.Ic[i][idir].Izy, fd_d_model.Ic[i][idir].Izy, 1e-6);
      CHECK_CLOSE(ad_d_model.Ic[i][idir].Izz, fd_d_model.Ic[i][idir].Izz, 1e-6);
    }
  }

  // nominal check
  CHECK_EQUAL(ad_model.X_lambda.size(), fd_model.X_lambda.size());
  for (unsigned i = 0; i < ad_model.X_lambda.size(); i++) {
    CHECK_ARRAY_CLOSE(ad_model.X_lambda[i].E.data(), fd_model.X_lambda[i].E.data(), 9, 1e-7);
    CHECK_ARRAY_CLOSE(ad_model.X_lambda[i].r.data(), fd_model.X_lambda[i].r.data(), 3, 1e-7);
  }
  // derivative check
  for (unsigned idir = 0; idir < ndirs; idir++) {
    for (unsigned i = 0; i < ad_model.X_lambda.size(); i++) {
      CHECK_ARRAY_CLOSE(ad_d_model.X_lambda[i][idir].E.data(), fd_d_model.X_lambda[i][idir].E.data(), 9, 1e-7);
      CHECK_ARRAY_CLOSE(ad_d_model.X_lambda[i][idir].r.data(), fd_d_model.X_lambda[i][idir].r.data(), 3, 1e-7);
    }
  }

  // nominal check
  CHECK_EQUAL(ad_model.U.size(), fd_model.U.size());
  for (unsigned i = 0; i < ad_model.U.size(); i++) {
    CHECK_ARRAY_CLOSE(ad_model.U[i].data(), fd_model.U[i].data(), 6, 1e-7);
  }
  // derivative check
  for (unsigned idir = 0; idir < ndirs; idir++) {
    for (unsigned i = 0; i < ad_model.U.size(); i++) {
      CHECK_ARRAY_CLOSE(ad_d_model.U[i][idir].data(), fd_d_model.U[i][idir].data(), 6, 1e-5);
    }
  }

  // nominal check
  CHECK_EQUAL(ad_model.d.rows(), fd_model.d.rows());
  CHECK_ARRAY_CLOSE(ad_model.d.data(), fd_model.d.data(), ad_model.d.rows(), 1e-7);
  // derivative check
  CHECK_ARRAY_CLOSE(ad_d_model.d.data(), fd_d_model.d.data(), ad_d_model.u.cols() * ad_model.d.rows(), 1e-5);

  // nominal check
  CHECK_EQUAL(ad_model.IA.size(), fd_model.IA.size());
  for (unsigned i = 0; i < ad_model.IA.size(); i++) {
    CHECK_ARRAY_CLOSE(ad_model.IA[i].data(), fd_model.IA[i].data(), 36, 1e-7);
  }
  // derivative check
  for (unsigned idir = 0; idir < ndirs; idir++) {
    for (unsigned i = 0; i < ad_model.IA.size(); i++) {
      CHECK_ARRAY_CLOSE(ad_d_model.IA[i][idir].data(), fd_d_model.IA[i][idir].data(), 36, 1e-5);
    }
  }

  // nominal check
  CHECK_EQUAL(ad_model.X_base.size(), fd_model.X_base.size());
  for (unsigned i = 0; i < ad_model.X_base.size(); i++) {
    CHECK_ARRAY_EQUAL(ad_model.X_base[i].E.data(), fd_model.X_base[i].E.data(), 9);
    CHECK_ARRAY_EQUAL(ad_model.X_base[i].r.data(), fd_model.X_base[i].r.data(), 3);
  }
  // derivative check
  for (unsigned idir = 0; idir < ndirs; idir++) {
    for (unsigned i = 0; i < ad_model.X_base.size(); i++) {
      CHECK_ARRAY_CLOSE(ad_d_model.X_base[i][idir].E.data(), fd_d_model.X_base[i][idir].E.data(), 9, 1e-7);
      CHECK_ARRAY_CLOSE(ad_d_model.X_base[i][idir].r.data(), fd_d_model.X_base[i][idir].r.data(), 3, 1e-7);
    }
  }

  // nominal check
  CHECK_EQUAL(ad_model.v.size(), fd_model.v.size());
  for (unsigned i = 0; i < ad_model.v.size(); i++) {
    CHECK_ARRAY_EQUAL(ad_model.v[i].data(), fd_model.v[i].data(), 6);
  }
  // derivative check
  for (unsigned idir = 0; idir < ndirs; idir++) {
    for (unsigned i = 0; i < ad_model.v.size(); i++) {
      CHECK_ARRAY_CLOSE(ad_d_model.v[i][idir].data(), fd_d_model.v[i][idir].data(), 6, 1e-7 );
    }
  }

  // nominal check
  CHECK_EQUAL(ad_model.c.size(), fd_model.c.size());
  for (unsigned i = 0; i < ad_model.c.size(); i++) {
    CHECK_ARRAY_EQUAL(ad_model.c[i].data(), fd_model.c[i].data(), 6);
  }
  // derivative check
  for (unsigned idir = 0; idir < ndirs; idir++) {
    for (unsigned i = 0; i < ad_model.c.size(); i++) {
      CHECK_ARRAY_CLOSE(ad_d_model.c[i][idir].data(), fd_d_model.c[i][idir].data(), 6, 1e-7 );
    }
  }

  // nominal check
  CHECK_EQUAL(ad_model.a.size(), fd_model.a.size());
  for (unsigned i = 0; i < ad_model.a.size(); i++) {
    CHECK_ARRAY_CLOSE(ad_model.a[i].data(), fd_model.a[i].data(), 6, 1e-7);
  }
  // derivative check
  for (unsigned idir = 0; idir < ndirs; idir++) {
    for (unsigned i = 0; i < ad_model.a.size(); i++) {
      CHECK_ARRAY_CLOSE(ad_d_model.a[i][idir].data(), fd_d_model.a[i][idir].data(), 6, 1e-5);
    }
  }

  // nominal check
  CHECK_EQUAL(ad_model.pA.size(), fd_model.pA.size());
  for (unsigned i = 0; i < ad_model.pA.size(); i++) {
    CHECK_ARRAY_CLOSE(ad_model.pA[i].data(), fd_model.pA[i].data(), 6, 1e-7 );
  }
  // derivative check
  for (unsigned idir = 0; idir < ndirs; idir++) {
    for (unsigned i = 0; i < ad_model.pA.size(); i++) {
      CHECK_ARRAY_CLOSE(ad_d_model.pA[i][idir].data(), fd_d_model.pA[i][idir].data(), 6, 1e-6);
    }
  }

  // nominal check
  CHECK_EQUAL(ad_model.u.rows(), fd_model.u.rows());
  CHECK_ARRAY_CLOSE(ad_model.u.data(), fd_model.u.data(), ad_model.u.rows(), 1e-7);
  // derivative check
  CHECK_ARRAY_CLOSE(ad_d_model.u.data(), fd_d_model.u.data(), ad_d_model.u.cols() * ad_model.u.rows(), 1e-6);

  // nominal check
  CHECK_EQUAL(ad_model.f.size(), fd_model.f.size());
  for (unsigned i = 0; i < ad_model.f.size(); i++) {
    CHECK_ARRAY_CLOSE(ad_model.f[i].data(), fd_model.f[i].data(), 6, 1e-7);
  }
  // derivative check
  for (unsigned idir = 0; idir < ndirs; idir++) {
    for (unsigned i = 0; i < ad_model.f.size(); i++) {
      CHECK_ARRAY_CLOSE(ad_d_model.f[i][idir].data(), fd_d_model.f[i][idir].data(), 6, 5e-5);
    }
  }
}

//cout << endl << i << "," << idir << endl;
//cout << "AD.f = " << ad_d_model.f[i][idir].transpose() << endl;
//cout << "FD.f = " << fd_d_model.f[i][idir].transpose() << endl;
//cout << "norm = " << (fd_d_model.f[i][idir] - ad_d_model.f[i][idir]).norm() << endl;

void checkConstraintSetsADvsFD (
    unsigned ndirs,
    ConstraintSet const & ad_cs,
    ADConstraintSet const & ad_d_cs,
    ConstraintSet const & fd_cs,
    ADConstraintSet const & fd_d_cs) {

  // nominal check
  CHECK_EQUAL(ad_cs.H.rows(), fd_cs.H.rows());
  CHECK_EQUAL(ad_cs.H.cols(), fd_cs.H.cols());
  CHECK_ARRAY_CLOSE(ad_cs.H.data(), fd_cs.H.data(),
                    fd_cs.H.rows() * fd_cs.H.cols(),
                    1e-6);
  // derivative check
  for (unsigned idir = 0; idir < ndirs; idir++) {
    CHECK_ARRAY_CLOSE(ad_d_cs.H[idir].data(), fd_d_cs.H[idir].data(),
                      fd_cs.H.rows() * fd_cs.H.cols(),
                      1e-6);
  }

  // nominal check
  CHECK_EQUAL(ad_cs.G.rows(), fd_cs.G.rows());
  CHECK_EQUAL(ad_cs.G.cols(), fd_cs.G.cols());
  CHECK_ARRAY_CLOSE(ad_cs.G.data(), fd_cs.G.data(),
                    fd_cs.G.rows() * fd_cs.G.cols(),
                    1e-6);
  // derivative check
  for (unsigned idir = 0; idir < ndirs; idir++) {
    CHECK_ARRAY_CLOSE(ad_d_cs.G[idir].data(), fd_d_cs.G[idir].data(),
                      fd_cs.G.rows() * fd_cs.G.cols(),
                      1e-6);
  }

  // nominal check
  CHECK_EQUAL(ad_cs.A.rows(), fd_cs.A.rows());
  CHECK_EQUAL(ad_cs.A.cols(), fd_cs.A.cols());
  CHECK_ARRAY_CLOSE(ad_cs.A.data(), fd_cs.A.data(),
                    fd_cs.A.rows() * fd_cs.A.cols(),
                    1e-6);
  // derivative check
  for (unsigned idir = 0; idir < ndirs; idir++) {
    CHECK_ARRAY_CLOSE(ad_d_cs.A[idir].data(), fd_d_cs.A[idir].data(),
                      fd_cs.A.rows() * fd_cs.A.cols(),
                      1e-6);
  }

  // nominal check
  CHECK_EQUAL(ad_cs.Gi.rows(), fd_cs.Gi.rows());
  CHECK_EQUAL(ad_cs.Gi.cols(), fd_cs.Gi.cols());
  CHECK_ARRAY_CLOSE(ad_cs.Gi.data(), fd_cs.Gi.data(),
                    fd_cs.Gi.rows() * fd_cs.Gi.cols(),
                    1e-6);
  // derivative check
  for (unsigned idir = 0; idir < ndirs; idir++) {
    CHECK_ARRAY_CLOSE(ad_d_cs.Gi[idir].data(), fd_d_cs.Gi[idir].data(),
                      fd_cs.Gi.rows() * fd_cs.Gi.cols(),
                      1e-6);
  }
}



// -----------------------------------------------------------------------------
} //namespace RigidBodyDynamics
// -----------------------------------------------------------------------------

