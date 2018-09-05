#include <unittest++/UnitTest++.h>

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
    ADModel const & fd_d_model
) {

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
//    std::cout << "--- AD --- " << i << std::endl;
//    std::cout << ad_model.IA[i] << std::endl;
//    std::cout << "--- FD --- " << i << std::endl;
//    std::cout << fd_model.IA[i] << std::endl;
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
      CHECK_ARRAY_CLOSE(ad_d_model.v[i][idir].data(), fd_d_model.v[i][idir].data(), 6, 1e-6 );
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
      CHECK_ARRAY_CLOSE(ad_d_model.a[i][idir].data(), fd_d_model.a[i][idir].data(), 6, 5e-5);
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
      CHECK_ARRAY_CLOSE(ad_d_model.pA[i][idir].data(), fd_d_model.pA[i][idir].data(), 6, 5e-5);
    }
  }

  // nominal check
  CHECK_EQUAL(ad_model.u.rows(), fd_model.u.rows());
  CHECK_ARRAY_CLOSE(ad_model.u.data(), fd_model.u.data(), ad_model.u.rows(), 1e-7);
  // derivative check
  CHECK_ARRAY_CLOSE(ad_d_model.u.data(), fd_d_model.u.data(), ad_d_model.u.cols() * ad_model.u.rows(), 5e-5);

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


void checkModelsADvsED(
    unsigned ndirs,
    Model const & ad_model,
    ADModel const & ad_d_model,
    Model const & ed_model,
    EDModel const & ed_d_model
) {

  // nominal check
  CHECK_EQUAL(ad_model.v_J.size(), ed_model.v_J.size());
  for (unsigned i = 0; i < ad_model.v_J.size(); i++) {
    CHECK_ARRAY_EQUAL(ad_model.v_J[i].data(), ed_model.v_J[i].data(), 6);
  }
  // // derivative check
  // for (unsigned idir = 0; idir < ndirs; idir++) {
  //   for (unsigned i = 0; i < ad_model.v_J.size(); i++) {
  //     CHECK_ARRAY_CLOSE(ad_d_model.v_J[i][idir].data(), ed_d_model.v_J[i].col(idir).data(), 6, 1e-7 );
  //   }
  // }

  // nominal check
  CHECK_EQUAL(ad_model.c_J.size(), ed_model.c_J.size());
  for (unsigned i = 0; i < ad_model.c_J.size(); i++) {
    CHECK_ARRAY_EQUAL(ad_model.c_J[i].data(), ed_model.c_J[i].data(), 6);
  }
  // derivative check
  // for (unsigned idir = 0; idir < ndirs; idir++) {
  //   for (unsigned i = 0; i < ad_model.c_J.size(); i++) {
  //     CHECK_ARRAY_CLOSE(ad_d_model.c_J[i][idir].data(), ed_d_model.c_J[i].col(idir).data(), 6, 1e-7 );
  //   }
  // }

  // nominal check
  CHECK_EQUAL(ad_model.X_J.size(), ed_model.X_J.size());
  for (unsigned i = 0; i < ad_model.X_J.size(); i++) {
    CHECK_ARRAY_CLOSE(ad_model.X_J[i].E.data(), ed_model.X_J[i].E.data(), 9, 1e-7);
    CHECK_ARRAY_CLOSE(ad_model.X_J[i].r.data(), ed_model.X_J[i].r.data(), 3, 1e-7);
  }
  // derivative check
  // for (unsigned idir = 0; idir < ndirs; idir++) {
  //   for (unsigned i = 0; i < ad_model.X_J.size(); i++) {
  //     CHECK_ARRAY_CLOSE(ad_d_model.X_J[i][idir].E.data(), ed_d_model.X_J[i][idir].E.data(), 9, 1e-7);
  //     CHECK_ARRAY_CLOSE(ad_d_model.X_J[i][idir].r.data(), ed_d_model.X_J[i][idir].r.data(), 3, 1e-7);
  //   }
  // }


  // nominal check
  CHECK_EQUAL(ad_model.Ic.size(), ed_model.Ic.size());
  for (unsigned i = 0; i < ad_model.Ic.size(); i++) {
    CHECK_CLOSE(ad_model.Ic[i].m, ed_model.Ic[i].m, 1e-12);
    CHECK_ARRAY_CLOSE(ad_model.Ic[i].h.data(), ed_model.Ic[i].h.data(), 3, 1e-12);
    CHECK_CLOSE(ad_model.Ic[i].Ixx, ed_model.Ic[i].Ixx, 1e-12);
    CHECK_CLOSE(ad_model.Ic[i].Iyx, ed_model.Ic[i].Iyx, 1e-12);
    CHECK_CLOSE(ad_model.Ic[i].Iyy, ed_model.Ic[i].Iyy, 1e-12);
    CHECK_CLOSE(ad_model.Ic[i].Izx, ed_model.Ic[i].Izx, 1e-12);
    CHECK_CLOSE(ad_model.Ic[i].Izy, ed_model.Ic[i].Izy, 1e-12);
    CHECK_CLOSE(ad_model.Ic[i].Izz, ed_model.Ic[i].Izz, 1e-12);
  }
  // derivative check
  for (unsigned i = 0; i < ad_model.Ic.size(); i++) {
    for (unsigned idir = 0; idir < ndirs; idir++) {
      // DEBUG OUTPUT
      RigidBodyDynamics::Math::SpatialMatrix error
        = (ad_d_model.Ic[i][idir].toMatrix() - ed_d_model.Ic[i][idir]).cwiseAbs();
      double max = error.maxCoeff();
      if (max > 1e-12) {
        std::cout << "error [" << i << "][" << idir << "] = \n" << error << std::endl;
        std::cout << "max   [" << i << "][" << idir << "] = " << max << std::endl;
        std::cout << "ad_d_model.Ic[" << i << "][" << idir << "] = \n"
          << ad_d_model.Ic[i][idir].toMatrix() << std::endl;
        std::cout << "ed_d_model.Ic[" << i << "][" << idir << "] = \n"
          << ed_d_model.Ic[i][idir] << std::endl;
        std::cout << endl;
      }
      CHECK_ARRAY_CLOSE(
        ad_d_model.Ic[i][idir].toMatrix().data(), ed_d_model.Ic[i][idir].data(),
        6*6, 1e-12
      );
    }
  }

  // nominal check
  CHECK_EQUAL(ad_model.X_lambda.size(), ed_model.X_lambda.size());
  for (unsigned i = 0; i < ad_model.X_lambda.size(); i++) {
    CHECK_ARRAY_CLOSE(ad_model.X_lambda[i].E.data(), ed_model.X_lambda[i].E.data(), 9, 1e-7);
    CHECK_ARRAY_CLOSE(ad_model.X_lambda[i].r.data(), ed_model.X_lambda[i].r.data(), 3, 1e-7);
  }
  // derivative check
  // for (unsigned idir = 0; idir < ndirs; idir++) {
  //   for (unsigned i = 0; i < ad_model.X_lambda.size(); i++) {
  //     CHECK_ARRAY_CLOSE(ad_d_model.X_lambda[i][idir].E.data(), ed_d_model.X_lambda[i][idir].E.data(), 9, 1e-7);
  //     CHECK_ARRAY_CLOSE(ad_d_model.X_lambda[i][idir].r.data(), ed_d_model.X_lambda[i][idir].r.data(), 3, 1e-7);
  //   }
  // }

  /*
  // nominal check
  CHECK_EQUAL(ad_model.U.size(), ed_model.U.size());
  for (unsigned i = 0; i < ad_model.U.size(); i++) {
    CHECK_ARRAY_CLOSE(ad_model.U[i].data(), ed_model.U[i].data(), 6, 1e-7);
  }
  // derivative check
  for (unsigned idir = 0; idir < ndirs; idir++) {
    for (unsigned i = 0; i < ad_model.U.size(); i++) {
      CHECK_ARRAY_CLOSE(ad_d_model.U[i][idir].data(), ed_d_model.U[i].col(idir).data(), 6, 1e-5);
    }
  }

  // nominal check
  CHECK_EQUAL(ad_model.d.rows(), ed_model.d.rows());
  CHECK_ARRAY_CLOSE(ad_model.d.data(), ed_model.d.data(), ad_model.d.rows(), 1e-7);
  // derivative check
  CHECK_ARRAY_CLOSE(ad_d_model.d.data(), ed_d_model.d.data(), ad_d_model.u.cols() * ad_model.d.rows(), 1e-5);

  // nominal check
  CHECK_EQUAL(ad_model.IA.size(), ed_model.IA.size());
  for (unsigned i = 0; i < ad_model.IA.size(); i++) {
//    std::cout << "--- AD --- " << i << std::endl;
//    std::cout << ad_model.IA[i] << std::endl;
//    std::cout << "--- FD --- " << i << std::endl;
//    std::cout << ed_model.IA[i] << std::endl;
    CHECK_ARRAY_CLOSE(ad_model.IA[i].data(), ed_model.IA[i].data(), 36, 1e-7);
  }
  // derivative check
  for (unsigned idir = 0; idir < ndirs; idir++) {
    for (unsigned i = 0; i < ad_model.IA.size(); i++) {
      CHECK_ARRAY_CLOSE(ad_d_model.IA[i][idir].data(), ed_d_model.IA[i][idir].data(), 36, 1e-5);
    }
  }

  // nominal check
  // std::cout << "Check X_base disabled!" << std::endl;
  // CHECK_EQUAL(ad_model.X_base.size(), ed_model.X_base.size());
  // for (unsigned i = 0; i < ad_model.X_base.size(); i++) {
  //   CHECK_ARRAY_EQUAL(ad_model.X_base[i].E.data(), ed_model.X_base[i].E.data(), 9);
  //   CHECK_ARRAY_EQUAL(ad_model.X_base[i].r.data(), ed_model.X_base[i].r.data(), 3);
  // }
  // // derivative check
  // for (unsigned idir = 0; idir < ndirs; idir++) {
  //   for (unsigned i = 0; i < ad_model.X_base.size(); i++) {
  //     CHECK_ARRAY_CLOSE(ad_d_model.X_base[i][idir].E.data(), ed_d_model.X_base[i][idir].E.data(), 9, 1e-7);
  //     CHECK_ARRAY_CLOSE(ad_d_model.X_base[i][idir].r.data(), ed_d_model.X_base[i][idir].r.data(), 3, 1e-7);
  //   }
  // }

  // nominal check
  CHECK_EQUAL(ad_model.v.size(), ed_model.v.size());
  for (unsigned i = 0; i < ad_model.v.size(); i++) {
    CHECK_ARRAY_EQUAL(ad_model.v[i].data(), ed_model.v[i].data(), 6);
  }
  // derivative check
  for (unsigned idir = 0; idir < ndirs; idir++) {
    for (unsigned i = 0; i < ad_model.v.size(); i++) {
      CHECK_ARRAY_CLOSE(ad_d_model.v[i][idir].data(), ed_d_model.v[i].col(idir).data(), 6, 1e-6 );
    }
  }

  // nominal check
  CHECK_EQUAL(ad_model.c.size(), ed_model.c.size());
  for (unsigned i = 0; i < ad_model.c.size(); i++) {
    CHECK_ARRAY_EQUAL(ad_model.c[i].data(), ed_model.c[i].data(), 6);
  }
  // derivative check
  for (unsigned idir = 0; idir < ndirs; idir++) {
    for (unsigned i = 0; i < ad_model.c.size(); i++) {
      CHECK_ARRAY_CLOSE(ad_d_model.c[i][idir].data(), ed_d_model.c[i].col(idir).data(), 6, 1e-7 );
    }
  }

  // nominal check
  CHECK_EQUAL(ad_model.a.size(), ed_model.a.size());
  for (unsigned i = 0; i < ad_model.a.size(); i++) {
    CHECK_ARRAY_CLOSE(ad_model.a[i].data(), ed_model.a[i].data(), 6, 1e-7);
  }
  // derivative check
  for (unsigned idir = 0; idir < ndirs; idir++) {
    for (unsigned i = 0; i < ad_model.a.size(); i++) {
      CHECK_ARRAY_CLOSE(ad_d_model.a[i][idir].data(), ed_d_model.a[i].col(idir).data(), 6, 5e-5);
    }
  }

  // nominal check
  CHECK_EQUAL(ad_model.pA.size(), ed_model.pA.size());
  for (unsigned i = 0; i < ad_model.pA.size(); i++) {
    CHECK_ARRAY_CLOSE(ad_model.pA[i].data(), ed_model.pA[i].data(), 6, 1e-7 );
  }
  // derivative check
  for (unsigned idir = 0; idir < ndirs; idir++) {
    for (unsigned i = 0; i < ad_model.pA.size(); i++) {
      CHECK_ARRAY_CLOSE(ad_d_model.pA[i][idir].data(), ed_d_model.pA[i].col(idir).data(), 6, 5e-5);
    }
  }

  // nominal check
  CHECK_EQUAL(ad_model.u.rows(), ed_model.u.rows());
  CHECK_ARRAY_CLOSE(ad_model.u.data(), ed_model.u.data(), ad_model.u.rows(), 1e-7);
  // derivative check
  CHECK_ARRAY_CLOSE(ad_d_model.u.data(), ed_d_model.u.data(), ad_d_model.u.cols() * ad_model.u.rows(), 5e-5);

  // nominal check
  CHECK_EQUAL(ad_model.f.size(), ed_model.f.size());
  for (unsigned i = 0; i < ad_model.f.size(); i++) {
    CHECK_ARRAY_CLOSE(ad_model.f[i].data(), ed_model.f[i].data(), 6, 1e-7);
  }
  // derivative check
  for (unsigned idir = 0; idir < ndirs; idir++) {
    for (unsigned i = 0; i < ad_model.f.size(); i++) {
      CHECK_ARRAY_CLOSE(ad_d_model.f[i][idir].data(), ed_d_model.f[i].col(idir).data(), 6, 5e-5);
    }
  }
*/
  return;
}


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

  // nominal check
  CHECK_EQUAL(ad_cs.GSpi.rows(), fd_cs.GSpi.rows());
  CHECK_EQUAL(ad_cs.GSpi.cols(), fd_cs.GSpi.cols());
  CHECK_ARRAY_CLOSE(ad_cs.GSpi.data(), fd_cs.GSpi.data(),
                    fd_cs.GSpi.rows() * fd_cs.GSpi.cols(),
                    1e-6);
  // derivative check
  for (unsigned idir = 0; idir < ndirs; idir++) {
    CHECK_ARRAY_CLOSE(ad_d_cs.GSpi[idir].data(), fd_d_cs.GSpi[idir].data(),
                      fd_cs.GSpi.rows() * fd_cs.GSpi.cols(),
                      1e-6);
  }

  // nominal check
  CHECK_EQUAL(ad_cs.GSsi.rows(), fd_cs.GSsi.rows());
  CHECK_EQUAL(ad_cs.GSsi.cols(), fd_cs.GSsi.cols());
  CHECK_ARRAY_CLOSE(ad_cs.GSsi.data(), fd_cs.GSsi.data(),
                    fd_cs.GSsi.rows() * fd_cs.GSsi.cols(),
                    1e-6);
  // derivative check
  for (unsigned idir = 0; idir < ndirs; idir++) {
    CHECK_ARRAY_CLOSE(ad_d_cs.GSsi[idir].data(), fd_d_cs.GSsi[idir].data(),
                      fd_cs.GSsi.rows() * fd_cs.GSsi.cols(),
                      1e-6);
  }

  // nominal check
  CHECK_EQUAL(ad_cs.C.rows(), fd_cs.C.rows());
  CHECK_ARRAY_CLOSE(ad_cs.C.data(), fd_cs.C.data(), ad_cs.C.rows(), 1e-6);
  // derivative check
  CHECK_EQUAL(ad_d_cs.C.rows(), fd_d_cs.C.rows());
  CHECK_EQUAL(ad_d_cs.C.cols(), fd_d_cs.C.cols());
  CHECK_ARRAY_CLOSE(ad_d_cs.C.data(), fd_d_cs.C.data(),
                    ad_d_cs.C.rows() * ad_d_cs.C.cols(), 1e-5);

  // nominal check
  CHECK_EQUAL(ad_cs.err.rows(), fd_cs.err.rows());
  CHECK_ARRAY_CLOSE(ad_cs.err.data(), fd_cs.err.data(), ad_cs.err.rows(), 1e-6);
  // derivative check
  CHECK_EQUAL(ad_d_cs.err.rows(), fd_d_cs.err.rows());
  CHECK_EQUAL(ad_d_cs.err.cols(), fd_d_cs.err.cols());
  CHECK_ARRAY_CLOSE(ad_d_cs.err.data(), fd_d_cs.err.data(),
                    ad_d_cs.err.rows() * ad_d_cs.err.cols(), 1e-6);

  // nominal check
  CHECK_EQUAL(ad_cs.errd.rows(), fd_cs.errd.rows());
  CHECK_ARRAY_CLOSE(ad_cs.errd.data(), fd_cs.errd.data(),
                    ad_cs.errd.rows(), 1e-6);
  // derivative check
  CHECK_EQUAL(ad_d_cs.errd.rows(), fd_d_cs.errd.rows());
  CHECK_EQUAL(ad_d_cs.errd.cols(), fd_d_cs.errd.cols());
  CHECK_ARRAY_CLOSE(ad_d_cs.errd.data(), fd_d_cs.errd.data(),
                    ad_d_cs.errd.rows() * ad_d_cs.errd.cols(), 1e-6);

  // nominal check
  CHECK_EQUAL(ad_cs.QDDot_0.rows(), fd_cs.QDDot_0.rows());
  CHECK_ARRAY_CLOSE(ad_cs.QDDot_0.data(), fd_cs.QDDot_0.data(),
                    ad_cs.QDDot_0.rows(), 1e-6);
  // derivative check
  CHECK_EQUAL(ad_d_cs.QDDot_0.rows(), fd_d_cs.QDDot_0.rows());
  CHECK_EQUAL(ad_d_cs.QDDot_0.cols(), fd_d_cs.QDDot_0.cols());
  CHECK_ARRAY_CLOSE(ad_d_cs.QDDot_0.data(), fd_d_cs.QDDot_0.data(),
                    ad_d_cs.QDDot_0.rows() * ad_d_cs.QDDot_0.cols(), 1e-6);

  // nominal check
  CHECK_EQUAL(ad_cs.gamma.rows(), fd_cs.gamma.rows());
  CHECK_ARRAY_CLOSE(ad_cs.gamma.data(), fd_cs.gamma.data(),
                    ad_cs.gamma.rows(), 1e-6);
  // derivative check
  CHECK_EQUAL(ad_d_cs.gamma.rows(), fd_d_cs.gamma.rows());
  CHECK_EQUAL(ad_d_cs.gamma.cols(), fd_d_cs.gamma.cols());
  CHECK_ARRAY_CLOSE(ad_d_cs.gamma.data(), fd_d_cs.gamma.data(),
                    ad_d_cs.gamma.rows() * ad_d_cs.gamma.cols(), 1e-6);

  // nominal check
  CHECK_EQUAL(ad_cs.d_u.rows(), fd_cs.d_u.rows());
  CHECK_ARRAY_CLOSE(ad_cs.d_u.data(), fd_cs.d_u.data(),
                    ad_cs.d_u.rows(), 1e-6);
  // derivative check
  CHECK_EQUAL(ad_d_cs.d_u.rows(), fd_d_cs.d_u.rows());
  CHECK_EQUAL(ad_d_cs.d_u.cols(), fd_d_cs.d_u.cols());
  CHECK_ARRAY_CLOSE(ad_d_cs.d_u.data(), fd_d_cs.d_u.data(),
                    ad_d_cs.d_u.rows() * ad_d_cs.d_u.cols(), 1e-6);

  // nominal check
  CHECK_EQUAL(ad_cs.d_a.size(), fd_cs.d_a.size());
  for (unsigned i = 0; i < ad_cs.d_a.size(); i++) {
    CHECK_ARRAY_CLOSE(ad_cs.d_a[i].data(), fd_cs.d_a[i].data(),
                      ad_cs.d_a[i].rows(), 1e-6);
    // derivative check
    CHECK_EQUAL(ad_d_cs.d_a[i].rows(), fd_d_cs.d_a[i].rows());
    CHECK_EQUAL(ad_d_cs.d_a[i].cols(), fd_d_cs.d_a[i].cols());
    CHECK_ARRAY_CLOSE(ad_d_cs.d_a[i].data(), fd_d_cs.d_a[i].data(),
                      ad_d_cs.d_a[i].rows() * ad_d_cs.d_a[i].cols(), 1e-6);
  }


}

// -----------------------------------------------------------------------------
} //namespace RigidBodyDynamics
// -----------------------------------------------------------------------------


