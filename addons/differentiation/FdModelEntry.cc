#include "FdModelEntry.h"

namespace RigidBodyDynamics {

void computeFDEntry(
    Model const & model,
    Model const & modelh,
    double h,
    int idir,
    ADModel & fd_model) {

  for (unsigned i = 0; i < model.X_lambda.size(); i++) {
    fd_model.X_lambda[i][idir].E = (modelh.X_lambda[i].E - model.X_lambda[i].E) / h;
    fd_model.X_lambda[i][idir].r = (modelh.X_lambda[i].r - model.X_lambda[i].r) / h;
  }
  for (unsigned i = 0; i < model.X_base.size(); i++) {
    fd_model.X_base[i][idir].E = (modelh.X_base[i].E - model.X_base[i].E) / h;
    fd_model.X_base[i][idir].r = (modelh.X_base[i].r - model.X_base[i].r) / h;
  }
  for (unsigned i = 0; i < model.X_J.size(); i++) {
    fd_model.X_J[i][idir].E = (modelh.X_J[i].E - model.X_J[i].E) / h;
    fd_model.X_J[i][idir].r = (modelh.X_J[i].r - model.X_J[i].r) / h;
  }
  for (unsigned i = 0; i < model.c.size(); i++) {
    fd_model.c[i][idir] = (modelh.c[i] - model.c[i]) / h;
  }
  for (unsigned i = 0; i < model.c_J.size(); i++) {
    fd_model.c_J[i][idir] = (modelh.c_J[i] - model.c_J[i]) / h;
  }
  for (unsigned i = 0; i < model.v.size(); i++) {
    fd_model.v[i][idir] = (modelh.v[i] - model.v[i]) / h;
  }
  for (unsigned i = 0; i < model.v_J.size(); i++) {
    fd_model.v_J[i][idir] = (modelh.v_J[i] - model.v_J[i]) / h;
  }
  for (unsigned i = 0; i < model.a.size(); i++) {
    fd_model.a[i][idir] = (modelh.a[i] - model.a[i]) / h;
  }
  for (unsigned i = 0; i < model.U.size(); i++) {
    fd_model.U[i][idir] = (modelh.U[i] - model.U[i]) / h;
  }
}

}
