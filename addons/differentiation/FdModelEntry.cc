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
  for (unsigned i = 0; i < model.Ic.size(); i++) {
    fd_model.Ic[i][idir].m = (modelh.Ic[i].m - model.Ic[i].m) / h;
    fd_model.Ic[i][idir].h = (modelh.Ic[i].h - model.Ic[i].h) / h;
    fd_model.Ic[i][idir].Ixx = (modelh.Ic[i].Ixx - model.Ic[i].Ixx) / h;
    fd_model.Ic[i][idir].Iyx = (modelh.Ic[i].Iyx - model.Ic[i].Iyx) / h;
    fd_model.Ic[i][idir].Iyy = (modelh.Ic[i].Iyy - model.Ic[i].Iyy) / h;
    fd_model.Ic[i][idir].Izx = (modelh.Ic[i].Izx - model.Ic[i].Izx) / h;
    fd_model.Ic[i][idir].Izy = (modelh.Ic[i].Izy - model.Ic[i].Izy) / h;
    fd_model.Ic[i][idir].Izz = (modelh.Ic[i].Izz - model.Ic[i].Izz) / h;
  }
  for (unsigned i = 0; i < model.IA.size(); i++) {
    fd_model.IA[i][idir] = (modelh.IA[i] - model.IA[i]) / h;
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
  for (unsigned i = 0; i < model.pA.size(); i++) {
    fd_model.pA[i][idir] = (modelh.pA[i] - model.pA[i]) / h;
  }
  for (unsigned i = 0; i < model.u.rows(); i++) {
    fd_model.u(i, idir) = (modelh.u(i) - model.u(i)) / h;
  }
  for (unsigned i = 0; i < model.d.rows(); i++) {
    fd_model.d(i, idir) = (modelh.d(i) - model.d(i)) / h;
  }
  for (unsigned i = 0; i < model.f.size(); i++) {
    fd_model.f[i][idir] = (modelh.f[i] - model.f[i]) / h;
  }

}

}
