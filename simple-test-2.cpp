/*  MUSCOD_SUITE/Apps/TEST/Src/SRC/simple-test-2.cpp
 *  (c) Tran, B Duy, 2022
 *  Simple test (using the convexification exactly 2 times)
 */

#include <cmath>
#include "def_usrmod.hpp"

#define  NMOS   1
#define  NP     0
#define  NRC    0
#define  NRCE   0

#define  NXD    2
#define  NXA    0
#define  NU     4
#define  NPR    0

#define  NRD_S  2
#define  NRDE_S 2

#define  NRD_I  10
#define  NRDE_I 0

static void mfcn(double *ts, double *sd, double *sa,
  double *p, double *pr, double *mval, long *dpnd, InfoPtr *info)
{
  if (*dpnd) {
    *dpnd = MFCN_DPND(0, *sd, 0, 0, 0);
    return;
  }

  *mval = sd[0]*sd[0] + sd[1]*sd[1];
}


static void ffcn(double *t, double *xd, double *xa, double *u,
  double *p, double *rhs, double *rwh, long *iwh, InfoPtr *info)
{
  rhs[0] = -u[2]*u[0];
  rhs[1] = u[3]*u[0];
}

static void rdfcn_s(double *ts, double *sd, double *sa, double *u,
  double *p, double *pr, double *res, long *dpnd, InfoPtr *info)
{
  if (*dpnd) {
    *dpnd = RFCN_DPND(0, *sd, 0, 0, 0, 0);
    return;
  }
  res[0] = sd[0] - 2;
  res[1] = sd[1] - 1;
}

static void rdfcn_i(double *ts, double *sd, double *sa, double *u, double *p, double *pr, double *res, long *dpnd, InfoPtr *info)
{
  if (*dpnd) { *dpnd = RFCN_DPND(0, *sd, 0, *u, 0, 0); return; }

  res[0] = u[0] + 1;
  res[1] = 1 - u[0];
  res[2] = u[1];
  res[3] = 1 - u[1];
  res[4] = u[1]*(sd[1] - sd[0]);
  res[5] = (1 - u[1])*(sd[0] - sd[1]);
  res[6] = u[1] - u[2];
  res[7] = u[1] + u[2];
  res[8] = 1 - u[1] - u[3];
  res[9] = 1 - u[1] + u[3];
}

extern "C" void def_model(void)
{
  def_mdims(NMOS, NP, NRC, NRCE);
  def_mstage(
    0,
    NXD, NXA, NU,
    mfcn, NULL,
    0, 0, 0, NULL, ffcn, NULL,
    NULL, NULL
  );
  def_mpc(0, "Start Point", NPR, NRD_S, NRDE_S, rdfcn_s, NULL);
  def_mpc(0, "Interior Point", NPR, NRD_I, NRDE_I, rdfcn_i, NULL);
}
