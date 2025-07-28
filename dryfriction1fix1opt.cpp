/*  MUSCOD_SUITE/Apps/TEST/Src/SRC/dryfriction1fix1opt.cpp
 *  (c) Tran, B Duy, May 2025
 *  Switched OCP with Dry Friction
 *  a system of 3 material points lie along horizontal straight line
 */

#include <cmath>
#include <stdio.h>
#include "def_usrmod.hpp"

#include <cstdlib>

#define  NMOS   3
#define  NP     0
#define  NRC    0
#define  NRCE   0

#define  NXD    6
#define  NXA    0
#define  NU     17
#define  NPR    0

#define  NRD_S  6
#define  NRDE_S 6

#define  NRD_E  6
#define  NRDE_E 6

#define  NRD_I  69
#define  NRDE_I 3

#define  sgn(v) ( ( (v) < 0 ) ? -1 : ( (v) > 0 ) ? 1 : 0)

#define  m      1
#define  k      0.5
#define  g      9.8
#define  eps    1e-7

static void mfcn(double *ts, double *sd, double *sa, double *p, double *pr, double *mval, long *dpnd, InfoPtr *info)
{
  if (*dpnd) {
    *dpnd = MFCN_DPND(0, *sd, 0, 0, 0);
    return;
  }

  *mval = -(sd[0] + sd[1] + sd[2])/3;
}

/*
   *  xd[i]  <-->  x_{i+1}, i=0,1,2  ... location
   *  xd[i]  <-->  v_{i-2}, i=3,4,5  ... algebraic projection of velocity at i-th point
   *  u[0]  <-->   \alpha_{1,1}      .  
   *  u[1]  <-->   \alpha_{1,2}      .
   *  u[2]  <-->   \alpha_{1,3}      .
   *  u[3]  <-->   \alpha_{2,1}      .
   *  u[4]  <-->   \alpha_{2,2}      .--> relaxed controls
   *  u[5]  <-->   \alpha_{2,3}      .
   *  u[6]  <-->   \alpha_{3,1}      .
   *  u[7]  <-->   \alpha_{3,2}      .
   *  u[8]  <-->   \alpha_{3,3}      .
   *  u[9]  <-->   \alpha_{4,1}      .
   *  u[10]  <-->   \alpha_{4,2}     .
   *  u[11]  <-->   \alpha_{4,3}     .
   *  u[12]  <-->   \alpha_{5,1}     .
   *  u[13]  <-->   \alpha_{5,2}     .
   *  u[14]  <-->   \alpha_{5,3}     .
   *  u[15]  <-->   f_{1}, i=1   ... algebraic projection of control force acting from 
   *  u[16] <-->   f_{2}            (i+1)th point to i-th point, where f_0 = f_3 = 0.  
   */
   
static void ffcn0(double *t, double *xd, double *xa, double *u,
  double *p, double *rhs, double *rwh, long *iwh, InfoPtr *info)
{
  rhs[0] = xd[3];
  rhs[1] = xd[4];
  rhs[2] = xd[5];
  rhs[3] = 0;
  rhs[4] = 0;
  rhs[5] = 0;
}

static void ffcn2(double *t, double *xd, double *xa, double *u,
  double *p, double *rhs, double *rwh, long *iwh, InfoPtr *info)
{
  rhs[0] = xd[3];
  rhs[1] = xd[4];
  rhs[2] = xd[5];
  rhs[3] = -u[15]/m + k*g*sgn(xd[3]);
  rhs[4] = (u[15]-u[16])/m + k*g*sgn(xd[4]);
  rhs[5] = u[16]/m + k*g*sgn(xd[5]);
}

static void ffcn1(double *t, double *xd, double *xa, double *u,
  double *p, double *rhs, double *rwh, long *iwh, InfoPtr *info)
{
  rhs[0] = xd[3];
  rhs[1] = xd[4];
  rhs[2] = xd[5];
  rhs[3] = -u[15]/m - k*g*sgn(-u[15]);
  rhs[4] = (u[15]-u[16])/m - k*g*sgn(u[15]-u[16]);
  rhs[5] = u[16]/m - k*g*sgn(u[16]);
}


static void rdfcn_s(double *ts, double *sd, double *sa, double *u,
  double *p, double *pr, double *res, long *dpnd, InfoPtr *info)
{
  if (*dpnd) {
    *dpnd = RFCN_DPND(0, *sd, 0, 0, 0, 0);
    return;
  }
  res[0] = sd[0] - 0.000000001;
  res[1] = sd[1] - 0.000000001;
  res[2] = sd[2] - 0.000000001;
  res[3] = sd[3] - 0.000000001;
  res[4] = sd[4] - 0.000000001;
  res[5] = sd[5] - 0.000000001;
}

static void rdfcn_e(double *ts, double *sd, double *sa, double *u,
  double *p, double *pr, double *res, long *dpnd, InfoPtr *info)
{
  if (*dpnd) {
    *dpnd = RFCN_DPND(0, *sd, 0, 0, 0, 0);
    return;
  }
  res[0] = sd[0] - 10;
  res[1] = sd[1] - 10;
  res[2] = sd[2] - 10;
  res[3] = sd[3] - 0.000000001;
  res[4] = sd[4] - 0.000000001;
  res[5] = sd[5] - 0.000000001;
}

static void rdfcn_i(double *ts, double *sd, double *sa, double *u, double *p, double *pr, double *res, long *dpnd, InfoPtr *info)
{
  if (*dpnd) { *dpnd = RFCN_DPND(0, *sd, 0, *u, 0, 0); return; }
  res[0] = 1 - u[0] - u[3] - u[6] - u[9] - u[12];
  res[1] = 1 - u[1] - u[4] - u[7] - u[10] - u[13];
  res[2] = 1 - u[2] - u[5] - u[8] - u[11] - u[14];
  res[3] = -u[6]*sd[3] + eps;
  res[4] = -u[7]*sd[4] + eps;
  res[5] = -u[8]*sd[5] + eps;
  res[6] = u[6]*sd[3] + eps;
  res[7] = u[7]*sd[4] + eps;
  res[8] = u[8]*sd[5] + eps;
  res[9] = -u[9]*sd[3] + eps;
  res[10] = -u[10]*sd[4] + eps;
  res[11] = -u[11]*sd[5] + eps;
  res[12] = u[9]*sd[3] + eps;
  res[13] = u[10]*sd[4] + eps;
  res[14] = u[11]*sd[5] + eps;
  res[15] = -u[12]*sd[3] + eps;
  res[16] = -u[13]*sd[4] + eps;
  res[17] = -u[14]*sd[5] + eps;  
  res[18] = u[12]*sd[3] + eps;
  res[19] = u[13]*sd[4] + eps;
  res[20] = u[14]*sd[5] + eps;
  res[21] = u[0];
  res[22] = 1 - u[0];
  res[23] = u[1];
  res[24] = 1 - u[1];
  res[25] = u[2];
  res[26] = 1 - u[2];
  res[27] = u[3];
  res[28] = 1 - u[3];
  res[29] = u[4];
  res[30] = 1 - u[4];
  res[31] = u[5];
  res[32] = 1 - u[5];
  res[33] = u[6];
  res[34] = 1 - u[6];
  res[35] = u[7];
  res[36] = 1 - u[7];
  res[37] = u[8];
  res[38] = 1 - u[8];
  res[39] = u[9];
  res[40] = 1 - u[9];
  res[41] = u[10];
  res[42] = 1 - u[10];
  res[43] = u[11];
  res[44] = 1 - u[11];
  res[45] = u[12];
  res[46] = 1 - u[12];
  res[47] = u[13];
  res[48] = 1 - u[13];
  res[49] = u[14];
  res[50] = 1 - u[14];
  res[51] = u[0]*sd[3] + eps;
  res[52] = u[1]*sd[4] + eps;
  res[53] = u[2]*sd[5] + eps;
  res[54] = -u[3]*sd[3] + eps;
  res[55] = -u[4]*sd[4] + eps;
  res[56] = -u[5]*sd[5] + eps;
  res[57] = -u[6]*(k*m*g + u[15]) + eps;
  res[58] = -u[7]*(k*m*g - u[15] + u[16]) + eps;
  res[59] = -u[8]*(k*m*g - u[16]) + eps;
  res[60] = u[6]*(-k*m*g + u[15]) + eps;
  res[61] = u[7]*(-k*m*g - u[15] + u[16]) + eps;
  res[62] = u[8]*(-k*m*g - u[16]) + eps;
  res[63] = u[9]*(-u[15]-k*m*g) + eps;
  res[64] = u[10]*(u[15]-u[16]-k*m*g) + eps;
  res[65] = u[11]*(u[16]-k*m*g) + eps;
  res[66] = -u[12]*(u[15]+k*m*g) + eps;
  res[67] = -u[13]*(u[15]-u[16]+k*m*g) + eps;
  res[68] = -u[14]*(-u[16]+k*m*g) + eps;
}

extern "C" void def_model(void)
{
  def_mdims(NMOS, NP, NRC, NRCE);
  
  def_mstage(0, NXD, NXA, NU, mfcn, NULL, 0, 0, 0, NULL, ffcn0, NULL, NULL, NULL);
  def_mstage(1, NXD, NXA, NU, mfcn, NULL, 0, 0, 0, NULL, ffcn1, NULL, NULL, NULL);
  def_mstage(2, NXD, NXA, NU, mfcn, NULL, 0, 0, 0, NULL, ffcn2, NULL, NULL, NULL);
  
  def_mpc(0, "Start Point", NPR, NRD_S, NRDE_S, rdfcn_s, NULL);
  def_mpc(1, "Interior Point", NPR, NRD_I, NRDE_I, rdfcn_i, NULL);
  def_mpc(2, "End Point", NPR, NRD_E, NRDE_E, rdfcn_e, NULL);
}
