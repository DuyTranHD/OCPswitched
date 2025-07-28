/*  MUSCOD_SUITE/Apps/TEST/Src/SRC/dryfrictionpmrpfixopt.cpp
 *  (c) Tran, Bao Duy, June 2025 (corrected)
 *  Switched OCP with Dry Friction
 *  a point mass on a rough plane
 */

#include <cmath>
#include <stdio.h>
#include "def_usrmod.hpp"

#include <cstdlib>

#define  NMOS   3
#define  NP     0
#define  NRC    0
#define  NRCE   0

#define  NXD    4
#define  NXA    0
#define  NU     6
#define  NPR    0

#define  NRD_S  4
#define  NRDE_S 4

#define  NRD_E  4
#define  NRDE_E 4

#define  NRD_I  16
#define  NRDE_I 1

#define  m      1
#define  k      0.5
#define  g      9.8
#define  eps    0.00000001
#define  Upa    10
#define  gamma  M_PI/3
#define  x_0    10
#define  y_0    1
#define  vx_0   10
#define  vy_0   5
#define  x_T    90
#define  y_T    9
#define  vx_T   0
#define  vy_T   0

/* -------------------------------------------- */

static void mfcn(double *ts, double *sd, double *sa,
  double *p, double *pr, double *mval, long *dpnd, InfoPtr *info)
{
  if (*dpnd) {
    *dpnd = MFCN_DPND(*ts, *sd, 0, *p, 0);
    return;
  }

  *mval = *ts;
}

/*
   *  xd[i]  <-->  x_{i+1}, i=0,1    ... location
   *  xd[i]  <-->  v_{i-2}, i=2,3,   ... algebraic projection of velocity at i-th point
   *  u[0]  <-->   f_x               .   
   *  u[1]  <-->   f_y               .--> control force
   *  u[2]  <-->   f_z               ...
   *  u[3]  <-->   \alpha_{+}        .
   *  u[4]  <-->   \alpha_{0}        .--> relaxed controls
   *  u[5]  <-->   \alpha_{-}        . 
   */

static void ffcn1(double *t, double *xd, double *xa, double *u,
  double *p, double *rhs, double *rwh, long *iwh, InfoPtr *info)
{
  rhs[0] = xd[2];
  rhs[1] = xd[3];
  rhs[2] = (u[0]-k*(m*g*cos(gamma)-u[2])*u[0]/sqrt(u[0]*u[0]+(u[1]-m*g*sin(gamma))*(u[1]-m*g*sin(gamma)) + eps))/m;
  rhs[3] = (u[1]-m*g*sin(gamma)-k*(m*g*cos(gamma)-u[2])*(u[1]-m*g*sin(gamma))/sqrt(u[0]*u[0]+(u[1]-m*g*sin(gamma))*(u[1]-m*g*sin(gamma)) +eps))/m;
}

static void ffcn2(double *t, double *xd, double *xa, double *u,
  double *p, double *rhs, double *rwh, long *iwh, InfoPtr *info)
{
  rhs[0] = xd[2];
  rhs[1] = xd[3];
  rhs[2] = 0;
  rhs[3] = 0;
}

static void ffcn3(double *t, double *xd, double *xa, double *u,
  double *p, double *rhs, double *rwh, long *iwh, InfoPtr *info)
{
  rhs[0] = xd[2];
  rhs[1] = xd[3];
  rhs[2] = (u[0]-k*(m*g*cos(gamma)-u[2])*xd[2]/sqrt(xd[2]*xd[2]+xd[3]*xd[3] +eps))/m;
  rhs[3] = (u[1]-m*g*sin(gamma)-k*(m*g*cos(gamma)-u[2])*xd[3]/sqrt(xd[2]*xd[2]+xd[3]*xd[3] + eps))/m;
}

/* -------------------------------------------- */

static void rdfcn_s(double *ts, double *sd, double *sa, double *u,
  double *p, double *pr, double *res, long *dpnd, InfoPtr *info)
{
  if (*dpnd) {
    *dpnd = RFCN_DPND(0, *sd, 0, 0, 0, 0);
    return;
  }
  res[0] = sd[0] - x_0;
  res[1] = sd[1] - y_0;
  res[2] = sd[2] - vx_0;
  res[3] = sd[3] - vy_0;
}

static void rdfcn_i(double *ts, double *sd, double *sa, double *u,
  double *p, double *pr, double *res, long *dpnd, InfoPtr *info)
{
  if (*dpnd) { *dpnd = RFCN_DPND(0, *sd, 0, *u, 0, 0); return; }

  res[0] = 1 - u[3] - u[4] - u[5];
  res[1] = u[3]*sqrt(sd[2]*sd[2]+sd[3]*sd[3]) - eps;
  res[2] = u[4]*sqrt(sd[2]*sd[2]+sd[3]*sd[3]) + eps;
  res[3] = -u[4]*sqrt(sd[2]*sd[2]+sd[3]*sd[3]) + eps;
  res[4] = -u[4]*(sqrt(u[0]*u[0]+(u[1]-m*g*sin(gamma))*(u[1]-m*g*sin(gamma))) - k*(m*g*cos(gamma)-u[2]));
  res[5] = -u[5]*sqrt(sd[2]*sd[2]+sd[3]*sd[3]) + eps;
  res[6] = u[5]*sqrt(sd[2]*sd[2]+sd[3]*sd[3]) + eps;
  res[7] = u[5]*(sqrt(u[0]*u[0]+(u[1]-m*g*sin(gamma))*(u[1]-m*g*sin(gamma))) - k*(m*g*cos(gamma)-u[2])) + eps;
  res[8] = u[3];
  res[9] = 1 - u[3];
  res[10] = u[4];
  res[11] = 1 - u[4];
  res[12] = u[5];
  res[13] = 1 - u[5];
  res[14] = -u[0]*u[0] - u[1]*u[1] - u[2]*u[2] + Upa*Upa;
  res[15] = -u[2] + m*g*cos(gamma);
}

static void rdfcn_e(double *ts, double *sd, double *sa, double *u,
  double *p, double *pr, double *res, long *dpnd, InfoPtr *info)
{
  if (*dpnd) {
    *dpnd = RFCN_DPND(0, *sd, 0, 0, 0, 0);
    return;
  }
  res[0] = sd[0] - x_T;
  res[1] = sd[1] - y_T;
  res[2] = sd[2] - vx_T;
  res[3] = sd[3] - vy_T;
}


extern "C" void def_model(void)
{
  def_mdims(NMOS, NP, NRC, NRCE);

  def_mstage   ( 0,  NXD, NXA, NU,  mfcn, NULL,  0, 0, 0, NULL, ffcn1, NULL,  NULL, NULL );
  def_mstage   ( 1,  NXD, NXA, NU,  mfcn, NULL,  0, 0, 0, NULL, ffcn2, NULL,  NULL, NULL );
  def_mstage   ( 2,  NXD, NXA, NU,  mfcn, NULL,  0, 0, 0, NULL, ffcn3, NULL,  NULL, NULL );
  
  def_mpc(0, "Start Point", NPR, NRD_S, NRDE_S, rdfcn_s, NULL);
  
  def_mpc(1, "Interior Point", NPR, NRD_I, NRDE_I, rdfcn_i, NULL);
  
  def_mpc(2, "End Point", NPR, NRD_E, NRDE_E, rdfcn_e, NULL);
}
