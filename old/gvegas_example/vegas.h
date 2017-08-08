#pragma once

void gVegas(double& avgi, double& sd, double& chi2a);

#ifndef __MAIN_LOGIC
#define EXTERN extern
#else
#define EXTERN
#endif

const int ndim_max = 20;
const double alph = 1.5;
EXTERN double dx[ndim_max];
EXTERN double randm[ndim_max];
const int nd_max = 50;
EXTERN double xin[nd_max];

EXTERN double xjac;

EXTERN double xl[ndim_max],xu[ndim_max];
EXTERN double acc;
EXTERN int ndim, ncall, itmx, nprn;

EXTERN double xi[ndim_max][nd_max];
EXTERN double si, si2, swgt, schi;
EXTERN int ndo, it;

//EXTERN double alph;
EXTERN int mds;

EXTERN double calls, ti, tsi;

EXTERN int npg, ng, nd;
EXTERN double dxg, xnd;

EXTERN unsigned nCubes;

#undef EXTERN

