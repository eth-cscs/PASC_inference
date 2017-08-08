#include "const.h"

#include "vegasconst.h"
#include "vegas.h"

__global__
void gVegasCallFunc(double* gFval, int* gIAval)
{
   //--------------------
   // Check the thread ID
   //--------------------
   const unsigned int tIdx  = threadIdx.x;
   const unsigned int bDimx = blockDim.x;

   const unsigned int bIdx  = blockIdx.x;
   const unsigned int gDimx = gridDim.x;
   const unsigned int bIdy  = blockIdx.y;

   unsigned int bid  = bIdy*gDimx+bIdx;
   const unsigned int tid = bid*bDimx+tIdx;

   int ig = tid/g_npg;

   unsigned nCubeNpg = g_nCubes*g_npg;

   if (tid<nCubeNpg) {

      unsigned ia[ndim_max];
      
      unsigned int tidRndm = tid;
      
      int kg[ndim_max];
      
      unsigned igg = ig;
      for (int j=0;j<g_ndim;j++) {
         kg[j] = igg%g_ng+1;
         igg /= g_ng;
      }
      
      double randm[ndim_max];
      fxorshift128(tidRndm, g_ndim, randm);
      
      double x[ndim_max];
      
      double wgt = g_xjac;
      for (int j=0;j<g_ndim;j++) {
         double xo,xn,rc;
         xn = (kg[j]-randm[j])*g_dxg+1.;
         ia[j] = (int)xn-1;
         if (ia[j]<=0) {
            xo = g_xi[j][ia[j]];
            rc = (xn-(double)(ia[j]+1))*xo;
         } else {
            xo = g_xi[j][ia[j]]-g_xi[j][ia[j]-1];
            rc = g_xi[j][ia[j]-1]+(xn-(double)(ia[j]+1))*xo;
         }
         x[j] = g_xl[j]+rc*g_dx[j];
         wgt *= xo*(double)g_nd;
      }
      
      double f = wgt * func(x,wgt);
      
      gFval[tid] = f;
      for (int idim=0;idim<g_ndim;idim++) {
         gIAval[idim*nCubeNpg+tid] = ia[idim];
      }
   }

}
