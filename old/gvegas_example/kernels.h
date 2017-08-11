#pragma once

#define Norum 0.0000000001*2.32830643653869628906
#define nLoopRndm (32)

__device__ __host__ __forceinline__
void fxorshift128(unsigned int seed,
                  int n,
                  double* a)
{
    unsigned int x,y,z,w,t;

    seed=seed*2357+123456789U;
    
    //1812433253 = 0x6C078965h;
    x=seed=1812433253U*(seed^(seed>>30))+1;
    y=seed=1812433253U*(seed^(seed>>30))+2;
    z=seed=1812433253U*(seed^(seed>>30))+3;
    w=seed=1812433253U*(seed^(seed>>30))+4;

    t=x^x<<11;x=y;y=z;z=w;w^=w>>19^t^t>>8;
    t=x^x<<11;x=y;y=z;z=w;w^=w>>19^t^t>>8;
    t=x^x<<11;x=y;y=z;z=w;w^=w>>19^t^t>>8;
    t=x^x<<11;x=y;y=z;z=w;w^=w>>19^t^t>>8;
    t=x^x<<11;x=y;y=z;z=w;w^=w>>19^t^t>>8;
    t=x^x<<11;x=y;y=z;z=w;w^=w>>19^t^t>>8;
    t=x^x<<11;x=y;y=z;z=w;w^=w>>19^t^t>>8;
    t=x^x<<11;x=y;y=z;z=w;w^=w>>19^t^t>>8;
    t=x^x<<11;x=y;y=z;z=w;w^=w>>19^t^t>>8;
    t=x^x<<11;x=y;y=z;z=w;w^=w>>19^t^t>>8;
    t=x^x<<11;x=y;y=z;z=w;w^=w>>19^t^t>>8;
    t=x^x<<11;x=y;y=z;z=w;w^=w>>19^t^t>>8;
    t=x^x<<11;x=y;y=z;z=w;w^=w>>19^t^t>>8;
    t=x^x<<11;x=y;y=z;z=w;w^=w>>19^t^t>>8;
    t=x^x<<11;x=y;y=z;z=w;w^=w>>19^t^t>>8;
    t=x^x<<11;x=y;y=z;z=w;w^=w>>19^t^t>>8;
    t=x^x<<11;x=y;y=z;z=w;w^=w>>19^t^t>>8;
    t=x^x<<11;x=y;y=z;z=w;w^=w>>19^t^t>>8;
    t=x^x<<11;x=y;y=z;z=w;w^=w>>19^t^t>>8;
    t=x^x<<11;x=y;y=z;z=w;w^=w>>19^t^t>>8;
    t=x^x<<11;x=y;y=z;z=w;w^=w>>19^t^t>>8;
    t=x^x<<11;x=y;y=z;z=w;w^=w>>19^t^t>>8;
    t=x^x<<11;x=y;y=z;z=w;w^=w>>19^t^t>>8;
    t=x^x<<11;x=y;y=z;z=w;w^=w>>19^t^t>>8;
    t=x^x<<11;x=y;y=z;z=w;w^=w>>19^t^t>>8;
    t=x^x<<11;x=y;y=z;z=w;w^=w>>19^t^t>>8;
    t=x^x<<11;x=y;y=z;z=w;w^=w>>19^t^t>>8;
    t=x^x<<11;x=y;y=z;z=w;w^=w>>19^t^t>>8;
    t=x^x<<11;x=y;y=z;z=w;w^=w>>19^t^t>>8;
    t=x^x<<11;x=y;y=z;z=w;w^=w>>19^t^t>>8;
    t=x^x<<11;x=y;y=z;z=w;w^=w>>19^t^t>>8;
    t=x^x<<11;x=y;y=z;z=w;w^=w>>19^t^t>>8;
    t=x^x<<11;x=y;y=z;z=w;w^=w>>19^t^t>>8;
    t=x^x<<11;x=y;y=z;z=w;w^=w>>19^t^t>>8;
    t=x^x<<11;x=y;y=z;z=w;w^=w>>19^t^t>>8;
    t=x^x<<11;x=y;y=z;z=w;w^=w>>19^t^t>>8;
    t=x^x<<11;x=y;y=z;z=w;w^=w>>19^t^t>>8;
    t=x^x<<11;x=y;y=z;z=w;w^=w>>19^t^t>>8;
    t=x^x<<11;x=y;y=z;z=w;w^=w>>19^t^t>>8;
    t=x^x<<11;x=y;y=z;z=w;w^=w>>19^t^t>>8;
    t=x^x<<11;x=y;y=z;z=w;w^=w>>19^t^t>>8;
    t=x^x<<11;x=y;y=z;z=w;w^=w>>19^t^t>>8;
    t=x^x<<11;x=y;y=z;z=w;w^=w>>19^t^t>>8;
    t=x^x<<11;x=y;y=z;z=w;w^=w>>19^t^t>>8;
    t=x^x<<11;x=y;y=z;z=w;w^=w>>19^t^t>>8;
    t=x^x<<11;x=y;y=z;z=w;w^=w>>19^t^t>>8;
    t=x^x<<11;x=y;y=z;z=w;w^=w>>19^t^t>>8;
    t=x^x<<11;x=y;y=z;z=w;w^=w>>19^t^t>>8;
    t=x^x<<11;x=y;y=z;z=w;w^=w>>19^t^t>>8;
    t=x^x<<11;x=y;y=z;z=w;w^=w>>19^t^t>>8;

    for(int i=0;i<n;++i){
        t=x^x<<11;
        x=y;
        y=z;
        z=w;
        w^=(w>>19)^(t^(t>>8));
        a[i]=w*Norum;
    }
    return;
}

__device__ __constant__ int g_ndim;
__device__ __constant__ int g_ng;
__device__ __constant__ int g_npg;
__device__ __constant__ int g_nd;
__device__ __constant__ double g_xjac;
__device__ __constant__ double g_dxg;
__device__ __constant__ double g_xl[ndim_max];
__device__ __constant__ double g_dx[ndim_max];
__device__ __constant__ double g_xi[ndim_max][nd_max];
__device__ __constant__ unsigned g_nCubes;

__device__
double func(double* rx, double wgt)
{
   double value = 1.;
   for (int i=0;i<g_ndim;i++) {
      value *= 2.*rx[i];
   }
   return value;

}

#include "gVegasCallFunc.cu"
#include "gVegas.cu"
