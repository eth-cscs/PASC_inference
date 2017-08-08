#include "vegasconst.h"

__device__
double func(double* rx, double wgt)
{
   double value = 1.;
   for (int i=0;i<g_ndim;i++) {
      value *= 2.*rx[i];
   }
   return value;

}


