#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <sstream>

#include <ctime>
#include <sys/time.h>
#include <sys/resource.h>

// includes, project
#include <cutil_inline.h>
// include initial files

#define __MAIN_LOGIC
#include "vegas.h"
#include "gvegas.h"
#undef __MAIN_LOGIC

#include "kernels.h"


int DecodeInt(std::string numstr)
{

   if (numstr.size()<5) return 0;
   if (numstr.substr(0,1)!="e") return 0;
   if (numstr.substr(3,1)!="n") return 0;

   std::string expstr = numstr.substr(1,2);
   std::string manstr = numstr.substr(4,1);

   std::istringstream iss;
   iss.clear();
   iss.str(expstr);
   int exp;
   iss>>exp;

   iss.clear();
   iss.str(manstr);
   int man;
   iss>>man;

   return man*(int)(pow(10.,(double)exp)+0.5);
}

int main(int argc, char** argv)
{

   //
   // program interface:
   //   program -ncall="ncall0" -itmx="itmx0" -acc="acc0" -b="nBlockSize0"
   //
   // parameters:
   //   ncall0 = "exxny"
   //   ncall = y*10^xx
   //   itmx  = itmx0
   //   acc   = 0.01*acc0
   //   nBlockSize = nBlockSize0
   //

   //------------------
   //  Initialization
   //------------------

   int itmx0 = 10;
   int nBlockSize0 = 256;
   int GPUdevice = 0;

   float acc0 = 0.0001f;

   char* nCallStr = "e06n1";
   cutGetCmdLineArgumentstr(argc, (const char**)argv, "ncall", &nCallStr);
   cutGetCmdLineArgumenti(argc, (const char**)argv, "itmx", &itmx0);
   cutGetCmdLineArgumentf(argc, (const char**)argv, "acc", &acc0);
   cutGetCmdLineArgumenti(argc, (const char**)argv, "blk", &nBlockSize0);
   cutGetCmdLineArgumenti(argc, (const char**)argv, "dev", &GPUdevice);

   ncall = DecodeInt((std::string)nCallStr);
   itmx = itmx0;
   acc = 0.01*acc0;
   nBlockSize = nBlockSize0;

   cutilSafeCallNoSync(cudaSetDevice(GPUdevice));

   mds = 1;
   ndim = 8;
   
   ng = 0;
   npg = 0;

   for (int i=0;i<ndim;i++) {
      xl[i] = 0.;
      xu[i] = 1.;
   }
   
   nprn = 1;
//   nprn = -1;

   double startTotal, endTotal, timeTotal;
   timeTotal = 0.;
   startTotal = getrusage_usec();

   timeVegasCall = 0.;
   timeVegasMove = 0.;
   timeVegasFill = 0.;
   timeVegasRefine = 0.;

   double avgi = 0.;
   double sd = 0.;
   double chi2a = 0.;

   gVegas(avgi, sd, chi2a);

   endTotal = getrusage_usec();
   timeTotal = endTotal - startTotal;

   //-------------------------
   //  Print out information
   //-------------------------
   std::cout.clear();
   std::cout<<std::setw(10)<<std::setprecision(6)<<std::endl;
   std::cout<<"#============================="<<std::endl;
   std::cout<<"# No. of Thread Block Size  : "<<nBlockSize<<std::endl;
   std::cout<<"#============================="<<std::endl;
   std::cout<<"# No. of dimensions         : "<<ndim<<std::endl;
   std::cout<<"# No. of func calls / iter  : "<<ncall<<std::endl;
   std::cout<<"# No. of max. iterations    : "<<itmx<<std::endl;
   std::cout<<"# Desired accuracy          : "<<acc<<std::endl;
   std::cout<<"#============================="<<std::endl;
   std::cout<<std::scientific;
   std::cout<<std::left<<std::setfill(' ');
   std::cout<<"# Result                    : "
            <<std::setw(12)<<std::setprecision(5)<<avgi<<" +- "
            <<std::setw(12)<<std::setprecision(5)<<sd<<" ( "
            <<std::setw(7)<<std::setprecision(4)
            <<std::fixed<<100.*sd/avgi<<"%)"<<std::endl;
   std::cout<<std::fixed;
   std::cout<<"# Chisquare                 : "<<std::setprecision(4)
            <<chi2a<<std::endl;
   std::cout<<"#============================="<<std::endl;
   std::cout<<std::right;
   std::cout<<"# Total Execution Time(sec) : "
            <<std::setw(10)<<std::setprecision(4)<<timeTotal<<std::endl;
   std::cout<<"#============================="<<std::endl;
   std::cout<<"# Time for func calls (sec) : "
            <<std::setw(10)<<std::setprecision(4)<<timeVegasCall
            <<" ( "<<std::setw(5)<<std::setprecision(2)
            <<100.*timeVegasCall/timeTotal<<"%)"<<std::endl;
   std::cout<<"# Time for data transf (sec): "
            <<std::setw(10)<<std::setprecision(4)<<timeVegasMove
            <<" ( "<<std::setw(5)<<std::setprecision(2)
            <<100.*timeVegasMove/timeTotal<<"%)"<<std::endl;
   std::cout<<"# Time for data fill (sec)  : "
            <<std::setw(10)<<std::setprecision(4)<<timeVegasFill
            <<" ( "<<std::setw(5)<<std::setprecision(2)
            <<100.*timeVegasFill/timeTotal<<"%)"<<std::endl;
   std::cout<<"# Time for grid refine (sec): "
            <<std::setw(10)<<std::setprecision(4)<<timeVegasRefine
            <<" ( "<<std::setw(5)<<std::setprecision(2)
            <<100.*timeVegasRefine/timeTotal<<"%)"<<std::endl;
   std::cout<<"#============================="<<std::endl;

   cudaThreadExit();

   return 0;
}
