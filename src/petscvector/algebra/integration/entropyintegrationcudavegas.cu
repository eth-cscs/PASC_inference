#include "external/petscvector/petscvector.cuh"
#include "external/petscvector/algebra/integration/entropyintegrationcudavegas.h"

namespace pascinference {
namespace algebra {

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


EntropyIntegrationCudaVegas<PetscVector>::ExternalContent::ExternalContent() {
	LOG_FUNC_BEGIN

	/* restart timers */
	this->timeVegasCall = 0.0;
	this->timeVegasMove = 0.0;
	this->timeVegasFill = 0.0;
	this->timeVegasRefine = 0.0;

	LOG_FUNC_END
}

void EntropyIntegrationCudaVegas<PetscVector>::ExternalContent::cuda_gVegas(double &avgi, double &sd, double &chi2a) {
	LOG_FUNC_BEGIN

	int nd_max = 50;
	int mds = 1;

	int it;
	int nd;
	int ng;
	int npg;
	int nCubes;
	double xi[this->ndim][nd_max];
	




	for (int j=0; j < this->ndim; j++) {
		xi[j][0] = 1.;
	}

	/* entry vegas1 */
	it = 0;

	/* entry vegas2 */
	nd = nd_max;
	ng = 1;
   
	npg = 0;
	if (mds!=0) {
		ng = (int)pow((0.5*(double)(this->ncall)),1./(double)(this->ndim));
		mds = 1;
		if (2*ng>=nd_max) {
			mds = -1;
			npg = ng/(double)nd_max+1;
			nd = ng/(double)npg;
			ng = npg*nd;
		}
	}

   gpuErrchk(cudaMemcpyToSymbol(g_ndim, &ndim, sizeof(int)));
   gpuErrchk(cudaMemcpyToSymbol(g_ng,   &ng,   sizeof(int)));
   gpuErrchk(cudaMemcpyToSymbol(g_nd,   &nd,   sizeof(int)));
   cudaThreadSynchronize(); /* wait for synchronize */

   nCubes = (unsigned)(pow(ng,this->ndim));
   gpuErrchk(cudaMemcpyToSymbol(g_nCubes, &nCubes, sizeof(nCubes)));
   cudaThreadSynchronize(); /* wait for synchronize */

   npg = ncall/(double)nCubes;
   if(npg < 2){
	   npg = 2;
   }
   calls = (double)(npg*nCubes);

   unsigned nCubeNpg = nCubes*npg;

   if (nprn!=0) {
      coutMaster << std::endl;
      coutMaster << " << vegas internal parameters >> " << std::endl;
      coutMaster << "            ng: " << std::setw(5) << ng << std::endl;
      coutMaster << "            nd: " << std::setw(5) << nd << std::endl;
      coutMaster << "           npg: " << std::setw(5) << npg << std::endl;
      coutMaster << "        nCubes: " << std::setw(12) << nCubes << std::endl;
      coutMaster << "    nCubes*npg: " << std::setw(12) << nCubeNpg << std::endl;
   }





	avgi = 11.1;
	sd = 22.2;
	chi2a = 33.33;



	LOG_FUNC_END
}



}
}

