#include "external/petscvector/petscvector.cuh"
#include "external/petscvector/algebra/integration/entropyintegrationcudavegas.h"

namespace pascinference {
namespace algebra {

const int ndim_max = 10;
const int nd_max = 50;

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

	int mds = 1;
	int nprn = 1;

	int it;
	int nd;
	int ng;
	int ndo;
	int npg;
	int nCubes;
	double calls;
	double dxg;
	double dnpg;
	double dv2g;
	
	int nGridSizeX, nGridSizeY;
	int nBlockTot;
	
	double xi[ndim_max][nd_max];
	double xl[ndim_max],xu[ndim_max];
	double dx[ndim_max];

	double xnd;
	double xjac;
	
	double si;
	double si2;
	double swgt;
	double schi;
	
	for (int i=0;i< this->ndim;i++) {
		xl[i] = 0.;
		xu[i] = 1.;
	}

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

	gpuErrchk(cudaMemcpyToSymbol(g_ndim, &(this->ndim), sizeof(int)));
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

	dxg = 1./(double)ng;
	dnpg = (double)npg;
	dv2g = calls*calls*pow(dxg,this->ndim)*pow(dxg,this->ndim)/(dnpg*dnpg*(dnpg-1.));
	xnd = (double)nd;
	dxg *= xnd;
	xjac = 1./(double)calls;
	for (int j=0;j<this->ndim;j++) {
		dx[j] = xu[j]-xl[j];
		xjac *= dx[j];
	}

	/* tranfer data to GPU */
	gpuErrchk(cudaMemcpyToSymbol(g_npg,  &npg,  sizeof(int)));
	gpuErrchk(cudaMemcpyToSymbol(g_xjac, &xjac, sizeof(double)));
	gpuErrchk(cudaMemcpyToSymbol(g_dxg,  &dxg,  sizeof(double)));
	cudaThreadSynchronize();

	ndo = 1;

	if (nd!=ndo) {
		double rc = (double)ndo/xnd;

		for (int j=0;j<this->ndim;j++) {
			int k = -1;
			double xn = 0.;
			double dr = 0.;
			int i = k;
			k++;
			dr += 1.;
			double xo = xn;
			xn = xi[j][k];

			while (i<nd-1) {
				while (dr<=rc) {
					k++;
					dr += 1.;
					xo = xn;
					xn = xi[j][k];
				}
				i++;
				dr -= rc;
				xin[i] = xn - (xn-xo)*dr;
			}

			for (int i=0;i<nd-1;i++) {
				xi[j][i] = (double)xin[i];
			}
			xi[j][nd-1] = 1.;
		}
		ndo = nd;
	}

	/* transfer data to GPU */
	gpuErrchk(cudaMemcpyToSymbol(g_xl, xl, sizeof(xl)));
	gpuErrchk(cudaMemcpyToSymbol(g_dx, dx, sizeof(dx)));
	gpuErrchk(cudaMemcpyToSymbol(g_xi, xi, sizeof(xi)));
	cudaThreadSynchronize();
	
	if (nprn!=0) {
		coutMaster << std::endl;
		coutMaster << " << input parameters for vegas >>" << std::endl;
		coutMaster << "     ndim =" << std::setw(3) << this->ndim
					<< "   ncall = " << std::setw(10) << this->calls <<std::endl;
		coutMaster << "     it   =  0"
					<< "   itmx = " << std::setw(5) << this->itmx << std::endl;
		coutMaster << "     acc  = " << std::fixed
					<< std::setw(9) << std::setprecision(3) << this->acc << std::endl;
		coutMaster << "     mds  = " << std::setw(3) << mds
					<< "   nd = " << std::setw(4) << nd <<std::endl;
		for(int j=0; j < this->ndim; j++){
			coutMaster << "    (xl,xu)= ( " << std::setw(6) << std::fixed
						<< xl[j] << ", " << xu[j] << " )" << std::endl;
		}
	}	

	/* entry vegas3 */
	it = 0;
	si = 0.;
	si2 = 0.;
	swgt = 0.;
	schi = 0.;

	/* --------------------------
	 * Set up kernel vaiables
     * --------------------------
     */
	const int nGridSizeMax =  65535;
   
	dim3 ThBk(nBlockSize);
	
	nBlockTot = (nCubeNpg-1)/nBlockSize+1;
	nGridSizeY = (nBlockTot-1)/nGridSizeMax+1;
	nGridSizeX = (nBlockTot-1)/nGridSizeY+1;
	dim3 BkGd(nGridSizeX, nGridSizeY);

	if (nprn!=0) {
		coutMaster << std::endl;
		coutMaster << " << kernel parameters for CUDA >> " << std::endl;
		coutMaster << "       Block size           = " << std::setw(7) << ThBk.x << std::endl;
		coutMaster << "       Grid size            = " << std::setw(7) << BkGd.x
					<< " x " << BkGd.y << std::endl;
		int nThreadsTot = ThBk.x*BkGd.x*BkGd.y;
		coutMaster << "     Actual Number of calls = " << std::setw(12)
					<< nThreadsTot << std::endl;
		coutMaster << "   Required Number of calls = " << std::setw(12)
					<< nCubeNpg << " ( " << std::setw(6) << std::setprecision(2)
					<< 100.*(double)nCubeNpg/(double)nThreadsTot << "%)" <<std::endl;
		coutMaster << std::endl;
	}

	int sizeFval;
	double* hFval;
	double* gFval;
	
	int sizeIAval;
	int* hIAval;
	int* gIAval;

	double startVegasCall, endVegasCall;
	double startVegasMove, endVegasMove;
	double startVegasFill, endVegasFill;
	double startVegasRefine, endVegasRefine;

	/* allocate Fval */
	sizeFval = nCubeNpg*sizeof(double);

	/* CPU */
	gpuErrchk(cudaMallocHost((void**)&hFval, sizeFval));
	memset(hFval, '\0', sizeFval);

	/* GPU */
	gpuErrchk(cudaMalloc((void**)&gFval, sizeFval));

	/* allocate IAval */
	sizeIAval = nCubeNpg*ndim*sizeof(int);

	/* CPU */
	gpuErrchk(cudaMallocHost((void**)&hIAval, sizeIAval));
	memset(hIAval, '\0', sizeIAval);

	/* GPU */
	gpuErrchk(cudaMalloc((void**)&gIAval, sizeIAval));

	/* perform main iterations */
	


	gpuErrchk(cudaFreeHost(hFval));
	gpuErrchk(cudaFree(gFval));

	gpuErrchk(cudaFreeHost(hIAval));
	gpuErrchk(cudaFree(gIAval));



	avgi = 11.1;
	sd = 22.2;
	chi2a = 33.33;



	LOG_FUNC_END
}



}
}

