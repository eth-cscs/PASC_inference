#include "external/petscvector/petscvector.cuh"
#include "external/petscvector/algebra/integration/entropyintegrationcudavegas.h"

#define Norum 0.0000000001*2.32830643653869628906

namespace pascinference {
namespace algebra {

const int xdim_max = 10;
const int nd_max = 50;

__device__ __constant__ int g_ng;
__device__ __constant__ int g_npg;
__device__ __constant__ int g_nd;
__device__ __constant__ double g_xjac;
__device__ __constant__ double g_dxg;
__device__ __constant__ double g_xl[xdim_max];
__device__ __constant__ double g_dx[xdim_max];
__device__ __constant__ double g_xi[xdim_max][nd_max];
__device__ __constant__ unsigned g_nCubes;

__global__ void gVegasCallFunc(double* gFval, int* gIAval, int xdim, int number_of_integrals, int number_of_moments, double *g_lambda, int *g_matrix_D_arr, double *g_fs);
__device__ void func_entropy(double *cvalues_out, double *cx, int xdim, int number_of_integrals, int number_of_moments, double *g_lambda, int *g_matrix_D_arr);
__device__ __host__ __forceinline__ void fxorshift128(unsigned int seed, int n, double* a);

__global__ void print_kernel(int xdim, int number_of_moments, int number_of_integrals, double *g_lambda, int *g_matrix_D_arr){
	printf("============ from GPU =========== \n");
	printf("xdim                : %d \n", xdim);
	printf("number_of_moments   : %d \n", number_of_moments);
	printf("number_of_integrals : %d \n", number_of_integrals);
	printf("g_lambda            : [");
	for(int i=0; i < number_of_moments-1; i++){
		printf("%f", g_lambda[i]);
		if(i<number_of_moments-2){
			printf(", ");
		}
	}
	printf("] \n");
	printf("matrix D: \n");
	for(int i_moment=0;i_moment < number_of_moments; i_moment++){
		for(int i_xdim=0;i_xdim < xdim; i_xdim++){
			printf("%d", g_matrix_D_arr[ i_moment*xdim + i_xdim]);
			
			if(i_xdim < xdim-1) printf(", ");
		}
		printf("\n");
	}
}

__global__ void kernel_set_zero(int size, double *arr1){
	
	//TODO: make it in a different way!!!
	for(int i=0;i<size;i++){
		arr1[i] = 0.0;
	}

}

EntropyIntegrationCudaVegas<PetscVector>::ExternalContent::ExternalContent(int xdim, int number_of_moments, int number_of_integrals, int *matrix_D_arr) {
	LOG_FUNC_BEGIN

	/* restart timers */
	this->timerVegasCall.restart();
	this->timerVegasMove.restart();
	this->timerVegasFill.restart();
	this->timerVegasRefine.restart();

	/* store host variables */
	this->xdim = xdim;
	this->number_of_moments = number_of_moments;
	this->number_of_integrals = number_of_integrals;
	this->matrix_D_arr = matrix_D_arr;

	/* allocate host arrays */
	this->sd = new double[this->get_number_of_integrals()];
	this->chi2a = new double[this->get_number_of_integrals()];

	/* allocate variables on cuda */
	gpuErrchk(cudaMalloc((void**)&(this->g_lambda), (number_of_moments-1)*sizeof(double)));
	gpuErrchk(cudaMalloc((void**)&(this->g_matrix_D_arr), number_of_moments*xdim*sizeof(int)));
	gpuErrchk(cudaMalloc((void**)&(this->g_fs), number_of_integrals*sizeof(int)));

	/* copy variables to CUDA */
	gpuErrchk( cudaMemcpy(this->g_matrix_D_arr, matrix_D_arr, number_of_moments*xdim*sizeof(int), cudaMemcpyHostToDevice ) );
	cudaThreadSynchronize(); /* wait for synchronize */


	LOG_FUNC_END
}

EntropyIntegrationCudaVegas<PetscVector>::ExternalContent::~ExternalContent(){
	LOG_FUNC_BEGIN

	/* free host arrays */
	free(this->sd);
	free(this->chi2a);

	/* free cuda variables */
	gpuErrchk(cudaFree(this->g_fs));
	gpuErrchk(cudaFree(this->g_lambda));
	gpuErrchk(cudaFree(this->g_matrix_D_arr));

	LOG_FUNC_END
}


void EntropyIntegrationCudaVegas<PetscVector>::ExternalContent::cuda_gVegas(double *avgi, double *lambda_arr) {
	LOG_FUNC_BEGIN

	/* copy given lambda to GPU */
	gpuErrchk( cudaMemcpy(this->g_lambda, lambda_arr, (number_of_moments-1)*sizeof(double), cudaMemcpyHostToDevice ) );
	cudaThreadSynchronize(); /* wait for synchronize */

	//TODO: temp
	print_kernel<<<1, 1>>>(xdim, number_of_moments, number_of_integrals, this->g_lambda, this->g_matrix_D_arr);
	cudaThreadSynchronize(); /* wait for synchronize */

	int mds = 1;
	int nprn = 1;
	const double alph = 1.5;

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
	double ti;
	double tsi;
	
	int nGridSizeX, nGridSizeY;
	int nBlockTot;
	
	double xi[xdim_max][nd_max];
	double xl[xdim_max],xu[xdim_max];
	double dx[xdim_max];
	double xin[nd_max];

	double xnd;
	double xjac;
	
	double si;
	double si2;
	double swgt;
	double schi;
	
	for (int i=0;i< this->xdim;i++) {
		xl[i] = 0.;
		xu[i] = 1.;
	}

	for (int j=0; j < this->xdim; j++) {
		xi[j][0] = 1.;
	}

	/* entry vegas1 */
	it = 0;

	/* entry vegas2 */
	nd = nd_max;
	ng = 1;
   
	npg = 0;
	if (mds!=0) {
		ng = (int)pow((0.5*(double)(this->ncall)),1./(double)(this->xdim));
		mds = 1;
		if (2*ng>=nd_max) {
			mds = -1;
			npg = ng/(double)nd_max+1;
			nd = ng/(double)npg;
			ng = npg*nd;
		}
	}

	gpuErrchk(cudaMemcpyToSymbol(g_ng,   &ng,   sizeof(int)));
	gpuErrchk(cudaMemcpyToSymbol(g_nd,   &nd,   sizeof(int)));
	cudaThreadSynchronize(); /* wait for synchronize */

	nCubes = (unsigned)(pow(ng,this->xdim));
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
		coutMaster << "            ng: " << std::setw(12) << ng << std::endl;
		coutMaster << "            nd: " << std::setw(12) << nd << std::endl;
		coutMaster << "           npg: " << std::setw(12) << npg << std::endl;
		coutMaster << "        nCubes: " << std::setw(12) << nCubes << std::endl;
		coutMaster << "    nCubes*npg: " << std::setw(12) << nCubeNpg << std::endl;
	}

	dxg = 1./(double)ng;
	dnpg = (double)npg;
	dv2g = calls*calls*pow(dxg,this->xdim)*pow(dxg,this->xdim)/(dnpg*dnpg*(dnpg-1.));
	xnd = (double)nd;
	dxg *= xnd;
	xjac = 1./(double)calls;
	for (int j=0;j<this->xdim;j++) {
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

		for (int j=0;j<this->xdim;j++) {
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
		coutMaster << "     xdim =" << std::setw(3) << this->xdim
					<< "   ncall = " << std::setw(10) << this->ncall <<std::endl;
		coutMaster << "     it   =  0"
					<< "   itmx = " << std::setw(5) << this->itmx << std::endl;
		coutMaster << "     acc  = " << std::fixed
					<< std::setw(9) << std::setprecision(3) << this->acc << std::endl;
		coutMaster << "     mds  = " << std::setw(3) << mds
					<< "   nd = " << std::setw(4) << nd <<std::endl;
		for(int j=0; j < this->xdim; j++){
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
		coutMaster << "       Block size           = " << std::setw(12) << ThBk.x << std::endl;
		coutMaster << "       Grid size            = " << std::setw(12) << BkGd.x
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

	/* allocate Fval */
	sizeFval = nCubeNpg*sizeof(double);

	/* CPU */
	gpuErrchk(cudaMallocHost((void**)&hFval, sizeFval));
	memset(hFval, '\0', sizeFval);

	/* GPU */
	gpuErrchk(cudaMalloc((void**)&gFval, sizeFval*this->number_of_integrals));

	/* allocate IAval */
	sizeIAval = nCubeNpg*xdim*sizeof(int);

	/* CPU */
	gpuErrchk(cudaMallocHost((void**)&hIAval, sizeIAval));
	memset(hIAval, '\0', sizeIAval);

	/* GPU */
	gpuErrchk(cudaMalloc((void**)&gIAval, sizeIAval));

	/* perform main iterations */
	do {
		it++;
		
		/* call integral function */
		timerVegasCall.start();
		 /* double* gFval, int* gIAval, int xdim, int number_of_integrals, int number_of_moments, int *g_lambda, int *g_matrix_D_arr */
		gVegasCallFunc<<<BkGd, ThBk>>>(gFval, gIAval, this->xdim, this->number_of_integrals, this->number_of_moments, this->g_lambda, this->g_matrix_D_arr, this->g_fs);
		cudaThreadSynchronize();
		timerVegasCall.stop();

		/* move computed results */
		timerVegasMove.start();
		gpuErrchk(cudaMemcpy(hFval, gFval,  sizeFval*this->number_of_integrals,
                               cudaMemcpyDeviceToHost));
		gpuErrchk(cudaMemcpy(hIAval, gIAval,  sizeIAval,
                               cudaMemcpyDeviceToHost));
		timerVegasMove.stop();

		/* fill */
		timerVegasFill.start();

		ti = 0.;
		tsi = 0.;

		double d[xdim_max][nd_max];

		for (int j=0;j<xdim;++j) {
			for (int i=0;i<nd;++i) {
				d[j][i] = 0.;
			}
		}

		for (unsigned ig=0;ig<nCubes;ig++) {
			double fb = 0.;
			double f2b = 0.;
			for(int ipg=0;ipg<npg;ipg++) {
				int idx = npg*ig+ipg;
				double f = hFval[idx];
				double f2 = f*f;
				fb += f;
				f2b += f2;
			}
			f2b = sqrt(f2b*npg);
			f2b = (f2b-fb)*(f2b+fb);
			ti += fb;
			tsi += f2b;
			if(mds<0){
				int idx = npg*ig;
				for(int idim=0;idim<xdim;idim++) {
					int iaj = hIAval[idim*nCubeNpg+idx];
					d[idim][iaj] += f2b;
				}
			}
		}

		if(mds>0){
			for (int idim=0;idim<xdim;idim++) {
				int idimCube = idim*nCubeNpg;
				for(int idx=0;idx<nCubeNpg;idx++) {
					double f = hFval[idx];
					int iaj = hIAval[idimCube+idx];
					d[idim][iaj] += f*f;
				}
			}
		}

		timerVegasFill.stop();

		tsi *= dv2g;
		double ti2 = ti*ti;
		double wgt = ti2/tsi;
		si += ti*wgt;
		si2 += ti2;
		swgt += wgt;
		schi += ti2*wgt;
		avgi[0] = si/swgt;
		sd[0] = swgt*it/si2;
		chi2a[0] = 0.;
		if(it>1) chi2a[0] = sd[0]*(schi/swgt-avgi[0]*avgi[0])/((double)it-1.);
		sd[0] = sqrt(1./sd[0]);

		if(nprn!=0) {
			tsi = sqrt(tsi);
			coutMaster << std::endl;
			coutMaster << " << integration by vegas >>" << std::endl;
			coutMaster << "     iteration no. " << std::setw(4) << it
						<< std::setw(10) << std::setprecision(6)
						<< "   integral=  " << ti << std::endl;
			coutMaster << "                          std dev  = " << tsi << std::endl;
			coutMaster << "     accumulated results: integral = " << avgi[0] << std::endl;
			coutMaster << "                          std dev  = " << sd[0] << std::endl;
			if(it > 1){
				coutMaster << "                          chi**2 per it'n = "
							<< std::setw(10) << std::setprecision(4) << chi2a << std::endl;
			}
			if(nprn<0){
				for (int j=0;j<xdim;j++) {
					coutMaster << "   == data for axis "
								<< std::setw(2) << j << " --" << std::endl;
					coutMaster << "    x    delt i   convce";
					coutMaster << "    x    delt i   convce";
					coutMaster << "    x    delt i   convce"<<std::endl;
				}
			}
		}

		/* refine grid */
		timerVegasRefine.start();

		double r[nd_max];
		double dt[xdim_max];
		for(int j=0;j<xdim;j++) {
			double xo = d[j][0];
			double xn = d[j][1];
			d[j][0] = 0.5*(xo+xn);
			dt[j] = d[j][0];
			for (int i=1;i<nd-1;i++) {
				d[j][i] = xo+xn;
				xo = xn;
				xn = d[j][i+1];
				d[j][i] = (d[j][i]+xn)/3.;
				dt[j] += d[j][i];
			}
			d[j][nd-1] = 0.5*(xn+xo);
			dt[j] += d[j][nd-1];
		}
      
		for(int j=0;j<xdim;j++) {
			double rc = 0.;
			for(int i=0;i<nd;i++) {
				r[i] = 0.;
				if(d[j][i]>0.) {
					double xo = dt[j]/d[j][i];
					if(!isinf(xo)){
						r[i] = pow(((xo-1.)/xo/log(xo)),alph);
					}
				}
				rc += r[i];
			}
			rc /= xnd;
			int k = -1;
			double xn = 0.;
			double dr = xn;
			int i = k;
			k++;
			dr += r[k];
			double xo = xn;
			xn = xi[j][k];

			do{
				while (dr<=rc) {
					k++;
					dr += r[k];
					xo = xn;
					xn = xi[j][k];
				}
				i++;
				dr -= rc;
				xin[i] = xn-(xn-xo)*dr/r[k];
			} while (i<nd-2);

			for (int i=0;i<nd-1;i++) {
				xi[j][i] = (double)xin[i];
			}
			
			xi[j][nd-1] = 1.;
		}

		gpuErrchk(cudaMemcpyToSymbol(g_xi, xi, sizeof(xi)));
		cudaThreadSynchronize();

		timerVegasRefine.stop();
   } while (it<itmx && acc*fabs(avgi[0])<sd[0]);

	gpuErrchk(cudaFreeHost(hFval));
	gpuErrchk(cudaFree(gFval));

	gpuErrchk(cudaFreeHost(hIAval));
	gpuErrchk(cudaFree(gIAval));

	LOG_FUNC_END
}

__global__
void gVegasCallFunc(double* gFval, int* gIAval, int xdim, int number_of_integrals, int number_of_moments, double *g_lambda, int *g_matrix_D_arr, double *g_fs){
	
	/* --------------------
	 * Check the thread ID
	 * -------------------- */
	 
	const unsigned int tIdx  = threadIdx.x;
	const unsigned int bDimx = blockDim.x;

	const unsigned int bIdx  = blockIdx.x;
	const unsigned int gDimx = gridDim.x;
	const unsigned int bIdy  = blockIdx.y;

	unsigned int bid  = bIdy*gDimx+bIdx;
	const unsigned int tid = bid*bDimx+tIdx;

	int ig = tid/g_npg;

	unsigned nCubeNpg = g_nCubes*g_npg;

	if(tid<nCubeNpg) {
		unsigned ia[xdim_max];
      
		unsigned int tidRndm = tid;
      
		int kg[xdim_max];
      
		unsigned igg = ig;
		for(int j=0;j<xdim;j++) {
			kg[j] = igg%g_ng+1;
			igg /= g_ng;
		}
      
		double randm[xdim_max];
		fxorshift128(tidRndm, xdim, randm);
      
		double x[xdim_max];
      
		double wgt = g_xjac;
		for(int j=0;j<xdim;j++){
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
      
		/* compute function value for this x */
		//func_entropy(double *g_values_out, double *xx, int xdim, int number_of_integrals, int number_of_moments, double *g_lambda, int *g_matrix_D_arr)
		func_entropy(g_fs, x, xdim, number_of_integrals, number_of_moments, g_lambda, g_matrix_D_arr);
		
		/* store computed function value */
		for(int i=0; i < number_of_integrals; i++){
			g_fs[i] = wgt*g_fs[i];
			
			//TODO: here
			gFval[i*nCubeNpg + tid] = g_fs[i];
		}

		for(int idim=0;idim<xdim;idim++) {
			gIAval[idim*nCubeNpg+tid] = ia[idim];
		}
	}
}

__device__ __host__ __forceinline__
void fxorshift128(unsigned int seed, int n, double* a){
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

__device__
void func_entropy(double *g_values_out, double *xx, int xdim, int number_of_integrals, int number_of_moments, double *g_lambda, int *g_matrix_D_arr){
    double V = 0.0;
    double p = 0.0;

    for (int i = 0; i < number_of_moments; i++){
        p = 1.0;
        for (int j = 0; j < xdim; j++){
            p = p*pow(xx[j], g_matrix_D_arr[(i+1)*xdim+j]);
		}

        V = V - p*g_lambda[i];
    }

	/* ff2[0] - type = 0 */
    g_values_out[0] = exp(V);

    /* ff2[1-n] - for gradient, type =2 */
	for(int order = 0; order < number_of_moments; order++){
		p = 1.0;
		for (int j = 0; j < xdim; j++){
			p = p*pow(xx[j], g_matrix_D_arr[(order+1)*xdim+j]);
		}
        g_values_out[1+order] = p*g_values_out[0];
    }

	/* ff2[n+1 - n+1+n*(n+1)/2] - for Hessian, type = 3 */
	int counter = 1+number_of_moments;
	for(int order = 0; order < number_of_moments; order++){
		for(int order2 = order; order2 < number_of_moments; order2++){
			p = 1.0;
			for(int j=0; j<xdim;j++){
				p = p*pow(xx[j], g_matrix_D_arr[(order+1)*xdim+j] + g_matrix_D_arr[(order2+1)*xdim+j]);
			}
			g_values_out[counter] = p*g_values_out[0];
			counter++;
		}
	}
}



}
}

