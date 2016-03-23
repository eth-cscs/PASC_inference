/** @file test_varx.cu
 *  @brief test the varx problem solver
 *
 *  Generate random problem and solve it.
 *
 *  @author Lukas Pospisil
 */

#include "pascinference.h"
#include "solver/tssolver.h"
#include "data/tsdata.h"
#include "model/varxh1fem.h"

#include "varx2D.h"

#ifndef USE_PETSCVECTOR
 #error 'This example is for PETSCVECTOR'
#endif
 
using namespace pascinference;

/* set what is what ( which type of vector to use where) */
typedef petscvector::PetscVector Global;
//typedef minlin::threx::HostVector<double> Global;

extern int pascinference::DEBUG_MODE;


int main( int argc, char *argv[] )
{
	
	Initialize(argc, argv); // TODO: load parameters from console input
	
	/* say hello */
	coutMaster << "- start program" << std::endl;

	/* dimension of the problem */
	int dimx = 2; /* data dimension */
	int T = 1000; /* length of time-series */
	int K = 2; /* number of clusters */
	int xmem = 3; /* coeficient of Var model */
	int umem = 3; /* coeficient of Var model */


	/* parameters of the model */
	double muK1[dimx] = {2.0, 6.0};
	double muK2[dimx] = {5.0, 4.0};
	double *mu_orig[K] = {muK1,muK2};

	double covarianceK1Q1[dimx*dimx] = {0.03, -0.07, 0.02, 0.07};
	double covarianceK1Q2[dimx*dimx] = {0.07, -0.03, 0.01, 0.06};
	double covarianceK1Q3[dimx*dimx] = {-0.4, -0.7, -0.2, -0.8};
	double *covarianceK1[xmem] = {covarianceK1Q1,covarianceK1Q2,covarianceK1Q3};

	double covarianceK2Q1[dimx*dimx] = {0.03, 0.07, -0.06, 0.04};
	double covarianceK2Q2[dimx*dimx] = {0.03, 0.08, 0.01, -0.09};
	double covarianceK2Q3[dimx*dimx] = {0.1, -0.2, -0.3, -0.4};
	double *covarianceK2[xmem] = {covarianceK2Q1,covarianceK2Q2,covarianceK2Q3};

	double **covariance_orig[K] = {covarianceK1,covarianceK2};

	double x0[dimx] = {0.3, -0.5};
	double x1[dimx] = {0.7, 0.1};
	double x2[dimx] = {0.1, -0.9};
	double *xbegin[xmem] = {x0,x1,x2};

	double ut[T+xmem];
	for(int i=0;i<T+xmem;i++){
		ut[i] = (-2.0/(double)(T+xmem)) * (i+1);
	}

	double gamma_orig[T*K];
	for(int t=0;t<T;t++){
		if((t>143 && t <=159) || (t>246 && t <= 303) || (t>346 && t <=382) || (t>433 && t <=475) || (t>577 && t <= 672) || (t>911 && t <=971) ){
			gamma_orig[t] = 0;
			gamma_orig[t+T] = 1;
		} else {
			gamma_orig[t] = 1;
			gamma_orig[t+T] = 0;
		} 
	}

	double epssqr = 10; /* penalty */
	double noise_sigma = 0.0005;
	
/* ----------- SOLUTION IN PETSC -----------*/
	/* prepare model */
	coutMaster << "--- PREPARING MODEL ---" << std::endl;
	VarxH1FEMModel<Global> mymodel(T, dimx, K, xmem, umem, epssqr);
	mymodel.print(coutMaster);

	/* prepare time-series data */
	coutMaster << "--- PREPARING DATA ---" << std::endl;
	TSData<Global> mydata(mymodel);
	mydata.print(coutMaster);

	/* generate some values to data */
	coutMaster << "--- GENERATING DATA ---" << std::endl;
//	example::Varx2D<Global>::generate(T,K,mu,dev,mydata.get_datavector());
//	mydata.printcontent(coutMaster);

	/* prepare time-series solver */
//	coutMaster << "--- PREPARING SOLVER ---" << std::endl;
//	TSSolver<Global> mysolver(mydata);

//	mysolver.setting.maxit = 50;
//	mysolver.setting.debug_mode = 0;

	/* solve the problem */
	/* gamma_solver = SOLVER_SPGQP, theta_solver = SOLVER_CG */
//	coutMaster << "--- SOLVING THE PROBLEM ---" << std::endl;
//	mysolver.solve(SOLVER_SPGQP, SOLVER_CG);

	/* save results into VTK file */
//	coutMaster << "--- SAVING VTK ---" << std::endl;
//	example::KMeans2D<Global>::saveVTK("output.vtk",T,K,mydata.get_datavector(),mydata.get_gammavector());

//	mysolver.printtimer(coutMaster);

	/* say bye */	
//	coutMaster << "- end program" << std::endl;
	
	Finalize();

	return 0;
}

