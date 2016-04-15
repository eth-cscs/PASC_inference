/** @file benchmark1.cu
 *  @brief test the varx problem solver
 *
 *  Generate random kmeans 2D problem and solve it using VarX solver.
 *
 *  @author Lukas Pospisil
 */

#include "pascinference.h"
#include "solver/tssolver.h"
#include "data/tsdata.h"
#include "model/varxh1fem.h"

#include "kmeans2D.h"

#ifndef USE_PETSCVECTOR
 #error 'This example is for PETSCVECTOR'
#endif
 
using namespace pascinference;

/* set what is what ( which type of vector to use where) */
typedef petscvector::PetscVector Global;
typedef minlin::threx::HostVector<double> Host;
typedef minlin::threx::DeviceVector<double> Device;


extern int pascinference::DEBUG_MODE;

int main( int argc, char *argv[] )
{

	Initialize(argc, argv); // TODO: load parameters from console input
	
	/* say hello */
	coutMaster << "- start program" << std::endl;

	/* parameters of the model */
	int dimx = 2; /* data dimension */
	int T = 50; /* length of time-series */
	int K = 3; /* number of clusters */
	int xmem = 0; /* coeficient of Var model - memory of x */
	int umem = 0; /* coeficient of Var model - memory of u */

	double mu[6] = {0.25, 0, 
					0.0, -0.5, 
					0.0, 0.5};

	double covariance[12] = {	0.001, 0.0, 0.0, 0.1,  
								0.005, 0.0, 0.0, 0.05,
								0.005, 0.0, 0.0, 0.05};
	
	double epssqr = 10; /* penalty */
	
/* ----------- SOLUTION IN PETSC -----------*/
	coutMaster << "-----------------------------------------------------------------" << std::endl;
	coutMaster << "------------------ SOLUTION IN PETSC ----------------------------" << std::endl;
	coutMaster << "-----------------------------------------------------------------" << std::endl;

	/* prepare model */
	coutMaster << "--- PREPARING MODEL ---" << std::endl;
	VarxH1FEMModel<Global> mymodel(T, dimx, K, xmem, umem, epssqr);
	mymodel.print(coutMaster);

	/* prepare time-series data */
	coutMaster << "--- PREPARING DATA ---" << std::endl;
	TSData<Global> mydata(mymodel);

	/* generate some values to data */
	coutMaster << "--- GENERATING DATA ---" << std::endl;
	example::KMeans2D<Global>::generate(T, K, mu, covariance, mydata.get_datavector());

	/* prepare time-series solver */
	coutMaster << "--- PREPARING SOLVER ---" << std::endl;
	TSSolver<Global> mysolver(mydata);

	/* store initial gamma - we will use it for other architectures */
	GeneralVector<Global> gamma0_global(*(mydata.get_gammavector()));
	gamma0_global = *(mydata.get_gammavector()); 

	/* solve the problem */
	coutMaster << "--- SOLVING THE PROBLEM ---" << std::endl;
	mysolver.setting.maxit = 1;
	mysolver.setting.debug_mode = 10;
	
	mysolver.solve();

	/* save results into VTK file */
	coutMaster << "--- SAVING VTK ---" << std::endl;
	example::KMeans2D<Global>::saveVTK("output_global.vtk",T,K,mydata.get_datavector(),mydata.get_gammavector());
	
	mysolver.printtimer(coutMaster);


/* ----------- SOLUTION IN MINLIN -----------*/
	coutMaster << "-----------------------------------------------------------------" << std::endl;
	coutMaster << "------------------ SOLUTION IN MINLIN HOST ----------------------" << std::endl;
	coutMaster << "-----------------------------------------------------------------" << std::endl;

	/* prepare model */
	coutMaster << "--- PREPARING MODEL ---" << std::endl;
	VarxH1FEMModel<Host> mymodel2(T, dimx, K, xmem, umem, epssqr);
	mymodel2.print(coutMaster);

	/* prepare time-series data */
	coutMaster << "--- PREPARING DATA ---" << std::endl;
	TSData<Host> mydata2(mymodel2);

	/* generate some values to data */
	coutMaster << "--- COPY DATA FROM PREVIOUS PROBLEM ---" << std::endl;
	GeneralVector<Global> *datavector_global = mydata.get_datavector();
	GeneralVector<Host> *datavector_host = mydata2.get_datavector();
	for(int i=0;i<datavector_global->size();i++){
		(*datavector_host)(i) = datavector_global->get(i);
	}

	/* prepare time-series solver */
	coutMaster << "--- PREPARING SOLVER ---" << std::endl;
	TSSolver<Host> mysolver2(mydata2);

	/* copy initial approximation from petsc solver */
	coutMaster << "--- SET INITIAL GAMMA FROM PREVIOUS PROBLEM ---" << std::endl;
	GeneralVector<Host> *gamma0_host = mydata2.get_gammavector();
	for(int i=0;i<gamma0_host->size();i++){
		(*gamma0_host)(i) = gamma0_global.get(i);
	}

	/* solve the problem */
	coutMaster << "--- SOLVING THE PROBLEM ---" << std::endl;
	mysolver2.setting.maxit = mysolver.setting.maxit;
	mysolver2.setting.debug_mode = mysolver.setting.debug_mode;
	
	mysolver2.solve();

	/* save results into VTK file */
	coutMaster << "--- SAVING VTK ---" << std::endl;
	example::KMeans2D<Host>::saveVTK("output_host.vtk",T,K,mydata2.get_datavector(),mydata2.get_gammavector());
	
	mysolver2.printtimer(coutMaster);


/* ----------- SOLUTION IN MINLIN -----------*/
	coutMaster << "-----------------------------------------------------------------" << std::endl;
	coutMaster << "------------------ SOLUTION IN MINLIN DEVICE --------------------" << std::endl;
	coutMaster << "-----------------------------------------------------------------" << std::endl;

	/* prepare model */
	coutMaster << "--- PREPARING MODEL ---" << std::endl;
	VarxH1FEMModel<Device> mymodel3(T, dimx, K, xmem, umem, epssqr);
	mymodel3.print(coutMaster);

	/* prepare time-series data */
	coutMaster << "--- PREPARING DATA ---" << std::endl;
	TSData<Device> mydata3(mymodel3);

	/* generate some values to data */
	coutMaster << "--- COPY DATA FROM PREVIOUS PROBLEM ---" << std::endl;
	GeneralVector<Device> *datavector_device = mydata3.get_datavector();
	for(int i=0;i<datavector_global->size();i++){
		(*datavector_device)(i) = datavector_global->get(i);
	}

	/* prepare time-series solver */
	coutMaster << "--- PREPARING SOLVER ---" << std::endl;
	TSSolver<Device> mysolver3(mydata3);

	/* copy initial approximation from petsc solver */
	coutMaster << "--- SET INITIAL GAMMA FROM PREVIOUS PROBLEM ---" << std::endl;
	GeneralVector<Device> *gamma0_device = mydata3.get_gammavector();
	for(int i=0;i<gamma0_device->size();i++){
		(*gamma0_device)(i) = gamma0_global.get(i);
	}

	/* solve the problem */
	coutMaster << "--- SOLVING THE PROBLEM ---" << std::endl;
	mysolver3.setting.maxit = mysolver.setting.maxit;
	mysolver3.setting.debug_mode = mysolver.setting.debug_mode;
	
	mysolver3.solve();

	/* save results into VTK file */
	coutMaster << "--- SAVING VTK ---" << std::endl;
	example::KMeans2D<Device>::saveVTK("output_device.vtk",T,K,mydata3.get_datavector(),mydata3.get_gammavector());
	
	mysolver3.printtimer(coutMaster);

	/* say bye */	
	coutMaster << "- end program" << std::endl;
	
	Finalize();

	return 0;
}

