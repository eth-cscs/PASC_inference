#include "pascinference.h"
#include "solver/tssolver.h"
#include "data/tsdata.h"
#include "result/tsresult.h"
#include "model/kmeansh1.h"

#include "matrix/blockdiaglaplace_explicit.h"
#include "feasibleset/simplexfeasibleset.h"


using namespace pascinference;

/* set what is what ( which type of vector to use where) */
typedef petscvector::PetscVector Global;
typedef minlin::threx::HostVector<double> Host;

extern bool petscvector::PETSC_INITIALIZED;
extern int pascinference::DEBUG_MODE;


int main( int argc, char *argv[] )
{
		
	Initialize(argc, argv); // TODO: load parameters from console input
	petscvector::PETSC_INITIALIZED = true;
	
	/* say hello */	
	Message("- start program");

	/* dimension of the problem */
	int dim = 2; /* data dimension */
	int T = 5; /* length of time-series (size of the block) */
	int K = 3; /* number of clusters (block) */

/* ----------- SOLUTION IN PETSC -----------*/
	/* prepare model */
	KmeansH1Model<Global> mymodel(T, dim, K);
	std::cout << mymodel << std::endl;

	/* prepare time-series data */
	TSData<Global> mydata(mymodel);
	std::cout << mydata << std::endl;

	/* prepare time-series results */
	TSResult<Global> myresult(mymodel);
	std::cout << myresult << std::endl;

	/* prepare time-series solver */
	TSSolver<Global> mysolver(mydata,myresult);
	std::cout << mysolver << std::endl;

	/* say bye */	
	Message("- end program");
	
	petscvector::PETSC_INITIALIZED = false;
	Finalize();

	return 0;
}

