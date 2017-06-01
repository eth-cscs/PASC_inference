#include "external/petscvector/common/initialize.h"

namespace pascinference {
namespace common {

bool PETSC_INITIALIZED;

char **argv_petsc;
int argc_petsc;

template<>
bool Initialize<PetscVector>(int argc, char *argv[]){
	/* console arguments */
	if(!consoleArg.init(argc,argv)){
		return false;
	}

	/* initialize random seed: */
	if(RANDOM_BY_TIME){
		srand(time(NULL));
	} else {
		srand(0);
	}

  	/* init Petsc */
	/* Petsc options are provided using "-petsc_options=string" */
	std::string petsc_options_string;
	consoleArg.set_option_value("petsc_options", &petsc_options_string, "");

	/* split string and create argc_petsc and argv_petsc */
	petsc_options_string.insert(0," "); /* add blank space in the beginning (?) */
	std::vector<std::string> petsc_options_vector = split(petsc_options_string, ' ');
	argc_petsc = petsc_options_vector.size();
	
	argv_petsc = new char*[argc_petsc]; /* the first parameter is ignored (?) */
	for(size_t i = 0; i < argc_petsc; ++i){
		argv_petsc[i] = new char[petsc_options_vector[i].size() + 1];
		std::strcpy(argv_petsc[i], petsc_options_vector[i].c_str());
	}

	#ifdef USE_PERMON
		FllopInitialize(&argc_petsc,&argv_petsc,PETSC_NULL);
//			FllopInitialize(PETSC_NULL,PETSC_NULL,PETSC_NULL);
	#else
		PetscInitialize(&argc_petsc,&argv_petsc,PETSC_NULL,PETSC_NULL);
//			PetscInitialize(PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);
	#endif

	petscvector::PETSC_INITIALIZED = true;
	
	/* cuda warm up */
	#ifdef USE_CUDA
		cuda_warmup();
	#endif
	
	return true;
}

template<>
void Finalize<PetscVector>(){
  	/* finalize Petsc */
	/* clean memory of arguments */
	for(size_t i = 0; i < argc_petsc; ++i)
		delete[] argv_petsc[i];

	#ifdef USE_PERMON
		FllopFinalize();
	#else
		PetscFinalize();
	#endif

	petscvector::PETSC_INITIALIZED = false;
}

template<>
void allbarrier<PetscVector>() {
	#ifdef USE_CUDA
		cuda_barrier();
	#endif

	#ifdef USE_PETSC
		TRYCXX(PetscBarrier(NULL));
	#endif
}

#ifdef USE_CUDA
void cuda_copytoGPU(Vec &x) {
	TRYCXX( VecCUDACopyToGPU(x) );
}
#endif

}
} /* end of namespace */
