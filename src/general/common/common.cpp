#include "common/common.h"


namespace pascinference {
namespace common {

/* for loading PETSc options */
char **argv_petsc;
int argc_petsc;

bool Initialize(int argc, char *argv[]){
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
  	#ifdef USE_PETSC
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
	#endif
	
	/* cuda warm up */
	#ifdef USE_CUDA
		kernel_warmup<<<1,1>>>();
		gpuErrchk( cudaDeviceSynchronize() );
	#endif
	
	return true;
}

void Finalize(){
  	/* finalize Petsc */
  	#ifdef USE_PETSC
		/* clean memory of arguments */
		for(size_t i = 0; i < argc_petsc; ++i)
			delete[] argv_petsc[i];

		#ifdef USE_PERMON
			FllopFinalize();
		#else
			PetscFinalize();
		#endif

		petscvector::PETSC_INITIALIZED = false;
	#endif

}

void allbarrier() {
	#ifdef USE_GPU
		gpuErrchk( cudaDeviceSynchronize() );
	#endif

	#ifdef USE_PETSC
		TRYCXX(PetscBarrier(NULL));
	#endif
}

void myround(double in, double *out){
	union myUnion {
		double dValue;
		uint64_t iValue;
	} myValue;
	myValue.dValue=in;

//	myValue.iValue = myValue.iValue*0.001;
//	myValue.iValue = myValue.iValue*1000;

	*out = myValue.dValue;
}

std::string printbool(bool input){
	std::string out;
	if(input){
		out = "true";
	} else {
		out = "false";
	}
	return out;
}

template<typename Out>
void split(const std::string &s, char delim, Out result) {
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        *(result++) = item;
    }
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, std::back_inserter(elems));
    return elems;
}



}
} /* end of namespace */
