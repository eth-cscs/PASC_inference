/** @file test_kmeans_generate.cpp
 *  @brief generate Kmeans problem and save it to separate files
 *
 *  @author Lukas Pospisil
 */

#include <iostream>
#include <list>
#include <algorithm>

#include "pascinference.h"
#include "data/kmeansdata.h"
#include "model/kmeansh1fem.h"

#ifndef USE_PETSCVECTOR
 #error 'This example is for PETSCVECTOR'
#endif
 
typedef petscvector::PetscVector PetscVector;

using namespace pascinference;

extern int pascinference::DEBUG_MODE;

int solution_get_cluster_id(int t, int T){
    int id_cluster;
    if((t >= ceil(0*T) && t < ceil(0.1*T)) || (t >= ceil(0.35*T) && t < ceil(0.4*T)) || (t >= ceil(0.5*T) && t < ceil(0.52*T)) || (t >= ceil(0.8*T) && t < ceil(0.9*T)) || (t >= ceil(0.2*T) && t < ceil(0.25*T))){
        id_cluster = 0;
    }

    if((t >= ceil(0.1*T) && t < ceil(0.2*T)) || (t >= ceil(0.4*T) && t < ceil(0.5*T)) || (t >= ceil(0.8*T) && t < ceil(0.8*T)) || (t >= ceil(0.95*T) && t <= ceil(1.0*T)) || (t >= ceil(0.6*T) && t < ceil(0.8*T))){
	id_cluster = 1;
    }
    
    if((t >= ceil(0.25*T) && t < ceil(0.35*T)) || (t >= ceil(0.52*T) && t < ceil(0.6*T)) || (t >= ceil(0.9*T) && t < ceil(0.95*T))){
	id_cluster = 2;
    }
    return id_cluster;
}

int main( int argc, char *argv[] )
{
	/* add local program options */
	boost::program_options::options_description opt_problem("PROBLEM EXAMPLE", consoleArg.get_console_nmb_cols());
	opt_problem.add_options()
		("test_T", boost::program_options::value<std::vector<int> >()->multitoken(), "dimensions of the problem")
		("test_K", boost::program_options::value<std::vector<int> >()->multitoken(), "numbers of clusters");
	consoleArg.get_description()->add(opt_problem);

	/* call initialize */
	if(!Initialize(argc, argv)){
		return 0;
	} 

	/* which times to generate */
	int T;
	std::vector<int> T_list;
	if(!consoleArg.set_option_value("test_T", &T_list)){
		std::cout << "test_T has to be set! Call application with parameter -h to see all parameters" << std::endl;
		return 0;
	}

	/* gamma0 to which K to generate */
	int K;
	std::vector<int> K_list;
	if(!consoleArg.set_option_value("test_K", &K_list)){
		std::cout << "test_K has to be set! Call application with parameter -h to see all parameters" << std::endl;
		return 0;
	}

	int i;
	int T_size = T_list.size();
	int K_size = K_list.size();
	std::ostringstream oss;
	int T_max = 0;

	/* print info about what we will compute */
	coutMaster << "- PROBLEM INFO --------------------------------------------------" << std::endl;
	coutMaster << " T      = ";
	for(i=0;i<T_size;i++){
		coutMaster << T_list[i];
		if(i < T_size-1){ 
				coutMaster << ", ";
		}
		/* store maximum value */
		if(T_list[i] > T_max){
			T_max = T_list[i];
		}
	}
	coutMaster << " (length of time-series)" << std::endl;

	coutMaster << " K      = ";
	for(i=0;i<K_size;i++){
		coutMaster << K_list[i];
		if(i < K_size-1){ 
				coutMaster << ", ";
		}
	}
	coutMaster << " (number of clusters)" << std::endl;

	coutMaster << "------------------------------------------------------------------" << std::endl;
	
	/* start logging */
	std::ostringstream oss_name_of_file_log;
	oss_name_of_file_log << "results/test_kmeans_generate_log_p" << GlobalManager.get_rank() << ".txt";
	logging.begin(oss_name_of_file_log.str());
		
	/* say hello */
	coutMaster << "- start program" << std::endl;

	/* parameters of the model */
	int xdim = 3; /* data dimension */

	/* solution - for generating the problem */
	int solution_K = 3;
	
	double solution_theta[9] = {
		 0.0, 0.0, 0.0,		/* K=1,n=1,2,3:  mu */
		 1.0, 0.0, 0.0,		/* K=2,n=1,2,3 */
		 0.0, 1.0, 0.0		/* K=3,n=1,2,3 */
	};

	double solution_xstart[3] = {
		0.0, 0.0, 0.0	/* n=1,2,3 */
	};

	double noise_covariance[9] = {
		0.05, 0.05, 0.01,  
		0.05, 0.05, 0.01,  
		0.05, 0.05, 0.01  
	};


	std::ostringstream oss_name_of_file;

	/* prepare random generator */
	PetscRandom rnd;
	TRY( PetscRandomCreate(PETSC_COMM_WORLD,&rnd) );
	TRY( PetscRandomSetType(rnd,PETSCRAND) );
	TRY( PetscRandomSetFromOptions(rnd) );
	TRY( PetscRandomSetSeed(rnd,13) );

	/* ---- GENERATE THE LARGEST PROBLEM ---- */
	coutMaster << "- generating largest problem: T = " << T_max << std::endl;
	KmeansData<PetscVector> mydata(T_max,xdim);
	KmeansH1FEMModel<PetscVector> mymodel(mydata, xdim, solution_K, 0);
	mydata.set_model(mymodel);

	mydata.generate(solution_K, solution_theta, &solution_get_cluster_id, false);
	mydata.add_noise(noise_covariance);

	/* we will generate gamma0 vector for all K */
	SimplexFeasibleSet_Local *feasibleset;
	Vec *gamma0s_Vec; /* array of vectors */
	GeneralVector<PetscVector> *gamma0; /* temp general vector for projection */
	gamma0s_Vec = (Vec *)(malloc(K_size*sizeof(Vec)));

	PetscViewer viewer_out;
	Vec data_Vec = mydata.get_datavector()->get_vector();

	int ki;
	for(ki = 0; ki < K_size; ki++){
		K = K_list[ki];

		coutMaster << "- computing gamma: K = " << K << std::endl;

		/* create feasible set - we will project */
		feasibleset = new SimplexFeasibleSet_Local(T_max,K); 

		///* create general vector */
		TRY( VecCreateSeq(PETSC_COMM_SELF, K*T_max, &(gamma0s_Vec[ki])) );
		gamma0 = new GeneralVector<PetscVector>(gamma0s_Vec[ki]);
		
		/* generate random gamma */
		TRY( VecSetRandom(gamma0s_Vec[ki], rnd) );

		///* project initial approximation to feasible set */
		feasibleset->project(*gamma0);
	
		///* destroy feasible set */
		free(feasibleset);
	}

	/* ---- GET SUB PROBLEMS ---- */
	Vec subdata_Vec;
	IS *subdata_ISs;
	subdata_ISs = (IS *)(malloc(xdim*sizeof(IS)));
	IS subdata_IS;

	Vec gamma0sum_Vec;
	IS *gamma0sub_ISs;
	IS gamma0sub_IS;

	int ti;
	for(ti=0;ti<T_size;ti++){
		T = T_list[ti];
		coutMaster << "- getting subproblem: T = " << T << std::endl;
	
		/* for every dimension of data create stride */
		for(i=0;i<xdim;i++){
			TRY( ISCreateStride(PETSC_COMM_WORLD, T, i, T_max/(double)(T), &(subdata_ISs[i])) );
		}
		TRY( ISConcatenate(PETSC_COMM_WORLD, xdim, subdata_ISs, &subdata_IS) );

		/* get subvector */
		TRY( VecGetSubVector(data_Vec, subdata_IS, &subdata_Vec) );
		
		/* save data */
		oss_name_of_file << "results/data_kmeans_T" << T << ".bin";
		coutMaster << "- saving: " << oss_name_of_file.str() << std::endl;
		TRY( PetscViewerBinaryOpen(PETSC_COMM_WORLD,oss_name_of_file.str().c_str(),FILE_MODE_WRITE,&viewer_out) );
		TRY( VecView( subdata_Vec, viewer_out) );
		TRY( PetscViewerDestroy(&viewer_out) );
		oss_name_of_file.str("");
		
		/* restore subvector */
		TRY( VecRestoreSubVector(data_Vec, subdata_IS, &subdata_Vec) );

		/* destroy indexsets */
		for(i=0;i<xdim;i++){
			TRY( ISDestroy(&(subdata_ISs[i])) );
		}
		TRY( ISDestroy(&subdata_IS) );

		/* subvectors from gamma0 */
		for(ki = 0; ki < K_size; ki++){
			K = K_list[ki];

			coutMaster << "- getting subgamma: K = " << K << std::endl;

			/* get subvectors from gamma0 */
			gamma0sub_ISs = (IS*)malloc(K*sizeof(IS));
			for(i=0;i<K;i++){
				TRY( ISCreateStride(PETSC_COMM_WORLD, T, i, T_max/(double)(T), &(gamma0sub_ISs[i])) );
			}
			TRY( ISConcatenate(PETSC_COMM_WORLD, K, gamma0sub_ISs, &gamma0sub_IS) );
			
			/* get subvector */
			TRY( VecGetSubVector(gamma0s_Vec[ki], gamma0sub_IS, &gamma0sum_Vec) );
		
			/* save data */
			oss_name_of_file << "results/gamma0_kmeans_T" << T << "K" << K << ".bin";
			coutMaster << "- saving: " << oss_name_of_file.str() << std::endl;
			TRY( PetscViewerBinaryOpen(PETSC_COMM_WORLD,oss_name_of_file.str().c_str(),FILE_MODE_WRITE,&viewer_out) );
			TRY( VecView( gamma0sum_Vec, viewer_out) );
			TRY( PetscViewerDestroy(&viewer_out) );
			oss_name_of_file.str("");
		
			/* restore subvector */
			TRY( VecRestoreSubVector(gamma0s_Vec[ki], gamma0sub_IS, &gamma0sum_Vec) );

			for(i=0;i<K;i++){
				TRY( ISDestroy(&(gamma0sub_ISs[i])) );
			}
			free(gamma0sub_ISs);
			TRY( ISDestroy(&gamma0sub_IS) );

		}

	}



	/* destroy the random generator */
	TRY( PetscRandomDestroy(&rnd) );

	/* say bye */	
	coutMaster << "- end program" << std::endl;

	logging.end();
	Finalize();

	return 0;
}

