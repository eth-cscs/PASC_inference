/** @file test_kmeans_generate.cpp
 *  @brief generate Kmeans problem and save it to separate files
 *
 *  @author Lukas Pospisil
 */

#include "pascinference.h"
#include "data/tsdata_global.h"
#include "model/varxh1fem_global.h"

#include "kmeans3D.h"

#ifndef USE_PETSCVECTOR
 #error 'This example is for PETSCVECTOR'
#endif
 
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
		("test_Tbegin", boost::program_options::value<int>(), "dimension of the problem")
		("test_Tstep", boost::program_options::value<int>(), "dimension of the problem")
		("test_Tend", boost::program_options::value<int>(), "dimension of the problem")
		("test_Kbegin", boost::program_options::value<int>(), "number of clusters")
		("test_Kstep", boost::program_options::value<int>(), "number of clusters")
		("test_Kend", boost::program_options::value<int>(), "number of clusters");
	consoleArg.get_description()->add(opt_problem);

	/* call initialize */
	if(!Initialize(argc, argv)){
		return 0;
	} 

	/* which times to generate */
	int T_begin,T_step, T_end;
	if(!consoleArg.set_option_value("test_Tbegin", &T_begin)){
		std::cout << "test_Tbegin has to be set! Call application with parameter -h to see all parameters" << std::endl;
		return 0;
	}
	consoleArg.set_option_value("test_Tend", &T_end, T_begin);
	consoleArg.set_option_value("test_Tstep", &T_step, 1);

	/* gamma0 to which K to generate */
	int K_begin,K_step, K_end;
	if(!consoleArg.set_option_value("test_Kbegin", &K_begin)){
		std::cout << "test_Kbegin has to be set! Call application with parameter -h to see all parameters" << std::endl;
		return 0;
	}
	consoleArg.set_option_value("test_Kend", &K_end, K_begin);
	consoleArg.set_option_value("test_Kstep", &K_step, 1);


	coutMaster << "- PROBLEM INFO --------------------------------------------------" << std::endl;
	coutMaster << " T_begin:T_step:T_end      = " << std::setw(7) << T_begin << std::setw(7) << T_step << std::setw(7) << T_end << " (length of time-series)" << std::endl;
	coutMaster << " K_begin:K_step:K_end      = " << std::setw(7) << K_begin << std::setw(7) << K_step << std::setw(7) << K_end << " (number of clusters)" << std::endl;
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
	int solution_xmem = 0;
	
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

	/* this is kmeans problem, but we use varx-model */
	int xmem = 0;

	std::ostringstream oss_name_of_file;

	/* prepare random generator */
	PetscRandom rnd;
	TRY( PetscRandomCreate(PETSC_COMM_WORLD,&rnd) );
	TRY( PetscRandomSetType(rnd,PETSCRAND) );
	TRY( PetscRandomSetFromOptions(rnd) );
	TRY( PetscRandomSetSeed(rnd,13) );

	/* ---- GENERATE THE LARGEST PROBLEM ---- */
	/* prepare data */
	coutMaster << "- generating largest problem: T = " << T_end << std::endl;
	VarxH1FEMModel_Global mymodel(T_end, xdim, 0, xmem, 0);
	TSData_Global mydata(mymodel);
	mymodel.generate_data(solution_K, solution_xmem, solution_theta, solution_xstart, &solution_get_cluster_id, &mydata, false);
	mymodel.generate_data_add_noise(&mydata, noise_covariance, &solution_get_cluster_id);

	Vec data_Vec = mydata.get_datavector()->get_vector();

	/* save data */
	oss_name_of_file << "results/data_kmeans_T" << T_end << ".bin";
	coutMaster << "- saving: " << oss_name_of_file.str() << std::endl;
	PetscViewer viewer_out;
	TRY( PetscViewerBinaryOpen(PETSC_COMM_WORLD,oss_name_of_file.str().c_str(),FILE_MODE_WRITE,&viewer_out) );
	TRY( VecView( data_Vec, viewer_out) );
	TRY( PetscViewerDestroy(&viewer_out) );
	oss_name_of_file.str("");

	/* generate gamma0 vector for all K */
	SimplexFeasibleSet_Local *feasibleset;
	Vec *gamma0s_Vec; /* array of vectors */
	GeneralVector<PetscVector> *gamma0; /* temp general vector for projection */
	int K_num = (K_end - K_begin)/(double)K_step;
	gamma0s_Vec = (Vec *)(malloc(K_num*sizeof(Vec)));

	int k, ki;
	for(ki = 0; ki <= K_num; ki++){
		k = K_begin + ki*K_step;

		coutMaster << "- getting subgamma: K = " << k << std::endl;

		/* create feasible set - we will project */
		feasibleset = new SimplexFeasibleSet_Local(T_end,k); 

		/* create general vector */
		TRY( VecCreateSeq(PETSC_COMM_SELF, k*T_end, &(gamma0s_Vec[ki])) );
		gamma0 = new GeneralVector<PetscVector>(gamma0s_Vec[ki]);
		
		/* generate random gamma */
		TRY( VecSetRandom(gamma0s_Vec[ki], rnd) );

		/* project initial approximation to feasible set */
		feasibleset->project(*gamma0);

		/* store initial approximation */
		oss_name_of_file << "results/gamma0_kmeans_T" << T_end << "K" << k << ".bin";
		coutMaster << "- saving: " << oss_name_of_file.str() << std::endl;
		TRY( PetscViewerBinaryOpen(PETSC_COMM_WORLD,oss_name_of_file.str().c_str(),FILE_MODE_WRITE,&viewer_out) );
		TRY( VecView( gamma0s_Vec[ki], viewer_out) );
		TRY( PetscViewerDestroy(&viewer_out) );
		oss_name_of_file.str("");
		
		/* destroy feasible set */
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

	int t,i;
	for(t=T_begin;t<T_end;t+=T_step){
		coutMaster << "- getting subproblem: T = " << t << std::endl;
	
		/* for every dimension of data create stride */
		for(i=0;i<xdim;i++){
			TRY( ISCreateStride(PETSC_COMM_WORLD, t, i*T_end, T_end/(double)(t), &(subdata_ISs[i])) );
		}
		TRY( ISConcatenate(PETSC_COMM_WORLD, xdim, subdata_ISs, &subdata_IS) );

		/* get subvector */
		TRY( VecGetSubVector(data_Vec, subdata_IS, &subdata_Vec) );
		
		/* save data */
		oss_name_of_file << "results/data_kmeans_T" << t << ".bin";
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
		for(ki = 0; ki <= K_num; ki++){
			k = K_begin + ki*K_step;

			coutMaster << "- getting subgamma: K = " << k << std::endl;

			/* get subvectors from gamma0 */
			gamma0sub_ISs = (IS*)malloc(k*sizeof(IS));
			for(i=0;i<k;i++){
				TRY( ISCreateStride(PETSC_COMM_WORLD, t, i*T_end, T_end/(double)(t), &(gamma0sub_ISs[i])) );
			}
			TRY( ISConcatenate(PETSC_COMM_WORLD, k, gamma0sub_ISs, &gamma0sub_IS) );
			
			/* get subvector */
			TRY( VecGetSubVector(gamma0s_Vec[ki], gamma0sub_IS, &gamma0sum_Vec) );
		
			/* save data */
			oss_name_of_file << "results/gamma0_kmeans_T" << t << "K" << k << ".bin";
			coutMaster << "- saving: " << oss_name_of_file.str() << std::endl;
			TRY( PetscViewerBinaryOpen(PETSC_COMM_WORLD,oss_name_of_file.str().c_str(),FILE_MODE_WRITE,&viewer_out) );
			TRY( VecView( gamma0sum_Vec, viewer_out) );
			TRY( PetscViewerDestroy(&viewer_out) );
			oss_name_of_file.str("");
		
			/* restore subvector */
			TRY( VecRestoreSubVector(gamma0s_Vec[ki], gamma0sub_IS, &gamma0sum_Vec) );

			for(i=0;i<k;i++){
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

