/** @file kmeansdata.h
 *  @brief
 * 
 *  @author Lukas Pospisil
 */

#ifndef PASC_KMEANSDATA_H
#define	PASC_KMEANSDATA_H

#ifndef USE_PETSCVECTOR
 #error 'KMEANSDATA is for PETSCVECTOR'
#endif

typedef petscvector::PetscVector PetscVector;

/* for debugging, if >= 100, then print info about ach called function */
extern int DEBUG_MODE;

#include <iostream>
#include "common/common.h"
#include "data/tsdata.h"
#include "model/tsmodel.h"

namespace pascinference {
namespace data {

template<class VectorBase>
class KmeansData: public TSData<VectorBase> {
	protected:

	public:
		KmeansData() : TSData<VectorBase>() {};
		KmeansData(GeneralVector<VectorBase> *datavector_new, GeneralVector<VectorBase> *gammavector_new, GeneralVector<VectorBase> *thetavector_new, int T) : TSData<VectorBase>(datavector_new, gammavector_new, thetavector_new, T) {};

		KmeansData(int T, int blocksize=1) : TSData<VectorBase>(T,blocksize) {};
		KmeansData(std::string filename , int xdim=1) : TSData<VectorBase>(filename, xdim) {};
		~KmeansData() {};

		virtual std::string get_name() const;

		void saveCSV(std::string filename) const;
		void saveVTK(std::string filename) const;
		
		void generate(int K_solution, double *theta_solution, int (*get_cluster_id)(int, int), bool scale_or_not=false);
		void add_noise(double *diag_covariance);
		
		void load_gammavector(std::string filename) const;
		void load_gammavector(VectorBase &gamma0) const;

};


}
} /* end of namespace */

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {
namespace data {

template<class VectorBase>
std::string KmeansData<VectorBase>::get_name() const {
	return "Kmeans Time-series Data";
}

template<>
void KmeansData<PetscVector>::saveCSV(std::string filename) const {
	LOG_FUNC_STATIC_BEGIN

	Timer timer_saveCSV; 
	timer_saveCSV.restart();
	timer_saveCSV.start();
	
	int nproc = GlobalManager.get_size();
	int my_rank = GlobalManager.get_rank();

	int n,t,k,proc_id;
	
	int K = get_K();
	int T = get_T();
	int Tlocal = get_Tlocal();
	int xdim = get_xdim();
				
	/* to manipulate with file */
	std::ofstream myfile;
						
	/* first processor write header */
	if(GlobalManager.get_rank() == 0){
		myfile.open(filename.c_str());

		/* write header to file */
		for(n=0; n<xdim; n++){
			myfile << "x" << n << "_orig,";
		}
		for(k=0; k<K; k++){
			myfile << "gamma" << k << ",";
		}
		for(n=0; n<xdim; n++){
			myfile << "x" << n << "_model";
			if(n+1 < xdim){
				myfile << ",";
			}
		}
		myfile << std::endl;
		myfile.close();
	}

	/* theta */
	double *theta_arr;
	TRY( VecGetArray(thetavector->get_vector(),&theta_arr) );

	/* gamma */
	double *gamma_arr;
	TRY( VecGetArray(gammavector->get_vector(),&gamma_arr) );

	/* data */
	double *data_arr;
	TRY( VecGetArray(datavector->get_vector(),&data_arr) );
	
	double xmodel_n;
	
	/* go through processors and write local portion of data */
	for(proc_id=0; proc_id < GlobalManager.get_size();proc_id++){
		if(GlobalManager.get_rank() == proc_id){
			/* open file - append to the end local data */
			myfile.open(filename.c_str(), std::fstream::in | std::fstream::out | std::fstream::app);
			/* write local portion of data */
			for(t=0;t<Tlocal;t++){
				/* original time_serie */
				for(n=0;n<xdim;n++){
					myfile << data_arr[t*xdim+n] << ",";
				}
				/* write gamma vectors */
				for(k=0;k<K;k++){
					myfile << gamma_arr[t*K+k] << ",";
				}
				
				/* compute model value and write it */
				for(n=0;n<xdim;n++){
					xmodel_n = 0;
					for(k=0;k<K;k++){
						xmodel_n += gamma_arr[t*K+k]*theta_arr[k*xdim+n];
					}
					myfile << xmodel_n;
					if(n+1 < xdim){
						myfile << ",";
					}
				}
				myfile << std::endl;
			}
			myfile.close();			
		}
		
		TRY(PetscBarrier(NULL));
	}

	TRY( VecRestoreArray(gammavector->get_vector(),&gamma_arr) );
	TRY( VecRestoreArray(thetavector->get_vector(),&theta_arr) );
	TRY( VecRestoreArray(datavector->get_vector(),&data_arr) );
			
	/* writing finished */
	timer_saveCSV.stop();
	coutMaster <<  " - problem saved to CSV in: " << timer_saveCSV.get_value_sum() << std::endl;

	LOG_FUNC_STATIC_END
}

template<>
void KmeansData<PetscVector>::generate(int K_solution, double *theta_solution, int (*get_cluster_id)(int, int), bool scale_or_not){
	LOG_FUNC_BEGIN
			
	/* get local data array */
	double *data_arr;
	TRY( VecGetArray(datavector->get_vector(),&data_arr) );

	int t, n, k;
	for(t=0;t<get_Tlocal();t++){
		k = get_cluster_id(get_Tbegin()+t, get_T());
		for(n=0;n<get_xdim();n++){
			data_arr[t*get_xdim()+n] = theta_solution[k*get_xdim()+n];
		}
	}

	/* restore local data array */
	TRY( VecRestoreArray(datavector->get_vector(),&data_arr) );

	/* scale data */
	if(scale_or_not){
		double max_value;
		TRY( VecMax(datavector->get_vector(), NULL, &max_value) );
		TRY( VecScale(datavector->get_vector(), 1.0/max_value) );
				
		coutAll << "--- scaling data with max value of x: " << max_value << std::endl;
		coutAll.synchronize();
	}

	LOG_FUNC_END
}

template<>
void KmeansData<PetscVector>::add_noise(double *diag_covariance){
	LOG_FUNC_BEGIN
	
	int n,t;
	
	/* prepare mu */
	double mu[get_xdim()];
	double values[get_xdim()];

	for(n = 0; n <get_xdim(); n++){
		mu[n] = 0.0;
	}

	/* get local data array */
	double *data_arr;
	TRY( VecGetArray(datavector->get_vector(),&data_arr) );
	
	/* add noise to generated data */
	for(t = 0; t < get_Tlocal(); t++){
		my_mvnrnd_Dn(get_xdim(), mu, diag_covariance, values);
		for(n = 0; n <get_xdim(); n++){
			data_arr[t*get_xdim()+n] += values[n];
		}
	}
	TRY( VecRestoreArray(datavector->get_vector(),&data_arr) );

	LOG_FUNC_END	
}

template<>
void KmeansData<PetscVector>::saveVTK(std::string filename) const{
	LOG_FUNC_BEGIN
	
	Timer timer_saveVTK; 
	timer_saveVTK.restart();
	timer_saveVTK.start();

	int t,k, id_proc;
	int T = get_T();
	int Tlocal = get_Tlocal();
	int K = get_K();
	int xdim = get_xdim();

	int prank = GlobalManager.get_rank();
	int psize = GlobalManager.get_size();

	/* to manipulate with file */
	std::ofstream myfile;

	/* master writes the header */
	if(prank == 0){
		myfile.open(filename.c_str());

		/* write header to file */
		myfile << "# vtk DataFile Version 3.1" << std::endl;
		myfile << "PASCInference: Kmeans solution" << std::endl;
		myfile << "ASCII" << std::endl;
		myfile << "DATASET UNSTRUCTURED_GRID" << std::endl;

		/* points - coordinates */
		myfile << "POINTS " << T << " FLOAT" << std::endl;
		
		myfile.close();
	}
	TRY( PetscBarrier(NULL) );

	double *data_arr;
	TRY( VecGetArray(datavector->get_vector(), &data_arr) );

	/* each processor writes its own portion of data */
	for(id_proc=0;id_proc<psize;id_proc++){
		if(id_proc == prank){
			myfile.open(filename.c_str(), std::fstream::in | std::fstream::out | std::fstream::app);

			for(t=0;t < Tlocal;t++){
				if(xdim == 1){ 
					/* 1D sample */
					myfile << data_arr[t] << " 0 0" << std::endl; /* x */
				}

				if(xdim == 2){ 
					/* 3D sample */
					myfile << data_arr[2*t] << " "; /* x */
					myfile << data_arr[2*t+1] << " 0" << std::endl; /* y */
				}

				if(xdim == 3){ 
					/* 3D sample */
					myfile << data_arr[3*t] << " "; /* x */
					myfile << data_arr[3*t+1] << " "; /* y */
					myfile << data_arr[3*t+2] << std::endl; /* z */
				}

				if(xdim > 3){
					//TODO ???
				}
			}
					
			myfile.close();
		}
		
		TRY( PetscBarrier(NULL) );
	}
	TRY( VecRestoreArray(datavector->get_vector(), &data_arr) );

	/* master writes the header */
	if(prank == 0){
		myfile.open(filename.c_str(), std::fstream::in | std::fstream::out | std::fstream::app);
	
		myfile << std::endl;
		myfile << "POINT_DATA " <<  T << std::endl;

		myfile.close();
	}
	TRY( PetscBarrier(NULL) );

	/* here I will store id of max gamma */
	int *gamma_max_id;
	gamma_max_id = (int *)malloc(Tlocal*sizeof(int));
	double *gamma_max;
	gamma_max = (double *)malloc(Tlocal*sizeof(double));
	for(t = 0; t<Tlocal;t++){
		gamma_max_id[t] = 0;
		gamma_max[t] = 0.0;
	}

	double *gamma_arr;
	TRY( VecGetArray(gammavector->get_vector(), &gamma_arr) );
	
	/* write gamma_k */
	for(k=0;k<K;k++){
		/* master writes header */
		if(prank == 0){
			myfile.open(filename.c_str(), std::fstream::in | std::fstream::out | std::fstream::app);

			myfile << std::endl;
			myfile << "SCALARS gamma_" << k << " float 1" << std::endl;
			myfile << "LOOKUP_TABLE default" << std::endl;

			myfile.close();
		}
		TRY( PetscBarrier(NULL) );
		
		/* each processor writes its own portion of data */
		for(id_proc=0;id_proc<psize;id_proc++){
			if(id_proc == prank){
				myfile.open(filename.c_str(), std::fstream::in | std::fstream::out | std::fstream::app);

				for(t=0;t < Tlocal;t++){
					myfile << gamma_arr[t*K+k] << std::endl;

					/* update maximum */
					if(gamma_arr[t*K+k] > gamma_max[t]){
						gamma_max[t] = gamma_arr[t*K+k];
						gamma_max_id[t] = k;
					}
				}

				myfile.close();
			}
		
			TRY( PetscBarrier(NULL) );
		}
	}
	
	TRY( VecRestoreArray(gammavector->get_vector(), &gamma_arr) );

	/* store gamma values */
	/* master writes the header */
	if(prank == 0){
		myfile.open(filename.c_str(), std::fstream::in | std::fstream::out | std::fstream::app);
	
		myfile << std::endl;
		myfile << "SCALARS gamma_max_id float 1" << std::endl;
		myfile << "LOOKUP_TABLE default" << std::endl;

		myfile.close();
	}
	TRY( PetscBarrier(NULL) );
	
	/* each processor writes its own portion of data */
	for(id_proc=0;id_proc<psize;id_proc++){
		if(id_proc == prank){
				myfile.open(filename.c_str(), std::fstream::in | std::fstream::out | std::fstream::app);

				for(t=0;t < Tlocal;t++){
					myfile << gamma_max_id[t] << std::endl;
				}

				myfile.close();
		}

		TRY( PetscBarrier(NULL) );
	}


	timer_saveVTK.stop();
	coutMaster <<  " - problem saved to VTK in: " << timer_saveVTK.get_value_sum() << std::endl;

	LOG_FUNC_END	
}

template<>
void KmeansData<PetscVector>::load_gammavector(PetscVector &gamma0) const {
	LOG_FUNC_BEGIN
	
	/* get petsc Vec from provided vector - this vector is stride */
	Vec gamma0_Vec = gamma0.get_vector();
	
	/* variables */
	int K = this->get_K();
	int T = this->get_T();
	int Tlocal = this->get_Tlocal();
	int Tbegin = this->get_Tbegin();
	int xdim = this->get_xdim();

	/* prepare IS with my indexes in provided vector */
	IS gamma_sub_IS;
	
	/* fill the index sets */
	TRY( ISCreateStride(PETSC_COMM_WORLD, Tlocal*K, Tbegin*K, 1, &(gamma_sub_IS)) );

	/* now get subvector with my local values from provided stride vector */
	Vec gamma_sub;
	TRY( VecGetSubVector(gamma0_Vec, gamma_sub_IS, &gamma_sub) );

	/* prepare local vector */
	Vec gamma_local;
	#ifndef USE_CUDA
		TRY( VecCreateSeq(PETSC_COMM_SELF, K*Tlocal, &gamma_local) );
	#else
		TRY( VecCreateSeqCUDA(PETSC_COMM_SELF, K*Tlocal, &gamma_local) );
	#endif

	/* get the vector where I will store my values */
	TRY( VecGetLocalVector(gammavector->get_vector(), gamma_local) );

	/* now copy values from subvector to local vector */
	TRY( VecCopy(gamma_sub, gamma_local) );

	/* restore subvector */
	TRY( VecRestoreLocalVector(gammavector->get_vector(), gamma_local) );
	TRY( VecRestoreSubVector(gamma0_Vec, gamma_sub_IS, &gamma_sub) );

	/* destroy auxiliary index sets */
	TRY( ISDestroy(&gamma_sub_IS) );

	LOG_FUNC_END
}

template<>
void KmeansData<PetscVector>::load_gammavector(std::string filename) const {
	LOG_FUNC_BEGIN
	
	//TODO: control existence of file

	/* aux vector, we first oad data and then distribute values to procs */
	Vec gamma_preload_Vec;
	TRY( VecCreate(PETSC_COMM_WORLD, &gamma_preload_Vec) );

	/* prepare viewer to load from file */
	PetscViewer mviewer;
	TRY( PetscViewerCreate(PETSC_COMM_WORLD, &mviewer) );
	TRY( PetscViewerBinaryOpen(PETSC_COMM_WORLD ,filename.c_str(), FILE_MODE_READ, &mviewer) );
	
	/* load vector from viewer */
	TRY( VecLoad(gamma_preload_Vec, mviewer) );

	/* destroy the viewer */
	TRY( PetscViewerDestroy(&mviewer) );	

	PetscVector gamma_preload(gamma_preload_Vec);
	this->load_gammavector(gamma_preload);

//	TRY( VecDestroy(&gamma_preload_Vec) );
	
	LOG_FUNC_END
}


}
} /* end namespace */

#endif
