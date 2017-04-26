#include "external/petscvector/algebra/graph/bgmgraph.h"

namespace pascinference {
namespace algebra {

template<>
BGMGraph<PetscVector>::BGMGraph(const double *coordinates_array, int n, int dim){
	/* prepare vector from values */
	Vec vec_arr;
	TRYCXX( VecCreateSeqWithArray(PETSC_COMM_SELF,1,n*dim,coordinates_array,&vec_arr) );

	coordinates = new GeneralVector<PetscVector>(vec_arr);

	this->dim = dim;
	this->n = n;

	m = 0;
	m_max = 0;
	threshold = -1;
	processed = false;

	DD_decomposed = false;
}

template<>
void BGMGraph<PetscVector>::process(double threshold) {
	LOG_FUNC_BEGIN
	
	this->threshold = threshold;
	
	/* prepare array for number of neighbors */
	neighbor_nmbs = (int*)malloc(n*sizeof(int));

	#pragma omp parallel for
	for(int i=0;i<n;i++){
		neighbor_nmbs[i] = 0;
	}
	
	/* get local array and work with it */
	const double *coordinates_arr;
	TRYCXX( VecGetArrayRead(coordinates->get_vector(),&coordinates_arr) );
	
	/* go throught graph - compute number of neighbors */
//	#pragma omp parallel for
	for(int i=0;i<n;i++){
		for(int j=i+1;j<n;j++){
			if(compute_normsqr(coordinates_arr, i, j) < threshold*threshold){
				neighbor_nmbs[i] += 1;
				neighbor_nmbs[j] += 1;
			}
		}
	}

	/* prepare storages for neightbors ids */
	neighbor_ids = (int**)malloc(n*sizeof(int*));

	#pragma omp parallel for
	for(int i=0;i<n;i++){
		neighbor_ids[i] = (int*)malloc(neighbor_nmbs[i]*sizeof(int));
	}

	/* go throught graph - fill indexes of neighbors */
	int *counters;
	counters = (int*)malloc(n*sizeof(int));

	#pragma omp parallel for
	for(int i=0;i<n;i++){
		counters[i] = 0;
	}

//	#pragma omp parallel for // TODO: here is a problem, cannot be used, maybe because of couter arrays?
	for(int i=0;i<n;i++){
		for(int j=i+1;j<n;j++){
			if(compute_normsqr(coordinates_arr, i, j) < threshold*threshold){
				neighbor_ids[i][counters[i]] = j;
				neighbor_ids[j][counters[j]] = i;

				counters[i] += 1;
				counters[j] += 1;
				
				this->m += 1;
			}
		}

		/* compute m_max (max degree of vertex) */
		if(neighbor_nmbs[i] > m_max){
			this->m_max = neighbor_nmbs[i];
		}
	}
	free(counters);
	
	/* restore array */
	TRYCXX( VecRestoreArrayRead(coordinates->get_vector(),&coordinates_arr) );

	#ifdef USE_CUDA
		cuda_BGMGraph_process(neighbor_nmbs_gpu, neighbor_ids_gpu, neighbor_ids_cpugpu, n);
	#endif
	
	processed = true;

	LOG_FUNC_END
}

template<>
void BGMGraph<PetscVector>::saveVTK(std::string filename) const {
	LOG_FUNC_BEGIN
	
	Timer timer_saveVTK; 
	timer_saveVTK.restart();
	timer_saveVTK.start();
	
	/* to manipulate with file */
	std::ofstream myfile;	
	
	/* master writes the file */
	if(GlobalManager.get_rank() == 0){
		myfile.open(filename.c_str());

		/* write header to file */
		myfile << "# vtk DataFile Version 3.1" << std::endl;
		myfile << "PASCInference: Graph" << std::endl;
		myfile << "ASCII" << std::endl;
		myfile << "DATASET UNSTRUCTURED_GRID" << std::endl;

		/* write points - coordinates */
		myfile << "POINTS " << n << " FLOAT" << std::endl;
		const double *coordinates_arr;
		TRYCXX( VecGetArrayRead(coordinates->get_vector(),&coordinates_arr) );
		for(int i=0;i<n;i++){
			if(dim == 1){ 
				/* 1D sample */
				myfile << coordinates_arr[i] << " 0 0" << std::endl; /* x */
			}

			if(dim == 2){ 
				/* 2D sample */
				myfile << coordinates_arr[i] << " "; /* x */
				myfile << coordinates_arr[n+i] << " 0" << std::endl; /* y */
			}

			if(dim == 3){ 
				/* 3D sample */
				myfile << coordinates_arr[i] << " "; /* x */
				myfile << coordinates_arr[n+i] << " "; /* y */
				myfile << coordinates_arr[2*n+i] << std::endl; /* z */
			}

			if(dim > 3){
				//TODO ???
			}
		}
		TRYCXX( VecRestoreArrayRead(coordinates->get_vector(),&coordinates_arr) );
		
		/* write edges */
		myfile << "\nCELLS " << 2*m << " " << 2*m*3 << std::endl; /* actually, edges are here twice */
		for(int i=0;i<n;i++){
			for(int j=0;j<neighbor_nmbs[i];j++){
				myfile << "2 " << i << " " << neighbor_ids[i][j] << std::endl;
			}
		}
		myfile << "\nCELL_TYPES " << 2*m << std::endl;
		for(int i=0;i<2*m;i++){
			myfile << "3" << std::endl;
		}
		
		/* write domain affiliation */
		myfile << "\nPOINT_DATA " << n << std::endl;
		myfile << "SCALARS domain float 1" << std::endl;
		myfile << "LOOKUP_TABLE default" << std::endl;
		for(int i=0;i<n;i++){
			if(DD_decomposed){
				myfile << DD_affiliation[i] << std::endl;
			} else {
				myfile << "-1" << std::endl;
			}
		}
		
		myfile.close();
	}
	TRYCXX( PetscBarrier(NULL) );
	
	timer_saveVTK.stop();
	coutMaster <<  " - graph saved to VTK in: " << timer_saveVTK.get_value_sum() << std::endl;

	LOG_FUNC_END
}


}
} /* end of namespace */

