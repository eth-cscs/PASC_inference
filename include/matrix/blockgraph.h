#ifndef BLOCKGRAPHMATRIX_H
#define	BLOCKGRAPHMATRIX_H

#include "pascinference.h"

#ifndef USE_PETSCVECTOR
 #error 'BLOCKGRAPH is for PETSCVECTOR only, sorry'
#endif

#ifdef USE_PETSCVECTOR
	typedef petscvector::PetscVector PetscVector;
#endif

#ifdef USE_MINLIN
	typedef minlin::threx::HostVector<double> MinlinHostVector;
	typedef minlin::threx::DeviceVector<double> MinlinDeviceVector;
#endif

#include "petscdmda.h"

namespace pascinference {

class BGM_Graph {
	private:
		GeneralVector<PetscVector> coordinates; /**< vector with coordinates [p1_x, ... pn_x, p1_y, ... pn_y] */
		int dim; /**< dimension of coordinates */	
		int n; /**< number of vertices */
		double threshold; /**< the value used for processing the graph */
		
		bool processed; /**< there was process with threshold value called? */
		int *neighbor_nmbs;
		int **neighbor_ids;

		double compute_norm(const double *values, int idx1, int idx2){
			int d;
			double mynorm = 0;
			for(d=0;d<dim;d++){
				mynorm += (values[idx1+d*n] - values[idx2+d*n])*(values[idx1+d*n] - values[idx2+d*n]);
			}
			return mynorm;
		}
		
	public:
		BGM_Graph(std::string filename, int dim=2){
			/* load vertices from file */
			coordinates.load_local(filename);

			this->dim = dim;
			n = coordinates.size()/(double)dim;

			threshold = -1;
			processed = false;
		};

		~BGM_Graph(){
//			TRY( VecDestroy(&coordinates.get_vector()));

			/* if the graph was processed, then free memory */
			if(processed){
				free(neighbor_nmbs);
			}
			
		};
		
		int get_n(){
			return this->n;
		}

		int get_dim(){
			return this->dim;
		}

		double get_threshold(){
			return this->threshold;
		}
		
		void process_cpu(double threshold) {
			this->threshold = threshold;
			
			int i,j,d;
			double mynorm;

			/* prepare array for number of neighbors */
			neighbor_nmbs = (int*)malloc(n*sizeof(int));
			for(i=0;i<n;i++){
				neighbor_nmbs[i] = 0;
			}
			
			/* get local array and work with it */
			const double *coordinates_arr;
			TRY( VecGetArrayRead(coordinates.get_vector(),&coordinates_arr) );
			
			/* go throught graph - compute number of neighbors */
			for(i=0;i<n;i++){
				for(j=i+1;j<n;j++){
					if(compute_norm(coordinates_arr, i, j) < threshold){
						neighbor_nmbs[i] += 1;
						neighbor_nmbs[j] += 1;
					}
				}
			}

			/* prepare storages for neightbors ids */
			neighbor_ids = (int**)malloc(n*sizeof(int*));
			for(i=0;i<n;i++){
				neighbor_ids[i] = (int*)malloc(neighbor_nmbs[i]*sizeof(int));
			}

			/* go throught graph - fill indexes of neighbors */
			int *counters;
			counters = (int*)malloc(n*sizeof(int));
			for(i=0;i<n;i++){
				counters[i] = 0;
			}
			for(i=0;i<n;i++){
				for(j=i+1;j<n;j++){
					if(compute_norm(coordinates_arr, i, j) < threshold){
						neighbor_ids[i][counters[i]] = j;
						neighbor_ids[j][counters[j]] = i;

						counters[i] += 1;
						counters[j] += 1;
					}
				}
			}
			free(counters);
			
			/* restore array */
			TRY( VecRestoreArrayRead(coordinates.get_vector(),&coordinates_arr) );
			
			processed = true;
		}
		
		void print(ConsoleOutput &output) const {
			output << "Graph" << std::endl;
			
			output.push();
			output << " - dim:       " << this->dim << std::endl;
			output << " - vertices:  " << this->n << std::endl;
			output << " - threshold: " << this->threshold << std::endl;
			output << " - processed: " << this->processed << std::endl;
			output.pop();
			
		};
		
		void print_content(ConsoleOutput &output) const {
			output << "Graph" << std::endl;
			
			output.push();
			output << " - dim:       " << this->dim << std::endl;
			output << " - vertices:  " << this->n << std::endl;
			output << " - threshold: " << this->threshold << std::endl;
			output << " - processed: " << this->processed << std::endl;

			output << " - coordinates: " << coordinates << std::endl;

			int i,j;
			if(this->processed){
				output << " - processed arrays: " << std::endl;
				output.push();
				for(i=0;i<n;i++){
					output << i << ": " << "(" << neighbor_nmbs[i] << "): ";
					for(j=0;j<neighbor_nmbs[i];j++){
						output << neighbor_ids[i][j];
						if(j < neighbor_nmbs[i]-1){
							output << ", ";
						}
					}					
					output << std::endl;
				}
				output.pop();
			}

			output.pop();
			
		};
		
};

template<class VectorBase>
class BlockGraphMatrix: public GeneralMatrix<VectorBase> {
	private:
		int T; /**< dimension of each block */
		int Tlocal; /**< local dimension of each block */
		int Tbegin; /**< ownership begin */
		int Tend; /**< ownership end */
		
		int R; /**< number of vertices = number of blocks in row,col */
		int K; /**< number of diagonal blocks */
		double alpha; /**< general matrix multiplicator */
		
		BGM_Graph *graph; /**< graph with stucture of the matrix */
		
		#ifdef USE_GPU
			int blockSize; /**< block size returned by the launch configurator */
			int minGridSize; /**< the minimum grid size needed to achieve the maximum occupancy for a full device launch */
			int gridSize; /**< the actual grid size needed, based on input size */
		#endif	

		void matmult_tridiag(VectorBase &y, const VectorBase &x) const; /* y = 3diag(1,1,1)*x */
				
	public:
		BlockGraphMatrix(const VectorBase &x, BGM_Graph &new_graph, int K, double alpha=1.0);
		~BlockGraphMatrix(); /* destructor - destroy inner matrix */

		void print(ConsoleOutput &output) const; /* print matrix */
		void print(ConsoleOutput &output_global, ConsoleOutput &output_local) const;
		
		void printcontent(ConsoleOutput &output) const;
		std::string get_name() const;
		
		void matmult(VectorBase &y, const VectorBase &x) const; /* y = A*x */

		int get_R() const;
		int get_K() const;
		int get_T() const;
		int get_Tlocal() const;
		double get_alpha() const;

};

template<class VectorBase>
std::string BlockGraphMatrix<VectorBase>::get_name() const {
	return "BlockGraphMatrix";
}


template<class VectorBase>
BlockGraphMatrix<VectorBase>::BlockGraphMatrix(const VectorBase &x, BGM_Graph &new_graph, int K, double alpha){
	LOG_FUNC_BEGIN

	this->graph = &new_graph;
//	this->R = this->graph->get_n();
	this->R = 1; // todo: hot-fix test 

	this->K = K;
	this->alpha = alpha;
	
	/* get informations from given vector */
	int size, size_local, low, high;
	TRY( VecGetSize(x.get_vector(), &size) );
	TRY( VecGetLocalSize(x.get_vector(), &size_local) );
	TRY( VecGetOwnershipRange(x.get_vector(), &low, &high) );

	this->T = size/(double)(R*K);
	this->Tlocal = size_local/(double)(R*K);
	this->Tbegin = low/(double)(R*K);
	this->Tend = high/(double)(R*K);

	#ifdef USE_GPU
		//gpuErrchk( cudaOccupancyMaxPotentialBlockSize( &minGridSize, &blockSize,kernel_mult, 0, 0) );
		//gridSize = (K*T + blockSize - 1)/ blockSize;
		minGridSize = 0;
		gridSize = 0;
		blockSize = 0;
	#endif

	LOG_FUNC_END
}	


template<class VectorBase>
BlockGraphMatrix<VectorBase>::~BlockGraphMatrix(){
	LOG_FUNC_BEGIN
	
	//TODO: what to do?
	
	LOG_FUNC_END	
}

/* print matrix */
template<class VectorBase>
void BlockGraphMatrix<VectorBase>::print(ConsoleOutput &output) const		
{
	LOG_FUNC_BEGIN

	output << this->get_name() << std::endl;
	output << " - T:     " << T << std::endl;
	output << " - R:     " << R << std::endl;
	output << " - K:     " << K << std::endl;
	output << " - size:  " << T*R*K << std::endl;

	output << " - alpha: " << alpha << std::endl;

	#ifdef USE_GPU
		output <<  " - blockSize:   " << blockSize << std::endl;
		output <<  " - gridSize:    " << gridSize << std::endl;
		output <<  " - minGridSize: " << minGridSize << std::endl;
	#endif
	
	LOG_FUNC_END	
}

template<class VectorBase>
void BlockGraphMatrix<VectorBase>::print(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN

	output_global << this->get_name() << std::endl;
	output_global << " - T:     " << T << std::endl;
	output_global.push();
		output_local << " - Tlocal:  " << Tlocal << " (" << Tbegin << "," << Tend << ")" << std::endl;
		output_local.synchronize();
	output_global.pop();
	
	output_global << " - R:     " << R << std::endl;
	output_global << " - K:     " << K << std::endl;
	output_global << " - size:  " << T*R*K << std::endl;
	output_global.push();
		output_local << " - sizelocal:  " << Tlocal*R*K << std::endl;
		output_local.synchronize();
	output_global.pop();

	output_global << " - alpha: " << alpha << std::endl;

	#ifdef USE_GPU
		output_global <<  " - blockSize:   " << blockSize << std::endl;
		output_global <<  " - gridSize:    " << gridSize << std::endl;
		output_global <<  " - minGridSize: " << minGridSize << std::endl;
	#endif

	LOG_FUNC_END	
}

/* print matrix */
template<class VectorBase>
void BlockGraphMatrix<VectorBase>::printcontent(ConsoleOutput &output) const		
{
	LOG_FUNC_BEGIN
	
	
	LOG_FUNC_END
}

/* matrix-vector multiplication */
template<class VectorBase>
void BlockGraphMatrix<VectorBase>::matmult(VectorBase &y, const VectorBase &x) const { 
	LOG_FUNC_BEGIN
	
	// TODO: maybe y is not initialized, who knows

	/* at first multiply y = 3diag(1,1,1)*x */
	matmult_tridiag(y, x);

	LOG_FUNC_END	
}

template<class VectorBase>
int BlockGraphMatrix<VectorBase>::get_R() const { 
	return this->R;
}

template<class VectorBase>
int BlockGraphMatrix<VectorBase>::get_K() const { 
	return this->K;
}

template<class VectorBase>
int BlockGraphMatrix<VectorBase>::get_T() const { 
	return this->T;
}

template<class VectorBase>
int BlockGraphMatrix<VectorBase>::get_Tlocal() const { 
	return this->Tlocal;
}

template<class VectorBase>
double BlockGraphMatrix<VectorBase>::get_alpha() const { 
	return this->alpha;
}

template<class VectorBase>
void BlockGraphMatrix<VectorBase>::matmult_tridiag(VectorBase &y, const VectorBase &x) const {
	double *y_arr;
	const double *x_arr;
	
	int coor_x, coor_y, coor_m, coor_n;
	
	DM da;
	TRY( DMDACreate2d(PETSC_COMM_WORLD, 
						DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, /* type of ghost nodes the array have */
						DMDA_STENCIL_STAR,	
                        K*R,T, /* global dimension of array */
                        1, PETSC_DECIDE, /* number of processors (I consider one row) */
                        1, /* dof ? */
                        1, /* stencil width */
                        PETSC_NULL,PETSC_NULL, // todo: ? 
                        &da) );

//	TRY( DMView(da,PETSC_VIEWER_STDOUT_WORLD) );
	TRY( DMDAGetCorners(da, &coor_x, &coor_y, NULL, &coor_m, &coor_n, NULL) );

	coutAll << "coors: " << coor_x << ", "  << coor_y << ", "  << coor_m << ", "  << coor_n << std::endl;
	coutAll.synchronize();

	Vec newvec;
	TRY( DMCreateGlobalVector(da, &newvec) );
	TRY( VecView(newvec, PETSC_VIEWER_STDOUT_WORLD));


	TRY( DMDAVecGetArray(da,newvec,&y_arr) );
//	TRY( DMDAVecGetArray(da,y.get_vector(),&y_arr) );
//	TRY( DMDAVecGetArrayRead(da,x.get_vector(),&x_arr) );

	for(int id_row=0;id_row<coor_m*coor_n;id_row++){
//		coutAll << id_row << ":" << x_arr[id_row] << std::endl;
//		y_arr[id_row] = 20;
	}
//	coutAll.synchronize();


	TRY( DMDAVecRestoreArray(da,newvec,&y_arr) );
//	TRY( DMDAVecRestoreArray(da,y.get_vector(),&y_arr) );
//	TRY( DMDAVecRestoreArrayRead(da,x.get_vector(),&x_arr) );

	
//	TRY( VecGetArray(y.get_vector(),&y_arr) );
//	TRY( VecGetArrayRead(x.get_vector(),&x_arr) );

	/* use openmp */
//	#pragma omp parallel for
//	for(int id_row=0;id_row<Tlocal*R*K;id_row++){
//		int k = (int)(id_row/(double)(Tlocal*R));
//		int r = (int)((id_row - k*Tlocal*R)/(double)Tlocal);
//		int t = id_row - k*Tlocal*R - r*Tlocal;
//		int t_global = Tbegin+t;

//		double value;

		/* first row */
//		if(t_global == 0){
//			value = 2;//alpha*x_arr[id_row] + alpha*x_arr[id_row+1];
//		}
		/* common row */
//		if(t_global > 0 && t_global < T-1){
//			value = 3;//alpha*x_arr[id_row-1] + alpha*x_arr[id_row] + alpha*x_arr[id_row+1];
//		}
		/* last row */
//		if(t_global == T-1){
//			value = 2;//alpha*x_arr[id_row-1] + alpha*x_arr[id_row];
//		}
		
//		y_arr[id_row] = value;
//	}

//	TRY( VecRestoreArray(y.get_vector(),&y_arr) );
//	TRY( VecRestoreArrayRead(x.get_vector(),&x_arr) );

	TRY( DMDestroy(&da) );

}



} /* end of namespace */


#endif
