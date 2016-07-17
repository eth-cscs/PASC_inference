#ifndef BLOCKGRAPHMATRIX_H
#define	BLOCKGRAPHMATRIX_H

#include "pascinference.h"

#ifndef USE_PETSCVECTOR
 #error 'BLOCKGRAPH is for PETSCVECTOR only, sorry'
#endif

#ifdef USE_CUDA
    #include <../src/vec/vec/impls/seq/seqcuda/cudavecimpl.h>
#endif

#ifdef USE_PETSCVECTOR
	typedef petscvector::PetscVector PetscVector;
#endif

namespace pascinference {

class BGM_Graph {
	private:
		GeneralVector<PetscVector> *coordinates; /**< vector with coordinates [p1_x, ... pn_x, p1_y, ... pn_y] */
		int dim; /**< dimension of coordinates */	
		int n; /**< number of vertices */
		double threshold; /**< the value used for processing the graph */
		
		bool processed; /**< there was process with threshold value called? */
		int *neighbor_nmbs; /**< number of neighbors for each node */
		int **neighbor_ids; /**< indexes of neighbors for each node */

		#ifdef USE_GPU
			int *neighbor_nmbs_gpu; /**< copy of values on GPU */
			int **neighbor_ids_gpu; /**< copy of values on GPU */
		#endif	

		double compute_norm(const double *values, int idx1, int idx2);
		
	public:
		BGM_Graph(std::string filename, int dim=2);

		BGM_Graph(const double *coordinates_array, int n, int dim);

		~BGM_Graph();
		
		int get_n();
		int get_dim();
		double get_threshold();
		int *get_neighbor_nmbs();
		int **get_neighbor_ids();

		int *get_neighbor_nmbs_gpu();
		int **get_neighbor_ids_gpu();
		GeneralVector<PetscVector> *get_coordinates();
		
		void process(double threshold);

		void print(ConsoleOutput &output) const;
		void print_content(ConsoleOutput &output) const;
		
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
			int blockSize1; /**< block size returned by the launch configurator for 3diag mult */
			int minGridSize1; /**< the minimum grid size needed to achieve the maximum occupancy for a full device launch for 3diag mult */
			int gridSize1; /**< the actual grid size needed, based on input size for 3diag mult */

			int blockSize2; /**< block size returned by the launch configurator for graph mult */
			int minGridSize2; /**< the minimum grid size needed to achieve the maximum occupancy for a full device launch for graph mult */
			int gridSize2; /**< the actual grid size needed, based on input size for graph mult */
		#endif

		void matmult_tridiag(const VectorBase &x) const; /* x_aux = 3diag(1,1,1)*x */
		void matmult_graph(VectorBase &y, const VectorBase &x) const; /* y = Wgraph *kr* x_aux */
		
		/* stuff for 3diag mult */
		#ifdef USE_PETSCVECTOR
			Vec x_aux;
			int left_overlap, right_overlap;
			IS left_overlap_is;
			IS right_overlap_is;
		#endif
		
		GeneralVector<VectorBase> *coeffs; /**< vector of coefficient for each block */
		
	public:
		BlockGraphMatrix(const VectorBase &x, BGM_Graph &new_graph, int K, double alpha=1.0, GeneralVector<VectorBase> *new_coeffs=NULL);
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

#ifdef USE_GPU
__global__ void kernel_BlockGraphMatrix_mult_tridiag(double* y_arr, double* x_arr, double *left_overlap_arr, double *right_overlap_arr, int T, int Tlocal, int Tbegin, int R, int K, double alpha);
__global__ void kernel_BlockGraphMatrix_mult_graph(double* y_arr, double* x_arr, double *x_aux_arr, int *neighbor_nmbs, int **neightbor_ids, int T, int Tlocal, int Tbegin, int R, int K, double alpha);
#endif


/* ---------------- IMPLEMENTATION -------------------- */

template<class VectorBase>
std::string BlockGraphMatrix<VectorBase>::get_name() const {
	return "BlockGraphMatrix";
}


template<class VectorBase>
BlockGraphMatrix<VectorBase>::BlockGraphMatrix(const VectorBase &x, BGM_Graph &new_graph, int K, double alpha, GeneralVector<VectorBase> *new_coeffs){
	LOG_FUNC_BEGIN

	this->graph = &new_graph;
	this->R = this->graph->get_n();

	this->K = K;
	this->alpha = alpha;
	
	this->coeffs = new_coeffs;
	
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
		gpuErrchk( cudaOccupancyMaxPotentialBlockSize( &minGridSize1, &blockSize1,kernel_BlockGraphMatrix_mult_tridiag, 0, 0) );
		gridSize1 = (Tlocal*K*R + blockSize1 - 1)/ blockSize1;

		gpuErrchk( cudaOccupancyMaxPotentialBlockSize( &minGridSize2, &blockSize2,kernel_BlockGraphMatrix_mult_graph, 0, 0) );
		gridSize2 = (Tlocal*K*R + blockSize2 - 1)/ blockSize2;
	#endif

	/* unfortunatelly, we need additional aux vector */
	TRY( VecDuplicate(x.get_vector(),&x_aux) );

	/* get ranges of all processors - necessary to compute overlaping indexes */
	const int *ranges;
	ranges = (int*)malloc((GlobalManager.get_size()+1)*sizeof(int));
    TRY( VecGetOwnershipRanges(x_aux,&ranges) );
	int myrank = GlobalManager.get_rank();

	/* create IS for 3diag mult */
	if(this->Tbegin > 0){
		left_overlap = 1;
		TRY( ISCreateStride(PETSC_COMM_SELF, R*K, ranges[myrank-1] + (ranges[myrank]-ranges[myrank-1])/(double)(R*K)-1 ,(ranges[myrank]-ranges[myrank-1])/(double)(R*K), &left_overlap_is) );
	} else {
		left_overlap = 0;
		TRY( ISCreateStride(PETSC_COMM_SELF, 0, 0, 0, &left_overlap_is) );
	}
	if(this->Tend < T){
		right_overlap = 1;
		TRY( ISCreateStride(PETSC_COMM_SELF, R*K, ranges[myrank+1], (ranges[myrank+2]-ranges[myrank+1])/(double)(R*K), &right_overlap_is) );
	} else {
		right_overlap = 0;
		TRY( ISCreateStride(PETSC_COMM_SELF, 0, 0, 0, &right_overlap_is) );
	}

	LOG_FUNC_END
}	


template<class VectorBase>
BlockGraphMatrix<VectorBase>::~BlockGraphMatrix(){
	LOG_FUNC_BEGIN
	
	/* destroy aux vector */
//	TRY( VecDestroy(&x_aux) );
	
	/* destroy IS for 3diag mult */
//	TRY( ISDestroy(&left_overlap_is) );
//	TRY( ISDestroy(&right_overlap_is) );
	
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

	if(coeffs){
		output << " - coeffs: " << *coeffs << std::endl;
	}

	#ifdef USE_GPU
		output <<  " - blockSize1:   " << blockSize1 << std::endl;
		output <<  " - gridSize1:    " << gridSize1 << std::endl;
		output <<  " - minGridSize1: " << minGridSize1 << std::endl;

		output <<  " - blockSize2:   " << blockSize2 << std::endl;
		output <<  " - gridSize2:    " << gridSize2 << std::endl;
		output <<  " - minGridSize2: " << minGridSize2 << std::endl;
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

	if(coeffs){
		output_local << " - coeffs: " << *coeffs << std::endl;
		output_local.synchronize();
	}

	#ifdef USE_GPU
		output_global <<  " - blockSize1:   " << blockSize1 << std::endl;
		output_global <<  " - gridSize1:    " << gridSize1 << std::endl;
		output_global <<  " - minGridSize1: " << minGridSize1 << std::endl;

		output_global <<  " - blockSize2:   " << blockSize2 << std::endl;
		output_global <<  " - gridSize2:    " << gridSize2 << std::endl;
		output_global <<  " - minGridSize2: " << minGridSize2 << std::endl;
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

	/* at first multiply x_aux = 3diag(1,1,1)*x */
	matmult_tridiag(x);

	/* now perform kronecker product with graph matrix */
	matmult_graph(y, x);

//	TRY(VecCopy(x_aux,y.get_vector()));

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

#ifndef USE_GPU
/* A*x using openmp */
template<>
void BlockGraphMatrix<PetscVector>::matmult_tridiag(const PetscVector &x) const {
	Vec x_sub;
	double *y_arr;
	const double *x_arr;
	
	Vec left_overlap_vec;
	Vec right_overlap_vec;
	const double *left_overlap_arr;
	const double *right_overlap_arr;

	/* get subvector and array */
	TRY( VecGetSubVector(x.get_vector(),left_overlap_is,&left_overlap_vec) );
	TRY( VecGetArrayRead(left_overlap_vec,&left_overlap_arr) );

	TRY( VecGetSubVector(x.get_vector(),right_overlap_is,&right_overlap_vec) );
	TRY( VecGetArrayRead(right_overlap_vec,&right_overlap_arr) );
	
	TRY( VecGetArrayRead(x.get_vector(),&x_arr) );
	TRY( VecGetArray(x_aux,&y_arr) );

	/* use openmp */
	#pragma omp parallel for
	for(int y_arr_idx=0;y_arr_idx<Tlocal*K*R;y_arr_idx++){
		int k = floor(y_arr_idx/(double)(Tlocal*R));
		int r = floor((y_arr_idx-k*Tlocal*R)/(double)(Tlocal));
		int tlocal = y_arr_idx - (k*R + r)*Tlocal;
		int tglobal = Tbegin+tlocal;
		int x_arr_idx = tlocal + (k*R+r)*Tlocal;
		int overlap_id = k*R+r;

		double value;
		double right_value;
		double left_value;

		/* first row */
		if(tglobal == 0){
			if(tlocal+1 >= Tlocal){
				right_value = right_overlap_arr[overlap_id];
			} else {
				right_value = x_arr[x_arr_idx+1];
			}
			value = alpha*x_arr[x_arr_idx] + alpha*right_value;
		}
		/* common row */
		if(tglobal > 0 && tglobal < T-1){
			if(tlocal+1 >= Tlocal){
				right_value = right_overlap_arr[overlap_id];
			} else {
				right_value = x_arr[x_arr_idx+1];
			}
			if(tlocal-1 < 0){
				left_value = left_overlap_arr[overlap_id];
			} else {
				left_value = x_arr[x_arr_idx-1];
			}
			value = alpha*left_value + alpha*x_arr[x_arr_idx] + alpha*right_value;
		}
		/* last row */
		if(tglobal == T-1){
			if(tlocal-1 < 0){
				left_value = left_overlap_arr[overlap_id];
			} else {
				left_value = x_arr[x_arr_idx-1];
			}
			value = alpha*left_value + alpha*x_arr[x_arr_idx];
		}
		
		y_arr[y_arr_idx] = value; /* = (k+1)*1000+(r+1)*100+(tglobal+1); test */
	}

	/* restore subvector and array */
	TRY( VecRestoreArrayRead(x.get_vector(),&x_arr) );
	TRY( VecRestoreArray(x_aux,&y_arr) );

	TRY( VecRestoreArrayRead(left_overlap_vec,&left_overlap_arr) );
	TRY( VecRestoreSubVector(x.get_vector(),left_overlap_is,&left_overlap_vec) );

	TRY( VecRestoreArrayRead(right_overlap_vec,&right_overlap_arr) );
	TRY( VecRestoreSubVector(x.get_vector(),right_overlap_is,&right_overlap_vec) );

}

template<>
void BlockGraphMatrix<PetscVector>::matmult_graph(PetscVector &y, const PetscVector &x) const {
	double *y_arr;
	const double *x_arr;
	const double *x_aux_arr;
	const double *coeffs_arr;
	
	int* neighbor_nmbs = graph->get_neighbor_nmbs();
	int **neightbor_ids = graph->get_neighbor_ids();
			
	/* get array */
	TRY( VecGetArrayRead(x.get_vector(),&x_arr) );
	TRY( VecGetArrayRead(x_aux,&x_aux_arr) );
	TRY( VecGetArray(y.get_vector(),&y_arr) );
	if(coeffs){
		TRY( VecGetArrayRead(coeffs->get_vector(),&coeffs_arr) );
	}

	/* use openmp */
	#pragma omp parallel for
	for(int y_arr_idx=0;y_arr_idx<Tlocal*K*R;y_arr_idx++){
		int k = floor(y_arr_idx/(double)(Tlocal*R));
		int r = floor((y_arr_idx-k*Tlocal*R)/(double)(Tlocal));
		int tlocal = y_arr_idx - (k*R + r)*Tlocal;
		int tglobal = Tbegin+tlocal;
		int x_arr_idx = tlocal + (k*R+r)*(Tlocal);

		double value;
		int Wsum;

		/* compute sum of W entries in row */
		if(tglobal == 0 || tglobal == T-1){
			Wsum = 2*neighbor_nmbs[r];
		} else {
			Wsum = 3*neighbor_nmbs[r];
		}

		/* diagonal entry */
		y_arr[y_arr_idx] = alpha*Wsum*x_arr[x_arr_idx];

		/* non-diagonal entries */
		for(int neighbor=0;neighbor<neighbor_nmbs[r];neighbor++){
			y_arr[y_arr_idx] -= x_aux_arr[k*Tlocal*R + (neightbor_ids[r][neighbor])*Tlocal + tlocal];
		}
		
		/* if coeffs are provided, then multiply with coefficient corresponding to this block */
		if(coeffs){
			y_arr[y_arr_idx] = coeffs_arr[k]*coeffs_arr[k]*y_arr[y_arr_idx];

			// temp: + Laplace
/*			if(tglobal == 0 || tglobal == T){
				y_arr[y_arr_idx] +=  alpha*(-x_aux_arr[x_arr_idx] + 2*x_arr[x_arr_idx]);
			} else {
				y_arr[y_arr_idx] +=  alpha*(-x_aux_arr[x_arr_idx] + 3*x_arr[x_arr_idx]);
			}
*/
		}

	}

	/* restore array */
	TRY( VecRestoreArrayRead(x.get_vector(),&x_arr) );
	TRY( VecRestoreArrayRead(x_aux,&x_aux_arr) );
	TRY( VecRestoreArray(y.get_vector(),&y_arr) );
	if(coeffs){
		TRY( VecRestoreArrayRead(coeffs->get_vector(),&coeffs_arr) );
	}

}

#else

/* A*x using CUDA kernel */
__global__ void kernel_BlockGraphMatrix_mult_tridiag(double* y_arr, double* x_arr, double *left_overlap_arr, double *right_overlap_arr, int T, int Tlocal, int Tbegin, int R, int K, double alpha)
{
	/* compute my id */
	int y_arr_idx = blockIdx.x*blockDim.x + threadIdx.x;

	if(y_arr_idx < K*Tlocal*R){
		int k = floor(y_arr_idx/(double)(Tlocal*R));
		int r = floor((y_arr_idx-k*Tlocal*R)/(double)(Tlocal));
		int tlocal = y_arr_idx - (k*R + r)*Tlocal;
		int tglobal = Tbegin+tlocal;
		int x_arr_idx = tlocal + (k*R+r)*Tlocal;
		int overlap_id = k*R+r;
		
		double value;
		double right_value;
		double left_value;

		/* first row */
		if(tglobal == 0){
			if(tlocal+1 >= Tlocal){
				right_value = right_overlap_arr[overlap_id];
			} else {
				right_value = x_arr[x_arr_idx+1];
			}
			value = alpha*x_arr[x_arr_idx] - alpha*right_value;
		}
		/* common row */
		if(tglobal > 0 && tglobal < T-1){
			if(tlocal+1 >= Tlocal){
				right_value = right_overlap_arr[overlap_id];
			} else {
				right_value = x_arr[x_arr_idx+1];
			}
			if(tlocal-1 < 0){
				left_value = left_overlap_arr[overlap_id];
			} else {
				left_value = x_arr[x_arr_idx-1];
			}
			value = -alpha*left_value + 2*alpha*x_arr[x_arr_idx] - alpha*right_value;
		}
		/* last row */
		if(tglobal == T-1){
			if(tlocal-1 < 0){
				left_value = left_overlap_arr[overlap_id];
			} else {
				left_value = x_arr[x_arr_idx-1];
			}
			value = -alpha*left_value + alpha*x_arr[x_arr_idx];
		}
		
		y_arr[y_arr_idx] = value; 
	}

	/* if id_row >= K*T*R then relax and do nothing */	

}

template<>
void BlockGraphMatrix<PetscVector>::matmult_tridiag(const PetscVector &x) const {
	double *y_arr;
	double *x_arr;
	
	Vec left_overlap_vec;
	Vec right_overlap_vec;
	double *left_overlap_arr;
	double *right_overlap_arr;

	/* get subvector and array */
	TRY( VecGetSubVector(x.get_vector(),left_overlap_is,&left_overlap_vec) );
	TRY( VecCUDAGetArrayReadWrite(left_overlap_vec,&left_overlap_arr) );

	TRY( VecGetSubVector(x.get_vector(),right_overlap_is,&right_overlap_vec) );
	TRY( VecCUDAGetArrayReadWrite(right_overlap_vec,&right_overlap_arr) );
	
	TRY( VecCUDAGetArrayReadWrite(x.get_vector(),&x_arr) );
	TRY( VecCUDAGetArrayReadWrite(x_aux,&y_arr) );

	/* call kernel */
	kernel_BlockGraphMatrix_mult_tridiag<<<gridSize1, blockSize1>>>(y_arr,x_arr,left_overlap_arr,right_overlap_arr,T,Tlocal,Tbegin,R,K,alpha);
	gpuErrchk( cudaDeviceSynchronize() );
	
	/* restore subvector and array */
	TRY( VecCUDARestoreArrayReadWrite(x.get_vector(),&x_arr) );
	TRY( VecCUDARestoreArrayReadWrite(x_aux,&y_arr) );

	TRY( VecCUDARestoreArrayReadWrite(left_overlap_vec,&left_overlap_arr) );
	TRY( VecRestoreSubVector(x.get_vector(),left_overlap_is,&left_overlap_vec) );

	TRY( VecCUDARestoreArrayReadWrite(right_overlap_vec,&right_overlap_arr) );
	TRY( VecRestoreSubVector(x.get_vector(),right_overlap_is,&right_overlap_vec) );
}

__global__ void kernel_BlockGraphMatrix_mult_graph(double* y_arr, double* x_arr, double *x_aux_arr, int *neighbor_nmbs, int **neightbor_ids, int T, int Tlocal, int Tbegin, int R, int K, double alpha)
{
	/* compute my id */
	int y_arr_idx = blockIdx.x*blockDim.x + threadIdx.x;

	if(y_arr_idx < K*Tlocal*R){
		int k = floor(y_arr_idx/(double)(Tlocal*R));
		int r = floor((y_arr_idx-k*Tlocal*R)/(double)(Tlocal));
		int tlocal = y_arr_idx - (k*R + r)*Tlocal;
		int tglobal = Tbegin+tlocal;
		int x_arr_idx = tlocal + (k*R+r)*(Tlocal);

		int Wsum;

		/* compute sum of W entries in row */
		if(tglobal == 0 || tglobal == T-1){
			Wsum = 2*neighbor_nmbs[r];
		} else {
			Wsum = 3*neighbor_nmbs[r];
		}

		/* diagonal entry */
		y_arr[y_arr_idx] = alpha*Wsum*x_arr[x_arr_idx];

		/* non-diagonal entries */
		for(int neighbor=0;neighbor<neighbor_nmbs[r];neighbor++){
			y_arr[y_arr_idx] -= x_aux_arr[k*Tlocal*R + (neightbor_ids[r][neighbor])*Tlocal + tlocal];
		}
	}

	/* if id_row >= K*T*R then relax and do nothing */	

}

template<>
void BlockGraphMatrix<PetscVector>::matmult_graph(PetscVector &y, const PetscVector &x) const {
	double *y_arr;
	double *x_arr;
	double *x_aux_arr;
	
	int* neighbor_nmbs = graph->get_neighbor_nmbs_gpu();
	int **neightbor_ids = graph->get_neighbor_ids_gpu();
			
	/* get array */
	TRY( VecCUDAGetArrayReadWrite(x.get_vector(),&x_arr) );
	TRY( VecCUDAGetArrayReadWrite(x_aux,&x_aux_arr) );
	TRY( VecCUDAGetArrayReadWrite(y.get_vector(),&y_arr) );

	/* call kernel */
	kernel_BlockGraphMatrix_mult_graph<<<gridSize2, blockSize2>>>(y_arr, x_arr, x_aux_arr, neighbor_nmbs, neightbor_ids, T, Tlocal, Tbegin, R, K, alpha);
	gpuErrchk( cudaDeviceSynchronize() );

	/* restore array */
	TRY( VecCUDARestoreArrayReadWrite(x.get_vector(),&x_arr) );
	TRY( VecCUDARestoreArrayReadWrite(x_aux,&x_aux_arr) );
	TRY( VecCUDARestoreArrayReadWrite(y.get_vector(),&y_arr) );

}

#endif

/* ----------------- BGM_Graph implementation ------------- */

double BGM_Graph::compute_norm(const double *values, int idx1, int idx2){
	int d;
	double mynorm = 0;
	for(d=0;d<dim;d++){
		mynorm += (values[idx1+d*n] - values[idx2+d*n])*(values[idx1+d*n] - values[idx2+d*n]);
	}
	return mynorm;
}

BGM_Graph::BGM_Graph(std::string filename, int dim){
	coordinates = new GeneralVector<PetscVector>();
	
	/* load vertices from file */
	coordinates->load_local(filename);

	this->dim = dim;
	n = coordinates->size()/(double)dim;

	threshold = -1;
	processed = false;
}

BGM_Graph::BGM_Graph(const double *coordinates_array, int n, int dim){
	/* prepare vector from values */
	Vec vec_arr;
	TRY( VecCreateSeqWithArray(PETSC_COMM_SELF,1,n*dim,coordinates_array,&vec_arr) );

	coordinates = new GeneralVector<PetscVector>(vec_arr);

	this->dim = dim;
	this->n = n;

	threshold = -1;
	processed = false;
}

BGM_Graph::~BGM_Graph(){
//	TRY( VecDestroy(&coordinates.get_vector()));
//	free(coordinates);

	/* if the graph was processed, then free memory */
	if(processed){
		free(neighbor_nmbs);
		int i;
		for(i=0;i<n;i++){
			free(neighbor_ids[i]);
		}
		free(neighbor_ids);

		#ifdef USE_GPU
			gpuErrchk( cudaFree(neighbor_nmbs_gpu) );
			for(i=0;i<n;i++){
				gpuErrchk( cudaFree(neighbor_ids_gpu[i]) );
			}
			gpuErrchk( cudaFree(neighbor_ids_gpu) );
		#endif

	}
	
}

int BGM_Graph::get_n(){
	return this->n;
}

int BGM_Graph::get_dim(){
	return this->dim;
}

double BGM_Graph::get_threshold(){
	return this->threshold;
}

int *BGM_Graph::get_neighbor_nmbs(){
	return neighbor_nmbs;
}

int **BGM_Graph::get_neighbor_ids(){
	return neighbor_ids;
}

int *BGM_Graph::get_neighbor_nmbs_gpu(){
	#ifdef USE_GPU
		return neighbor_nmbs_gpu;
	#else
		return neighbor_nmbs;
	#endif
}

int **BGM_Graph::get_neighbor_ids_gpu(){
	#ifdef USE_GPU
		return neighbor_ids_gpu;
	#else
		return neighbor_ids;
	#endif
}

GeneralVector<PetscVector> *BGM_Graph::get_coordinates(){
	return coordinates;
}

void BGM_Graph::process(double threshold) {
	this->threshold = threshold;
	
	int i,j;

	/* prepare array for number of neighbors */
	neighbor_nmbs = (int*)malloc(n*sizeof(int));
	for(i=0;i<n;i++){
		neighbor_nmbs[i] = 0;
	}
	
	/* get local array and work with it */
	const double *coordinates_arr;
	TRY( VecGetArrayRead(coordinates->get_vector(),&coordinates_arr) );
	
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
	TRY( VecRestoreArrayRead(coordinates->get_vector(),&coordinates_arr) );
	
	#ifdef USE_GPU
		/* copy data to gpu */
		gpuErrchk( cudaMalloc((void **)&neighbor_nmbs_gpu, n*sizeof(int)) );	
		gpuErrchk( cudaMemcpy( neighbor_nmbs_gpu, neighbor_nmbs, n*sizeof(int), cudaMemcpyHostToDevice) );
		
		gpuErrchk( cudaMalloc((void **)&neighbor_ids_gpu, n*sizeof(int)) );	
		for(i=0;i<n;i++){
			gpuErrchk( cudaMalloc((void **)&(neighbor_ids_gpu[i]), neighbor_nmbs[i]*sizeof(int)) );
			gpuErrchk( cudaMemcpy( neighbor_ids_gpu[i], neighbor_ids[i], n*sizeof(int), cudaMemcpyHostToDevice) );
		}

		gpuErrchk( cudaDeviceSynchronize() );
	#endif
	
	processed = true;
}

void BGM_Graph::print(ConsoleOutput &output) const {
	output << "Graph" << std::endl;
	
	output.push();
	output << " - dim:       " << this->dim << std::endl;
	output << " - vertices:  " << this->n << std::endl;
	output << " - threshold: " << this->threshold << std::endl;
	output << " - processed: " << this->processed << std::endl;
	output.pop();
	
}

void BGM_Graph::print_content(ConsoleOutput &output) const {
	output << "Graph" << std::endl;
	
	output.push();
	output << " - dim:       " << this->dim << std::endl;
	output << " - vertices:  " << this->n << std::endl;
	output << " - threshold: " << this->threshold << std::endl;
	output << " - processed: " << this->processed << std::endl;

	output << " - coordinates: " << *coordinates << std::endl;

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
}



} /* end of namespace */


#endif
