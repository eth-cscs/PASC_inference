/** @file bgmgraph.h
 *  @brief class for manipulation with graphs
 *
 *  Defines some basic functions for manipulaton with graphs, especially for BLOCKGRAPHMATRIX.
 *
 *  @author Lukas Pospisil
 */
 
#ifndef PASC_COMMON_BGMGRAPH_H
#define	PASC_COMMON_BGMGRAPH_H

namespace pascinference {

/** \class BGMGraph
 *  \brief General class for manipulation with graphs.
 *
*/
class BGMGraph {
	protected:
		GeneralVector<PetscVector> *coordinates; /**< vector with coordinates [p1_x, ... pn_x, p1_y, ... pn_y] */
		int dim; /**< dimension of coordinates */	
		int n; /**< number of nodes */
		int m; /**< number of edges */
		int m_max; /**< maximum degree of vertices */
		double threshold; /**< the value used for processing the graph */
		
		bool processed; /**< there was process with threshold value called? */
		int *neighbor_nmbs; /**< number of neighbors for each node */
		int **neighbor_ids; /**< indexes of neighbors for each node */

		#ifdef USE_GPU
			int *neighbor_nmbs_gpu; /**< copy of values on GPU */
			int **neighbor_ids_gpu; /**< copy of values on GPU */
		#endif	

		/** @brief compute distance between two vertices
		*
		*  @param values array with coordinates of vertices
		*  @param idx1 index of vertex
		*  @param idx2 index of vertex
		*/
		double compute_norm(const double *values, int idx1, int idx2);
		
	public:
		BGMGraph();
		BGMGraph(std::string filename, int dim=2);
		BGMGraph(const double *coordinates_array, int n, int dim);

		~BGMGraph();
		
		/** @brief return number of vertices
		*/
		int get_n();

		/** @brief return number of edges
		*/
		int get_m();

		/** @brief return degree of vertex with max degree 
		*/
		int get_m_max();

		/** @brief return dimension of coordinates of vertices
		*/
		int get_dim();

		/** @brief return value used for filling graph with edges
		*/
		double get_threshold();

		/** @brief return array containing number of neighbors of vertices
		*/
		int *get_neighbor_nmbs();

		/** @brief return array of arrays containing indexes of neighbors of vertices
		*/
		int **get_neighbor_ids();

		/** @brief return array containing number of neighbors of vertices stored on GPU
		*/
		int *get_neighbor_nmbs_gpu();

		/** @brief return array of arrays containing indexes of neighbors of vertices stored on GPU
		*/
		int **get_neighbor_ids_gpu();

		/** @brief return vector with coordinates of vertices
		*/
		GeneralVector<PetscVector> *get_coordinates();
		
		/** @brief fill graph with edges based on given length of edge
		*/
		virtual void process(double threshold);

		/** @brief print basic informations of graph
		*/
		void print(ConsoleOutput &output) const;

		/** @brief print all vertices, edges, neighbors
		*/
		void print_content(ConsoleOutput &output) const;
		
};

/* ----------------- BGMGraph implementation ------------- */

double BGMGraph::compute_norm(const double *values, int idx1, int idx2){
	int d;
	double mynorm = 0;
	for(d=0;d<dim;d++){
		mynorm += (values[idx1+d*n] - values[idx2+d*n])*(values[idx1+d*n] - values[idx2+d*n]);
	}
	return mynorm;
}

BGMGraph::BGMGraph(std::string filename, int dim){
	coordinates = new GeneralVector<PetscVector>();
	
	/* load nodes from file */
	coordinates->load_local(filename);

	this->dim = dim;
	n = coordinates->size()/(double)dim;

	m = 0;
	m_max = 0;
	threshold = -1;
	processed = false;
}

BGMGraph::BGMGraph(const double *coordinates_array, int n, int dim){
	/* prepare vector from values */
	Vec vec_arr;
	TRY( VecCreateSeqWithArray(PETSC_COMM_SELF,1,n*dim,coordinates_array,&vec_arr) );

	coordinates = new GeneralVector<PetscVector>(vec_arr);

	this->dim = dim;
	this->n = n;

	m = 0;
	m_max = 0;
	threshold = -1;
	processed = false;
}

BGMGraph::BGMGraph(){
	this->n = 0;
	this->m = 0;
	this->m_max = 0;
	threshold = -1;
	processed = false;
}

BGMGraph::~BGMGraph(){
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

		#ifdef USE_CUDA
			gpuErrchk( cudaFree(neighbor_nmbs_gpu) );
			for(i=0;i<n;i++){
				gpuErrchk( cudaFree(neighbor_ids_gpu[i]) );
			}
			gpuErrchk( cudaFree(neighbor_ids_gpu) );
		#endif

	}
	
}

int BGMGraph::get_n(){
	return this->n;
}

int BGMGraph::get_m(){
	return this->m;
}

int BGMGraph::get_m_max(){
	return this->m_max;
}

int BGMGraph::get_dim(){
	return this->dim;
}

double BGMGraph::get_threshold(){
	return this->threshold;
}

int *BGMGraph::get_neighbor_nmbs(){
	return neighbor_nmbs;
}

int **BGMGraph::get_neighbor_ids(){
	return neighbor_ids;
}

int *BGMGraph::get_neighbor_nmbs_gpu(){
	#ifdef USE_GPU
		return neighbor_nmbs_gpu;
	#else
		return neighbor_nmbs;
	#endif
}

int **BGMGraph::get_neighbor_ids_gpu(){
	#ifdef USE_GPU
		return neighbor_ids_gpu;
	#else
		return neighbor_ids;
	#endif
}

GeneralVector<PetscVector> *BGMGraph::get_coordinates(){
	return coordinates;
}

void BGMGraph::process(double threshold) {
	this->threshold = threshold;
	
	/* prepare array for number of neighbors */
	neighbor_nmbs = (int*)malloc(n*sizeof(int));
	for(int i=0;i<n;i++){
		neighbor_nmbs[i] = 0;
	}
	
	/* get local array and work with it */
	const double *coordinates_arr;
	TRY( VecGetArrayRead(coordinates->get_vector(),&coordinates_arr) );
	
	/* go throught graph - compute number of neighbors */
//	#pragma omp parallel for
	for(int i=0;i<n;i++){
		for(int j=i+1;j<n;j++){
			if(compute_norm(coordinates_arr, i, j) < threshold){
				neighbor_nmbs[i] += 1;
				neighbor_nmbs[j] += 1;
			}
		}
	}

	/* prepare storages for neightbors ids */
	neighbor_ids = (int**)malloc(n*sizeof(int*));
//	#pragma omp parallel for
	for(int i=0;i<n;i++){
		neighbor_ids[i] = (int*)malloc(neighbor_nmbs[i]*sizeof(int));
	}

	/* go throught graph - fill indexes of neighbors */
	int *counters;
	counters = (int*)malloc(n*sizeof(int));

//	#pragma omp parallel for
	for(int i=0;i<n;i++){
		counters[i] = 0;
	}

//	#pragma omp parallel for // TODO: here is a problem, cannot be used, maybe because of couter arrays?
	for(int i=0;i<n;i++){
		for(int j=i+1;j<n;j++){
			if(compute_norm(coordinates_arr, i, j) < threshold){
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
	TRY( VecRestoreArrayRead(coordinates->get_vector(),&coordinates_arr) );

	#ifdef USE_CUDA
		/* copy data to gpu */
		gpuErrchk( cudaMalloc((void **)&neighbor_nmbs_gpu, n*sizeof(int)) );	
		gpuErrchk( cudaMemcpy( neighbor_nmbs_gpu, neighbor_nmbs, n*sizeof(int), cudaMemcpyHostToDevice) );
		
		gpuErrchk( cudaMalloc((void **)&neighbor_ids_gpu, n*sizeof(int)) );	
		for(int i=0;i<n;i++){
			gpuErrchk( cudaMalloc((void **)&(neighbor_ids_gpu[i]), neighbor_nmbs[i]*sizeof(int)) );
			gpuErrchk( cudaMemcpy( neighbor_ids_gpu[i], neighbor_ids[i], n*sizeof(int), cudaMemcpyHostToDevice) );
		}

		gpuErrchk( cudaDeviceSynchronize() );
	#endif
	
	processed = true;
}

void BGMGraph::print(ConsoleOutput &output) const {
	output << "Graph" << std::endl;
	
	output.push();
	output << " - dim:        " << this->dim << std::endl;
	output << " - vertices:   " << this->n << std::endl;
	output << " - edges:      " << this->m << std::endl;
	output << " - max_degree: " << this->m_max << std::endl;
	output << " - threshold:  " << this->threshold << std::endl;
	output << " - processed:  " << this->processed << std::endl;
	output.pop();
	
}

void BGMGraph::print_content(ConsoleOutput &output) const {
	output << "Graph" << std::endl;
	
	output.push();
	output << " - dim:        " << this->dim << std::endl;
	output << " - vertices:   " << this->n << std::endl;
	output << " - edges:      " << this->m << std::endl;
	output << " - max_degree: " << this->m_max << std::endl;
	output << " - threshold:  " << this->threshold << std::endl;
	output << " - processed:  " << this->processed << std::endl;

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


}

#endif
