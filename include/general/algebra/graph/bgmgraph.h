/** @file bgmgraph.h
 *  @brief class for manipulation with graphs
 *
 *  Defines some basic functions for manipulaton with graphs, especially for BLOCKGRAPHMATRIX.
 *
 *  @author Lukas Pospisil
 */
 
#ifndef PASC_COMMON_BGMGRAPH_H
#define	PASC_COMMON_BGMGRAPH_H

/* include metis */
#ifdef USE_METIS
	#include "metis.h"
#endif

#include "general/common/timer.h"
#include "general/algebra/vector/generalvector.h"
#include "general/common/logging.h"

namespace pascinference {
using namespace common;

namespace algebra {

/** \class BGMGraph
 *  \brief General class for manipulation with graphs.
 *
*/
template<class VectorBase>
class BGMGraph {
	public:
		class ExternalContent;
	
	protected:
		friend class ExternalContent;
		ExternalContent *externalcontent;			/**< for manipulation with external-specific stuff */

		GeneralVector<VectorBase> *coordinates; /**< vector with coordinates [p1_x, ... pn_x, p1_y, ... pn_y] */
		int dim; /**< dimension of coordinates */	
		int n; /**< number of nodes */
		int m; /**< number of edges */
		int m_max; /**< maximum degree of vertices */
		double threshold; /**< the value used for processing the graph */
		
		bool processed; /**< there was process with threshold value called? */
		int *neighbor_nmbs; /**< number of neighbors for each node */
		int **neighbor_ids; /**< indexes of neighbors for each node */

		bool DD_decomposed; /**< the decomposition was computed? */
		int DD_size; /**< number of domains for decomposition */
		int *DD_affiliation; /**< domain affiliation of vertices */
		int *DD_permutation; /**< permutation of global indexes Rorig = P(Rnew) */
		int *DD_invpermutation; /**< inverse permutation of global indexes Rnew = invP(Rorig) */
		int *DD_lengths; /**< array of local lengths */
		int *DD_ranges; /**< ranges in permutation array */

		/** @brief compute distance between two vertices
		*
		*  @param values array with coordinates of vertices
		*  @param idx1 index of vertex
		*  @param idx2 index of vertex
		*/
		double compute_normsqr(const double *values, int idx1, int idx2);
		
	public:
		BGMGraph();
		BGMGraph(std::string filename, int dim=2);
		BGMGraph(const double *coordinates_array, int n, int dim);

		~BGMGraph();
		
		/** @brief print the name of graph
		 * 
		 */
		virtual std::string get_name() const;
		
		/** @brief return number of vertices
		*/
		int get_n() const;

		/** @brief return number of edges
		*/
		int get_m() const;

		/** @brief return degree of vertex with max degree 
		*/
		int get_m_max() const;

		/** @brief return dimension of coordinates of vertices
		*/
		int get_dim() const;

		/** @brief return value used for filling graph with edges
		*/
		double get_threshold() const;

		/** @brief return array containing number of neighbors of vertices
		*/
		int *get_neighbor_nmbs() const;

		/** @brief return array of arrays containing indexes of neighbors of vertices
		*/
		int **get_neighbor_ids() const;

		/** @brief return array containing number of neighbors of vertices stored on GPU
		*/
		int *get_neighbor_nmbs_gpu() const;

		/** @brief return array of arrays containing indexes of neighbors of vertices stored on GPU
		*/
		int **get_neighbor_ids_gpu() const;

		/** @brief return vector with coordinates of vertices
		*/
		GeneralVector<VectorBase> *get_coordinates() const;

		bool get_DD_decomposed() const;
		
		/** @brief return number of domains for decomposition
		*/
		int get_DD_size() const;


		int *get_DD_affiliation() const;
		int *get_DD_permutation() const;
		int *get_DD_invpermutation() const;
		int *get_DD_lengths() const;
		int *get_DD_ranges() const;
		
		/** @brief fill graph with edges based on given length of edge
		*/
		virtual void process(double threshold);

		/** @brief decompose graph to given number of domains using METIS
		*/
		virtual void decompose(int nmb_domains);

		/** @brief print basic informations of graph
		*/
		void print(ConsoleOutput &output) const;

		/** @brief print all vertices, edges, neighbors
		*/
		void print_content(ConsoleOutput &output) const;
		
		/** @brief save content of graph to VTK
		*/
		void saveVTK(std::string filename) const;
		
		ExternalContent *get_externalcontent() const;
};


}
} /* end of namespace */


/* --------- IMPLEMENTATION ------ */
namespace pascinference {
namespace algebra {

template<class VectorBase>
double BGMGraph<VectorBase>::compute_normsqr(const double *values, int idx1, int idx2){
	int d;
	double mynorm = 0;
	for(d=0;d<dim;d++){
		mynorm += (values[idx1+d*n] - values[idx2+d*n])*(values[idx1+d*n] - values[idx2+d*n]);
	}
	return mynorm;
}

template<class VectorBase>
BGMGraph<VectorBase>::BGMGraph(std::string filename, int dim){
	coordinates = new GeneralVector<VectorBase>();
	
	/* load nodes from file */
	coordinates->load_local(filename);

	this->dim = dim;
	n = coordinates->size()/(double)dim;

	m = 0;
	m_max = 0;
	threshold = -1;
	processed = false;

	DD_decomposed = false;
}

template<class VectorBase>
BGMGraph<VectorBase>::BGMGraph(const double *coordinates_array, int n, int dim){
	//TODO

}

template<class VectorBase>
BGMGraph<VectorBase>::BGMGraph(){
	this->n = 0;
	this->m = 0;
	this->m_max = 0;
	threshold = -1;
	processed = false;
	DD_decomposed = false;
}

template<class VectorBase>
BGMGraph<VectorBase>::~BGMGraph(){
//	TRYCXX( VecDestroy(&coordinates.get_vector()));
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
			cuda_BGMGraph_destroy(neighbor_nmbs_gpu, neighbor_ids_gpu, neighbor_ids_cpugpu, int n);
		#endif

	}
	
	if(DD_decomposed){
		free(DD_affiliation);
		free(DD_permutation);
		free(DD_invpermutation);
		free(DD_lengths);
		free(DD_ranges);
	}
}

template<class VectorBase>
std::string BGMGraph<VectorBase>::get_name() const {
	return "BGMGraph";
}

template<class VectorBase>
int BGMGraph<VectorBase>::get_n() const {
	return this->n;
}

template<class VectorBase>
int BGMGraph<VectorBase>::get_m() const {
	return this->m;
}

template<class VectorBase>
int BGMGraph<VectorBase>::get_m_max() const {
	return this->m_max;
}

template<class VectorBase>
int BGMGraph<VectorBase>::get_dim() const {
	return this->dim;
}

template<class VectorBase>
double BGMGraph<VectorBase>::get_threshold() const {
	return this->threshold;
}

template<class VectorBase>
int *BGMGraph<VectorBase>::get_neighbor_nmbs() const {
	return this->neighbor_nmbs;
}

template<class VectorBase>
int **BGMGraph<VectorBase>::get_neighbor_ids() const {
	return this->neighbor_ids;
}

template<class VectorBase>
int *BGMGraph<VectorBase>::get_neighbor_nmbs_gpu() const {
	#ifdef USE_CUDA
		return this->neighbor_nmbs_gpu;
	#else
		return this->neighbor_nmbs;
	#endif
}

template<class VectorBase>
int **BGMGraph<VectorBase>::get_neighbor_ids_gpu() const {
	#ifdef USE_CUDA
		return this->neighbor_ids_gpu;
	#else
		return this->neighbor_ids;
	#endif
}

template<class VectorBase>
GeneralVector<VectorBase> *BGMGraph<VectorBase>::get_coordinates() const {
	return this->coordinates;
}

template<class VectorBase>
bool BGMGraph<VectorBase>::get_DD_decomposed() const {
	return this->DD_decomposed;
}

template<class VectorBase>
int BGMGraph<VectorBase>::get_DD_size() const {
	return this->DD_size;
}

template<class VectorBase>
int *BGMGraph<VectorBase>::get_DD_affiliation() const {
	return this->DD_affiliation;
}

template<class VectorBase>
int *BGMGraph<VectorBase>::get_DD_permutation() const {
	return this->DD_permutation;
}

template<class VectorBase>
int *BGMGraph<VectorBase>::get_DD_invpermutation() const {
	return this->DD_invpermutation;
}

template<class VectorBase>
int *BGMGraph<VectorBase>::get_DD_lengths() const {
	return this->DD_lengths;
}

template<class VectorBase>
int *BGMGraph<VectorBase>::get_DD_ranges() const {
	return this->DD_ranges;
}

template<class VectorBase>
void BGMGraph<VectorBase>::decompose(int nmb_domains){
	LOG_FUNC_BEGIN
	
	if(!this->DD_decomposed){ /* there wasn't decompose called yet */
		this->DD_decomposed = true;
		this->DD_size = nmb_domains;

		/* allocate arrays */
		DD_affiliation = (int*)malloc(n*sizeof(int));
		DD_permutation = (int*)malloc(n*sizeof(int));
		DD_invpermutation = (int*)malloc(n*sizeof(int));
		DD_lengths = (int*)malloc(DD_size*sizeof(int));
		DD_ranges = (int*)malloc((DD_size+1)*sizeof(int));

		if(nmb_domains > 1){
			/* ---- METIS STUFF ---- */
			int *xadj; /* Indexes of starting points in adjacent array */
			xadj = (int*)malloc((n+1)*sizeof(int));
		
			int *adjncy; /* Adjacent vertices in consecutive index order */
			adjncy = (int*)malloc(2*m*sizeof(int));

			/* fill aux metis stuff */
			int counter = 0;
			for(int i=0;i<n;i++){
				xadj[i] = counter;
				for(int j=0;j<neighbor_nmbs[i];j++){
					adjncy[counter] = neighbor_ids[i][j];
					counter++;
				}
			}	
			xadj[n] = counter;

			int objval;
			int nWeights = 1; /* something with weights of graph, I really don't know, sorry */

			/* run decomposition */
			int metis_ret = METIS_PartGraphKway(&n,&nWeights, xadj, adjncy,
							   NULL, NULL, NULL, &DD_size, NULL,
							   NULL, NULL, &objval, DD_affiliation);

			/* free aux stuff */
			free(xadj);
			free(adjncy);
			/* --------------------- */

			/* compute local lengths and permutation of global indexes */
			for(int i=0;i<DD_size;i++){ /* use DD_length as counters */
				DD_lengths[i] = 0;
			}
			for(int i=0;i<n;i++){
				DD_permutation[i] = DD_lengths[DD_affiliation[i]]; /* set index as a value of counter */
				DD_lengths[DD_affiliation[i]] += 1;
			}

			/* compute ranges */
			DD_ranges[0] = 0;
			for(int i=0;i<DD_size;i++){
				DD_ranges[i+1] = DD_ranges[i] + DD_lengths[i];
			}
			for(int i=0;i<n;i++){
				DD_permutation[i] += DD_ranges[DD_affiliation[i]]; /* shift local indexes to global */

				/* compute inverse permutation */
				DD_invpermutation[DD_permutation[i]] = i;
			}

			
		} else {
			
			/* nmb_domains <= 1 */
			/* fill arrays manually, METIS is not able to do it with number of domains equal to 1 */
			for(int i=0;i<n;i++){
				DD_affiliation[i] = 0;
			}

			for(int i=0;i<n;i++){
				DD_permutation[i] = i;
				DD_invpermutation[i] = i;
			}

			DD_lengths[0] = n;

			DD_ranges[0] = 0;
			DD_ranges[1] = n;
		}
	
	} else {
		// TODO: give error that decompose was already called, or clean stuff and make it again?
	}
	
	LOG_FUNC_END
}

template<class VectorBase>
void BGMGraph<VectorBase>::process(double threshold) {
	LOG_FUNC_BEGIN
	
	//TODO

	LOG_FUNC_END
}

template<class VectorBase>
void BGMGraph<VectorBase>::print(ConsoleOutput &output) const {
	output << this->get_name() << std::endl;
	
	output.push();
	output << " - dim:        " << this->dim << std::endl;
	output << " - vertices:   " << this->n << std::endl;
	output << " - edges:      " << this->m << std::endl;
	output << " - max_degree: " << this->m_max << std::endl;
	output << " - threshold:  " << this->threshold << std::endl;
	output << " - processed:  " << this->processed << std::endl;

	output << " - decomposed: " << this->DD_decomposed << std::endl;
	output.push();
	if(DD_decomposed){
		output << " - nmb of domains: " << this->DD_size << std::endl;
		output << " - lengths:        ";
		for(int i=0;i<DD_size;i++){ /* use DD_length as counters */
			output << DD_lengths[i];
			if(i < DD_size-1){
				output << ", ";
			}
		}
		output << std::endl;
	}
	output.pop();
	
	output.pop();
	
}

template<class VectorBase>
void BGMGraph<VectorBase>::print_content(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN
	
	output << "Graph" << std::endl;
	
	output.push();
	output << " - dim         : " << this->dim << std::endl;
	output << " - vertices    : " << this->n << std::endl;
	output << " - edges       : " << this->m << std::endl;
	output << " - max_degree  : " << this->m_max << std::endl;
	output << " - threshold   : " << this->threshold << std::endl;
	output << " - coordinates : " << *coordinates << std::endl;

	output << " - decomposed  : " << this->DD_decomposed << std::endl;
	output.push();
	if(DD_decomposed){
		output << " - DD_size: " << this->DD_size << std::endl;
		output << " - DD_lengths:        ";
		for(int i=0;i<DD_size;i++){
			output << DD_lengths[i];
			if(i < DD_size-1){
				output << ", ";
			}
		}
		output << std::endl;	
		output << " - DD_ranges:        ";
		for(int i=0;i<DD_size+1;i++){
			output << DD_ranges[i];
			if(i < DD_size){
				output << ", ";
			}
		}
		output << std::endl;	
		output << " - DD_affiliation: ";
		for(int i=0;i<n;i++){
			output << DD_affiliation[i];
			if(i < n-1){
				output << ", ";
			}
		}
		output << std::endl;	
		output << " - DD_permutation: ";
		for(int i=0;i<n;i++){
			output << DD_permutation[i];
			if(i < n-1){
				output << ", ";
			}
		}
		output << std::endl;	
		output << " - DD_invpermutation: ";
		for(int i=0;i<n;i++){
			output << DD_invpermutation[i];
			if(i < n-1){
				output << ", ";
			}
		}
		output << std::endl;	
			
	}
	output.pop();

	output << " - processed:  " << this->processed << std::endl;
	if(this->processed){
		output << " - arrays of neighbors: " << std::endl;
		output.push();
		for(int i=0;i<n;i++){
			output << i << ": " << "(" << neighbor_nmbs[i] << "): ";
			for(int j=0;j<neighbor_nmbs[i];j++){
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
	
	LOG_FUNC_END
}

template<class VectorBase>
void BGMGraph<VectorBase>::saveVTK(std::string filename) const {
	LOG_FUNC_BEGIN

	//TODO
	
	LOG_FUNC_END
}


}
} /* end of namespace */

#endif
