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

#include "common/timer.h"
#include "algebra/vector/generalvector.h"

namespace pascinference {
using namespace common;

namespace algebra {

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

		bool DD_decomposed; /**< the decomposition was computed? */
		int DD_size; /**< number of domains for decomposition */
		int *DD_affiliation; /**< domain affiliation of vertices */
		int *DD_permutation; /**< permutation of global indexes Rorig = P(Rnew) */
		int *DD_invpermutation; /**< inverse permutation of global indexes Rnew = invP(Rorig) */
		int *DD_lengths; /**< array of local lengths */
		int *DD_ranges; /**< ranges in permutation array */
		
		#ifdef USE_CUDA
			int *neighbor_nmbs_gpu; /**< copy of values on GPU */
			int **neighbor_ids_cpugpu; /**< pointers to GPU arrays on CPU */
			int **neighbor_ids_gpu; /**< copy of values on GPU */
		#endif

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
		GeneralVector<PetscVector> *get_coordinates() const;

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
};


}
} /* end of namespace */

#endif
