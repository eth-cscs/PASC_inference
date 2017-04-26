/** @file decomposition.h
 *  @brief class for storing and manipulation with decomposition informations
 *
 *  @author Lukas Pospisil
 */
 
#ifndef PASC_COMMON_DECOMPOSITION_H
#define	PASC_COMMON_DECOMPOSITION_H

#include "general/algebra/graph/bgmgraph.h"

namespace pascinference {
namespace algebra {

/** \class Decomposition
 *  \brief Manipulation with problem layout.
 *
*/
template<class VectorBase>
class Decomposition {
	protected:
		/* T decomposition */
		int T; /**< time - global length */
		int DDT_size; /**< time - number of domains for decomposition */
		int *DDT_ranges; /**< time - local ownerships */
		bool destroy_DDT_arrays; /**< if decomposition of T was not provided, I have to create and destroy arrays here */

		/* R decomposition */
		BGMGraph<VectorBase> *graph; /**< graph with stucture of the matrix */

		int R; /**< space - global length */
		int DDR_size; /**< space - number of domains for decomposition */
		int *DDR_affiliation; /**< space - domain affiliation of vertices */
		int *DDR_permutation; /**< space - permutation of global indexes */
		int *DDR_invpermutation; /**< space - permutation of global indexes */
		int *DDR_lengths; /**< space - array of local lengths */
		int *DDR_ranges; /**< space - local ownership */
		bool destroy_DDR_arrays; /**< if graph was not provided, I have to create and destroy arrays here */
	
		/* other stuff */
		int K; /**< number of clusters */
		int xdim; /**< data dimension */
	
		/* my coordinates */
		int DDT_rank; /**< my coordinate in time decomposition */
		int DDR_rank; /**< my coordinate in space decomposition */

		/** @brief compute coordinates in decomposition based on rank of the processor
		*/		
		void compute_rank();
		
	public:
		/** @brief decomposition only in time
		*/
		Decomposition(int T, int R, int K, int xdim, int DDT_size);

		/** @brief decomposition only in space
		*/
		Decomposition(int T, BGMGraph<VectorBase> &new_graph, int K, int xdim, int DDR_size);

		/** @brief decomposition in time and space
		 * 
		 * @todo has to be tested, probably is not working
		*/
		Decomposition(int T, BGMGraph<VectorBase> &new_graph, int K, int xdim, int DDT_size, int DDR_size);

		/** @brief destructor
		*/
		~Decomposition();
		
		/** @brief get global time length
		 * 
		 * @return global length of time
		 */
		int get_T() const;
		
		/** @brief get local portion of time
		 * 
		 * @return local number of time steps
		 */
		int get_Tlocal() const;

		/** @brief get first local index of time in global scope
		 * 
		 * @return starting index of local portion of time steps
		 */
		int get_Tbegin() const;

		/** @brief get last local index (+1) of time in global scope
		 *
		 * @return last index of local portion of time steps plus one
		 */
		int get_Tend() const;

		/** @brief get number of parts of time decomposition
		 * 
		 * @return number of parts of time decomposition
		 */
		int get_DDT_size() const;

		/** @brief get my coordinate in time decomposition
		 * 
		 * @return the coortinate of this MPI process in time decomposition
		 */
		int get_DDT_rank() const;
		
		/** @brief get the array with first indexes of local portion of time
		 * 
		 * @return the array (of size DDT_size+1) with global starting indexes
		 */
		const int *get_DDT_ranges() const;

		/** @brief get global space length
		 * 
		 * @return number of nodes in spatial graph
		 */
		int get_R() const;

		/** @brief get first local index of space in global scope
		 * 
		 * @return starting index of local portion of graph nodes
		 */
		int get_Rbegin() const;

		/** @brief get last local index (+1) of space decomposition in global scope
		 *
		 * @return last index of local portion of graph nodes plus one
		 */
		int get_Rend() const;

		/** @brief get local portion of decomposed spatial graph
		 * 
		 * @return local number of graph nodes
		 */
		int get_Rlocal() const;

		/** @brief get number of parts of space decomposition
		 * 
		 * @return the number of parts of graph decomposition
		 */
		int get_DDR_size() const;

		/** @brief get coordinate of this MPI process in space decomposition
		 * 
		 * @return coordinate of local part in graph decomposition
		 */
		int get_DDR_rank() const;


		int *get_DDR_affiliation() const;
		int *get_DDR_permutation() const;
		int *get_DDR_invpermutation() const;
		int *get_DDR_lengths() const;
		int *get_DDR_ranges() const;
		
		/** @brief get graph of space decomposition
		 * 
		 * @return pointer to graph
		 */
		BGMGraph<VectorBase> *get_graph() const;
		
		/** @brief set graph of space decomposition
		 * 
		 * @param graph the new graph of decomposition
		 * @param DDR_size number of parts of graph decomposition
		 */
		void set_graph(BGMGraph<VectorBase> &graph, int DDR_size=1);

		/** @brief get number of clusters
		 * 
		 * @return number of clusters
		 */
		int get_K() const;

		/** @brief get data dimension (number of components)
		 * 
		 * @return number of data components
		 */
		int get_xdim() const;

		/** @brief print informations of decomposition
		*/
		void print(ConsoleOutput &output) const;

		/** @brief print informations of decomposition
		*/
		void print_content(ConsoleOutput &output_master, ConsoleOutput &output_local, bool print_details=true) const;

		/** @brief create global PETSc gamma vector with respect to this decomposition
		 * 
		 * @param x_Vec pointer to new gamma vector
		*/
		void createGlobalVec_gamma(Vec *x_Vec) const;

		/** @brief create global PETSc data vector with respect to this decomposition
		 * 
		 * @param x_Vec pointer to new vector
		*/
		void createGlobalVec_data(Vec *x_Vec) const;

		/** @brief get local index in gamma vector from global indexes
		 * 
		 * @param t_global global time index
		 * @param r_global global space index
		 * @param k index of cluster
		 * @return local index in gamma vector
		 */
		int get_idxglobal(int t_global, int r_global, int k) const;

		/** @brief get the index of node in original graph from index in permutated graph
		 * 
		 * @param r_global global node index in permutated graph
		 * @return node index in original graph
		 * @todo has to be tested
		 */
		int get_invPr(int r_global) const;

		/** @brief get the index of node in permutated graph from index in original graph
		 * 
		 * @param r_global global node index in original graph
		 * @return node index in permutated graph
		 * @todo has to be tested
		 */
		int get_Pr(int r_global) const;

#ifdef USE_PETSC
		void permute_TRxdim(Vec orig_Vec, Vec new_Vec, bool invert=false) const;
		void permute_TRK(Vec orig_Vec, Vec new_Vec, bool invert=false) const;
		void permute_TRblocksize(Vec orig_Vec, Vec new_Vec, int blocksize, bool invert) const;
		
		/** @brief create PETSC index set with local gamma indexes which correspond to given cluster index
		 * 
		 * @param is pointer to new index set
		 * @param k index of cluster
		 */
		void createIS_gammaK(IS *is, int k) const;
#endif
};

/* ----------------- Decomposition implementation ------------- */
template<class VectorBase>
Decomposition<VectorBase>::Decomposition(int T, int R, int K, int xdim, int DDT_size){
	LOG_FUNC_BEGIN

	//TODO: write something for general case

	LOG_FUNC_END
}

template<class VectorBase>
Decomposition<VectorBase>::Decomposition(int T, BGMGraph<VectorBase> &new_graph, int K, int xdim, int DDT_size, int DDR_size){
	LOG_FUNC_BEGIN

	//TODO: write something for general case

	LOG_FUNC_END
}

template<class VectorBase>
Decomposition<VectorBase>::~Decomposition(){
	LOG_FUNC_BEGIN

	//TODO: write something for general case

	LOG_FUNC_END
}

template<class VectorBase>
void Decomposition<VectorBase>::compute_rank(){
	LOG_FUNC_BEGIN

	//TODO: write something for general case
	
	LOG_FUNC_END
}

template<class VectorBase>
int Decomposition<VectorBase>::get_T() const{
	return T;
}

template<class VectorBase>
int Decomposition<VectorBase>::get_Tlocal() const{
	return DDT_ranges[DDT_rank+1] - DDT_ranges[DDT_rank];
}

template<class VectorBase>
int Decomposition<VectorBase>::get_Tbegin() const{
	return DDT_ranges[DDT_rank];
}

template<class VectorBase>
int Decomposition<VectorBase>::get_Tend() const{
	return DDT_ranges[DDT_rank+1];
}

template<class VectorBase>
int Decomposition<VectorBase>::get_DDT_size() const{
	return DDT_size;
}

template<class VectorBase>
int Decomposition<VectorBase>::get_DDT_rank() const{
	return DDT_rank;
}

template<class VectorBase>
const int *Decomposition<VectorBase>::get_DDT_ranges() const{
	return DDT_ranges;
}

template<class VectorBase>
int Decomposition<VectorBase>::get_R() const{
	return R;
}

template<class VectorBase>
int Decomposition<VectorBase>::get_Rlocal() const{
	return DDR_ranges[DDR_rank+1] - DDR_ranges[DDR_rank];
}

template<class VectorBase>
int Decomposition<VectorBase>::get_Rbegin() const{
	return DDR_ranges[DDR_rank];
}

template<class VectorBase>
int Decomposition<VectorBase>::get_Rend() const{
	int rank = GlobalManager.get_rank();
	return DDR_ranges[DDR_rank+1];
}

template<class VectorBase>
int Decomposition<VectorBase>::get_DDR_size() const{
	return DDR_size;
}

template<class VectorBase>
int Decomposition<VectorBase>::get_DDR_rank() const{
	return DDR_rank;
}

template<class VectorBase>
int *Decomposition<VectorBase>::get_DDR_affiliation() const{
	return DDR_affiliation;
}

template<class VectorBase>
int *Decomposition<VectorBase>::get_DDR_permutation() const{
	return DDR_permutation;
}

template<class VectorBase>
int *Decomposition<VectorBase>::get_DDR_invpermutation() const{
	return DDR_invpermutation;
}

template<class VectorBase>
int *Decomposition<VectorBase>::get_DDR_lengths() const{
	return DDR_lengths;
}

template<class VectorBase>
int *Decomposition<VectorBase>::get_DDR_ranges() const{
	return DDR_ranges;
}

template<class VectorBase>
BGMGraph<VectorBase> *Decomposition<VectorBase>::get_graph() const{
	return graph;
}

template<class VectorBase>
void Decomposition<VectorBase>::set_graph(BGMGraph<VectorBase> &new_graph, int DDR_size) {

}

template<class VectorBase>
int Decomposition<VectorBase>::get_K() const{
	return K;
}

template<class VectorBase>
int Decomposition<VectorBase>::get_xdim() const{
	return xdim;
}

template<class VectorBase>
void Decomposition<VectorBase>::print_content(ConsoleOutput &output_master, ConsoleOutput &output_local, bool print_details) const {
	output_master << "Decomposition" << std::endl;
	
	output_master.push();
	output_master << " K                     : " << this->K << std::endl;
	output_master << " Data dimension        : " << this->xdim << std::endl;
	output_master << " Time                  : " << this->T << std::endl;
	output_master << " - DDT_size            : " << this->DDT_size << std::endl;
	output_master << " - DDT_ranges          : ";
	for(int i=0;i< this->DDT_size+1;i++){
		output_master << this->DDT_ranges[i];
		if(i< this->DDT_size){
			output_master << ", ";
		}
	}
	output_master << std::endl;
	output_master.pop();

	output_master.push();
	output_master << " Space                 : " << this->R << std::endl;
	output_master << " - DDR_size            : " << this->DDR_size << std::endl;
	if(print_details){
		output_master << " - DDR_affiliation     : ";
		for(int i=0;i< this->R;i++){
			output_master << this->DDR_affiliation[i];
			if(i< this->R-1){
				output_master << ", ";
			}
		}
		output_master << std::endl;
		output_master << " - DDR_permutation     : ";
		for(int i=0;i< this->R;i++){
			output_master << this->DDR_permutation[i];
			if(i< this->R-1){
				output_master << ", ";
			}
		}
		output_master << std::endl;
	}
	output_master << " - DDR_ranges          : ";
	for(int i=0;i< this->DDR_size+1;i++){
		output_master << this->DDR_ranges[i];
		if(i< this->DDR_size){
			output_master << ", ";
		}
	}
	output_master << std::endl;
	output_master << " - DDR_lengths         : ";
	for(int i=0;i< this->DDR_size;i++){
		output_master << this->DDR_lengths[i];
		if(i< this->DDR_size-1){
			output_master << ", ";
		}
	}
	output_master << std::endl;
	output_master.pop();

	output_master.push();
	output_master << " - coordinates [T,R]: " << std::endl;
	output_local  << "   [ " << this->DDT_rank << ", " << this->DDR_rank << " ]" << std::endl;
	output_local.synchronize();
	output_master.pop();
	
}

template<class VectorBase>
void Decomposition<VectorBase>::print(ConsoleOutput &output) const {
	output << "Decomposition" << std::endl;
	
	output.push();
	output << " Clusters              : " << this->K << std::endl;
	output << " Data dimension        : " << this->xdim << std::endl;
	output << " Time                  : " << this->T << std::endl;
	output << " - DDT_size            : " << this->DDT_size << std::endl;
	output << " - DDT_ranges          : ";
	for(int i=0;i< this->DDT_size+1;i++){
		output << this->DDT_ranges[i];
		if(i< this->DDT_size){
			output << ", ";
		}
	}
	output << std::endl;
	output.pop();

	output.push();
	output << " Space                 : " << this->R << std::endl;
	output << " - DDR_size            : " << this->DDR_size << std::endl;
	output << " - DDR_ranges          : ";
	for(int i=0;i< this->DDR_size+1;i++){
		output << this->DDR_ranges[i];
		if(i< this->DDR_size){
			output << ", ";
		}
	}
	output << std::endl;
	output.pop();
	
}

template<class VectorBase>
int Decomposition<VectorBase>::get_idxglobal(int t_global, int r_global, int k) const {
	int Pr = DDR_permutation[r_global];
	return t_global*R*K + Pr*K + k;
}

template<class VectorBase>
int Decomposition<VectorBase>::get_invPr(int r_global) const {
	return DDR_invpermutation[r_global];
}

template<class VectorBase>
int Decomposition<VectorBase>::get_Pr(int r_global) const {
	return DDR_permutation[r_global];
}


}
} /* end of namespace */


#endif
