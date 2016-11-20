/** @file decomposition.h
 *  @brief class for storing and manipulation with decomposition informations
 *
 *  @author Lukas Pospisil
 */
 
#ifndef PASC_COMMON_DECOMPOSITION_H
#define	PASC_COMMON_DECOMPOSITION_H

#include "algebra/bgmgraph.h"

/* this class is for petscvector */
typedef petscvector::PetscVector PetscVector;

namespace pascinference {
namespace algebra {

/** \class Decomposition
 *  \brief Manipulation with problem layout.
 *
*/
class Decomposition {
	protected:
		/* T decomposition */
		int T; /**< time - global length */
		int DDT_size; /**< time - number of domains for decomposition */
		int *DDT_ranges; /**< time - local ownerships */
		bool destroy_DDT_arrays; /**< if decomposition of T was not provided, I have to create and destroy arrays here */

		/* R decomposition */
		BGMGraph *graph; /**< graph with stucture of the matrix */

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
		Decomposition(int T, BGMGraph &new_graph, int K, int xdim, int DDR_size);

		/** @brief decomposition in time and space
		 * 
		 * @todo has to be tested, probably is not working
		*/
		Decomposition(int T, BGMGraph &new_graph, int K, int xdim, int DDT_size, int DDR_size);

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
		
		/** @
		BGMGraph *get_graph() const;
		
		/** @brief set graph of space decomposition
		 * 
		 * @param graph the new graph of decomposition
		 * @param DDR_size number of parts of graph decomposition
		 */
		void set_graph(BGMGraph &graph, int DDR_size=1);

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

		/** @brief create global PETSc gamma vector on GPU with respect to this decomposition
		 * 
		 * VECMPICUDA
		 * 
		 * @param x_Vec pointer to new gamma vector
		*/
		void createGlobalCudaVec_gamma(Vec *x_Vec) const;

		/** @brief create global PETSc data vector with respect to this decomposition
		 * 
		 * @param x_Vec pointer to new vector
		*/
		void createGlobalVec_data(Vec *x_Vec) const;

		/** @brief create global PETSc data vector on GPU with respect to this decomposition
		 * 
		 * VECMPICUDA
		 * 
		 * @param x_Vec pointer to new vector
		*/
		void createGlobalCudaVec_data(Vec *x_Vec) const;

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

		void permute_TRxdim(Vec orig_Vec, Vec new_Vec, bool invert=false) const;
		void permute_TRK(Vec orig_Vec, Vec new_Vec, bool invert=false) const;
		void permute_TRblocksize(Vec orig_Vec, Vec new_Vec, int blocksize, bool invert) const;
		
		/** @brief create PETSC index set with local gamma indexes which correspond to given cluster index
		 * 
		 * @param is pointer to new index set
		 * @param k index of cluster
		 */
		void createIS_gammaK(IS *is, int k) const;
};

/* ----------------- Decomposition implementation ------------- */

Decomposition::Decomposition(int T, int R, int K, int xdim, int DDT_size){
	LOG_FUNC_BEGIN

	this->T = T;
	this->R = R;
	this->K = K;
	this->xdim = xdim;
		
	this->DDT_size = DDT_size;
	this->DDR_size = 1;
	
	/* prepare new layout for T */
	destroy_DDT_arrays = true;	
	DDT_ranges = (int *)malloc((this->DDT_size+1)*sizeof(int));
	
	Vec DDT_layout;
	TRYCXX( VecCreate(PETSC_COMM_WORLD,&DDT_layout) );
	TRYCXX( VecSetSizes(DDT_layout, PETSC_DECIDE, T ));
	TRYCXX( VecSetFromOptions(DDT_layout) );
	
	/* get properties of this layout */
	const int *DDT_ranges_const;
	TRYCXX( VecGetOwnershipRanges(DDT_layout, &DDT_ranges_const) );
	for(int i=0;i<DDT_size+1;i++){
		DDT_ranges[i] = DDT_ranges_const[i]; // TODO: how to deal with const int in ranges form PETSc?
	}
	
	/* destroy temp vector */
	TRYCXX( VecDestroy(&DDT_layout) );

	/* prepare new layout for R */
	/* no graph provided - allocate arrays */
	destroy_DDR_arrays = true;

	DDR_affiliation = (int *)malloc(R*sizeof(int));
	DDR_permutation = (int *)malloc(R*sizeof(int));
	DDR_invpermutation = (int *)malloc(R*sizeof(int));
	for(int i=0;i<R;i++){
		DDR_affiliation[i] = 0;
		DDR_permutation[i] = i;
		DDR_invpermutation[i] = i;
	}

	DDR_lengths = (int *)malloc((this->DDR_size)*sizeof(int));
	DDR_lengths[0] = R;

	DDR_ranges = (int *)malloc((this->DDR_size+1)*sizeof(int));
	DDR_ranges[0] = 0;
	DDR_ranges[1] = R;

	/* no graph provided */
	graph = NULL;

	compute_rank();

	LOG_FUNC_END
}

Decomposition::Decomposition(int T, BGMGraph &new_graph, int K, int xdim, int DDR_size){
	LOG_FUNC_BEGIN

	this->T = T;
	this->R = new_graph.get_n();
	this->K = K;
	this->xdim = xdim;

	/* prepare new layout for R */
	set_graph(new_graph,DDR_size);

	this->DDT_size = 1;
	this->DDR_size = new_graph.get_DD_size();

	/* prepare new layout for T */
	destroy_DDT_arrays = true;	
	DDT_ranges = (int *)malloc((this->DDT_size+1)*sizeof(int));
	DDT_ranges[0] = 0;
	DDT_ranges[1] = T;

	compute_rank();

	LOG_FUNC_END
}

Decomposition::Decomposition(int T, BGMGraph &new_graph, int K, int xdim, int DDT_size, int DDR_size){
	LOG_FUNC_BEGIN

	this->T = T;
	this->R = new_graph.get_n();
	this->K = K;
	this->xdim = xdim;

	/* prepare new layout for R */
	destroy_DDR_arrays = false;
	set_graph(new_graph, DDR_size);

	this->DDT_size = DDT_size;
	this->DDR_size = graph->get_DD_size();
	
	/* prepare new layout for T */
	destroy_DDT_arrays = true;	
	/* unfortunatelly, we have to compute distribution of T manually */
	int DDT_optimal_local_size = T/(double)DDT_size;
	int DDT_optimal_local_size_residue = T - DDT_optimal_local_size*DDT_size;
	DDT_ranges = (int *)malloc((this->DDT_size+1)*sizeof(int));
	DDT_ranges[0] = 0;
	for(int i=0;i<DDT_size;i++){
		DDT_ranges[i+1] = DDT_ranges[i] + DDT_optimal_local_size;
		if(i < DDT_optimal_local_size_residue){
			DDT_ranges[i+1] += 1;
		}
	}

	compute_rank();

	LOG_FUNC_END
}

Decomposition::~Decomposition(){
	LOG_FUNC_BEGIN

	if(destroy_DDT_arrays){
		free(DDT_ranges);
	}

	if(destroy_DDR_arrays){
		free(DDR_affiliation);
		free(DDR_permutation);
		free(DDR_invpermutation);
		free(DDR_lengths);
		free(DDR_ranges);
	}

	LOG_FUNC_END
}

void Decomposition::compute_rank(){
	LOG_FUNC_BEGIN

	/* get rank of this processor */
	int rank = GlobalManager.get_rank();

	this->DDT_rank = rank/(double)this->DDR_size;
	this->DDR_rank = rank - (this->DDT_rank)*(this->DDR_size);

	/* control the decomposition */
//	if(this->DDT_size*this->DDR_size != GlobalManager.get_size()){
//		coutMaster << "Sorry, DDT_size*DDR_size != nproc" << std::endl;
//		coutMaster << " DDT_size = " << this->DDT_size << std::endl;
//		coutMaster << " DDR_size = " << this->DDR_size << std::endl;
//		coutMaster << " nproc    = " << GlobalManager.get_size() << std::endl;

		// TODO: throw error
//	}
	
	LOG_FUNC_END
}

int Decomposition::get_T() const{
	return T;
}

int Decomposition::get_Tlocal() const{
	return DDT_ranges[DDT_rank+1] - DDT_ranges[DDT_rank];
}

int Decomposition::get_Tbegin() const{
	return DDT_ranges[DDT_rank];
}

int Decomposition::get_Tend() const{
	return DDT_ranges[DDT_rank+1];
}

int Decomposition::get_DDT_size() const{
	return DDT_size;
}

int Decomposition::get_DDT_rank() const{
	return DDT_rank;
}

const int *Decomposition::get_DDT_ranges() const{
	return DDT_ranges;
}

int Decomposition::get_R() const{
	return R;
}

int Decomposition::get_Rlocal() const{
	return DDR_ranges[DDR_rank+1] - DDR_ranges[DDR_rank];
}

int Decomposition::get_Rbegin() const{
	return DDR_ranges[DDR_rank];
}

int Decomposition::get_Rend() const{
	int rank = GlobalManager.get_rank();
	return DDR_ranges[DDR_rank+1];
}

int Decomposition::get_DDR_size() const{
	return DDR_size;
}

int Decomposition::get_DDR_rank() const{
	return DDR_rank;
}

int *Decomposition::get_DDR_affiliation() const{
	return DDR_affiliation;
}

int *Decomposition::get_DDR_permutation() const{
	return DDR_permutation;
}

int *Decomposition::get_DDR_invpermutation() const{
	return DDR_invpermutation;
}

int *Decomposition::get_DDR_lengths() const{
	return DDR_lengths;
}

int *Decomposition::get_DDR_ranges() const{
	return DDR_ranges;
}

BGMGraph *Decomposition::get_graph() const{
	return graph;
}

void Decomposition::set_graph(BGMGraph &new_graph, int DDR_size) {

	if(destroy_DDR_arrays){
		free(DDR_affiliation);
		free(DDR_permutation);
		free(DDR_invpermutation);
		free(DDR_lengths);
		free(DDR_ranges);
	}

	/* decompose graph */
	this->DDR_size = DDR_size;
	new_graph.decompose(DDR_size);

	this->graph = &new_graph;
	destroy_DDR_arrays = false;
	DDR_affiliation = new_graph.get_DD_affiliation();
	DDR_permutation = new_graph.get_DD_permutation();
	DDR_invpermutation = new_graph.get_DD_invpermutation();
	DDR_lengths = new_graph.get_DD_lengths();
	DDR_ranges = new_graph.get_DD_ranges();	
	
	compute_rank();
}


int Decomposition::get_K() const{
	return K;
}

int Decomposition::get_xdim() const{
	return xdim;
}

void Decomposition::print_content(ConsoleOutput &output_master, ConsoleOutput &output_local, bool print_details) const {
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

void Decomposition::print(ConsoleOutput &output) const {
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

void Decomposition::createGlobalVec_gamma(Vec *x_Vec) const {
	LOG_FUNC_BEGIN

	int T = this->get_T();
	int R = this->get_R();
	int K = this->get_K();
	int Tlocal = this->get_Tlocal();
	int Rlocal = this->get_Rlocal();

	TRYCXX( VecCreate(PETSC_COMM_WORLD,x_Vec) );
	TRYCXX( VecSetSizes(*x_Vec,Tlocal*Rlocal*K,T*R*K) );
	TRYCXX( VecSetFromOptions(*x_Vec) );

	LOG_FUNC_END
}

void Decomposition::createGlobalCudaVec_gamma(Vec *x_Vec) const { //TODO: how about call it with GeneralVector<PetscVector> ?
	LOG_FUNC_BEGIN

	int T = this->get_T();
	int R = this->get_R();
	int K = this->get_K();
	int Tlocal = this->get_Tlocal();
	int Rlocal = this->get_Rlocal();

	TRYCXX( VecCreate(PETSC_COMM_WORLD,x_Vec) );
	TRYCXX( VecSetType(*x_Vec, VECMPICUDA) );
	TRYCXX( VecSetSizes(*x_Vec,Tlocal*Rlocal*K,T*R*K) );
	TRYCXX( VecSetFromOptions(*x_Vec) );

	LOG_FUNC_END
}

void Decomposition::createGlobalVec_data(Vec *x_Vec) const { //TODO: how about call it with GeneralVector<PetscVector> ?
	LOG_FUNC_BEGIN

	int T = this->get_T();
	int R = this->get_R();
	int xdim = this->get_xdim();
	int Tlocal = this->get_Tlocal();
	int Rlocal = this->get_Rlocal();

	TRYCXX( VecCreate(PETSC_COMM_WORLD,x_Vec) );
	TRYCXX( VecSetSizes(*x_Vec,Tlocal*Rlocal*xdim,T*R*xdim) );
	TRYCXX( VecSetFromOptions(*x_Vec) );

	LOG_FUNC_END
}

void Decomposition::createGlobalCudaVec_data(Vec *x_Vec) const { //TODO: how about call it with GeneralVector<PetscVector> ?
	LOG_FUNC_BEGIN

	int T = this->get_T();
	int R = this->get_R();
	int xdim = this->get_xdim();
	int Tlocal = this->get_Tlocal();
	int Rlocal = this->get_Rlocal();

	TRYCXX( VecCreate(PETSC_COMM_WORLD,x_Vec) );
	TRYCXX( VecSetType(*x_Vec, VECMPICUDA) );
	TRYCXX( VecSetSizes(*x_Vec,Tlocal*Rlocal*xdim,T*R*xdim) );
	TRYCXX( VecSetFromOptions(*x_Vec) );

//	TRYCXX( VecAssemblyBegin(*x_Vec));
//	TRYCXX( VecAssemblyEnd(*x_Vec));

	LOG_FUNC_END
}

int Decomposition::get_idxglobal(int t_global, int r_global, int k) const {
	int Pr = DDR_permutation[r_global];
	return t_global*R*K + Pr*K + k;
}

int Decomposition::get_invPr(int r_global) const {
	return DDR_invpermutation[r_global];
}

int Decomposition::get_Pr(int r_global) const {
	return DDR_permutation[r_global];
}

void Decomposition::permute_TRxdim(Vec orig_Vec, Vec new_Vec, bool invert) const {
	permute_TRblocksize(orig_Vec, new_Vec, xdim, invert);
}

void Decomposition::permute_TRK(Vec orig_Vec, Vec new_Vec, bool invert) const {
	permute_TRblocksize(orig_Vec, new_Vec, K, invert);
}

void Decomposition::permute_TRblocksize(Vec orig_Vec, Vec new_Vec, int blocksize, bool invert) const {
	LOG_FUNC_BEGIN

	int Tbegin = get_Tbegin();
	int Tend = get_Tend();
	int Tlocal = get_Tlocal();
	int Rbegin = get_Rbegin();
	int Rlocal = get_Rlocal();
	
	int local_size = Tlocal*Rlocal*blocksize;

	IS orig_local_is;
	IS new_local_is;

	Vec orig_local_Vec;
	Vec new_local_Vec;

	/* prepare index set with local data */
	int orig_local_arr[local_size];
	int j = 0;
	for(int t=Tbegin;t<Tend;t++){
		for(int i=0;i<R;i++){
			if(DDR_affiliation[i] == DDR_rank){
				for(int k=0;k<blocksize;k++){
					orig_local_arr[j*blocksize+k] = t*R*blocksize + i*blocksize + k;
				}
				j++;
			}
		}
	}
	TRYCXX( ISCreateGeneral(PETSC_COMM_WORLD, local_size, orig_local_arr, PETSC_COPY_VALUES,&orig_local_is) );
//	TRYCXX( ISCreateGeneral(PETSC_COMM_WORLD, local_size, orig_local_arr, PETSC_OWN_POINTER,&orig_local_is) );

	TRYCXX( ISCreateStride(PETSC_COMM_WORLD, local_size, Tbegin*R*blocksize + Rbegin*blocksize, 1, &new_local_is) );

	/* get subvector with local values from original data */
	TRYCXX( VecGetSubVector(new_Vec, new_local_is, &new_local_Vec) );
	TRYCXX( VecGetSubVector(orig_Vec, orig_local_is, &orig_local_Vec) );

	/* copy values */
	if(!invert){
		TRYCXX( VecCopy(orig_local_Vec, new_local_Vec) );
	} else {
		TRYCXX( VecCopy(new_local_Vec, orig_local_Vec) );
	}

	/* restore subvector with local values from original data */
	TRYCXX( VecRestoreSubVector(new_Vec, new_local_is, &new_local_Vec) );
	TRYCXX( VecRestoreSubVector(orig_Vec, orig_local_is, &orig_local_Vec) );
	
	/* destroy used stuff */
	TRYCXX( ISDestroy(&orig_local_is) );
	TRYCXX( ISDestroy(&new_local_is) );

	TRYCXX( PetscBarrier(NULL));

	LOG_FUNC_END
}

void Decomposition::createIS_gammaK(IS *is, int k) const {
	LOG_FUNC_BEGIN
	
	TRYCXX( ISCreateStride(PETSC_COMM_WORLD, get_Tlocal()*get_Rlocal(), get_Tbegin()*get_R()*get_K() + get_Rbegin()*get_K() + k, get_K(), is) );

	LOG_FUNC_END
}


}
} /* end of namespace */

#endif
