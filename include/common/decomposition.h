/** @file decomposition.h
 *  @brief class for storing and manipulation with decomposition informations
 *
 *  @author Lukas Pospisil
 */
 
#ifndef PASC_COMMON_DECOMPOSITION_H
#define	PASC_COMMON_DECOMPOSITION_H

#include "common/bgmgraph.h"

namespace pascinference {

class Decomposition {
	protected:
		/* T decomposition */
		int T;
		int DDT_size; /**< time - number of domains for decomposition */
		int *DDT_ranges; /**< time - local ownerships */
		bool destroy_DDT_arrays; /**< if decomposition of T was not provided, I have to create and destroy arrays here */

		/* R decomposition */
		int R;
		int DDR_size; /**< space - number of domains for decomposition */
		int *DDR_affiliation; /**< space - domain affiliation of vertices */
		int *DDR_permutation; /**< space - permutation of global indexes */
		int *DDR_lengths; /**< space - array of local lengths */
		int *DDR_ranges; /**< space - local ownership */
		bool destroy_DDR_arrays; /**< if graph was not provided, I have to create and destroy arrays here */
	
		/* my coordinates */
		int DDT_rank; /**< my coordinate in time decomposition */
		int DDR_rank; /**< my coordinate in space decomposition */

		/** @brief compute coordinates in decomposition based on rank of the processor
		*/		
		void compute_rank();
		
	public:
		/** @brief decomposition only in time
		*/
		Decomposition(int T, int DDT_size, int R=1);

		/** @brief decomposition only in space
		*/
		Decomposition(int T, BGMGraph &new_graph, int DDR_size);

		/** @brief decomposition in time and space
		*/
		Decomposition(int T, int DDT_size, BGMGraph &new_graph, int DDR_size);

		/** @brief destructor
		*/
		~Decomposition();
		
		int get_T() const;
		int get_DDT_size() const;
		const int *get_DDT_ranges() const;

		int get_R() const;
		int get_DDR_size() const;
		int *get_DDR_affiliation() const;
		int *get_DDR_permutation() const;
		int *get_DDR_lengths() const;
		int *get_DDR_ranges() const;
	
		int get_DDT_rank() const;
		int get_DDR_rank() const;

		/** @brief print informations of decomposition
		*/
		void print(ConsoleOutput &output) const;

		/** @brief print informations of decomposition
		*/
		void print_content(ConsoleOutput &output_master, ConsoleOutput &output_local, bool print_details=true) const;

};

/* ----------------- Decomposition implementation ------------- */

Decomposition::Decomposition(int T, int DDT_size, int R){
	LOG_FUNC_BEGIN

	this->T = T;
	this->R = R;
	
	this->DDT_size = DDT_size;
	this->DDR_size = 1;
	
	/* prepare new layout for T */
	destroy_DDT_arrays = true;	
	DDT_ranges = (int *)malloc((this->DDT_size+1)*sizeof(int));
	
	Vec DDT_layout;
	TRY( VecCreate(PETSC_COMM_WORLD,&DDT_layout) );
	TRY( VecSetSizes(DDT_layout, PETSC_DECIDE, T ));
	TRY( VecSetFromOptions(DDT_layout) );
	
	/* get properties of this layout */
	const int *DDT_ranges_const;
	TRY( VecGetOwnershipRanges(DDT_layout, &DDT_ranges_const) );
	for(int i=0;i<DDT_size+1;i++){
		DDT_ranges[i] = DDT_ranges_const[i]; // TODO: how to deal with const int in ranges form PETSc?
	}
	
	/* destroy temp vector */
	TRY( VecDestroy(&DDT_layout) );

	/* prepare new layout for R */
	/* no graph provided - allocate arrays */
	destroy_DDR_arrays = true;
	DDR_affiliation = (int *)malloc(R*sizeof(int));
	for(int i=0;i<R;i++){
		DDR_affiliation[i] = 0;
	}

	DDR_permutation = (int *)malloc(R*sizeof(int));
	for(int i=0;i<R;i++){
		DDR_permutation[i] = i;
	}

	DDR_lengths = (int *)malloc((this->DDR_size)*sizeof(int));
	DDR_lengths[0] = R;

	DDR_ranges = (int *)malloc((this->DDR_size+1)*sizeof(int));
	DDR_ranges[0] = 0;
	DDR_ranges[1] = R;

	compute_rank();

	LOG_FUNC_END
}

Decomposition::Decomposition(int T, BGMGraph &new_graph, int DDR_size){
	LOG_FUNC_BEGIN

	this->T = T;
	this->R = new_graph.get_n();

	/* decompose graph */
	new_graph.decompose(DDR_size);
	
	this->DDT_size = 1;
	this->DDR_size = new_graph.get_DD_size();
	
	/* prepare new layout for T */
	destroy_DDT_arrays = true;	
	DDT_ranges = (int *)malloc((this->DDT_size+1)*sizeof(int));
	DDT_ranges[0] = 0;
	DDT_ranges[1] = T;

	/* prepare new layout for R */
	/* no graph provided - allocate arrays */
	destroy_DDR_arrays = false;
	DDR_affiliation = new_graph.get_DD_affiliation();
	DDR_permutation = new_graph.get_DD_permutation();
	DDR_lengths = new_graph.get_DD_lengths();
	DDR_ranges = new_graph.get_DD_ranges();

	compute_rank();

	LOG_FUNC_END
}

Decomposition::Decomposition(int T, int DDT_size, BGMGraph &new_graph, int DDR_size){
	LOG_FUNC_BEGIN

	this->T = T;
	this->R = new_graph.get_n();

	/* decompose graph */
	new_graph.decompose(DDR_size);
	
	this->DDT_size = DDT_size;
	this->DDR_size = new_graph.get_DD_size();
	
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

	/* prepare new layout for R */
	/* no graph provided - allocate arrays */
	destroy_DDR_arrays = false;
	DDR_affiliation = new_graph.get_DD_affiliation();
	DDR_permutation = new_graph.get_DD_permutation();
	DDR_lengths = new_graph.get_DD_lengths();
	DDR_ranges = new_graph.get_DD_ranges();

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
		free(DDR_lengths);
		free(DDR_ranges);
	}

	LOG_FUNC_END
}

void Decomposition::compute_rank(){
	LOG_FUNC_BEGIN

	/* control the decomposition */
	if(DDT_size*DDR_size != GlobalManager.get_size()){
		coutMaster << "Sorry, DDT_size*DDR_size != nproc" << std::endl;
		coutMaster << " DDT_size = " << DDT_size << std::endl;
		coutMaster << " DDR_size = " << DDR_size << std::endl;
		coutMaster << " nproc    = " << GlobalManager.get_size() << std::endl;
		
		// TODO: throw error
	}

	/* get rank of this processor */
	int rank = GlobalManager.get_rank();

	this->DDT_rank = rank/(double)this->DDR_size;
	this->DDR_rank = rank - (this->DDT_rank)*(this->DDR_size);
	
	LOG_FUNC_END
}

int Decomposition::get_T() const{
	return T;
}

int Decomposition::get_DDT_size() const{
	return DDT_size;
}

const int *Decomposition::get_DDT_ranges() const{
	return DDT_ranges;
}

int Decomposition::get_R() const{
	return R;
}

int Decomposition::get_DDR_size() const{
	return DDR_size;
}

int *Decomposition::get_DDR_affiliation() const{
	return DDR_affiliation;
}

int *Decomposition::get_DDR_permutation() const{
	return DDR_permutation;
}

int *Decomposition::get_DDR_lengths() const{
	return DDR_lengths;
}

int *Decomposition::get_DDR_ranges() const{
	return DDR_ranges;
}

int Decomposition::get_DDT_rank() const{
	return DDT_rank;
}

int Decomposition::get_DDR_rank() const{
	return DDR_rank;
}

void Decomposition::print_content(ConsoleOutput &output_master, ConsoleOutput &output_local, bool print_details) const {
	output_master << "Decomposition" << std::endl;
	
	output_master.push();
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
	output.pop();
	
}


}

#endif
