#include "external/seqarrayvector/common/decomposition.h"

namespace pascinference {
namespace algebra {

template<>
Decomposition<SeqArrayVector>::Decomposition(int T, int R, int K, int xdim, int DDT_size){
	LOG_FUNC_BEGIN

	this->T = T;
	this->R = R;
	this->K = K;
	this->xdim = xdim;
		
	this->DDT_size = 1;
	this->DDR_size = 1;
	
	/* no graph provided */
	graph = NULL;

	LOG_FUNC_END
}

template<>
Decomposition<SeqArrayVector>::Decomposition(int T, BGMGraph<SeqArrayVector> &new_graph, int K, int xdim, int DDR_size){
	LOG_FUNC_BEGIN

	this->T = T;
	this->R = new_graph.get_n();
	this->K = K;
	this->xdim = xdim;

	this->DDT_size = 1;
	this->DDR_size = 1;

	/* prepare new layout for T */
	destroy_DDT_arrays = true;	
	DDT_ranges = (int *)malloc((this->DDT_size+1)*sizeof(int));
	DDT_ranges[0] = 0;
	DDT_ranges[1] = T;

	compute_rank();

	LOG_FUNC_END
}

template<>
Decomposition<SeqArrayVector>::Decomposition(int T, BGMGraph<SeqArrayVector> &new_graph, int K, int xdim, int DDT_size, int DDR_size){
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

template<>
Decomposition<SeqArrayVector>::~Decomposition(){
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

template<>
void Decomposition<SeqArrayVector>::compute_rank(){
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

template<>
void Decomposition<SeqArrayVector>::set_graph(BGMGraph<SeqArrayVector> &new_graph, int DDR_size) {

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


}
} /* end of namespace */

