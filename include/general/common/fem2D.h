/** @file femhat.h
 *  @brief class for reduction and prolongation on fem meshes using hat functions
 *
 *  @author Lukas Pospisil
 */
 
#ifndef PASC_FEM2D_H
#define	PASC_FEM2D_H

#include "general/common/fem.h"
#include "general/algebra/graph/bgmgraphgrid2D.h"

namespace pascinference {
namespace common {

/** \class FEM 2D manipulation
 *  \brief Reduction/prolongation between FEM meshes.
 *
*/
template<class VectorBase>
class Fem2D : public Fem<VectorBase> {
	public:
		class ExternalContent;

	protected:
		friend class ExternalContent;
		ExternalContent *externalcontent;			/**< for manipulation with external-specific stuff */

		bool left_overlap;			/**< is there overlap to the left side of space? */
		bool right_overlap;			/**< is there overlap to the right side of space? */
		bool top_overlap;			/**< is there overlap to the top side of space? */
		bool bottom_overlap;		/**< is there overlap to the bottom side of space? */
		
		void compute_overlaps();
		int overlap1_idx_size;		/**< number of elements in overlap1_idx */
		int *overlap1_idx;			/**< permutated indexes of overlap part in grid1 */
		int overlap2_idx_size;		/**< number of elements in overlap2_idx */
		int *overlap2_idx;			/**< permutated indexes of overlap part in grid2 */

		double diff_x;
		double diff_y;

		BGMGraphGrid2D<VectorBase> *grid1;
		BGMGraphGrid2D<VectorBase> *grid2;

		int *bounding_box1;		/**< bounds of local domain [x1_min,x1_max,y1_min,y1_max] of grid1 */
		int *bounding_box2;		/**< bounds of local domain [x2_min,x2_max,y2_min,y2_max] of grid2 */

	public:
		/** @brief create FEM mapping between two decompositions
		*/
		Fem2D(Decomposition<VectorBase> *decomposition1, Decomposition<VectorBase> *decomposition2, double fem_reduce);

		/** @brief create general FEM mapping
		 * 
		 * do not forget to call Fem::set_decomposition_original() and afterwards Fem::compute_decomposition_reduced() to compute decomposition2 internally
		 * 
		 */
		Fem2D(double fem_reduce = 1.0);

		/** @brief destructor
		*/
		~Fem2D();

		void print(ConsoleOutput &output_global, ConsoleOutput &output_local) const;
		std::string get_name() const;
		
		void reduce_gamma(GeneralVector<VectorBase> *gamma1, GeneralVector<VectorBase> *gamma2) const;
		void prolongate_gamma(GeneralVector<VectorBase> *gamma2, GeneralVector<VectorBase> *gamma1) const;

		void compute_decomposition_reduced();

		ExternalContent *get_externalcontent() const;	
};

/* ----------------- Fem implementation ------------- */
template<class VectorBase>
Fem2D<VectorBase>::Fem2D(double fem_reduce) : Fem<VectorBase>(fem_reduce){
	LOG_FUNC_BEGIN
	
	/* I don't have this information without decompositions */
	this->diff_x = 0;
	this->diff_y = 0;
	
	this->grid1 = NULL;
	this->grid2 = NULL;
	
	this->bounding_box1 = new int[4];
	set_value_array(4, this->bounding_box1, 0); /* initial values */

	this->bounding_box2 = new int[4];
	set_value_array(4, this->bounding_box2, 0); /* initial values */

	
	LOG_FUNC_END
}

template<class VectorBase>
Fem2D<VectorBase>::Fem2D(Decomposition<VectorBase> *decomposition1, Decomposition<VectorBase> *decomposition2, double fem_reduce) : Fem<VectorBase>(decomposition1, decomposition2, fem_reduce){
	LOG_FUNC_BEGIN

	this->grid1 = (BGMGraphGrid2D<VectorBase>*)(this->decomposition1->get_graph());
	this->grid2 = (BGMGraphGrid2D<VectorBase>*)(this->decomposition2->get_graph());

	this->diff = 1; /* time */
	this->diff_x = (grid1->get_width()-1)/(double)(grid2->get_width()-1);
	this->diff_y = (grid1->get_height()-1)/(double)(grid2->get_height()-1);

	if(this->is_reduced()){
		this->bounding_box1 = new int[4];
		set_value_array(4, this->bounding_box1, 0); /* initial values */
		this->bounding_box2 = new int[4];
		set_value_array(4, this->bounding_box2, 0); /* initial values */

		compute_overlaps();
	}

	LOG_FUNC_END
}

template<class VectorBase>
Fem2D<VectorBase>::~Fem2D(){
	LOG_FUNC_BEGIN

	LOG_FUNC_END
}

template<class VectorBase>
std::string Fem2D<VectorBase>::get_name() const {
	return "FEM2D";
}

template<class VectorBase>
void Fem2D<VectorBase>::print(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN

	output_global << this->get_name() << std::endl;
	
	/* information of reduced problem */
	output_global <<  " - is reduced       : " << printbool(this->is_reduced()) << std::endl;
	output_global <<  " - diff             : " << this->diff << std::endl;
	output_global <<  " - diff_x           : " << this->diff_x << std::endl;
	output_global <<  " - diff_y           : " << this->diff_y << std::endl;

	output_global <<  " - bounding_box" << std::endl;
	output_local <<   "   - bounding_box1    : " << print_array(this->bounding_box1, 4) << std::endl;
	output_local <<   "   - bounding_box2    : " << print_array(this->bounding_box2, 4) << std::endl;
	output_local.synchronize();

	output_global <<  " - fem_reduce       : " << this->fem_reduce << std::endl;
	output_global <<  " - fem_type         : " << get_name() << std::endl;
	
	if(this->decomposition1 == NULL){
		output_global <<  " - decomposition1   : NO" << std::endl;
	} else {
		output_global <<  " - decomposition1   : YES" << std::endl;
		output_global.push();
		this->decomposition1->print(output_global);
		output_global.pop();
	}
	if(grid1 == NULL){
		output_global <<  " - grid1            : NO" << std::endl;
	} else {
		output_global <<  " - grid1            : YES [" << grid1->get_width() << ", " << grid1->get_height() << "]" << std::endl;
		output_global.push();
		grid1->print(output_global);
		output_global.pop();
	}

	if(this->decomposition2 == NULL){
		output_global <<  " - decomposition2   : NO" << std::endl;
	} else {
		output_global <<  " - decomposition2   : YES" << std::endl;
		output_global.push();
		this->decomposition2->print(output_global);
		output_global.pop();
	}
	if(grid2 == NULL){
		output_global <<  " - grid2            : NO" << std::endl;
	} else {
		output_global <<  " - grid2            : YES [" << grid2->get_width() << ", " << grid2->get_height() << "]" << std::endl;
		output_global.push();
		grid2->print(output_global);
		output_global.pop();
	}
	
	output_global.synchronize();	

	LOG_FUNC_END
}

template<class VectorBase>
void Fem2D<VectorBase>::compute_overlaps() {
	LOG_FUNC_BEGIN
	
	if(this->is_reduced()){
		int width1 = grid1->get_width();
		int height1 = grid1->get_height();
		int width2 = grid2->get_width();
		int height2 = grid2->get_height();
		
		/* get arrays of grids */
		int *DD_affiliation1 = grid1->get_DD_affiliation(); 
		int *DD_permutation1 = grid1->get_DD_permutation(); 
		int *DD_invpermutation1 = grid1->get_DD_invpermutation(); 

		int *DD_affiliation2 = grid2->get_DD_affiliation(); 
		int *DD_permutation2 = grid2->get_DD_permutation(); 
		int *DD_invpermutation2 = grid2->get_DD_invpermutation(); 

		/* prepare overlap indexes */
		overlap1_idx_size = (bounding_box1[1]-bounding_box1[0]+1)*(bounding_box1[3]-bounding_box1[2]+1);
		overlap1_idx = new int[overlap1_idx_size];
		overlap2_idx_size = (bounding_box2[1]-bounding_box2[0]+1)*(bounding_box2[3]-bounding_box2[2]+1);
		overlap2_idx = new int[overlap2_idx_size];
		
		/* fill overlapping indexes with.. indexes */
		for(int id_x1 = bounding_box1[0]; id_x1 <= bounding_box1[1]; id_x1++){
			for(int id_y1 = bounding_box1[2]; id_y1 <= bounding_box1[3]; id_y1++){
				overlap1_idx[(id_y1-bounding_box1[2])*(bounding_box1[1]-bounding_box1[0]+1) + (id_x1-bounding_box1[0])] = DD_permutation1[id_y1*width1 + id_x1];
			}
		}
		for(int id_x2 = bounding_box2[0]; id_x2 <= bounding_box2[1]; id_x2++){
			for(int id_y2 = bounding_box2[2]; id_y2 <= bounding_box2[3]; id_y2++){
				overlap2_idx[(id_y2-bounding_box2[2])*(bounding_box2[1]-bounding_box2[0]+1) + (id_x2-bounding_box2[0])] = DD_permutation2[id_y2*width2 + id_x2];
			}
		}
		
	}
	
	LOG_FUNC_END
}

template<class VectorBase>
void Fem2D<VectorBase>::reduce_gamma(GeneralVector<VectorBase> *gamma1, GeneralVector<VectorBase> *gamma2) const {
	LOG_FUNC_BEGIN

	//TODO

	LOG_FUNC_END
}

template<class VectorBase>
void Fem2D<VectorBase>::prolongate_gamma(GeneralVector<VectorBase> *gamma2, GeneralVector<VectorBase> *gamma1) const {
	LOG_FUNC_BEGIN

	//TODO

	LOG_FUNC_END
}

template<class VectorBase>
void Fem2D<VectorBase>::compute_decomposition_reduced() {
	LOG_FUNC_BEGIN
	
	/* decomposition1 has to be set */
	this->grid1 = (BGMGraphGrid2D<VectorBase>*)(this->decomposition1->get_graph());

	if(this->is_reduced()){
		int T_reduced = 1;
		int width_reduced = ceil(grid1->get_width()*this->fem_reduce);
		int height_reduced = ceil(grid1->get_height()*this->fem_reduce);
		
		this->grid2 = new BGMGraphGrid2D<VectorBase>(width_reduced, height_reduced);
		this->grid2->process_grid();
		
		/* decompose second grid based on the decomposition of the first grid */
		this->grid2->decompose(this->grid1, this->bounding_box1, this->bounding_box2);
		this->grid2->print(coutMaster);
		
		/* compute new decomposition */
		this->decomposition2 = new Decomposition<VectorBase>(T_reduced, 
				*(this->grid2), 
				this->decomposition1->get_K(), 
				this->decomposition1->get_xdim(), 
				this->decomposition1->get_DDT_size(), 
				this->decomposition1->get_DDR_size());

		compute_overlaps();
	} else {
		/* there is not reduction of the data, we can reuse the decomposition */
		this->set_decomposition_reduced(this->decomposition1);
		this->grid2 = this->grid1;
	}

	#ifdef USE_CUDA
		/* compute optimal kernel calls */
		gpuErrchk( cudaOccupancyMaxPotentialBlockSize( &minGridSize_reduce, &blockSize_reduce, kernel_femhat_reduce_data, 0, 0) );
		gridSize_reduce = (this->decomposition2->get_Tlocal() + blockSize_reduce - 1)/ blockSize_reduce;

		gpuErrchk( cudaOccupancyMaxPotentialBlockSize( &minGridSize_prolongate, &blockSize_prolongate, kernel_femhat_prolongate_data, 0, 0) );
		gridSize_prolongate = (this->decomposition1->get_Tlocal() + blockSize_prolongate - 1)/ blockSize_prolongate;
	#endif

	this->diff = 1; /* time */
	this->diff_x = (grid1->get_width()-1)/(double)(grid2->get_width()-1);
	this->diff_y = (grid1->get_height()-1)/(double)(grid2->get_height()-1);
			
	LOG_FUNC_END
}

#ifdef USE_CUDA
__global__ void kernel_femhat_reduce_data(double *data1, double *data2, int T1, int T2, int Tbegin1, int Tbegin2, int T1local, int T2local, int left_t1_idx, int left_t2_idx, double diff) {
	int t2 = blockIdx.x*blockDim.x + threadIdx.x;

	if(t2 < T2local){
		double center_t1 = (Tbegin2+t2)*diff;
		double left_t1 = (Tbegin2+t2-1)*diff;
		double right_t1 = (Tbegin2+t2+1)*diff;
				
		int id_counter = floor(left_t1) - left_t1_idx; /* first index in provided local t1 array */

		double phi_value; /* value of basis function */

		/* left part of hat function */
		double mysum = 0.0;
		int t1 = floor(left_t1);

		/* compute linear combination with coefficients given by basis functions */
		while(t1 <= center_t1){
			phi_value = (t1 - left_t1)/(center_t1 - left_t1);
			mysum += phi_value*data1[id_counter];
			t1 += 1;
			id_counter += 1;
		}

		/* right part of hat function */
		while(t1 < right_t1){
			phi_value = (t1 - right_t1)/(center_t1 - right_t1);
			mysum += phi_value*data1[id_counter];
			t1 += 1;
			id_counter += 1;
		}

		data2[t2] = mysum;
	}
}


__global__ void kernel_femhat_prolongate_data(double *data1, double *data2, int T1, int T2, int Tbegin1, int Tbegin2, int T1local, int T2local, int left_t1_idx, int left_t2_idx, double diff) {
	int t1 = blockIdx.x*blockDim.x + threadIdx.x;

	if(t1 < T1local){
		int t2_left_id_orig = floor((t1 + Tbegin1)/diff);
		int t2_right_id_orig = floor((t1 + Tbegin1)/diff) + 1;

		double t1_left = t2_left_id_orig*diff;
		double t1_right = t2_right_id_orig*diff;

		int t2_left_id = t2_left_id_orig - left_t2_idx;
		int t2_right_id = t2_right_id_orig - left_t2_idx;

		/* value of basis functions */
		double t1_value = 0.0;
		double phi_value_left = (t1 + Tbegin1 - t1_left)/(t1_right - t1_left); 
		t1_value += phi_value_left*data2[t2_right_id];
				
		double phi_value_right = (t1 + Tbegin1 - t1_right)/(t1_left - t1_right); 
		t1_value += phi_value_right*data2[t2_left_id];

		data1[t1] = t1_value;
	}
}

#endif




}
} /* end of namespace */

#endif
