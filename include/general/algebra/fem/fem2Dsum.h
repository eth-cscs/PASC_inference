/** @file fem2Dsum.h
 *  @brief class for reduction and prolongation on fem meshes using constant functions in 2D
 *
 *  @author Lukas Pospisil
 */

#ifndef PASC_FEM2DSUM_H
#define	PASC_FEM2DSUM_H

#include "general/algebra/fem/fem.h"
#include "general/algebra/graph/bgmgraphgrid2D.h"

namespace pascinference {
namespace algebra {

/** \class Fem2DSum
 *  \brief Reduction/prolongation between FEM meshes using averaging.
 *
*/
template<class VectorBase>
class Fem2DSum : public Fem<VectorBase> {
	public:
		class ExternalContent;

	protected:
		friend class ExternalContent;
		ExternalContent *externalcontent;			/**< for manipulation with external-specific stuff */

		void compute_overlaps();
		void compute_bounding_box();

		int overlap1_idx_size;		/**< number of elements in overlap1_idx */
		int *overlap1_idx;			/**< permutated indexes of overlap part in grid1 */
		int overlap2_idx_size;		/**< number of elements in overlap2_idx */
		int *overlap2_idx;			/**< permutated indexes of overlap part in grid2 */

		double diff_x;
		double diff_y;

		BGMGraphGrid2D<VectorBase> *grid1;
		BGMGraphGrid2D<VectorBase> *grid2;

		int *bounding_box1;		/**< bounds of local domain [x1_min,x1_max,y1_min,y1_max] of grid1 in R format */
		int *bounding_box2;		/**< bounds of local domain [x2_min,x2_max,y2_min,y2_max] of grid2 in R format */

		virtual double get_overlap_extension() const;
	public:
		/** @brief create FEM mapping between two decompositions
		*/
		Fem2DSum(Decomposition<VectorBase> *decomposition1, Decomposition<VectorBase> *decomposition2, double fem_reduce);

		/** @brief create general FEM mapping
		 *
		 * do not forget to call Fem::set_decomposition_original() and afterwards Fem::compute_decomposition_reduced() to compute decomposition2 internally
		 *
		 */
		Fem2DSum(double fem_reduce = 1.0);

		/** @brief destructor
		*/
		~Fem2DSum();

		void print(ConsoleOutput &output_global, ConsoleOutput &output_local) const;
		virtual std::string get_name() const;

		virtual void reduce_gamma(GeneralVector<VectorBase> *gamma1, GeneralVector<VectorBase> *gamma2) const;
		virtual void prolongate_gamma(GeneralVector<VectorBase> *gamma2, GeneralVector<VectorBase> *gamma1) const;

		void compute_decomposition_reduced();

		void set_decomposition_original(Decomposition<VectorBase> *decomposition1);
		void set_decomposition_reduced(Decomposition<VectorBase> *decomposition2);

		BGMGraphGrid2D<VectorBase>* get_grid_original() const;
		BGMGraphGrid2D<VectorBase>* get_grid_reduced() const;

        void saveVTK_bounding_box(std::string filename_bounding_box1, std::string filename_bounding_box2) const;

		ExternalContent *get_externalcontent() const;
};

/* ----------------- Fem implementation ------------- */
template<class VectorBase>
Fem2DSum<VectorBase>::Fem2DSum(double fem_reduce) : Fem<VectorBase>(fem_reduce){
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
Fem2DSum<VectorBase>::Fem2DSum(Decomposition<VectorBase> *decomposition1, Decomposition<VectorBase> *decomposition2, double fem_reduce) : Fem<VectorBase>(decomposition1, decomposition2, fem_reduce){
	LOG_FUNC_BEGIN

	this->grid1 = (BGMGraphGrid2D<VectorBase>*)(this->decomposition1->get_graph());
	this->grid2 = (BGMGraphGrid2D<VectorBase>*)(this->decomposition2->get_graph());

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
Fem2DSum<VectorBase>::~Fem2DSum(){
	LOG_FUNC_BEGIN

	LOG_FUNC_END
}

template<class VectorBase>
std::string Fem2DSum<VectorBase>::get_name() const {
	return "FEM-2D-SUM";
}

template<class VectorBase>
void Fem2DSum<VectorBase>::print(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN

	output_global << this->get_name() << std::endl;

	/* information of reduced problem */
	output_global <<  " - is reduced        : " << printbool(this->is_reduced()) << std::endl;
	output_global <<  " - diff_x            : " << this->diff_x << std::endl;
	output_global <<  " - diff_y            : " << this->diff_y << std::endl;
	output_global <<  " - overlap_extension : " << this->get_overlap_extension() << std::endl;

	output_global <<  " - bounding_box" << std::endl;
	output_local <<   "   - bounding_box1     : " << print_array(this->bounding_box1, 4) << std::endl;
	output_local <<   "   - bounding_box2     : " << print_array(this->bounding_box2, 4) << std::endl;
	output_local.synchronize();

	output_global <<  " - fem_reduce        : " << this->get_fem_reduce() << std::endl;
	output_global <<  " - fem_type          : " << get_name() << std::endl;

	if(this->decomposition1 == NULL){
		output_global <<  " - decomposition1    : NO" << std::endl;
	} else {
		output_global <<  " - decomposition1    : YES" << std::endl;
		output_global.push();
		this->decomposition1->print(output_global);
		output_global.pop();
	}
	if(grid1 == NULL){
		output_global <<  " - grid1             : NO" << std::endl;
	} else {
		output_global <<  " - grid1             : YES [" << grid1->get_width() << ", " << grid1->get_height() << "]" << std::endl;
		output_global.push();
		grid1->print(output_global);
		output_global.pop();
	}

	if(this->decomposition2 == NULL){
		output_global <<  " - decomposition2    : NO" << std::endl;
	} else {
		output_global <<  " - decomposition2    : YES" << std::endl;
		output_global.push();
		this->decomposition2->print(output_global);
		output_global.pop();
	}
	if(grid2 == NULL){
		output_global <<  " - grid2             : NO" << std::endl;
	} else {
		output_global <<  " - grid2             : YES [" << grid2->get_width() << ", " << grid2->get_height() << "]" << std::endl;
		output_global.push();
		grid2->print(output_global);
		output_global.pop();
	}

	output_global.synchronize();

	LOG_FUNC_END
}

template<class VectorBase>
void Fem2DSum<VectorBase>::compute_overlaps() {
	LOG_FUNC_BEGIN

	if(this->is_reduced()){
        compute_bounding_box();

        int TRbegin1 = this->decomposition1->get_TRbegin();
        int TRbegin2 = this->decomposition1->get_TRbegin();
        int Tlocal1 = this->decomposition1->get_Tlocal();
        int Tlocal2 = this->decomposition2->get_Tlocal(); /* = Tlocal1 */

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
		overlap1_idx = new int[Tlocal1*overlap1_idx_size];
		overlap2_idx_size = (bounding_box2[1]-bounding_box2[0]+1)*(bounding_box2[3]-bounding_box2[2]+1);
		overlap2_idx = new int[Tlocal2*overlap2_idx_size];

		/* fill overlapping indexes with.. indexes  */
		for(int t=0; t < Tlocal1; t++){
            for(int x = bounding_box1[0]; x <= bounding_box1[1]; x++){
                for(int y = bounding_box1[2]; y <= bounding_box1[3]; y++){
                    int r = y*width1 + x; /* in original R format */
                    overlap1_idx[t*overlap1_idx_size + (y-bounding_box1[2])*(bounding_box1[1]-bounding_box1[0]+1) + (x-bounding_box1[0])] = DD_permutation1[r];
                }
            }
		}

		for(int t=0; t < Tlocal2; t++){
            for(int x = bounding_box2[0]; x <= bounding_box2[1]; x++){
                for(int y = bounding_box2[2]; y <= bounding_box2[3]; y++){
                    int r = y*width2 + x; /* in original R format */
                    overlap2_idx[t*overlap2_idx_size + (y-bounding_box2[2])*(bounding_box2[1]-bounding_box2[0]+1) + (x-bounding_box2[0])] = DD_permutation2[r];
                }
            }
        }
	}

	LOG_FUNC_END
}

template<class VectorBase>
void Fem2DSum<VectorBase>::compute_bounding_box() {
	LOG_FUNC_BEGIN

	if(this->is_reduced()){
		/* here is the trick: bounding boxes are swaped because of mapping !!! */
        this->grid1->compute_local_bounding_box(bounding_box2, this->decomposition1->get_DDR_invpermutation(), this->decomposition1->get_Rbegin(), this->decomposition1->get_Rlocal());
        this->grid2->compute_local_bounding_box(bounding_box1, this->decomposition2->get_DDR_invpermutation(), this->decomposition2->get_Rbegin(), this->decomposition2->get_Rlocal());

		int width1 = grid1->get_width();
		int height1 = grid1->get_height();
		int width2 = grid2->get_width();
		int height2 = grid2->get_height();

		this->diff_x = (width1-1)/(double)(width2-1);
		this->diff_y = (height1-1)/(double)(height2-1);

        /* scale bounding boxes */
        double overlap_extension = this->get_overlap_extension(); /* extension of bounding box inside the neighbour */

        bounding_box1[0] = floor((bounding_box1[0]-overlap_extension)*this->diff_x);
        bounding_box1[1] = ceil((bounding_box1[1]+overlap_extension)*this->diff_x);
        bounding_box1[2] = floor((bounding_box1[2]-overlap_extension)*this->diff_y);
        bounding_box1[3] = ceil((bounding_box1[3]+overlap_extension)*this->diff_y);

        bounding_box2[0] = floor((bounding_box2[0]-overlap_extension)/this->diff_x);
        bounding_box2[1] = ceil((bounding_box2[1]+overlap_extension)/this->diff_x);
        bounding_box2[2] = floor((bounding_box2[2]-overlap_extension)/this->diff_y);
        bounding_box2[3] = ceil((bounding_box2[3]+overlap_extension)/this->diff_y);

		if(bounding_box1[0] < 0) {
			bounding_box1[0] = 0;
		}
		if(bounding_box1[1] >= width1) {
			bounding_box1[1] = width1-1;
		}
		if(bounding_box1[2] < 0) {
			bounding_box1[2] = 0;
		}
		if(bounding_box1[3] >= height1) {
			bounding_box1[3] = height1-1;
		}

		if(bounding_box2[0] < 0) {
			bounding_box2[0] = 0;
		}
		if(bounding_box2[1] >= width2) {
			bounding_box2[1] = width2-1;
		}
		if(bounding_box2[2] < 0) {
			bounding_box2[2] = 0;
		}
		if(bounding_box2[3] >= height2) {
			bounding_box2[3] = height2-1;
		}

	}

	LOG_FUNC_END
}


template<class VectorBase>
void Fem2DSum<VectorBase>::reduce_gamma(GeneralVector<VectorBase> *gamma1, GeneralVector<VectorBase> *gamma2) const {
	LOG_FUNC_BEGIN

	//TODO

	LOG_FUNC_END
}

template<class VectorBase>
void Fem2DSum<VectorBase>::prolongate_gamma(GeneralVector<VectorBase> *gamma2, GeneralVector<VectorBase> *gamma1) const {
	LOG_FUNC_BEGIN

	//TODO

	LOG_FUNC_END
}

template<class VectorBase>
void Fem2DSum<VectorBase>::compute_decomposition_reduced() {
	LOG_FUNC_BEGIN

	if(this->is_reduced()){

		int T_reduced = this->decomposition1->get_T(); /* this is the reduction only in space, not in time */
		int width_reduced = ceil(grid1->get_width()*this->fem_reduce);
		int height_reduced = ceil(grid1->get_height()*this->fem_reduce);

		this->grid2 = new BGMGraphGrid2D<VectorBase>(width_reduced, height_reduced);
		this->grid2->process_grid();

		/* compute new decomposition */
		/* this automatically decompose new given grid based on this new decomposition */
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

	this->diff_x = (grid1->get_width()-1)/(double)(grid2->get_width()-1);
	this->diff_y = (grid1->get_height()-1)/(double)(grid2->get_height()-1);

	LOG_FUNC_END
}

template<class VectorBase>
void Fem2DSum<VectorBase>::set_decomposition_original(Decomposition<VectorBase> *decomposition1) {
	LOG_FUNC_BEGIN

	this->decomposition1 = decomposition1;

    /* this works if and only if on regular 2D mesh (BGMGraphGrid2D) !! as whole FEM2D */
	this->grid1 = (BGMGraphGrid2D<VectorBase>*)(this->decomposition1->get_graph());

	LOG_FUNC_END
}

template<class VectorBase>
void Fem2DSum<VectorBase>::set_decomposition_reduced(Decomposition<VectorBase> *decomposition2) {
	LOG_FUNC_BEGIN

	this->decomposition2 = decomposition2;

    /* this works if and only if on regular 2D mesh (BGMGraphGrid2D) !! as whole FEM2D */
	this->grid2 = (BGMGraphGrid2D<VectorBase>*)(this->decomposition2->get_graph());

	LOG_FUNC_END
}

template<class VectorBase>
void Fem2DSum<VectorBase>::saveVTK_bounding_box(std::string filename_bounding_box1, std::string filename_bounding_box2) const{
	LOG_FUNC_BEGIN

    if(grid1 != NULL){
        grid1->saveVTK_bounding_box(filename_bounding_box1, bounding_box1);
    }

    if(grid2 != NULL){
        grid2->saveVTK_bounding_box(filename_bounding_box2, bounding_box2);
    }

	LOG_FUNC_END
}

template<class VectorBase>
BGMGraphGrid2D<VectorBase>* Fem2DSum<VectorBase>::get_grid_original() const {
    return this->grid1;
}

template<class VectorBase>
BGMGraphGrid2D<VectorBase>* Fem2DSum<VectorBase>::get_grid_reduced() const {
    return this->grid2;
}

template<class VectorBase>
double Fem2DSum<VectorBase>::get_overlap_extension() const {
    return 1.0;
}


}
} /* end of namespace */

#endif
