/** @file bgmgraphgrid2D.h
 *  @brief class for manipulation with 2D grid
 *
 *  Defines some basic functions for manipulaton with grids, especially for BLOCKGRAPHMATRIX.
 *
 *  @author Lukas Pospisil
 */

#ifndef PASC_COMMON_BGMGRAPHGRID2D_H
#define	PASC_COMMON_BGMGRAPHGRID2D_H

#include "general/algebra/graph/bgmgraph.h"

namespace pascinference {
namespace algebra {

/** \class BGMGraphGrid2D
 *  \brief Graph of two dimensional grid.
 *
 *  Could be used for faster and simplier manipulation with image graph.
 *
*/
template<class VectorBase>
class BGMGraphGrid2D: public BGMGraph<VectorBase> {
	public:
		class ExternalContent;

	protected:
		friend class ExternalContent;
		ExternalContent *externalcontent;			/**< for manipulation with external-specific stuff */

		int width; /**< dimension of grid */
		int height; /**< dimension of grid */

		void process_grid_cuda();
	public:

		BGMGraphGrid2D(int width, int height);
		BGMGraphGrid2D(std::string filename, int dim=2) : BGMGraph<VectorBase>(filename, dim) {};
		BGMGraphGrid2D(const double *coordinates_array, int n, int dim) : BGMGraph<VectorBase>(coordinates_array, n, dim) {};

		~BGMGraphGrid2D();

		virtual std::string get_name() const;
		virtual void process_grid();

		int get_width() const;
		int get_height() const;

        void compute_local_bounding_box(int *bounding_box_local, int *permutation, int start_idx, int local_nmb) const;
        void saveVTK_bounding_box(std::string filename, int *bounding_box_local) const;

		ExternalContent *get_externalcontent() const;
};


}
} /* end of namespace */



/* -------- IMPLEMENTATION ------- */
namespace pascinference {
namespace algebra {

template<class VectorBase>
BGMGraphGrid2D<VectorBase>::BGMGraphGrid2D(int width, int height) : BGMGraph<VectorBase>(){
	LOG_FUNC_BEGIN

	this->width = width;
	this->height = height;

	this->dim = 2;
	this->n = width*height;

	//TODO

	LOG_FUNC_END
}

template<class VectorBase>
BGMGraphGrid2D<VectorBase>::~BGMGraphGrid2D(){
	LOG_FUNC_BEGIN

	LOG_FUNC_END
}

template<class VectorBase>
void BGMGraphGrid2D<VectorBase>::process_grid(){
	LOG_FUNC_BEGIN

	this->threshold = 1.1;
	this->m = height*(width-1) + width*(height-1);
	this->m_max = 4;

	/* prepare array for number of neighbors */
	this->neighbor_nmbs = (int*)malloc(this->n*sizeof(int));
	this->neighbor_ids = (int**)malloc(this->n*sizeof(int*));

//	#pragma omp parallel for
	for(int idx=0;idx<width*height;idx++){
		int i = idx/(double)width; /* index of row */
		int j = idx - i*width; /* index of column */

		/* compute number of neighbors */
		int nmb = 0;
		if(j>0){
			nmb+=1;
		}
		if(j<width-1){
			nmb+=1;
		}
		if(i>0){
			nmb+=1;
		}
		if(i<height-1){
			nmb+=1;
		}
		this->neighbor_nmbs[idx] = nmb;
		this->neighbor_ids[idx] = (int*)malloc(this->neighbor_nmbs[idx]*sizeof(int));

		/* fill neighbors */
		nmb = 0;
		if(j>0){ /* left */
			this->neighbor_ids[idx][nmb] = idx-1;
			nmb+=1;
		}
		if(j<width-1){ /* right */
			this->neighbor_ids[idx][nmb] = idx+1;
			nmb+=1;
		}
		if(i>0){ /* down */
			this->neighbor_ids[idx][nmb] = idx-width;
			nmb+=1;
		}
		if(i<height-1){ /* up */
			this->neighbor_ids[idx][nmb] = idx+width;
			nmb+=1;
		}
	}

	this->processed = true;

	LOG_FUNC_END
}

template<class VectorBase>
std::string BGMGraphGrid2D<VectorBase>::get_name() const {
	return "BGMGraphGrid2D";
}

template<class VectorBase>
int BGMGraphGrid2D<VectorBase>::get_width() const {
	return this->width;
}

template<class VectorBase>
int BGMGraphGrid2D<VectorBase>::get_height() const {
	return this->height;
}

template<class VectorBase>
void BGMGraphGrid2D<VectorBase>::compute_local_bounding_box(int *bounding_box, int *permutation, int start_idx, int local_nmb) const {
	LOG_FUNC_BEGIN

    bounding_box[0] = this->width-1; /* inf */
    bounding_box[1] = 0; /* -inf */
    bounding_box[2] = this->height-1; /* inf */
    bounding_box[3] = 0; /* -inf */

    int Ridx; /* index in R */
    int x,y; /* coordinates in grid */
    for(int r=start_idx;r<start_idx+local_nmb;r++){
        /* convert dR to R */
        Ridx = permutation[r];

        /* compute coordinates in grid */
        y = (int)(Ridx/(double)this->width);
        x = Ridx - y*this->width;

        /* update bounding box values */
        if(x < bounding_box[0]){
            bounding_box[0] = x;
        }
        if(x > bounding_box[1]){
            bounding_box[1] = x;
        }
        if(y < bounding_box[2]){
            bounding_box[2] = y;
        }
        if(y > bounding_box[3]){
            bounding_box[3] = y;
        }
    }

	LOG_FUNC_END
}

template<class VectorBase>
void BGMGraphGrid2D<VectorBase>::saveVTK_bounding_box(std::string filename, int *bounding_box) const {
	LOG_FUNC_BEGIN

    //TODO

	LOG_FUNC_END
}


}
} /* end of namespace */

#endif
