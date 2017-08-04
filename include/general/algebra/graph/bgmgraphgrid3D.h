/** @file bgmgraphgrid3D.h
 *  @brief class for manipulation with 3D grid
 *
 *  Defines some basic functions for manipulaton with grids, especially for BLOCKGRAPHMATRIX.
 *
 *  @author Lukas Pospisil
 */

#ifndef PASC_COMMON_BGMGRAPHGRID3D_H
#define	PASC_COMMON_BGMGRAPHGRID3D_H

#include "general/algebra/graph/bgmgraph.h"

namespace pascinference {
namespace algebra {

/** \class BGMGraphGrid3D
 *  \brief Graph of three dimensional grid.
 *
 *  Could be used for faster and simplier manipulation with 3D image graph.
 *
*/
template<class VectorBase>
class BGMGraphGrid3D: public BGMGraph<VectorBase> {
	public:
		class ExternalContent;

	protected:
		friend class ExternalContent;
		ExternalContent *externalcontent;			/**< for manipulation with external-specific stuff */

		int x_size; /**< dimension of grid */
		int y_size; /**< dimension of grid */
		int z_size; /**< dimension of grid */

		void process_grid_cuda();
	public:

		BGMGraphGrid3D(int x_size, int y_size, int z_size);
		BGMGraphGrid3D(std::string filename) : BGMGraph<VectorBase>(filename, 3) {};
		BGMGraphGrid3D(const double *coordinates_array, int n) : BGMGraph<VectorBase>(coordinates_array, n, 3) {};

		~BGMGraphGrid3D();

		virtual std::string get_name() const;
		virtual void process_grid();

		int get_x_size() const;
		int get_y_size() const;
		int get_z_size() const;

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
BGMGraphGrid3D<VectorBase>::BGMGraphGrid3D(int x_size, int y_size, int z_size) : BGMGraph<VectorBase>(){
	LOG_FUNC_BEGIN

	this->x_size = x_size;
	this->y_size = y_size;
	this->z_size = z_size;

	this->dim = 3;
	this->n = x_size*y_size*z_size;

	//TODO

	LOG_FUNC_END
}

template<class VectorBase>
BGMGraphGrid3D<VectorBase>::~BGMGraphGrid3D(){
	LOG_FUNC_BEGIN

	LOG_FUNC_END
}

template<class VectorBase>
void BGMGraphGrid3D<VectorBase>::process_grid(){
	LOG_FUNC_BEGIN

	this->threshold = 1.1;
	this->m = (2*x_size - 1)*(2*y_size - 1)*(2*z_size - 1) - x_size*y_size*z_size; /* number of edges */ 
	this->m_max = 6; /* max number of neighbours */

	/* prepare array for number of neighbors */
	this->neighbor_nmbs = (int*)malloc(this->n*sizeof(int));
	this->neighbor_ids = (int**)malloc(this->n*sizeof(int*));

//	#pragma omp parallel for
	for(int idx=0;idx<x_size*y_size*z_size;idx++){
		int z = idx/(double)(x_size*y_size); /* index of z */
		int y = (idx-z*x_size*y_size)/(double)(x_size); /* index of z */
		int x = idx- z*x_size*y_size - y*x_size; /* index of x */

		/* compute number of neighbors */
		int nmb = 0;
		if(x > 0)	nmb+=1;
		if(x < x_size-1) nmb+=1;
		if(y > 0)	nmb+=1;
		if(y < y_size-1) nmb+=1;
		if(z > 0)	nmb+=1;
		if(z < z_size-1) nmb+=1;

		this->neighbor_nmbs[idx] = nmb;
		this->neighbor_ids[idx] = (int*)malloc(this->neighbor_nmbs[idx]*sizeof(int));

		/* fill neighbors */
		nmb = 0;
		if(x > 0){ /* left */
			this->neighbor_ids[idx][nmb] = idx-1;
			nmb+=1;
		}
		if(x < x_size-1){ /* right */
			this->neighbor_ids[idx][nmb] = idx+1;
			nmb+=1;
		}
		if(y > 0){ /* bottom */
			this->neighbor_ids[idx][nmb] = idx-x_size;
			nmb+=1;
		}
		if(y < y_size-1){ /* top */
			this->neighbor_ids[idx][nmb] = idx+x_size;
			nmb+=1;
		}
		if(z > 0){ /* closer */
			this->neighbor_ids[idx][nmb] = idx-x_size*y_size;
			nmb+=1;
		}
		if(z < z_size-1){ /* distant */
			this->neighbor_ids[idx][nmb] = idx+x_size*y_size;
			nmb+=1;
		}
	}

	this->processed = true;

	LOG_FUNC_END
}

template<class VectorBase>
std::string BGMGraphGrid3D<VectorBase>::get_name() const {
	return "BGMGraphGrid3D";
}

template<class VectorBase>
int BGMGraphGrid3D<VectorBase>::get_x_size() const {
	return this->x_size;
}

template<class VectorBase>
int BGMGraphGrid3D<VectorBase>::get_y_size() const {
	return this->y_size;
}

template<class VectorBase>
int BGMGraphGrid3D<VectorBase>::get_z_size() const {
	return this->z_size;
}

template<class VectorBase>
void BGMGraphGrid3D<VectorBase>::compute_local_bounding_box(int *bounding_box, int *permutation, int start_idx, int local_nmb) const {
	LOG_FUNC_BEGIN

    bounding_box[0] = this->x_size-1; /* inf */
    bounding_box[1] = 0; /* -inf */
    bounding_box[2] = this->y_size-1; /* inf */
    bounding_box[3] = 0; /* -inf */
    bounding_box[4] = this->z_size-1; /* inf */
    bounding_box[5] = 0; /* -inf */

    for(int r=start_idx;r<start_idx+local_nmb;r++){
        /* convert dR to R */
		int Ridx = permutation[r];

        /* compute coordinates in grid */
		int z = Ridx/(double)(x_size*y_size); /* index of z */
		int y = (Ridx-z*x_size*y_size)/(double)(x_size); /* index of z */
		int x = Ridx- z*x_size*y_size - y*x_size; /* index of x */

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
        if(z < bounding_box[2]){
            bounding_box[4] = z;
        }
        if(z > bounding_box[3]){
            bounding_box[5] = z;
        }
    }

	LOG_FUNC_END
}

template<class VectorBase>
void BGMGraphGrid3D<VectorBase>::saveVTK_bounding_box(std::string filename, int *bounding_box) const {
	LOG_FUNC_BEGIN

    //TODO

	LOG_FUNC_END
}


}
} /* end of namespace */

#endif
