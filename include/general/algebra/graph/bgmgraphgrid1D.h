/** @file bgmgraphgrid1D.h
 *  @brief class for manipulation with 1D grid
 *
 *  Defines some basic functions for manipulaton with grids, especially for BLOCKGRAPHMATRIX.
 *
 *  @author Lukas Pospisil
 */
 
#ifndef PASC_COMMON_BGMGRAPHGRID1D_H
#define	PASC_COMMON_BGMGRAPHGRID1D_H

#include "general/algebra/graph/bgmgraphgrid1D.h"
#include "external/petscvector/algebra/graph/bgmgraph.h"

namespace pascinference {
namespace algebra {

/** \class BGMGraphGrid1D
 *  \brief Graph of one dimensional grid.
 *
*/
template<class VectorBase>
class BGMGraphGrid1D: public BGMGraph<VectorBase> {
	protected:
		int width;
		void process_grid_cuda();
		
	public:
		BGMGraphGrid1D(int width);
		BGMGraphGrid1D(std::string filename, int dim=2) : BGMGraph<VectorBase>(filename, dim) {};
		BGMGraphGrid1D(const double *coordinates_array, int n, int dim) : BGMGraph<VectorBase>(coordinates_array, n, dim) {};

		~BGMGraphGrid1D();
		
		virtual std::string get_name() const;
		virtual void process_grid();

		int get_width() const;
};

}
} /* end of namespace */

/* -------- IMPLEMENTATION ------ */
namespace pascinference {
namespace algebra {

template<class VectorBase>
BGMGraphGrid1D<VectorBase>::BGMGraphGrid1D(int width) : BGMGraph<VectorBase>(){
	LOG_FUNC_BEGIN

	this->width = width;

	this->dim = 2;
	this->n = width;
	
	//TODO

	LOG_FUNC_END
}

template<class VectorBase>
BGMGraphGrid1D<VectorBase>::~BGMGraphGrid1D(){
	
}

template<class VectorBase>
std::string BGMGraphGrid1D<VectorBase>::get_name() const {
	return "BGMGraphGrid1D";
}

template<class VectorBase>
void BGMGraphGrid1D<VectorBase>::process_grid(){
	this->threshold = 1.1;
	this->m = width-1;
	this->m_max = 2;

	/* prepare array for number of neighbors */
	this->neighbor_nmbs = (int*)malloc(this->n*sizeof(int));
	this->neighbor_ids = (int**)malloc(this->n*sizeof(int*));

//	#pragma omp parallel for
	for(int i=0;i<width;i++){
		int idx = i;

		/* compute number of neighbors */
		int nmb = 0;
		if(i>0){
			nmb+=1;				
		}
		if(i<width-1){
			nmb+=1;				
		}
		this->neighbor_nmbs[idx] = nmb;
		this->neighbor_ids[idx] = (int*)malloc(this->neighbor_nmbs[idx]*sizeof(int));
			
		/* fill neighbors */
		nmb = 0;
		if(i>0){ /* left */
			this->neighbor_ids[idx][nmb] = idx-1;
			nmb++;
		}
		if(i<width-1){ /* right */
			this->neighbor_ids[idx][nmb] = idx+1;
			nmb++;
		}
	}

	#ifdef USE_CUDA
		process_grid_cuda();
	#endif
	
	this->processed = true;
}

template<class VectorBase>
int BGMGraphGrid1D<VectorBase>::get_width() const {
	return this->width;
}



}
} /* end of namespace */


#endif
