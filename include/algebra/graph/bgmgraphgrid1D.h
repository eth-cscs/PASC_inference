/** @file bgmgraphgrid1D.h
 *  @brief class for manipulation with 1D grid
 *
 *  Defines some basic functions for manipulaton with grids, especially for BLOCKGRAPHMATRIX.
 *
 *  @author Lukas Pospisil
 */
 
#ifndef PASC_COMMON_BGMGRAPHGRID1D_H
#define	PASC_COMMON_BGMGRAPHGRID1D_H

#include "algebra/graph/bgmgraph.h"

namespace pascinference {
namespace algebra {

/** \class BGMGraphGrid1D
 *  \brief Graph of one dimensional grid.
 *
*/
class BGMGraphGrid1D: public BGMGraph {
	protected:
		int width;
		void process_grid_cuda();
		
	public:
		BGMGraphGrid1D(int width);
		BGMGraphGrid1D(std::string filename, int dim=2) : BGMGraph(filename, dim) {};
		BGMGraphGrid1D(const double *coordinates_array, int n, int dim) : BGMGraph(coordinates_array, n, dim) {};

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

BGMGraphGrid1D::BGMGraphGrid1D(int width) : BGMGraph(){
	LOG_FUNC_BEGIN

	this->width = width;

	this->dim = 2;
	this->n = width;
	
	/* fill coordinates */
	Vec coordinates_Vec;
	TRYCXX( VecCreateSeq(PETSC_COMM_SELF, this->n*this->dim, &coordinates_Vec) );
	
	double *coordinates_arr;
	TRYCXX( VecGetArray(coordinates_Vec, &coordinates_arr) );
	for(int i=0;i<width;i++){
		coordinates_arr[i] = i;
		coordinates_arr[i + this->n] = 0;
	}
	TRYCXX( VecRestoreArray(coordinates_Vec, &coordinates_arr) );
	
	this->coordinates = new GeneralVector<PetscVector>(coordinates_Vec);

	this->threshold = -1;
	processed = false;

	LOG_FUNC_END
}

BGMGraphGrid1D::~BGMGraphGrid1D(){
	
}

std::string BGMGraphGrid1D::get_name() const {
	return "BGMGraphGrid1D";
}

void BGMGraphGrid1D::process_grid(){
	this->threshold = 1.1;
	this->m = width-1;
	this->m_max = 2;

	/* prepare array for number of neighbors */
	neighbor_nmbs = (int*)malloc(n*sizeof(int));
	neighbor_ids = (int**)malloc(n*sizeof(int*));

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
		neighbor_nmbs[idx] = nmb;
		neighbor_ids[idx] = (int*)malloc(neighbor_nmbs[idx]*sizeof(int));
			
		/* fill neighbors */
		nmb = 0;
		if(i>0){ /* left */
			neighbor_ids[idx][nmb] = idx-1;
			nmb++;
		}
		if(i<width-1){ /* right */
			neighbor_ids[idx][nmb] = idx+1;
			nmb++;
		}
	}

	#ifdef USE_CUDA
		process_grid_cuda();
	#endif
	
	processed = true;
}

int BGMGraphGrid1D::get_width() const {
	return this->width;
}



}
} /* end of namespace */


#endif
