/** @file bgmgraphgrid2D.h
 *  @brief class for manipulation with 2D grid
 *
 *  Defines some basic functions for manipulaton with grids, especially for BLOCKGRAPHMATRIX.
 *
 *  @author Lukas Pospisil
 */
 
#ifndef PASC_COMMON_BGMGRAPHGRID2D_H
#define	PASC_COMMON_BGMGRAPHGRID2D_H

#include "algebra/graph/bgmgraph.h"

namespace pascinference {
namespace algebra {

/** \class BGMGraphGrid2D
 *  \brief Graph of two dimensional grid.
 *
 *  Could be used for faster and simplier manipulation with image graph.
 * 
*/
class BGMGraphGrid2D: public BGMGraph {
	protected:
		int width; /**< dimension of grid */
		int height; /**< dimension of grid */
		
		void process_grid_cuda();
	public:
	
		BGMGraphGrid2D(int width, int height);
		BGMGraphGrid2D(std::string filename, int dim=2) : BGMGraph(filename, dim) {};
		BGMGraphGrid2D(const double *coordinates_array, int n, int dim) : BGMGraph(coordinates_array, n, dim) {};

		~BGMGraphGrid2D();
		
		virtual std::string get_name() const;
		virtual void process_grid();
		
		int get_width() const;
		int get_height() const;

		void decompose(BGMGraphGrid2D *finer_grid, int *bounding_box1, int *bounding_box2);
};


}
} /* end of namespace */

#endif
