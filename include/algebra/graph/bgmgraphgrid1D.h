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

#endif
