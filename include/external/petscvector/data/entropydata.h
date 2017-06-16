#ifndef PASC_PETSCVECTOR_ENTROPYDATA_H
#define	PASC_PETSCVECTOR_ENTROPYDATA_H

#include "external/petscvector/algebra/vector/generalvector.h"
#include "general/data/entropydata.h"

namespace pascinference {
namespace data {

/* external-specific stuff */
template<> class EntropyData<PetscVector>::ExternalContent {
	public:
		void prepare_matrix_D_recursion(double *values, int *idx, int top_i, int level, int Km, int xdim);

		void set_D_value(double *values, double value, int row, int col, int ncols);
		double get_D_value(double *values, int row, int col, int ncols) const;

};

template<> EntropyData<PetscVector>::EntropyData(int T, int xdim, int K, int Km);
template<> void EntropyData<PetscVector>::prepare_matrix_D();
template<> void EntropyData<PetscVector>::print_matrix_D(ConsoleOutput &output) const;

template<> EntropyData<PetscVector>::ExternalContent * EntropyData<PetscVector>::get_externalcontent() const;

}
} /* end namespace */

#endif
