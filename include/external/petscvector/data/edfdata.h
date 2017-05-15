#ifndef PASC_PETSCVECTOR_EDFDATA_H
#define	PASC_PETSCVECTOR_EDFDATA_H

#ifdef USE_BOOST
	#include <boost/filesystem.hpp>
#endif

//TODO: if USE_BOOST is not defined?


#include "external/petscvector/algebra/vector/generalvector.h"
#include "general/data/edfdata.h"
#include "external/petscvector/data/tsdata.h"

namespace pascinference {
namespace data {

template<> void EdfData<PetscVector>::edfRead(std::string filename, int max_record_nmb);
template<> void EdfData<PetscVector>::set_decomposition(Decomposition<PetscVector> &new_decomposition);
template<> EdfData<PetscVector>::EdfData(std::string filename_data, int max_record_nmb);
template<> EdfData<PetscVector>::~EdfData();
template<> void EdfData<PetscVector>::saveVTK(std::string filename) const;
template<> void EdfData<PetscVector>::saveVector(std::string filename, bool save_original) const;

}
} /* end namespace */

#endif
