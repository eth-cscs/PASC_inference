#ifndef PASC_EDFDATA_PETSC_H
#define	PASC_EDFDATA_PETSC_H

#include "algebra/vector/generalvector.h"
#include "data/edfdata.h"

namespace pascinference {
namespace data {

template<> void EdfData<PetscVector>::edfRead(std::string filename, int max_record_nmb);
template<> void EdfData<PetscVector>::set_decomposition(Decomposition &new_decomposition);
template<> EdfData<PetscVector>::EdfData(std::string filename_data, int max_record_nmb);
template<> EdfData<PetscVector>::~EdfData();
template<> void EdfData<PetscVector>::saveVTK(std::string filename) const;
template<> void EdfData<PetscVector>::saveVector(std::string filename, bool save_original) const;

}
} /* end namespace */

#endif
