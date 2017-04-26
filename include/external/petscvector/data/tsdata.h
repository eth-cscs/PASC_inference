#ifndef PASC_PETSCVECTOR_TSDATA_H
#define	PASC_PETSCVECTOR_TSDATA_H

#include "external/petscvector/algebra/vector/generalvector.h"
#include "general/data/tsdata.h"

namespace pascinference {

/* Maybe these classes are not defined yet */ 
/*namespace model {
	template<class VectorBase>
	class TSModel;
}
*/
using namespace model;

namespace data {

template<> TSData<PetscVector>::TSData(Decomposition<PetscVector> &new_decomposition, GeneralVector<PetscVector> *datavector_new, GeneralVector<PetscVector> *gammavector_new, GeneralVector<PetscVector> *thetavector_new);
template<> TSData<PetscVector>::TSData(Decomposition<PetscVector> &new_decomposition);
template<> TSData<PetscVector>::TSData(Decomposition<PetscVector> &new_decomposition, std::string filename);
template<> void TSData<PetscVector>::set_model(TSModel<PetscVector> &tsmodel);
template<> void TSData<PetscVector>::cutgamma() const;
template<> void TSData<PetscVector>::save_datavector(std::string filename) const;
template<> void TSData<PetscVector>::save_thetavector(std::string filename) const;
template<> void TSData<PetscVector>::print_thetavector(ConsoleOutput &output) const;
template<> std::string TSData<PetscVector>::print_thetavector() const;
template<> void TSData<PetscVector>::save_gammavector(std::string filename, int blocksize) const;
template<> void TSData<PetscVector>::printstats(ConsoleOutput &output, bool printdetails) const;
template<> void TSData<PetscVector>::scaledata(double a, double b);
template<> void TSData<PetscVector>::unscaledata(double a, double b);
template<> void TSData<PetscVector>::cutdata(double a, double b);
template<> void TSData<PetscVector>::shiftdata(double a);
template<> void TSData<PetscVector>::scaledata(double a, double b, double scale_min, double scale_max);
template<> void TSData<PetscVector>::load_gammavector(PetscVector &gamma0) const;
template<> void TSData<PetscVector>::load_gammavector(std::string filename) const;


}
} /* end namespace */

#endif
