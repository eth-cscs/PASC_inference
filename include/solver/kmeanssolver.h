#ifndef PASC_KMEANSSOLVER_H
#define	PASC_KMEANSSOLVER_H

#ifndef USE_PETSCVECTOR
 #error 'KMEANSSOLVER is for PETSCVECTOR'
 typedef petscvector::PetscVector PetscVector;
#endif


#include "pascinference.h"

#include "solver/tssolver.h"
#include "data/kmeansdata.h"
#include "model/tsmodel.h"

//temp
#include "solver/diagsolver.h"
#include "solver/qpsolver.h"

#define KMEANSSOLVER_DEFAULT_MAXIT 1000
#define KMEANSSOLVER_DEFAULT_EPS 0.001
#define KMEANSSOLVER_DEFAULT_DEBUG_MODE 0

namespace pascinference {

/* KmeansSolver */ 
template<class VectorBase>
class KmeansSolver: public TSSolver<VectorBase> {
	protected:
		KmeansData<VectorBase> *tsdata; /**< tsdata on which the solver operates */

	public:
		KmeansSolver() : TSSolver<VectorBase> () {};
		KmeansSolver(KmeansData<VectorBase> &new_tsdata) : TSSolver<VectorBase> (new_tsdata) {}; 

		std::string get_name() const;

		KmeansData<VectorBase> *get_data() const;

};

} // end of namespace

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {

template<class VectorBase>
std::string KmeansSolver<VectorBase>::get_name() const {
	return "KmeansSolver<VectorBase>";
}


template<class VectorBase>
KmeansData<VectorBase> *KmeansSolver<VectorBase>::get_data() const {
	return tsdata;
}

} /* end namespace */

#endif
