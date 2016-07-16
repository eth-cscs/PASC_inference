#ifndef PASC_EDFSOLVER_H
#define	PASC_EDFSOLVER_H

#ifndef USE_PETSCVECTOR
 #error 'EDF is for PETSCVECTOR'
 typedef petscvector::PetscVector PetscVector;
#endif


#include "pascinference.h"

#include "solver/tssolver.h"
#include "data/edfdata.h"
#include "model/tsmodel.h"

#define EDFSOLVER_DEFAULT_MAXIT 1000
#define EDFSOLVER_DEFAULT_EPS 0.001
#define EDFSOLVER_DEFAULT_DEBUG_MODE 0

namespace pascinference {

/* EdfSolver */ 
template<class VectorBase>
class EdfSolver: public TSSolver<VectorBase> {
	protected:
		EdfData<VectorBase> *tsdata; /**< tsdata on which the solver operates */

	public:
		EdfSolver() : TSSolver<VectorBase> () {};
		EdfSolver(EdfData<VectorBase> &new_tsdata, int annealing = 1) : TSSolver<VectorBase> (new_tsdata, annealing) {}; 

		std::string get_name() const;

		EdfData<VectorBase> *get_data() const;

};

} // end of namespace

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {

template<class VectorBase>
std::string EdfSolver<VectorBase>::get_name() const {
	return "EdfSolver<VectorBase>";
}


template<class VectorBase>
EdfData<VectorBase> *EdfSolver<VectorBase>::get_data() const {
	return tsdata;
}

} /* end namespace */

#endif
