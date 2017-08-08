#ifndef PASC_PETSCVECTOR_ENTROPYSOLVERDLIB_H
#define	PASC_PETSCVECTOR_ENTROPYSOLVERDLIB_H

#include "general/solver/entropysolverdlib.h"

/* include integration algorithms */
#include "external/petscvector/algebra/integration/entropyintegration.h"
#include "external/petscvector/algebra/integration/entropyintegrationdlib.h"
//#include "external/petscvector/algebra/integration/entropyintegrationcuba.h"

//#include "external/petscvector/solver/generalsolver.h"
#include "external/petscvector/data/entropydata.h"

#ifdef USE_DLIB

/* include Dlib stuff */
#include "dlib/matrix.h"
#include "dlib/numeric_constants.h"
#include "dlib/numerical_integration.h"
#include "dlib/optimization.h"

/* Dlib column vector */
typedef dlib::matrix<double,0,1> column_vector;
typedef dlib::matrix<double,1,0> row_vector;

namespace pascinference {
namespace solver {

/* external-specific stuff */
template<> class EntropySolverDlib<PetscVector>::ExternalContent {
	private:
		column_vector cLM;
		column_vector cgrad;
		dlib::matrix<double> chess;

		double cF;

		bool debug_print_content;
		bool debug_print_integration;

		EntropyIntegration<PetscVector> *entropyintegration;

	public:

		ExternalContent(EntropyIntegration<PetscVector> *entropyintegration, bool debug_print_content = false, bool debug_print_integration = false);
		double get_functions_obj(const column_vector& _LM, const column_vector& _Mom, double eps, int k, const dlib::matrix<double>& mom_powers);
		column_vector get_functions_grad(const column_vector& _LM, const column_vector& _Mom, int k, const dlib::matrix<double>& mom_powers);
		dlib::matrix<double> get_functions_hess(const column_vector& _LM, const column_vector& _Mom, int k, const dlib::matrix<double>& mom_powers);

		double get_F() const;

		double integration_time; /**< total integration time */
		double *computed_integrals; /**< array for storing computed integrals */

};

template<> EntropySolverDlib<PetscVector>::EntropySolverDlib(EntropyData<PetscVector> &new_entropydata);
template<> EntropySolverDlib<PetscVector>::~EntropySolverDlib();

template<> void EntropySolverDlib<PetscVector>::solve();

template<> double EntropySolverDlib<PetscVector>::get_integration_time() const;

template<> EntropySolverDlib<PetscVector>::ExternalContent * EntropySolverDlib<PetscVector>::get_externalcontent() const;

}
} /* end namespace */

#endif /* Dlib */
#endif
