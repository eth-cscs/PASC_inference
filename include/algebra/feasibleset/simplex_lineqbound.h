/** @file simplex_lineqbound.h
 *  @brief simplex feasible set defined by linear equality constraints and bound constraints
 *
 *  @author Lukas Pospisil
 */

#ifndef PASC_SIMPLEX_LINEQBOUNDFEASIBLESET_LOCAL_H
#define	PASC_SIMPLEX_LINEQBOUNDFEASIBLESET_LOCAL_H

#include "pascinference.h"

#ifndef USE_PETSC
 #error 'SIMPLEX_LINEQBOUNDFEASIBLESET is for PETSC'
#endif

typedef petscvector::PetscVector PetscVector;

namespace pascinference {
namespace algebra {

/** \class SimplexFeasibleSet_LinEqBound
 *  \brief simplex feasible set defined by linear equality constraints and bound constraints
 *
 *  Provides structure for manipulation with feasible set in form
 * \f[
 * \Omega = 
 *  \left\lbrace x \in R^{KT}: 
 *    \forall t = 0,\dots,T-1: \sum\limits_{k=0}^{K-1} x_{tK+k} = 1,  
 *    x \geq 0
 *  \right\rbrace
 *	\f] =
 *  \left\lbrace x \in R^{KT}: 
 *    Bx = c \wedge x \geq 0
 *  \right\rbrace
 * 
*/
template<class VectorBase>
class SimplexFeasibleSet_LinEqBound: public GeneralFeasibleSet<VectorBase> {
	private:
		
		int T; /**< number of global disjoint simplex subsets */
		int Tlocal; /**< number of local disjoint simplex subsets */
		int K; /**< size of each simplex subset */

		Mat B; /**< matrix of equality constraints Bx=c */
		Vec c; /**< vector of equality constraints Bx=c */
		Vec lb; /**< vector of lower bounds */
		
	public:
		/** @brief default constructor
		*/ 	
		SimplexFeasibleSet_LinEqBound(int T, int Tlocal, int K);
		
		/** @brief default destructor
		 */ 
		~SimplexFeasibleSet_LinEqBound();

		/** @brief print properties of this feasible set
		 * 
		 * @param output where to print
		 */ 
		void print(ConsoleOutput &output) const;

		/** @brief get name of this feasible set
		 */
		virtual std::string get_name() const;

		/** @brief compute projection onto feasible set
		 * 
		 * @param x point which will be projected
		 */		
		void project(GeneralVector<PetscVector> &x);
		
		Mat get_B() const;
		Vec get_c() const;
		Vec get_lb() const;
		
};


}
} /* end of namespace */

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {
namespace algebra {

/* constructor */
template<class VectorBase>
SimplexFeasibleSet_LinEqBound<VectorBase>::SimplexFeasibleSet_LinEqBound(int T, int Tlocal, int K){
	LOG_FUNC_BEGIN

	this->T = T;
	this->Tlocal = Tlocal;
	this->K = K;

	/* create global matrix from local blocks */
	TRYCXX( MatCreate(PETSC_COMM_WORLD, &B) );
	TRYCXX( MatSetSizes(B,this->Tlocal,this->K*this->Tlocal,this->T,this->K*this->T) );
	TRYCXX( MatSetFromOptions(B) );
	TRYCXX( MatMPIAIJSetPreallocation(B,this->K,NULL,this->K,NULL) );
	TRYCXX( MatSeqAIJSetPreallocation(B,this->K,NULL) );

	TRYCXX( MatAssemblyBegin(B,MAT_FLUSH_ASSEMBLY) );
	TRYCXX( MatAssemblyEnd(B,MAT_FLUSH_ASSEMBLY) );

	int row_begin;
	int col_begin;
	TRYCXX( MatGetOwnershipRange(B, &row_begin, NULL) );
	TRYCXX( MatGetOwnershipRangeColumn(B, &col_begin, NULL) );

	/* set local content */
	double value = 1.0;///(double)(this->K*this->T);
	for(int t=0; t<this->Tlocal;t++){
		for(int k=0;k<this->K;k++){
			TRYCXX( MatSetValue(B, row_begin + t, col_begin + t*this->K+k, value, INSERT_VALUES) );
		}
	}

	TRYCXX( MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY) );
	TRYCXX( MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY) );
	TRYCXX( PetscObjectSetName((PetscObject)B,"equality constraint") );

	/* create RHS vector of equality constraints */
	TRYCXX( MatCreateVecs(B, &lb, &c) );
	TRYCXX( VecSet(c,value) );
	TRYCXX( VecAssemblyBegin(c) );
	TRYCXX( VecAssemblyEnd(c) );
	TRYCXX( PetscObjectSetName((PetscObject)c,"RHS eq") );

	/* create RHS vector of lower bounds */
	TRYCXX( VecSet(lb,0.0) );
	TRYCXX( VecAssemblyBegin(lb) );
	TRYCXX( VecAssemblyEnd(lb) );
	TRYCXX( PetscObjectSetName((PetscObject)lb,"lower bounds") );

	LOG_FUNC_END
}

/* general destructor */
template<class VectorBase>
SimplexFeasibleSet_LinEqBound<VectorBase>::~SimplexFeasibleSet_LinEqBound(){
	LOG_FUNC_BEGIN
	
	
	LOG_FUNC_END	
}

/* print info about feasible set */
template<class VectorBase>
void SimplexFeasibleSet_LinEqBound<VectorBase>::print(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN
	
	output <<  this->get_name() << std::endl;
	
	/* give information about presence of the data */
	output <<  " - nmb of subsets:     " << T << std::endl;
	output <<  " - size of subset:     " << K << std::endl;


	LOG_FUNC_END
}

template<class VectorBase>
std::string SimplexFeasibleSet_LinEqBound<VectorBase>::get_name() const {
	return "SimplexFeasibleSet_LinEqBound";
}

template<>
void SimplexFeasibleSet_LinEqBound<PetscVector>::project(GeneralVector<PetscVector> &x) {
	LOG_FUNC_BEGIN
	
	//TODO: give error, the projection is not defined for this type of feasible set
	TRYCXX( PetscBarrier(NULL) );

	LOG_FUNC_END
}

template<>
Mat SimplexFeasibleSet_LinEqBound<PetscVector>::get_B() const {
	return this->B;
}

template<>
Vec SimplexFeasibleSet_LinEqBound<PetscVector>::get_c() const {
	return this->c;
}

template<>
Vec SimplexFeasibleSet_LinEqBound<PetscVector>::get_lb() const {
	return this->lb;
}



}
} /* end namespace */

#endif
