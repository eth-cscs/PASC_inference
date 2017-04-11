/** @file generalmatrixrhs.h
 *  @brief class for manipulation with right hand side vector A*x
 *
 *  @author Lukas Pospisil
 */

#ifndef PASC_GENERALMATRIXRHS_H
#define	PASC_GENERALMATRIXRHS_H

#include "common/logging.h"
#include "common/consoleoutput.h"

namespace pascinference {
using namespace common;

namespace algebra {

template<class VectorBase> class GeneralVector;
template<class VectorBase> class GeneralMatrix;


/** \class GeneralMatrixRHS
 *  \brief General class for manipulation with matrices.
 *
 *  Right hand-side vector of y=A*x - will be provided into overloaded operator y = RHS
 *	
*/
template<class VectorBase>
class GeneralMatrixRHS{
	private:
		const GeneralMatrix<VectorBase> *matrix; /**< pointer to general matrix */
		const GeneralVector<VectorBase> *x; /**< pointer to input vector */
	public:

		/** @brief constructor from matrix and vector
		 * 
		 * create RHS from given pointers to matrix & vector, rhs = A*x
		 * 
		 * @param newmatrix pointer to input matrix on right side of y=A*x
		 * @param newx pointer to input vector on right side of y=A*x
		 */ 
		GeneralMatrixRHS(const GeneralMatrix<VectorBase> *newmatrix, const GeneralVector<VectorBase> *newx){
			matrix = newmatrix;
			x = newx;
		}	

		/** @brief multiplicate and store result
		 * 
		 * Call multiplication function from matrix class to perform y = A*x 
		 * 
		 * @param y result vector y=A*x
		 */ 
		void matmult(GeneralVector<VectorBase> &y){ 
//			LOG_FUNC_BEGIN
			
//			(*matrix).matmult(y, *x); /* call virtual function (of Matrix) for multiplication */

	//		LOG_FUNC_END
		}	

};



}
} /* end of namespace */


#endif
