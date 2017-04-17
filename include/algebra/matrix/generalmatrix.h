/** @file generalmatrix.h
 *  @brief class for manipulation with matrices
 *
 *  @author Lukas Pospisil
 */

#ifndef PASC_GENERALMATRIX_H
#define	PASC_GENERALMATRIX_H

#include "common/logging.h"
#include "common/consoleoutput.h"
#include "algebra/vector/generalvector.h"
#include "algebra/matrix/generalmatrixrhs.h"



namespace pascinference {
using namespace common;

namespace algebra {

/* there will be a class with vectors inside namespace vector... later defined, but I need it now */
template<class VectorBase> class GeneralVector;

/* class for manipulation with A*x as one object, will be defined later */
template<class VectorBase> class GeneralMatrixRHS;

/** \class GeneralMatrix
 *  \brief General class for manipulation with matrices.
 *
 *  Parent class for manipulation with matrices.
 *  All specific matrix implementations should be defined as inherited classes from this class.
 *	
 *  @todo add timers and other general stuff
*/
template<class VectorBase>
class GeneralMatrix {
	protected:
		// TODO: timers & other general funny stuff 

	public:

		/** @brief print informations of matrix
		 * 
		 *  Print the basic informations of matrix.
		 * 
		 * @param output where to print
		 */ 
		virtual void print(ConsoleOutput &output) const {};

		/** @brief print the content of matrix
		 * 
		 *  Print the content of matrix.
		 * 
		 * @param output where to print
		 */ 
		virtual void printcontent(ConsoleOutput &output) const {};

		/** @brief get the type name of matrix
		 */ 
		virtual std::string get_name() const {
			return "GeneralMatrix";
		}

		/** @brief perform matrix-vector multiplication
		 * 
		 *  Perform matrix-vector multiplication
		 *  \f[ y = Ax \f]
		 * 
		 * @param y output vector
		 * @param x input vector
		 */ 
		virtual void matmult(VectorBase &y, const VectorBase &x) const {};

		/** @brief print name of matrix to stream
		 * 
		 * @param output where to print
		 * @param matrix 
		 */ 
		template<class VectorBase2>
		friend ConsoleOutput &operator<<(ConsoleOutput &output, const GeneralMatrix<VectorBase2> &matrix);

		virtual double get_coeff() const {
			return 1.0;
		}

		virtual void set_coeff(double coeff){
			// TODO: write here something really funny
		}

};

/* print general matrix, call virtual print() */
template<class VectorBase>
ConsoleOutput &operator<<(ConsoleOutput &output, const GeneralMatrix<VectorBase> &matrix){
	output << matrix.get_name();
	return output;
}

/* operator A*x (creates RHS to be proceeded into overloaded operator Vector = RHS */
template<class VectorBase>
GeneralMatrixRHS<VectorBase> operator*(const GeneralMatrix<VectorBase> &matrix, const GeneralVector<VectorBase> &x){
	return GeneralMatrixRHS<VectorBase>(&matrix,&x);	
}


}
} /* end of namespace */

#ifdef USE_PETSC
 #include "external/petsc/algebra/matrix/generalmatrix.h"
#endif


#endif
