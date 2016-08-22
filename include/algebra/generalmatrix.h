/** @file generalmatrix.h
 *  @brief class for manipulation with matrices
 *
 *  @author Lukas Pospisil
 */

#ifndef PASC_GENERALMATRIX_H
#define	PASC_GENERALMATRIX_H

/* for debugging, if >= 100, then print info about ach called function */
extern int DEBUG_MODE; 

namespace pascinference {

/* there will be a class with vectors... later */
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
		virtual std::string get_name() const;

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

		virtual double get_coeff() const;
		virtual void set_coeff(double coeff);

};

/* print general matrix, call virtual print() */
template<class VectorBase>
ConsoleOutput &operator<<(ConsoleOutput &output, const GeneralMatrix<VectorBase> &matrix){
	if(DEBUG_MODE >= 100) coutMaster << "(GeneralMatrixRHS)OPERATOR: <<" << std::endl;
	output << matrix.get_name();
	return output;
}

template<class VectorBase>
std::string GeneralMatrix<VectorBase>::get_name() const {
	return "GeneralMatrix";
}


/* operator A*x (creates RHS to be proceeded into overloaded operator Vector = RHS */
template<class VectorBase>
GeneralMatrixRHS<VectorBase> operator*(const GeneralMatrix<VectorBase> &matrix, const GeneralVector<VectorBase> &x){
	if(DEBUG_MODE >= 100) coutMaster << "(GeneralMatrixRHS)OPERATOR: *" << std::endl;
	return GeneralMatrixRHS<VectorBase>(&matrix,&x);	
}

template<class VectorBase>
double GeneralMatrix<VectorBase>::get_coeff() const {
	return 1.0;
}

template<class VectorBase>
void GeneralMatrix<VectorBase>::set_coeff(double coeff) {
	// TODO: write here something really funny
}


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
		 * @param newx pointer to input vector on right side
		 */ 
		void matmult(GeneralVector<VectorBase> &y){ 
			if(DEBUG_MODE >= 100) coutMaster << "(GeneralMatrixRHS)FUNCTION: matmult" << std::endl;
			
			(*matrix).matmult(y, *x); /* call virtual function (of Matrix) for multiplication */
		}	

};




} /* end of namespace */


#endif
