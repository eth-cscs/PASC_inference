#ifndef GENERALMATRIX_H
#define	GENERALMATRIX_H

#include <iostream>

/* for debugging, if >= 100, then print info about ach called function */
extern int DEBUG_MODE; 

namespace pascinference {

/* there will be a class with vectors... later */
template<class VectorBase> class Vector;
	
/* class for manipulation with A*x as one object, will be defined later */
template<class VectorBase> class GeneralMatrixRHS;
//template<class VectorBase, typename enable_if<!is_base_of<A, VectorBase>::value>::type* = nullptr> class GeneralMatrixRHS;

/* class for manipulation with general matrix */
template<class VectorBase>
class GeneralMatrix {
	protected:
		// TODO: timers & other general funny stuff 

	public:

		virtual void print(std::ostream &output) const {}; /* print matrix */
		virtual void matmult(VectorBase &y, const VectorBase &x) const {}; /* y = A*x */

		template<class VectorBase2>
		friend std::ostream &operator<<(std::ostream &output, const GeneralMatrix<VectorBase2> &matrix); /* cannot be virtual, therefore it call virtual print() */

};

/* print general matrix, call virtual print() */
template<class VectorBase>
std::ostream &operator<<(std::ostream &output, const GeneralMatrix<VectorBase> &matrix){
	if(DEBUG_MODE >= 100) std::cout << "(GeneralMatrixRHS)OPERATOR: <<" << std::endl;
	matrix.print(output);
	return output;
}

/* operator A*x (creates RHS to be proceeded into overloaded operator Vector = RHS */
template<class VectorBase>
GeneralMatrixRHS<VectorBase> operator*(const GeneralMatrix<VectorBase> &matrix, const Vector<VectorBase> &x){
	if(DEBUG_MODE >= 100) std::cout << "(GeneralMatrixRHS)OPERATOR: *" << std::endl;
	return GeneralMatrixRHS<VectorBase>(&matrix,&x);	
}


/* right hand-side vector of y=Ax - will be provided into overloaded operator  Vector = RHS */
template<class VectorBase>
class GeneralMatrixRHS{
	private:
		const GeneralMatrix<VectorBase> *matrix; /* pointer to general matrix */
		const Vector<VectorBase> *x; /* pointer to vector */
	public:

		/* constructor: create RHS from given pointers to matrix & vector */
		GeneralMatrixRHS(const GeneralMatrix<VectorBase> *newmatrix, const Vector<VectorBase> *newx){
			matrix = newmatrix;
			x = newx;
		}	

		void matmult(Vector<VectorBase> &y){ /* call multiplication function from matrix class to perform y = A*x */
			if(DEBUG_MODE >= 100) std::cout << "(GeneralMatrixRHS)FUNCTION: matmult" << std::endl;
			
			(*matrix).matmult(y, *x); /* call virtual function (of Matrix) for multiplication */
		}	

};




} /* end of namespace */


#endif
