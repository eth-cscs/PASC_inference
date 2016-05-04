#ifndef BLOCKDIAGLAPLACEVECTORMATRIX_H
#define	BLOCKDIAGLAPLACEVECTORMATRIX_H

/* for debugging, if >= 100, then print info about ach called function */
extern int DEBUG_MODE;

#include <iostream>
#include "algebra.h" /* parent GeneralMatrix class */

#ifdef USE_PETSCVECTOR
	typedef petscvector::PetscVector PetscVector;
#endif

#ifdef USE_MINLIN
	typedef minlin::threx::HostVector<double> MinlinHostVector;
	typedef minlin::threx::DeviceVector<double> MinlinDeviceVector;
#endif

namespace pascinference {

/* laplace matrix */ 
template<class VectorBase>
class BlockDiagLaplaceVectorMatrix: public GeneralMatrix<VectorBase> {
	private:
		int K; /* number of block */
		int N; /* size of the matrix */
		int T; /* size of each block N = K*T */
		double alpha; /* scale of whole matrix alpha*A */
	
	public:
		BlockDiagLaplaceVectorMatrix(const VectorBase &x, int K, double alpha = 1.0); /* constructor from vector and number of blocks */

		~BlockDiagLaplaceVectorMatrix(); /* destructor - destroy inner matrix */

		void print(std::ostream &output) const; /* print matrix */
		std::string get_name() const;

		void matmult(VectorBase &y, const VectorBase &x) const; /* y = A*x */

};

template<class VectorBase>
std::string BlockDiagLaplaceVectorMatrix<VectorBase>::get_name() const {
	return "BlockDiagLaplaceVectorMatrix";
}

template<class VectorBase>
BlockDiagLaplaceVectorMatrix<VectorBase>::BlockDiagLaplaceVectorMatrix(const VectorBase &x, int newK, double newalpha){
	if(DEBUG_MODE >= 100) coutMaster << "(BlockDiagLaplaceVectorMatrix)CONSTRUCTOR: from PetscVector" << std::endl;

	/* get informations from given vector */
	alpha = newalpha;
	K = newK; /* number of blocks */
	N = x.size(); /* length of whole matrix N = K*T */
	T = N/(double)K; /* size of each block */
	
}

template<class VectorBase>
BlockDiagLaplaceVectorMatrix<VectorBase>::~BlockDiagLaplaceVectorMatrix(){
	if(DEBUG_MODE >= 100) coutMaster << "(BlockDiagLaplaceVectorMatrix)DESTRUCTOR" << std::endl;

}

template<class VectorBase>
void BlockDiagLaplaceVectorMatrix<VectorBase>::print(std::ostream &output) const		
{
	if(DEBUG_MODE >= 100) coutMaster << "(BlockDiagLaplaceVectorMatrix)OPERATOR: << print" << std::endl;

	output << "Laplace matrix" << std::endl;
	output << " - K = " << K << std::endl;
	output << " - T = " << T << std::endl;
	output << " - N = " << N << std::endl;

}

template<class VectorBase>
void BlockDiagLaplaceVectorMatrix<VectorBase>::matmult(VectorBase &y, const VectorBase &x) const { 
	if(DEBUG_MODE >= 100) coutMaster << "(BlockDiagLaplaceVectorMatrix)FUNCTION: matmult" << std::endl;

	// TODO: maybe y is not initialized, who knows

	y(1,N-2) = alpha*2*x(1,N-2) - alpha*x(0,N-3) - alpha*x(2,N-1);
	
	/* begin and end of each block */
	for(int k=0;k<K;k++){
		y(k*T) = alpha*x(k*T) - alpha*x(k*T+1);
		y((k+1)*T-1) = alpha*x((k+1)*T-1) - alpha*x((k+1)*T-2);
	}
	
}

#ifdef USE_PETSCVECTOR

typedef petscvector::PetscVector PetscVector;

template<>
void BlockDiagLaplaceVectorMatrix<PetscVector>::matmult(PetscVector &y, const PetscVector &x) const { 
	if(DEBUG_MODE >= 100) coutMaster << "(BlockDiagLaplaceVectorMatrix)FUNCTION: matmult" << std::endl;

	// TODO: maybe y is not initialized, who knows
	double *y_arr;
	const double *x_arr;
	TRY( VecGetArray(y.get_vector(),&y_arr) );
	TRY( VecGetArrayRead(x.get_vector(),&x_arr) );

	int k,t,id_row;
	for(k=0;k<K;k++){
		for(t=0;t<T;t++){
			id_row = k*T+t;

			/* first row */
			if(t == 0){
				y_arr[id_row] = alpha*x_arr[id_row] - alpha*x_arr[id_row+1];
			}
			/* common row */
			if(t > 0 && t < T-1){
				y_arr[id_row] = -alpha*x_arr[id_row-1] + 2.0*alpha*x_arr[id_row] - alpha*x_arr[id_row+1];
			}
			/* last row */
			if(t == T-1){
				y_arr[id_row] = -alpha*x_arr[id_row-1] + alpha*x_arr[id_row];
			}
		}
	}

	TRY( VecRestoreArray(y.get_vector(),&y_arr) );
	TRY( VecRestoreArrayRead(x.get_vector(),&x_arr) );
	
}

#endif



} /* end of namespace */


#endif
