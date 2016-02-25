#ifndef LAPLACEEXPLICITMATRIX_H
#define	LAPLACEEXPLICITMATRIX_H

/* for debugging, if >= 100, then print info about ach called function */
extern int DEBUG_MODE;

#include <iostream>
#include "algebra/generalmatrix.h" /* parent GeneralMatrix class */

typedef petscvector::PetscVector PetscVector;
typedef Mat PetscMatrix;

typedef minlin::threx::HostMatrix<double> MinlinHostMatrix;
typedef minlin::threx::HostVector<double> MinlinHostVector;

typedef minlin::threx::DeviceMatrix<double> MinlinDeviceMatrix;
typedef minlin::threx::DeviceVector<double> MinlinDeviceVector;



namespace pascinference {

/* laplace matrix */ 
template<class VectorType>
class LaplaceExplicitMatrix: public GeneralMatrix<VectorType> {
	private:
		/* Petsc stuff */ // TODO: if USE_PETSC
		PetscMatrix A_petsc;

		/* MINLIN stuff */ // TODO: if USE_MINLIN
		MinlinHostMatrix A_minlinhost;
		MinlinDeviceMatrix A_minlindevice;
		
	
	public:
		LaplaceExplicitMatrix(const VectorType &x); /* constructor from vector */
		~LaplaceExplicitMatrix(); /* destructor - destroy inner matrix */

		void print(std::ostream &output) const; /* print matrix */
		void matmult(VectorType &y, const VectorType &x) const; /* y = A*x */

};



/* -------------------------------- PETSC VECTOR -------------------------*/

/* Petsc: constructor from given right PetscVector */
template<>
LaplaceExplicitMatrix<PetscVector>::LaplaceExplicitMatrix(const PetscVector &x){
	/* init Petsc Vector */
	if(DEBUG_MODE >= 100) std::cout << "(LaplaceExplicitMatrix)CONSTRUCTOR: from PetscVector" << std::endl;

	int N, n;

	/* get informations from given vector */
	N = x.size();
	n = x.local_size();

	TRY( MatCreate(PETSC_COMM_WORLD,&A_petsc) );
	TRY( MatSetSizes(A_petsc,n,n,N,N) );
	TRY( MatSetFromOptions(A_petsc) ); 
//	TRY( MatSetType(A,MATMPIAIJ) ); 
	TRY( MatMPIAIJSetPreallocation(A_petsc,5,NULL,5,NULL) ); 
	TRY( MatSeqAIJSetPreallocation(A_petsc,5,NULL) );

	/* fill matrix */
	int row,col;
	double new_value;
	for(row=0;row<N;row++){
		for(col=row-1;col<=row+1;col++){
			/* first row */
			if(row == 0){
				new_value = 1;
				if(col > row){
					new_value = -1;
				}
			}
				
			/* last row */
			if(row == N-1){
				new_value = 1;
				if(col < row){
					new_value = -1;
				}
			}

			/* ordinary row */
			if(row > 0 && row < N-1){
				new_value = 2;
				if(col > row || col < row){
					new_value = -1;
				}
			}

			// TODO: only for testing string problem - regularization - remove this hotfix 
			if(true){
				if((row == 0 && col == 1) || (row == 1 && col == 0) || (row == N-2 && col == N-1) || (row == N-1 && col == N-2)){
					new_value = 0;
				}
			}

			/* set value */
			if(row >= 0 && row <= N-1 && col >=0 && col <= N-1){
				TRY( MatSetValue(A_petsc,row,col,new_value,INSERT_VALUES) );
			}
		}
	}
		
	/* assemble matrix */
	TRY( MatAssemblyBegin(A_petsc,MAT_FINAL_ASSEMBLY) );
	TRY( MatAssemblyEnd(A_petsc,MAT_FINAL_ASSEMBLY) );
	
}

/* Petsc: destructor - destroy the matrix */
template<>
LaplaceExplicitMatrix<PetscVector>::~LaplaceExplicitMatrix(){
	/* init Petsc Vector */
	if(DEBUG_MODE >= 100) std::cout << "(LaplaceExplicitMatrix)DESTRUCTOR" << std::endl;

	if(petscvector::PETSC_INITIALIZED){ /* maybe Petsc was already finalized and there is nothing to destroy */
		TRY( MatDestroy(&A_petsc) );
	}
}

/* print matrix */
template<>
void LaplaceExplicitMatrix<PetscVector>::print(std::ostream &output) const		
{
	if(DEBUG_MODE >= 100) std::cout << "(LaplaceExplicitMatrix)OPERATOR: << print" << std::endl;

	output << "Laplace matrix (sorry, 'only' MatView from Petsc follows):" << std::endl;
	output << "----------------------------------------------------------" << std::endl;
	
	TRY( MatView(A_petsc, PETSC_VIEWER_STDOUT_WORLD) );

	output << "----------------------------------------------------------" << std::endl;
}

/* Petsc: matrix-vector multiplication */
template<>
void LaplaceExplicitMatrix<PetscVector>::matmult(PetscVector &y, const PetscVector &x) const { 
	if(DEBUG_MODE >= 100) std::cout << "(LaplaceExplicitMatrix)FUNCTION: matmult" << std::endl;

	// TODO: maybe y is not initialized, who knows
	
	TRY( MatMult(A_petsc, x.get_vector(), y.get_vector()) ); // TODO: I dont want to use get_vector :( friend in PetscVector? and in MinLin?
}



/* -------------------------------- MINLIN HOST -------------------------*/

/* MinLinHost: constructor from given right HostVector<double> */
template<>
LaplaceExplicitMatrix<MinlinHostVector>::LaplaceExplicitMatrix(const MinlinHostVector &x){
	/* init Petsc Vector */
	if(DEBUG_MODE >= 100) std::cout << "(LaplaceExplicitMatrix)CONSTRUCTOR: from MinLin host" << std::endl;

	int N = x.size();

	/* get informations from given vector */
	MinlinHostMatrix A_new(N,N);

	int row,col;
	double new_value;
	for(row=0;row<N;row++){
		for(col=row-1;col<=row+1;col++){
			/* first row */
			if(row == 0){
				new_value = 1;
				if(col > row){
					new_value = -1;
				}
			}
				
			/* last row */
			if(row == N-1){
				new_value = 1;
				if(col < row){
					new_value = -1;
				}
			}

			/* ordinary row */
			if(row > 0 && row < N-1){
				new_value = 2;
				if(col > row || col < row){
					new_value = -1;
				}
			}

			// TODO: only for testing string problem - regularization - remove this hotfix 
			if(true){
				if((row == 0 && col == 1) || (row == 1 && col == 0) || (row == N-2 && col == N-1) || (row == N-1 && col == N-2)){
					new_value = 0;
				}
			}

			/* set value */
			if(row >= 0 && row <= N-1 && col >=0 && col <= N-1){
				A_new(row,col) = new_value;
			}
		}
	}
		
	A_minlinhost = A_new;

}


/* MinLinHost: destructor - destroy the matrix */
template<>
LaplaceExplicitMatrix<MinlinHostVector>::~LaplaceExplicitMatrix(){
	/* init Petsc Vector */
	if(DEBUG_MODE >= 100) std::cout << "(LaplaceExplicitMatrix)DESTRUCTOR" << std::endl;

	// TODO: how to destroy minlin matrix?
}

/* MinLinHost: print matrix */
template<>
void LaplaceExplicitMatrix<MinlinHostVector>::print(std::ostream &output) const		
{
	if(DEBUG_MODE >= 100) std::cout << "(LaplaceExplicitMatrix)OPERATOR: << print" << std::endl;
	output << A_minlinhost << std::endl;
}

/* MinLinHost: matrix-vector multiplication */
template<>
void LaplaceExplicitMatrix<MinlinHostVector>::matmult(MinlinHostVector &y, const MinlinHostVector &x) const { 
	if(DEBUG_MODE >= 100) std::cout << "(LaplaceExplicitMatrix)FUNCTION: matmult" << std::endl;

	y = A_minlinhost*x;	

}




} /* end of namespace */


#endif
