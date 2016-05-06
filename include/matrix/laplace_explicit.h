#ifndef LAPLACEEXPLICITMATRIX_H
#define	LAPLACEEXPLICITMATRIX_H

/* for debugging, if >= 100, then print info about ach called function */
extern int DEBUG_MODE;

#include <iostream>
#include "algebra.h" /* parent GeneralMatrix class */

#ifdef USE_PETSCVECTOR
typedef petscvector::PetscVector PetscVector;
typedef Mat PetscMatrix;
#endif

#ifdef USE_MINLIN
typedef minlin::threx::HostMatrix<double> MinlinHostMatrix;
typedef minlin::threx::HostVector<double> MinlinHostVector;

typedef minlin::threx::DeviceMatrix<double> MinlinDeviceMatrix;
typedef minlin::threx::DeviceVector<double> MinlinDeviceVector;
#endif


namespace pascinference {

/* laplace matrix */ 
template<class VectorBase>
class LaplaceExplicitMatrix: public GeneralMatrix<VectorBase> {
	private:
		#ifdef USE_PETSC
		/* Petsc stuff */
		PetscMatrix A_petsc;
		#endif

		#ifdef USE_MINLIN
		/* MINLIN stuff */ 
		MinlinHostMatrix A_minlinhost;
		MinlinDeviceMatrix A_minlindevice;
		#endif
	
	public:
		LaplaceExplicitMatrix(const VectorBase &x); /* constructor from vector */
		~LaplaceExplicitMatrix(); /* destructor - destroy inner matrix */

		void print(ConsoleOutput &output) const; /* print matrix */
		std::string get_name() const;
		
		void matmult(VectorBase &y, const VectorBase &x) const; /* y = A*x */

};

template<class VectorBase>
std::string LaplaceExplicitMatrix<VectorBase>::get_name() const {
	return "LaplaceExplicitMatrix";
}


/* -------------------------------- PETSC VECTOR -------------------------*/

#ifdef USE_PETSCVECTOR

/* Petsc: constructor from given right PetscVector */
template<>
LaplaceExplicitMatrix<PetscVector>::LaplaceExplicitMatrix(const PetscVector &x){
	if(DEBUG_MODE >= 100) coutMaster << "(LaplaceExplicitMatrix)CONSTRUCTOR: from PetscVector" << std::endl;

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
	if(DEBUG_MODE >= 100) coutMaster << "(LaplaceExplicitMatrix)DESTRUCTOR" << std::endl;

	if(petscvector::PETSC_INITIALIZED){ /* maybe Petsc was already finalized and there is nothing to destroy */
		TRY( MatDestroy(&A_petsc) );
	}
}

/* print matrix */
template<>
void LaplaceExplicitMatrix<PetscVector>::print(ConsoleOutput &output) const		
{
	if(DEBUG_MODE >= 100) coutMaster << "(LaplaceExplicitMatrix)OPERATOR: << print" << std::endl;

	output << "Laplace matrix (sorry, 'only' MatView from Petsc follows):" << std::endl;
	output << "----------------------------------------------------------" << std::endl;
	
	TRY( MatView(A_petsc, PETSC_VIEWER_STDOUT_WORLD) );

	output << "----------------------------------------------------------" << std::endl;
}

/* Petsc: matrix-vector multiplication */
template<>
void LaplaceExplicitMatrix<PetscVector>::matmult(PetscVector &y, const PetscVector &x) const { 
	if(DEBUG_MODE >= 100) coutMaster << "(LaplaceExplicitMatrix)FUNCTION: matmult" << std::endl;

	// TODO: maybe y is not initialized, who knows
	
	TRY( MatMult(A_petsc, x.get_vector(), y.get_vector()) ); // TODO: I dont want to use get_vector :( friend in PetscVector? and in MinLin?
}


#endif /* ifdef USE_PETSC */

/* -------------------------------- MINLIN HOST -------------------------*/

#ifdef USE_MINLIN

/* MinLinHost: constructor from given right HostVector<double> */
template<>
LaplaceExplicitMatrix<MinlinHostVector>::LaplaceExplicitMatrix(const MinlinHostVector &x){
	/* init Petsc Vector */
	if(DEBUG_MODE >= 100) coutMaster << "(LaplaceExplicitMatrix)CONSTRUCTOR: from MinLin host" << std::endl;

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
	if(DEBUG_MODE >= 100) coutMaster << "(LaplaceExplicitMatrix)DESTRUCTOR" << std::endl;

	// TODO: how to destroy minlin matrix?
}

/* MinLinHost: print matrix */
template<>
void LaplaceExplicitMatrix<MinlinHostVector>::print(ConsoleOutput &output) const		
{
	if(DEBUG_MODE >= 100) coutMaster << "(LaplaceExplicitMatrix)OPERATOR: << print" << std::endl;
	output << A_minlinhost << std::endl;
}

/* MinLinHost: matrix-vector multiplication */
template<>
void LaplaceExplicitMatrix<MinlinHostVector>::matmult(MinlinHostVector &y, const MinlinHostVector &x) const { 
	if(DEBUG_MODE >= 100) coutMaster << "(LaplaceExplicitMatrix)FUNCTION: matmult" << std::endl;

	y = A_minlinhost*x;	

}

#endif /* ifdef USE_MINLIN */

/* -------------------------------- MINLIN DEVICE -------------------------*/

#ifdef USE_MINLIN

/* MinLinDevice: constructor from given right DeviceVector<double> */
template<>
LaplaceExplicitMatrix<MinlinDeviceVector>::LaplaceExplicitMatrix(const MinlinDeviceVector &x){
	/* init Petsc Vector */
	if(DEBUG_MODE >= 100) coutMaster << "(LaplaceExplicitMatrix)CONSTRUCTOR: from MinLin host" << std::endl;

	int N = x.size();

	/* get informations from given vector */
	MinlinDeviceMatrix A_new(N,N);

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

			/* set value */
			if(row >= 0 && row <= N-1 && col >=0 && col <= N-1){
				A_new(row,col) = new_value;
			}
		}
	}
		
	A_minlindevice = A_new;

}


/* MinLinDevice: destructor - destroy the matrix */
template<>
LaplaceExplicitMatrix<MinlinDeviceVector>::~LaplaceExplicitMatrix(){
	/* init Petsc Vector */
	if(DEBUG_MODE >= 100) coutMaster << "(LaplaceExplicitMatrix)DESTRUCTOR" << std::endl;

	// TODO: how to destroy minlin matrix?
}

/* MinLinDevice: print matrix */
template<>
void LaplaceExplicitMatrix<MinlinDeviceVector>::print(ConsoleOutput &output) const		
{
	if(DEBUG_MODE >= 100) coutMaster << "(LaplaceExplicitMatrix)OPERATOR: << print" << std::endl;
	output << A_minlindevice << std::endl;
}

/* MinLinDevice: matrix-vector multiplication */
template<>
void LaplaceExplicitMatrix<MinlinDeviceVector>::matmult(MinlinDeviceVector &y, const MinlinDeviceVector &x) const { 
	if(DEBUG_MODE >= 100) coutMaster << "(LaplaceExplicitMatrix)FUNCTION: matmult" << std::endl;

	y = A_minlindevice*x;	

}

#endif /* ifdef USE_MINLIN */ 

} /* end of namespace */


#endif
