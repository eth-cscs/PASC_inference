#ifndef BLOCKDIAGLAPLACEEXPLICITMATRIX_H
#define	BLOCKDIAGLAPLACEEXPLICITMATRIX_H

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
class BlockDiagLaplaceExplicitMatrix: public GeneralMatrix<VectorBase> {
	private:
		#ifdef USE_PETSCVECTOR
			/* Petsc stuff */ 
			PetscMatrix A_petsc;
		#endif

		#ifdef USE_MINLIN
			/* MINLIN stuff */ 
			MinlinHostMatrix A_minlinhost;
			MinlinDeviceMatrix A_minlindevice;
		#endif
		
		int K; /* number of block */
		double alpha; /* scale of whole matrix alpha*A */
	
	public:
		BlockDiagLaplaceExplicitMatrix(const VectorBase &x, int K, double alpha = 1.0); /* constructor from vector and number of blocks */

		~BlockDiagLaplaceExplicitMatrix(); /* destructor - destroy inner matrix */

		void print(std::ostream &output) const; /* print matrix */
		void matmult(VectorBase &y, const VectorBase &x) const; /* y = A*x */

};



/* -------------------------------- PETSC VECTOR -------------------------*/

#ifdef USE_PETSCVECTOR

/* Petsc: constructor from given right PetscVector */
template<>
BlockDiagLaplaceExplicitMatrix<PetscVector>::BlockDiagLaplaceExplicitMatrix(const PetscVector &x, int newK, double newalpha){
	if(DEBUG_MODE >= 100) std::cout << "(BlockDiagLaplaceExplicitMatrix)CONSTRUCTOR: from PetscVector" << std::endl;

	int N, n, T;

	/* get informations from given vector */
	alpha = newalpha;
	K = newK; /* number of blocks */
	N = x.size(); /* length of whole matrix N = K*T */
	n = x.local_size();
	T = N/(double)K; /* size of each block */

	std::cout << "N = " << N << std::endl;
	std::cout << "K = " << K << std::endl;
	std::cout << "T = " << T << std::endl;

	TRY( MatCreate(PETSC_COMM_WORLD,&A_petsc) );
	TRY( MatSetSizes(A_petsc,n,n,N,N) );
	TRY( MatSetFromOptions(A_petsc) ); 
//	TRY( MatSetType(A,MATMPIAIJ) ); 
	TRY( MatMPIAIJSetPreallocation(A_petsc,5,NULL,5,NULL) ); 
	TRY( MatSeqAIJSetPreallocation(A_petsc,5,NULL) );

	/* fill matrix */
	int row,col; /* coordinates in block */
	int k, row_global, col_global; /* coordinates in global matrix */
	double new_value;
	
	for(k=0;k<K;k++){
		for(row=0;row<T;row++){
			for(col=row-1;col<=row+1;col++){
				/* first row */
				if(row == 0){
					new_value = 1;
					if(col > row){
						new_value = -1;
					}
				}
				
				/* last row */
				if(row == T-1){
					new_value = 1;
					if(col < row){
						new_value = -1;
					}
				}

				/* ordinary row */
				if(row > 0 && row < T-1){
					new_value = 2;
					if(col > row || col < row){
						new_value = -1;
					}
				}

				/* set value */
				row_global = k*T + row;
				col_global = k*T + col;
				if(row >= 0 && row <= T-1 && col >=0 && col <= T-1){ /* I am in my block */
					if(row_global >= 0 && row_global <= N-1 && col_global >=0 && col_global <= N-1){ /* I am in global matrix */
						TRY( MatSetValue(A_petsc,row_global,col_global,alpha*new_value,INSERT_VALUES) );
					}
				}
				
			} /* endfor col */
		} /* endfor row */
	} /* endfor k */

	/* assemble matrix */
	TRY( MatAssemblyBegin(A_petsc,MAT_FINAL_ASSEMBLY) );
	TRY( MatAssemblyEnd(A_petsc,MAT_FINAL_ASSEMBLY) );
	
}

/* Petsc: destructor - destroy the matrix */
template<>
BlockDiagLaplaceExplicitMatrix<PetscVector>::~BlockDiagLaplaceExplicitMatrix(){
	if(DEBUG_MODE >= 100) std::cout << "(BlockDiagLaplaceExplicitMatrix)DESTRUCTOR" << std::endl;

	if(petscvector::PETSC_INITIALIZED){ /* maybe Petsc was already finalized and there is nothing to destroy */
		TRY( MatDestroy(&A_petsc) );
	}
}

/* print matrix */
template<>
void BlockDiagLaplaceExplicitMatrix<PetscVector>::print(std::ostream &output) const		
{
	if(DEBUG_MODE >= 100) std::cout << "(BlockDiagLaplaceExplicitMatrix)OPERATOR: << print" << std::endl;

	output << "Laplace matrix (sorry, 'only' MatView from Petsc follows):" << std::endl;
	output << "----------------------------------------------------------" << std::endl;
	
	TRY( MatView(A_petsc, PETSC_VIEWER_STDOUT_WORLD) );

	output << "----------------------------------------------------------" << std::endl;
}

/* Petsc: matrix-vector multiplication */
template<>
void BlockDiagLaplaceExplicitMatrix<PetscVector>::matmult(PetscVector &y, const PetscVector &x) const { 
	if(DEBUG_MODE >= 100) std::cout << "(BlockDiagLaplaceExplicitMatrix)FUNCTION: matmult" << std::endl;

	// TODO: maybe y is not initialized, who knows
	
	TRY( MatMult(A_petsc, x.get_vector(), y.get_vector()) ); // TODO: I dont want to use get_vector :( friend in GlobalVector?
}

#endif

/* -------------------------------- MINLIN HOST -------------------------*/

#ifdef USE_MINLIN

/* MinLinHost: constructor from given right HostVector<double> */
template<>
BlockDiagLaplaceExplicitMatrix<MinlinHostVector>::BlockDiagLaplaceExplicitMatrix(const MinlinHostVector &x, int newK, double newalpha){
	/* init Petsc Vector */
	if(DEBUG_MODE >= 100) std::cout << "(BlockDiagLaplaceExplicitMatrix)CONSTRUCTOR: from MinLin host" << std::endl;

	int N = x.size();
	alpha = newalpha;
	K = newK; /* number of blocks */
	int T = N/(double)K; /* size of each block */

	/* get informations from given vector */
	MinlinHostMatrix A_new(N,N);

	int row,col; /* coordinates in block */
	int k, row_global, col_global; /* coordinates in global matrix */
	double new_value;

	for(k=0;k<K;k++){
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
				if(row == T-1){
					new_value = 1;
					if(col < row){
						new_value = -1;
					}
				}

				/* ordinary row */
				if(row > 0 && row < T-1){
					new_value = 2;
					if(col > row || col < row){
						new_value = -1;
					}
				}

				/* set value */
				row_global = k*T + row;
				col_global = k*T + col;
				if(row >= 0 && row <= T-1 && col >=0 && col <= T-1){ /* I am in my block */
					if(row_global >= 0 && row_global <= N-1 && col_global >=0 && col_global <= N-1){ /* I am in global matrix */
						A_new(row_global,col_global) = alpha*new_value;
					}
				}
				
			} /* endfor col */
		} /* endfor row */
	} /* endfor k */
		
	A_minlinhost = A_new;

}


/* MinLinHost: destructor - destroy the matrix */
template<>
BlockDiagLaplaceExplicitMatrix<MinlinHostVector>::~BlockDiagLaplaceExplicitMatrix(){
	/* init Petsc Vector */
	if(DEBUG_MODE >= 100) std::cout << "(BlockDiagLaplaceExplicitMatrix)DESTRUCTOR" << std::endl;

	// TODO: how to destroy minlin matrix?
}

/* MinLinHost: print matrix */
template<>
void BlockDiagLaplaceExplicitMatrix<MinlinHostVector>::print(std::ostream &output) const		
{
	if(DEBUG_MODE >= 100) std::cout << "(BlockDiagLaplaceExplicitMatrix)OPERATOR: << print" << std::endl;
	output << A_minlinhost << std::endl;
}

/* MinLinHost: matrix-vector multiplication */
template<>
void BlockDiagLaplaceExplicitMatrix<MinlinHostVector>::matmult(MinlinHostVector &y, const MinlinHostVector &x) const { 
	if(DEBUG_MODE >= 100) std::cout << "(BlockDiagLaplaceExplicitMatrix)FUNCTION: matmult" << std::endl;

	y = A_minlinhost*x;	

}

#endif

/* -------------------------------- MINLIN DEVICE -------------------------*/

#ifdef USE_MINLIN

/* MinLinDevice: constructor from given right DeviceVector<double> */
template<>
BlockDiagLaplaceExplicitMatrix<MinlinDeviceVector>::BlockDiagLaplaceExplicitMatrix(const MinlinDeviceVector &x, int newK, double newalpha){
	/* init Petsc Vector */
	if(DEBUG_MODE >= 100) std::cout << "(BlockDiagLaplaceExplicitMatrix)CONSTRUCTOR: from MinLinDevice" << std::endl;

	int N = x.size();
	alpha = newalpha;
	K = newK; /* number of blocks */
	int T = N/(double)K; /* size of each block */

	/* get informations from given vector */
	MinlinDeviceMatrix A_new(N,N);

	int row,col; /* coordinates in block */
	int k, row_global, col_global; /* coordinates in global matrix */
	double new_value;

	for(k=0;k<K;k++){
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
				if(row == T-1){
					new_value = 1;
					if(col < row){
						new_value = -1;
					}
				}

				/* ordinary row */
				if(row > 0 && row < T-1){
					new_value = 2;
					if(col > row || col < row){
						new_value = -1;
					}
				}

				/* set value */
				row_global = k*T + row;
				col_global = k*T + col;
				if(row >= 0 && row <= T-1 && col >=0 && col <= T-1){ /* I am in my block */
					if(row_global >= 0 && row_global <= N-1 && col_global >=0 && col_global <= N-1){ /* I am in global matrix */
						A_new(row_global,col_global) = alpha*new_value;
					}
				}
				
			} /* endfor col */
		} /* endfor row */
	} /* endfor k */
			
	A_minlindevice = A_new;

}


/* MinLinDevice: destructor - destroy the matrix */
template<>
BlockDiagLaplaceExplicitMatrix<MinlinDeviceVector>::~BlockDiagLaplaceExplicitMatrix(){
	/* init Petsc Vector */
	if(DEBUG_MODE >= 100) std::cout << "(BlockDiagLaplaceExplicitMatrix)DESTRUCTOR" << std::endl;

	// TODO: how to destroy minlin matrix?
}

/* MinLinDevice: print matrix */
template<>
void BlockDiagLaplaceExplicitMatrix<MinlinDeviceVector>::print(std::ostream &output) const		
{
	if(DEBUG_MODE >= 100) std::cout << "(BlockDiagLaplaceExplicitMatrix)OPERATOR: << print" << std::endl;
	output << A_minlindevice << std::endl;
}

/* MinLinDevice: matrix-vector multiplication */
template<>
void BlockDiagLaplaceExplicitMatrix<MinlinDeviceVector>::matmult(MinlinDeviceVector &y, const MinlinDeviceVector &x) const { 
	if(DEBUG_MODE >= 100) std::cout << "(BlockDiagLaplaceExplicitMatrix)FUNCTION: matmult" << std::endl;

	y = A_minlindevice*x;	

}

#endif

} /* end of namespace */


#endif
