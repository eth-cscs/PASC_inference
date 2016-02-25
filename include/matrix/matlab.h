#ifndef MATLABMATRIX_H
#define	MATLABMATRIX_H

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

/* matlab matrix */ 
template<class VectorType>
class MatlabMatrix: public GeneralMatrix<VectorType> {
	private:
		/* Petsc stuff */ // TODO: if USE_PETSC
		PetscMatrix A_petsc;

		/* MINLIN stuff */ // TODO: if USE_MINLIN
		MinlinHostMatrix A_minlinhost;
		MinlinDeviceMatrix A_minlindevice;
		
	
	public:
		MatlabMatrix(const VectorType &x, std::string filename); /* constructor */
		~MatlabMatrix(); /* destructor - destroy inner matrix */

		void print(std::ostream &output) const; /* print matrix */
		void matmult(VectorType &y, const VectorType &x) const; /* y = A*x */

};



/* -------------------------------- PETSC VECTOR -------------------------*/

/* Petsc: constructor from given right PetscVector */
template<>
MatlabMatrix<PetscVector>::MatlabMatrix(const PetscVector &x, std::string filename){
	/* init Petsc Vector */
	if(DEBUG_MODE >= 100){
		std::cout << "(MatlabMatrix)CONSTRUCTOR: from filename" << std::endl;
		std::cout << " - read matrix in petsc format from: " << filename << std::endl;
	}

	/* get informations from given vector */
	int N, n;

	N = x.size();
	n = x.local_size();

	/* prepare matrix */
	TRY( MatCreate(PETSC_COMM_WORLD,&A_petsc) );
	TRY( MatSetSizes(A_petsc,n,n,N,N) );
	TRY( MatSetFromOptions(A_petsc) ); 

	/* prepare viewer to load from file */
	PetscViewer mviewer;
	TRY( PetscViewerCreate(PETSC_COMM_SELF, &mviewer) );
	TRY( PetscViewerBinaryOpen(PETSC_COMM_WORLD , filename.c_str(), FILE_MODE_READ, &mviewer) );
	
	/* load matrix from viewer */
	TRY( MatLoad(A_petsc, mviewer) );

	/* destroy the viewer */
	TRY( PetscViewerDestroy(&mviewer) );
	
}

/* Petsc: destructor - destroy the matrix */
template<>
MatlabMatrix<PetscVector>::~MatlabMatrix(){
	/* init Petsc Vector */
	if(DEBUG_MODE >= 100) std::cout << "(MatlabMatrix)DESTRUCTOR" << std::endl;

	if(petscvector::PETSC_INITIALIZED){ /* maybe Petsc was already finalized and there is nothing to destroy */
		TRY( MatDestroy(&A_petsc) );
	}
}

/* print matrix */
template<>
void MatlabMatrix<PetscVector>::print(std::ostream &output) const		
{
	if(DEBUG_MODE >= 100) std::cout << "(MatlabMatrix)OPERATOR: << print" << std::endl;

	output << "Matlab matrix (sorry, 'only' MatView from Petsc follows):" << std::endl;
	output << "----------------------------------------------------------" << std::endl;
	
	TRY( MatView(A_petsc, PETSC_VIEWER_STDOUT_WORLD) );

	output << "----------------------------------------------------------" << std::endl;
}

/* Petsc: matrix-vector multiplication */
template<>
void MatlabMatrix<PetscVector>::matmult(PetscVector &y, const PetscVector &x) const { 
	if(DEBUG_MODE >= 100) std::cout << "(MatlabMatrix)FUNCTION: matmult" << std::endl;

	// TODO: maybe y is not initialized, who knows
	
	TRY( MatMult(A_petsc, x.get_vector(), y.get_vector()) ); // TODO: I dont want to use get_vector :( friend in PetscVector? and in MinLin?
}



/* -------------------------------- MINLIN HOST -------------------------*/
template<>
MatlabMatrix<MinlinHostVector>::MatlabMatrix(const MinlinHostVector &x, std::string filename){
	if(DEBUG_MODE >= 100){
		std::cout << "(MatlabMatrix)CONSTRUCTOR: from filename" << std::endl;
		std::cout << " - read matrix in petsc format from: " << filename << std::endl;
	}

	/* get informations from given vector */
	int N = x.size();

	/* prepare matrix */
	MinlinHostMatrix A_new(N,N);

	/* set initial value */
	A_new(all) = 0.0;

	/* open file */

	
	/* set new matrix */	
	A_minlinhost = A_new;
}

/* Petsc: destructor - destroy the matrix */
template<>
MatlabMatrix<MinlinHostVector>::~MatlabMatrix(){
	/* init Petsc Vector */
	if(DEBUG_MODE >= 100) std::cout << "(MatlabMatrix)DESTRUCTOR" << std::endl;

	// TODO: destroy minlin matrix
}

/* print matrix */
template<>
void MatlabMatrix<MinlinHostVector>::print(std::ostream &output) const {
	if(DEBUG_MODE >= 100) std::cout << "(MatlabMatrix)OPERATOR: << print" << std::endl;

	output << A_minlinhost << std::endl;
}

/* Petsc: matrix-vector multiplication */
template<>
void MatlabMatrix<MinlinHostVector>::matmult(MinlinHostVector &y, const MinlinHostVector &x) const { 
	if(DEBUG_MODE >= 100) std::cout << "(MatlabMatrix)FUNCTION: matmult" << std::endl;

	y = A_minlinhost*x;	
		
}



/* -------------------------------- MINLIN DEVICE -------------------------*/
// TODO: same as for host




} /* end of namespace */


#endif
