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
		MatlabMatrix(const VectorType &x); /* constructor */
		~MatlabMatrix(); /* destructor - destroy inner matrix */

		void print(std::ostream &output) const; /* print matrix */
		void matmult(VectorType &y, const VectorType &x) const; /* y = A*x */

};



/* -------------------------------- PETSC VECTOR -------------------------*/

/* Petsc: constructor from given right PetscVector */
template<>
MatlabMatrix<PetscVector>::MatlabMatrix(const PetscVector &x){
	/* init Petsc Vector */
	if(DEBUG_MODE >= 100) std::cout << "(MatlabMatrix)CONSTRUCTOR: from filename" << std::endl;

	int N, n;

	/* get informations from given vector */
	N = x.size();
	n = x.local_size();

	if(DEBUG_MODE >= 100) std::cout << " - read matrix in petsc format" << std::endl;

	/* prepare viewer to load from file */
	PetscViewer mviewer;
	TRY( PetscViewerCreate(PETSC_COMM_WORLD, &mviewer) );
	TRY( PetscViewerSetType(mviewer,PETSCVIEWERMATLAB) );
	TRY( PetscViewerBinaryOpen(PETSC_COMM_WORLD , "mymatrix.mat", FILE_MODE_READ, &mviewer) );
	
//	TRY( PetscViewerMatlabOpen(PETSC_COMM_WORLD, "mymatrix.mat", FILE_MODE_READ, &mviewer) );
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
// TODO: find a way how to load matlab matrix without Petsc

/* -------------------------------- MINLIN DEVICE -------------------------*/
// TODO: find a way how to load matlab matrix without Petsc




} /* end of namespace */


#endif
