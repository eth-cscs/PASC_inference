#ifndef MATLABMATRIX_H
#define	MATLABMATRIX_H

/* for debugging, if >= 100, then print info about ach called function */
extern int DEBUG_MODE;

#include <iostream>
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
	std::ifstream myfile(filename.c_str(), std::ios::in | std::ios::binary);

	/* petsc binary file content: (http://www.mcs.anl.gov/petsc/petsc-dev/docs/manualpages/Mat/MatLoad.html)
	 * int    MAT_FILE_CLASSID
	 * int    number of rows
	 * int    number of columns
	 * int    total number of nonzeros
     * int    *number nonzeros in each row
     * int    *column indices of all nonzeros (starting index is zero)
     * PetscScalar *values of all nonzeros
	 */ 

	/* variables for reading the file */
//	char buffer[sizeof(int32_t)];
	
	int32_t mat_file_classid[4];
//	int nmb_of_rows;
/*	int nmb_of_cols;
	int nmb_of_nz;
*/

	/* set reader to the begining of the file */
//	myfile.seekg(0, std::ios::beg);
	
	myfile.clear();
	
	/* read MAT_FILE_CLASSID */
	myfile.read((char*)&mat_file_classid, sizeof(mat_file_classid)); /* read block of memory */ 
	
	std::cout<< "sizeof(mat_file_class_id): " << sizeof(mat_file_classid) << std::endl;
	std::cout<< "content0: " << (int)mat_file_classid[0] << std::endl;
	std::cout<< "content1: " << (int)mat_file_classid[1] << std::endl;
	std::cout<< "content2: " << (int)mat_file_classid[2] << std::endl;
	std::cout<< "content3: " << (int)mat_file_classid[3] << std::endl;
	
//	myfile.seekg(sizeof(mat_file_classid),std::ios::cur); /* change location of reader */

	/* read nmb_of_rows */
//	myfile.read(memblock_int, sizeof(int)); /* read block of memory */
//	nmb_of_rows = atoi(memblock_int); /* convert to int */
//	myfile.seekg(sizeof(int),std::ios::cur); /* change location of reader */


//	std::cout << "MAT_FILE_CLASSID: " << mat_file_classid << std::endl;
//	std::cout << "nmb of rows:      " << nmb_of_rows << std::endl;

	/* close file */
    myfile.close();
	
//	std::cout<< "kurveliak: " << sizeof(int) << ", picus: " << memblock_int << ", kokotar: " << atoi(memblock_int) << std::endl;
	
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
