#ifndef FILECRSMATRIX_H
#define	FILECRSMATRIX_H

/* for debugging, if >= 100, then print info about ach called function */
extern int DEBUG_MODE;

#include <iostream>
#include <fstream>
#include "algebra.h" /* parent GeneralMatrix class */

typedef petscvector::PetscVector PetscVector;
typedef Mat PetscMatrix;

typedef minlin::threx::HostMatrix<double> MinlinHostMatrix;
typedef minlin::threx::HostVector<double> MinlinHostVector;

typedef minlin::threx::DeviceMatrix<double> MinlinDeviceMatrix;
typedef minlin::threx::DeviceVector<double> MinlinDeviceVector;



namespace pascinference {

/* matlab matrix */ 
template<class VectorBase>
class FileCRSMatrix: public GeneralMatrix<VectorBase> {
	private:
		/* Petsc stuff */ // TODO: if USE_PETSC
		PetscMatrix A_petsc;

		/* MINLIN stuff */ // TODO: if USE_MINLIN
		MinlinHostMatrix A_minlinhost;
		MinlinDeviceMatrix A_minlindevice;
		
		template<class MatrixType>
		void write_aij(std::ifstream &myfile, MatrixType &matrix); /* write to given matrix using A(i,j) */
 		
 		/* reading from file */
 		int read_int_from_file(std::ifstream &myfile);
 		double read_double_from_file(std::ifstream &myfile);
 		int *read_int_array_from_file(std::ifstream &myfile, const int size);
 		double *read_double_array_from_file(std::ifstream &myfile, const int size);
 		int read_filesize(std::ifstream &myfile);
 		
	public:
		FileCRSMatrix(const VectorBase &x, std::string filename); /* constructor */
		~FileCRSMatrix(); /* destructor - destroy inner matrix */

		void print(std::ostream &output) const; /* print matrix */
		void matmult(VectorBase &y, const VectorBase &x) const; /* y = A*x */

};



/* -------------------------------- PETSC VECTOR -------------------------*/

/* Petsc: constructor from given right PetscVector */
template<>
FileCRSMatrix<PetscVector>::FileCRSMatrix(const PetscVector &x, std::string filename){
	/* init Petsc Vector */
	if(DEBUG_MODE >= 100){
		std::cout << "(FileCRSMatrix)CONSTRUCTOR: from filename" << std::endl;
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
FileCRSMatrix<PetscVector>::~FileCRSMatrix(){
	/* init Petsc Vector */
	if(DEBUG_MODE >= 100) std::cout << "(FileCRSMatrix)DESTRUCTOR" << std::endl;

	if(petscvector::PETSC_INITIALIZED){ /* maybe Petsc was already finalized and there is nothing to destroy */
		TRY( MatDestroy(&A_petsc) );
	}
}

/* print matrix */
template<>
void FileCRSMatrix<PetscVector>::print(std::ostream &output) const		
{
	if(DEBUG_MODE >= 100) std::cout << "(FileCRSMatrix)OPERATOR: << print" << std::endl;

	output << "FileCRS matrix (sorry, 'only' MatView from Petsc follows):" << std::endl;
	output << "----------------------------------------------------------" << std::endl;
	
	TRY( MatView(A_petsc, PETSC_VIEWER_STDOUT_WORLD) );

	output << "----------------------------------------------------------" << std::endl;
}

/* Petsc: matrix-vector multiplication */
template<>
void FileCRSMatrix<PetscVector>::matmult(PetscVector &y, const PetscVector &x) const { 
	if(DEBUG_MODE >= 100) std::cout << "(FileCRSMatrix)FUNCTION: matmult" << std::endl;

	// TODO: maybe y is not initialized, who knows
	
	TRY( MatMult(A_petsc, x.get_vector(), y.get_vector()) ); // TODO: I dont want to use get_vector :( friend in PetscVector? and in MinLin?
}



/* -------------------------------- MINLIN HOST -------------------------*/
template<>
FileCRSMatrix<MinlinHostVector>::FileCRSMatrix(const MinlinHostVector &x, std::string filename){
	if(DEBUG_MODE >= 100){
		std::cout << "(FileCRSMatrix)CONSTRUCTOR: from filename" << std::endl;
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

	/* write as Aij */
	write_aij(myfile, A_new);
	
	/* close file */
    myfile.close();
	
	/* set new matrix */	
	A_minlinhost = A_new;
}

/* Petsc: destructor - destroy the matrix */
template<>
FileCRSMatrix<MinlinHostVector>::~FileCRSMatrix(){
	/* init Petsc Vector */
	if(DEBUG_MODE >= 100) std::cout << "(FileCRSMatrix)DESTRUCTOR" << std::endl;

	// TODO: destroy minlin matrix
}

/* print matrix */
template<>
void FileCRSMatrix<MinlinHostVector>::print(std::ostream &output) const {
	if(DEBUG_MODE >= 100) std::cout << "(FileCRSMatrix)OPERATOR: << print" << std::endl;

	output << A_minlinhost << std::endl;
}

/* Petsc: matrix-vector multiplication */
template<>
void FileCRSMatrix<MinlinHostVector>::matmult(MinlinHostVector &y, const MinlinHostVector &x) const { 
	if(DEBUG_MODE >= 100) std::cout << "(FileCRSMatrix)FUNCTION: matmult" << std::endl;

	y = A_minlinhost*x;	
		
}



/* -------------------------------- MINLIN DEVICE -------------------------*/
template<>
FileCRSMatrix<MinlinDeviceVector>::FileCRSMatrix(const MinlinDeviceVector &x, std::string filename){
	if(DEBUG_MODE >= 100){
		std::cout << "(FileCRSMatrix)CONSTRUCTOR: from filename" << std::endl;
		std::cout << " - read matrix in petsc format from: " << filename << std::endl;
	}

	/* get informations from given vector */
	int N = x.size();

	/* prepare matrix */
	MinlinDeviceMatrix A_new(N,N);

	/* set initial value */
	A_new(all) = 0.0;

	/* open file */
	std::ifstream myfile(filename.c_str(), std::ios::in | std::ios::binary);

	/* write as Aij */
	write_aij(myfile, A_new);
	
	/* close file */
    myfile.close();
	
	/* set new matrix */	
	A_minlinhost = A_new;
}

/* Petsc: destructor - destroy the matrix */
template<>
FileCRSMatrix<MinlinDeviceVector>::~FileCRSMatrix(){
	/* init Petsc Vector */
	if(DEBUG_MODE >= 100) std::cout << "(FileCRSMatrix)DESTRUCTOR" << std::endl;

	// TODO: destroy minlin matrix
}

/* print matrix */
template<>
void FileCRSMatrix<MinlinDeviceVector>::print(std::ostream &output) const {
	if(DEBUG_MODE >= 100) std::cout << "(FileCRSMatrix)OPERATOR: << print" << std::endl;

	output << A_minlinhost << std::endl;
}

/* Petsc: matrix-vector multiplication */
template<>
void FileCRSMatrix<MinlinDeviceVector>::matmult(MinlinDeviceVector &y, const MinlinDeviceVector &x) const { 
	if(DEBUG_MODE >= 100) std::cout << "(FileCRSMatrix)FUNCTION: matmult" << std::endl;

	y = A_minlinhost*x;	
		
}




/* -------------------------------- GENERAL FUNCTIONS ---------------------*/
template<class VectorBase>
template<class MatrixType>
void FileCRSMatrix<VectorBase>::write_aij(std::ifstream &myfile, MatrixType &A){

	/* petsc binary file content: (http://www.mcs.anl.gov/petsc/petsc-dev/docs/manualpages/Mat/MatLoad.html)
	 * int    MAT_FILE_CLASSID
	 * int    number of rows
	 * int    number of columns
	 * int    total number of nonzeros
     * int    *number nonzeros in each row
     * int    *column indices of all nonzeros (starting index is zero)
     * PetscScalar *values of all nonzeros
	 */ 

	int i;

	/* read values from file */
	int mat_file_classid = read_int_from_file(myfile);
	int nmb_of_rows = read_int_from_file(myfile);
	int nmb_of_cols = read_int_from_file(myfile);
	int nmb_of_nz = read_int_from_file(myfile);

	int *nmb_of_row_nz = read_int_array_from_file(myfile,nmb_of_rows);
	int *col_idx_nz = read_int_array_from_file(myfile,nmb_of_nz);
	double *values = read_double_array_from_file(myfile,nmb_of_nz);

	if( mat_file_classid != 1211216){
		// TODO: give error, this is not PetscBin file
	}

	/* print info */
	if(DEBUG_MODE >= 100){
		/* get the size of the file */
		int file_size = read_filesize(myfile);
		
		std::cout << "size of file:     " << file_size << std::endl;
		std::cout << "mat_file_classid: " << mat_file_classid << std::endl;
		std::cout << "nmb_of_rows:      " << nmb_of_rows << std::endl;
		std::cout << "nmb_of_cols:      " << nmb_of_cols << std::endl;
		std::cout << "nmb_of_nz:        " << nmb_of_nz << std::endl;

		std::cout << "nmb_of_row_nz:    ";
		for(i=0;i<nmb_of_rows;i++) std::cout << nmb_of_row_nz[i] << ", ";
		std::cout << std::endl;

		std::cout << "col_idx_nz:       ";
		for(i=0;i<nmb_of_nz;i++) std::cout << col_idx_nz[i] << ", ";
		std::cout << std::endl;

		std::cout << "values:           ";
		for(i=0;i<nmb_of_nz;i++) std::cout << values[i] << ", ";
		std::cout << std::endl;
	}

	int id_row, id_in_row;
	int col_idid=0;
	for(id_row=0; id_row<nmb_of_rows; id_row++){
		for(id_in_row=0; id_in_row<nmb_of_row_nz[id_row]; id_in_row++){
			A(id_row,col_idx_nz[col_idid]) = values[col_idid];
			col_idid++;
		}
	}


}

template<class VectorBase>
int FileCRSMatrix<VectorBase>::read_int_from_file(std::ifstream &myfile) {
	int value;
	myfile.read((char *)&value, sizeof(int)); /* read block of memory */
	value = __builtin_bswap32(value);
	return value;
}

template<class VectorBase>
double FileCRSMatrix<VectorBase>::read_double_from_file(std::ifstream &myfile) {
	long int int_value;
	myfile.read((char *)&int_value, sizeof(long int)); /* read block of memory */

	int_value = __builtin_bswap64(int_value);

	return *((double*)&int_value); 
}

template<class VectorBase>
int *FileCRSMatrix<VectorBase>::read_int_array_from_file(std::ifstream &myfile, const int size) {
	int *return_array = new int[size];
	int i;
	for(i = 0; i < size; i++){
		return_array[i] = read_int_from_file(myfile);
	}
	
	return return_array; /* deal with MATLAB 'ieee-be' */
}

template<class VectorBase>
double *FileCRSMatrix<VectorBase>::read_double_array_from_file(std::ifstream &myfile, const int size) {
	double *return_array = new double[size];
	int i;
	for(i = 0; i < size; i++){
		return_array[i] = read_double_from_file(myfile);
	}
	
	return return_array; /* deal with MATLAB 'ieee-be' */
}

template<class VectorBase>
int FileCRSMatrix<VectorBase>::read_filesize(std::ifstream &myfile) {
	std::streampos begin,end;

	/* get begin */
	myfile.seekg(0, std::ios::beg);
	begin = myfile.tellg();

	/* get end */
	myfile.seekg(0, std::ios::end);
	end = myfile.tellg();

	/* restart the file seeker */
	myfile.seekg(0, std::ios::beg);

	/* return difference, retype to int */
	return (int)(end-begin);
}




} /* end of namespace */


#endif
