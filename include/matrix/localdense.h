#ifndef LOCALDENSEMATRIX_H
#define	LOCALDENSEMATRIX_H

/* for debugging, if >= 100, then print info about ach called function */
extern int DEBUG_MODE;

#include <iostream>
#include <fstream>
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

/* matlab matrix */ 
template<class VectorBase>
class LocalDenseMatrix: public GeneralMatrix<VectorBase> {
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
		
		template<class MatrixType>
		void write_aij(std::ifstream &myfile, MatrixType &matrix); /* write to given matrix using A(i,j) */
 		
 		/* reading from file */
 		int read_int_from_file(std::ifstream &myfile);
 		double read_double_from_file(std::ifstream &myfile);
 		int *read_int_array_from_file(std::ifstream &myfile, const int size);
 		double *read_double_array_from_file(std::ifstream &myfile, const int size);
 		int read_filesize(std::ifstream &myfile);
 		
 		int nmb_rows, nmb_cols;
 		
	public:
		LocalDenseMatrix(std::string filename); /* constructor from filename */
		LocalDenseMatrix(int nmb_rows, int nmb_cols); /* constructor with dimension */
		LocalDenseMatrix(double *values, int nmb_rows, int nmb_cols); /* constructor from array values */
		~LocalDenseMatrix(); /* destructor - destroy inner matrix */

		void print(ConsoleOutput &output) const; /* print matrix */
		void printcontent(ConsoleOutput &output) const;

		std::string get_name() const;
		
		void matmult(VectorBase &y, const VectorBase &x) const; /* y = A*x */

		void set_value(int row, int col, double value);
		void add_value(int row, int col, double value);

		void assemble();
};

template<class VectorBase>
std::string LocalDenseMatrix<VectorBase>::get_name() const {
	return "LocalDenseMatrix";
}

/* print matrix */
template<class VectorBase>
void LocalDenseMatrix<VectorBase>::print(ConsoleOutput &output) const		
{
	LOG_FUNC_BEGIN

	output << "LocalDense matrix (rows = " << this->nmb_rows << ", nmb_cols = "<< this->nmb_cols << ")";
	
//	TRY( MatView(A_petsc, PETSC_VIEWER_STDOUT_SELF) );

	LOG_FUNC_END
}

/* -------------------------------- PETSC VECTOR -------------------------*/

#ifdef USE_PETSCVECTOR

#include <petscmat.h>

/* Petsc: constructor from given filename */
template<>
LocalDenseMatrix<PetscVector>::LocalDenseMatrix(std::string filename){
	LOG_FUNC_BEGIN

	/* init Petsc Vector */
	if(DEBUG_MODE >= 100){
		coutMaster << "(LocalDenseMatrix)CONSTRUCTOR: from filename" << std::endl;
		coutMaster << " - read matrix in petsc format from: " << filename << std::endl;
	}

	/* prepare matrix */
	TRY( MatCreate(PETSC_COMM_SELF,&A_petsc) ); /* this is always local matrix! */
	TRY( MatSetType(A_petsc, MATSEQDENSE) );

	/* prepare viewer to load from file */
	PetscViewer mviewer;
	TRY( PetscViewerCreate(PETSC_COMM_SELF, &mviewer) );
	TRY( PetscViewerBinaryOpen(PETSC_COMM_WORLD , filename.c_str(), FILE_MODE_READ, &mviewer) );
	
	/* load matrix from viewer */
	TRY( MatLoad(A_petsc, mviewer) );

	/* destroy the viewer */
	TRY( PetscViewerDestroy(&mviewer) );
	
	/* assembly matrix */ //TODO: is it necessary?
	TRY( MatAssemblyBegin(A_petsc, MAT_FINAL_ASSEMBLY) );
	TRY( MatAssemblyEnd(A_petsc, MAT_FINAL_ASSEMBLY) );
	
	TRY( MatGetSize(A_petsc, &nmb_rows, &nmb_cols) );

	LOG_FUNC_END
}

/* Petsc: constructor from given array of values and size */
template<>
LocalDenseMatrix<PetscVector>::LocalDenseMatrix(double *values, int nmb_rows, int nmb_cols){
	LOG_FUNC_BEGIN

	/* prepare matrix, values in column major order! */
	TRY( MatCreateSeqDense(PETSC_COMM_SELF, nmb_rows, nmb_cols, values, &A_petsc) );

	/* assembly matrix */ //TODO: is it necessary?
	TRY( MatAssemblyBegin(A_petsc, MAT_FINAL_ASSEMBLY) );
	TRY( MatAssemblyEnd(A_petsc, MAT_FINAL_ASSEMBLY) );
	
	this->nmb_rows = nmb_rows;
	this->nmb_cols = nmb_cols;

	LOG_FUNC_END
}

/* Petsc: constructor from given size */
template<>
LocalDenseMatrix<PetscVector>::LocalDenseMatrix(int nmb_rows, int nmb_cols){
	LOG_FUNC_BEGIN

	/* prepare matrix, values in column major order! */
	TRY( MatCreateSeqDense(PETSC_COMM_SELF, nmb_rows, nmb_cols, NULL, &A_petsc) );

	/* assembly matrix */ //TODO: is it necessary?
	TRY( MatAssemblyBegin(A_petsc, MAT_FINAL_ASSEMBLY) );
	TRY( MatAssemblyEnd(A_petsc, MAT_FINAL_ASSEMBLY) );

	this->nmb_rows = nmb_rows;
	this->nmb_cols = nmb_cols;
	
	LOG_FUNC_END
}

/* Petsc: destructor - destroy the matrix */
template<>
LocalDenseMatrix<PetscVector>::~LocalDenseMatrix(){
	LOG_FUNC_BEGIN

	if(petscvector::PETSC_INITIALIZED){ /* maybe Petsc was already finalized and there is nothing to destroy */
		TRY( MatDestroy(&A_petsc) );
	}

	LOG_FUNC_END
}

/* Petsc: matrix-vector multiplication */
template<>
void LocalDenseMatrix<PetscVector>::matmult(PetscVector &y, const PetscVector &x) const { 
	LOG_FUNC_BEGIN

	// TODO: maybe y is not initialized, who knows
	
	TRY( MatMult(A_petsc, x.get_vector(), y.get_vector()) );

	LOG_FUNC_END
}

/* Petsc: set value of matrix */
template<>
void LocalDenseMatrix<PetscVector>::set_value(int row, int col, double value) { 
	if(DEBUG_MODE >= 100) coutMaster << "(LocalDenseMatrix)FUNCTION: set value" << std::endl;

	TRY( MatSetValue(A_petsc, row, col, value, INSERT_VALUES) );

}

/* Petsc: set value of matrix */
template<>
void LocalDenseMatrix<PetscVector>::add_value(int row, int col, double value) { 
	if(DEBUG_MODE >= 100) coutMaster << "(LocalDenseMatrix)FUNCTION: set value" << std::endl;

	TRY( MatSetValue(A_petsc, row, col, value, ADD_VALUES) );

}

template<>
void LocalDenseMatrix<PetscVector>::assemble() { 
	/* assembly matrix */ 
	TRY( MatAssemblyBegin(A_petsc, MAT_FINAL_ASSEMBLY) );
	TRY( MatAssemblyEnd(A_petsc, MAT_FINAL_ASSEMBLY) );
}

template<>
void LocalDenseMatrix<PetscVector>::printcontent(ConsoleOutput &output) const		
{
	LOG_FUNC_BEGIN

	const double *values_row;
	double value;

	std::ostringstream temp;

	int i,j;
	for(i=0;i<nmb_rows;i++){
		/* get one row */
		TRY( MatGetRow(A_petsc,i,NULL,NULL,&values_row) );

		for(j=0;j<nmb_cols;j++){
			temp << values_row[j];
			output << std::setw(8) << temp.str();
			temp.str("");
			if(j+1 < nmb_cols){
				output << ",";
			}
		}
		output << std::endl;

		TRY( MatRestoreRow(A_petsc,i,NULL,NULL,&values_row) );

	}

	LOG_FUNC_END
}

#endif

/* -------------------------------- MINLIN HOST -------------------------*/

#ifdef USE_MINLIN

template<>
LocalDenseMatrix<MinlinHostVector>::LocalDenseMatrix(std::string filename){
	LOG_FUNC_BEGIN

	if(DEBUG_MODE >= 100){
		coutMaster << "(LocalDenseMatrix)CONSTRUCTOR: from filename" << std::endl;
		coutMaster << " - read matrix in petsc format from: " << filename << std::endl;
	}

	/* open file */
	std::ifstream myfile(filename.c_str(), std::ios::in | std::ios::binary);

	/* load header from file */
	int mat_file_classid = read_int_from_file(myfile);
	if( mat_file_classid != 1211216){
		// TODO: give error, this is not PetscBin file
	}

	int nmb_of_rows = read_int_from_file(myfile);
	int nmb_of_cols = read_int_from_file(myfile);

	this->nmb_rows = nmb_of_rows;
	this->nmb_cols = nmb_of_cols;

	/* prepare matrix */
	MinlinHostMatrix A_new(nmb_of_rows,nmb_of_cols);

	/* set initial value */
	A_new(all) = 0.0;

	/* reset reader to the begining of file */
	myfile.clear();

	/* write as Aij */
	write_aij(myfile, A_new);
	
	/* close file */
    myfile.close();
	
	/* set new matrix */	
	A_minlinhost = A_new;

	LOG_FUNC_END
}

/* Minlin: constructor from given array of values and size */
template<>
LocalDenseMatrix<MinlinHostVector>::LocalDenseMatrix(double *values, int nmb_rows, int nmb_cols){
	LOG_FUNC_BEGIN

	MinlinHostMatrix A_new(nmb_rows,nmb_cols);

	int row,col;
	for(col=0;col<nmb_cols;col++){
		for(row=0;row<nmb_rows;row++){
			A_new(row,col) = values[col*nmb_rows+row];
		}
	}

	/* set new matrix */	
	A_minlinhost = A_new;
	
	this->nmb_rows = nmb_rows;
	this->nmb_cols = nmb_cols;

	LOG_FUNC_END
}

/* MinLin: constructor from given size */
template<>
LocalDenseMatrix<MinlinHostVector>::LocalDenseMatrix(int nmb_rows, int nmb_cols){
	LOG_FUNC_BEGIN

	MinlinHostMatrix A_new(nmb_rows,nmb_cols);
	A_minlinhost = A_new;
	
	LOG_FUNC_END
}

/* MinLin: destructor - destroy the matrix */
template<>
LocalDenseMatrix<MinlinHostVector>::~LocalDenseMatrix(){
	LOG_FUNC_BEGIN

	// TODO: destroy minlin matrix

	LOG_FUNC_END
}

/* Minlin: matrix-vector multiplication */
template<>
void LocalDenseMatrix<MinlinHostVector>::matmult(MinlinHostVector &y, const MinlinHostVector &x) const { 
	LOG_FUNC_BEGIN

	y = A_minlinhost*x;	
		
	LOG_FUNC_END
}

/* MinLin: set value of matrix */
template<>
void LocalDenseMatrix<MinlinHostVector>::set_value(int row, int col, double value) { 

	A_minlinhost(row,col) = value;	

}

#endif

/* -------------------------------- MINLIN DEVICE -------------------------*/

#ifdef USE_MINLIN

template<>
LocalDenseMatrix<MinlinDeviceVector>::LocalDenseMatrix(std::string filename){
	LOG_FUNC_BEGIN

	if(DEBUG_MODE >= 100){
		coutMaster << "(LocalDenseMatrix)CONSTRUCTOR: from filename" << std::endl;
		coutMaster << " - read matrix in petsc format from: " << filename << std::endl;
	}


	/* open file */
	std::ifstream myfile(filename.c_str(), std::ios::in | std::ios::binary);

	/* load header from file */
	int mat_file_classid = read_int_from_file(myfile);
	if( mat_file_classid != 1211216){
		// TODO: give error, this is not PetscBin file
	}

	int nmb_of_rows = read_int_from_file(myfile);
	int nmb_of_cols = read_int_from_file(myfile);

	this->nmb_rows = nmb_of_rows;
	this->nmb_cols = nmb_of_cols;
	
	/* prepare matrix */
	MinlinDeviceMatrix A_new(nmb_of_rows,nmb_of_cols);

	/* reset reader to the begining of file */
	myfile.clear();

	/* write as Aij */
	write_aij(myfile, A_new);
	
	/* close file */
    myfile.close();
	
	/* set new matrix */	
	A_minlindevice = A_new;

	LOG_FUNC_END
}

/* MinLin: constructor from given size */
template<>
LocalDenseMatrix<MinlinDeviceVector>::LocalDenseMatrix(int nmb_rows, int nmb_cols){
	LOG_FUNC_BEGIN

	/* init Petsc Vector */
	if(DEBUG_MODE >= 100){
		coutMaster << "(LocalDenseMatrix)CONSTRUCTOR: from size" << std::endl;
	}

	MinlinDeviceMatrix A_new(nmb_rows,nmb_cols);
	A_minlindevice = A_new;
	
	this->nmb_rows = nmb_rows;
	this->nmb_cols = nmb_cols;	
		
	LOG_FUNC_END
}

/* Minlin: destructor - destroy the matrix */
template<>
LocalDenseMatrix<MinlinDeviceVector>::~LocalDenseMatrix(){
	LOG_FUNC_BEGIN

	// TODO: destroy minlin matrix

	LOG_FUNC_END
}

/* MinLin: matrix-vector multiplication */
template<>
void LocalDenseMatrix<MinlinDeviceVector>::matmult(MinlinDeviceVector &y, const MinlinDeviceVector &x) const { 
	LOG_FUNC_BEGIN

	y = A_minlindevice*x;	

	LOG_FUNC_END
}

/* MinLin: set value of matrix */
template<>
void LocalDenseMatrix<MinlinDeviceVector>::set_value(int row, int col, double value) { 
	if(DEBUG_MODE >= 100) coutMaster << "(LocalDenseMatrix)FUNCTION: set value" << std::endl;

	A_minlindevice(row,col) = value;	
		
}


#endif


/* -------------------------------- GENERAL FUNCTIONS ---------------------*/
template<class VectorBase>
template<class MatrixType>
void LocalDenseMatrix<VectorBase>::write_aij(std::ifstream &myfile, MatrixType &A){
	LOG_FUNC_BEGIN

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

	/* print info */
	if(DEBUG_MODE >= 100){
		/* get the size of the file */
		int file_size = read_filesize(myfile);
		
		coutMaster << "size of file:     " << file_size << std::endl;
		coutMaster << "mat_file_classid: " << mat_file_classid << std::endl;
		coutMaster << "nmb_of_rows:      " << nmb_of_rows << std::endl;
		coutMaster << "nmb_of_cols:      " << nmb_of_cols << std::endl;
		coutMaster << "nmb_of_nz:        " << nmb_of_nz << std::endl;

		coutMaster << "nmb_of_row_nz:    ";
		for(i=0;i<nmb_of_rows;i++) coutMaster << nmb_of_row_nz[i] << ", ";
		coutMaster << std::endl;

		coutMaster << "col_idx_nz:       ";
		for(i=0;i<nmb_of_nz;i++) coutMaster << col_idx_nz[i] << ", ";
		coutMaster << std::endl;

		coutMaster << "values:           ";
		for(i=0;i<nmb_of_nz;i++) coutMaster << values[i] << ", ";
		coutMaster << std::endl;
	}

	int id_row, id_in_row;
	int col_idid=0;
	for(id_row=0; id_row<nmb_of_rows; id_row++){
		for(id_in_row=0; id_in_row<nmb_of_row_nz[id_row]; id_in_row++){
			A(id_row,col_idx_nz[col_idid]) = values[col_idid];
			col_idid++;
		}
	}

	LOG_FUNC_END
}

template<class VectorBase>
int LocalDenseMatrix<VectorBase>::read_int_from_file(std::ifstream &myfile) {
	int value;
	myfile.read((char *)&value, sizeof(int)); /* read block of memory */
	value = __builtin_bswap32(value);
	return value;
}

template<class VectorBase>
double LocalDenseMatrix<VectorBase>::read_double_from_file(std::ifstream &myfile) {
	long int int_value;
	myfile.read((char *)&int_value, sizeof(long int)); /* read block of memory */

	int_value = __builtin_bswap64(int_value);

	return *((double*)&int_value); 
}

template<class VectorBase>
int *LocalDenseMatrix<VectorBase>::read_int_array_from_file(std::ifstream &myfile, const int size) {
	int *return_array = new int[size];
	int i;
	for(i = 0; i < size; i++){
		return_array[i] = read_int_from_file(myfile);
	}
	
	return return_array; /* deal with MATLAB 'ieee-be' */
}

template<class VectorBase>
double *LocalDenseMatrix<VectorBase>::read_double_array_from_file(std::ifstream &myfile, const int size) {
	double *return_array = new double[size];
	int i;
	for(i = 0; i < size; i++){
		return_array[i] = read_double_from_file(myfile);
	}
	
	return return_array; /* deal with MATLAB 'ieee-be' */
}

template<class VectorBase>
int LocalDenseMatrix<VectorBase>::read_filesize(std::ifstream &myfile) {
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
