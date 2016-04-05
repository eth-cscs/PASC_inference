#ifndef BLOCKDIAGMATRIX_H
#define	BLOCKDIAGMATRIX_H

/* for debugging, if >= 100, then print info about ach called function */
extern int DEBUG_MODE;

#include <iostream>
#include <fstream>
#include "algebra.h" /* parent GeneralMatrix class */

#ifdef USE_PETSCVECTOR
	typedef petscvector::PetscVector PetscVector;
#endif

#ifdef USE_MINLIN
	typedef minlin::threx::HostVector<double> MinlinHostVector;
	typedef minlin::threx::DeviceVector<double> MinlinDeviceVector;
#endif


namespace pascinference {

template<class VectorBase, class MatrixBase>
class BlockDiagMatrix: public GeneralMatrix<VectorBase> {
	private:
		MatrixBase **blocks; /**< array of pointers to matrices */
		int nmb_blocks; /**< number of diagonal blocks */
		
	public:
		BlockDiagMatrix(int nmb_block, MatrixBase **new_blocks); /* constructor with number_of_blocks and blocks */
		~BlockDiagMatrix(); /* destructor - destroy inner matrix */

		void print(std::ostream &output) const; /* print matrix */
		std::string get_name() const;
		
		void matmult(VectorBase &y, const VectorBase &x) const; /* y = A*x */
};

template<class VectorBase, class MatrixBase>
std::string BlockDiagMatrix<VectorBase,MatrixBase>::get_name() const {
	return "BlockDiagMatrix";
}


template<class VectorBase, class MatrixBase>
BlockDiagMatrix<VectorBase,MatrixBase>::BlockDiagMatrix(int new_nmb_blocks, MatrixBase **new_blocks){
	if(DEBUG_MODE >= 100){
		coutMaster << "(BlockDiagMatrix)CONSTRUCTOR: from blocks" << std::endl;
	}

	this->blocks = new_blocks;
	this->nmb_blocks = new_nmb_blocks;
}


template<class VectorBase, class MatrixBase>
BlockDiagMatrix<VectorBase,MatrixBase>::~BlockDiagMatrix(){
	/* init Petsc Vector */
	if(DEBUG_MODE >= 100) coutMaster << "(BlockDiagMatrix)DESTRUCTOR" << std::endl;

	//TODO: what to do?
}

/* print matrix */
template<class VectorBase, class MatrixBase>
void BlockDiagMatrix<VectorBase,MatrixBase>::print(std::ostream &output) const		
{
	if(DEBUG_MODE >= 100) coutMaster << "(BlockDiagMatrix)OPERATOR: << print" << std::endl;

	int i;
	for(i=0;i<this->nmb_blocks;i++){ /* print each block */
		output << "block (" << i << ")" << std::endl;
		this->blocks[i]->print(output);
	}
}

/* Petsc: matrix-vector multiplication */
template<class VectorBase, class MatrixBase>
void BlockDiagMatrix<VectorBase,MatrixBase>::matmult(VectorBase &y, const VectorBase &x) const { 
	if(DEBUG_MODE >= 100) coutMaster << "(BlockDiagMatrix)FUNCTION: matmult" << std::endl;

	// TODO: maybe y is not initialized, who knows

	// TODO: write block matmult
}






} /* end of namespace */


#endif
