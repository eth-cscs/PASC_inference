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
		int blocksize; /**< constant size of all blocks */
		
	public:
		BlockDiagMatrix(int nmb_block, MatrixBase **new_blocks, int blocksize); /* constructor with number_of_blocks and blocks */
		~BlockDiagMatrix(); /* destructor - destroy inner matrix */

		void print(std::ostream &output) const; /* print matrix */
		std::string get_name() const;
		
		void matmult(VectorBase &y, const VectorBase &x) const; /* y = A*x */

		int get_nmb_blocks() const;
		MatrixBase **get_blocks() const;
		int get_blocksize() const;

};

template<class VectorBase, class MatrixBase>
std::string BlockDiagMatrix<VectorBase,MatrixBase>::get_name() const {
	return "BlockDiagMatrix";
}


template<class VectorBase, class MatrixBase>
BlockDiagMatrix<VectorBase,MatrixBase>::BlockDiagMatrix(int new_nmb_blocks, MatrixBase **new_blocks, int new_blocksize){
	if(DEBUG_MODE >= 100){
		coutMaster << "(BlockDiagMatrix)CONSTRUCTOR: from blocks" << std::endl;
	}

	this->blocks = new_blocks;
	this->nmb_blocks = new_nmb_blocks;
	this->blocksize = new_blocksize;
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

	output << "blocksize:     " << this->blocksize << std::endl;
	output << "nmb of blocks: " << this->nmb_blocks << std::endl;

// TODO: following does not work if output is coutAll because of barrier
//	int i;
//	for(i=0;i<this->nmb_blocks;i++){ /* print each block */
//		output << "block (" << i << ")" << std::endl;
//		this->blocks[i]->print(output);
//	}
}

/* matrix-vector multiplication */
template<class VectorBase, class MatrixBase>
void BlockDiagMatrix<VectorBase,MatrixBase>::matmult(VectorBase &y, const VectorBase &x) const { 
	if(DEBUG_MODE >= 100) coutMaster << "(BlockDiagMatrix)FUNCTION: matmult" << std::endl;

	// TODO: maybe y is not initialized, who knows

	// TODO: write block matmult
}

/* get number of blocks */
template<class VectorBase, class MatrixBase>
int BlockDiagMatrix<VectorBase,MatrixBase>::get_nmb_blocks() const { 
	if(DEBUG_MODE >= 100) coutMaster << "(BlockDiagMatrix)FUNCTION: get_nmb_blocks" << std::endl;

	return this->nmb_blocks;
}


/* get number of blocks */
template<class VectorBase, class MatrixBase>
MatrixBase **BlockDiagMatrix<VectorBase,MatrixBase>::get_blocks() const { 
	if(DEBUG_MODE >= 100) coutMaster << "(BlockDiagMatrix)FUNCTION: get_blocks" << std::endl;

	return this->blocks;
}

/* get blocksize */
template<class VectorBase, class MatrixBase>
int BlockDiagMatrix<VectorBase,MatrixBase>::get_blocksize() const { 
	if(DEBUG_MODE >= 100) coutMaster << "(BlockDiagMatrix)FUNCTION: get_blocksize" << std::endl;

	return this->blocksize;
}



} /* end of namespace */


#endif
