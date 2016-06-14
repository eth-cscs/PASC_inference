#ifndef BLOCKGRAPHMATRIX_H
#define	BLOCKGRAPHMATRIX_H

#include "pascinference.h"

#ifdef USE_PETSCVECTOR
	typedef petscvector::PetscVector PetscVector;
#endif

#ifdef USE_MINLIN
	typedef minlin::threx::HostVector<double> MinlinHostVector;
	typedef minlin::threx::DeviceVector<double> MinlinDeviceVector;
#endif


namespace pascinference {

template<class VectorBase, class MatrixBase>
class BlockGraphMatrix: public GeneralMatrix<VectorBase> {
	private:
		MatrixBase **blocks; /**< array of pointers to matrices */
		int nmb_blocks; /**< number of diagonal blocks */
		int blocksize; /**< constant size of all blocks */
		
	public:
		BlockGraphMatrix(int nmb_block, MatrixBase **new_blocks, int blocksize); /* constructor with number_of_blocks and blocks */
		~BlockGraphMatrix(); /* destructor - destroy inner matrix */

		void print(ConsoleOutput &output) const; /* print matrix */
		void printcontent(ConsoleOutput &output) const;
		std::string get_name() const;
		
		void matmult(VectorBase &y, const VectorBase &x) const; /* y = A*x */

		int get_nmb_blocks() const;
		MatrixBase **get_blocks() const;
		int get_blocksize() const;

};

template<class VectorBase, class MatrixBase>
std::string BlockGraphMatrix<VectorBase,MatrixBase>::get_name() const {
	return "BlockGraphMatrix";
}


template<class VectorBase, class MatrixBase>
BlockGraphMatrix<VectorBase,MatrixBase>::BlockGraphMatrix(int new_nmb_blocks, MatrixBase **new_blocks, int new_blocksize){
	LOG_FUNC_BEGIN
	
	this->blocks = new_blocks;
	this->nmb_blocks = new_nmb_blocks;
	this->blocksize = new_blocksize;

	LOG_FUNC_END
}	


template<class VectorBase, class MatrixBase>
BlockGraphMatrix<VectorBase,MatrixBase>::~BlockGraphMatrix(){
	LOG_FUNC_BEGIN
	
	//TODO: what to do?
	
	LOG_FUNC_END	
}

/* print matrix */
template<class VectorBase, class MatrixBase>
void BlockGraphMatrix<VectorBase,MatrixBase>::print(ConsoleOutput &output) const		
{
	LOG_FUNC_BEGIN
	
	output << "Blocks: (blocksize = " << this->blocksize << ", ";
	output << "nmb of blocks = " << this->nmb_blocks << ")" << std::endl;

	int i;
	output.push();
	for(i=0;i<this->nmb_blocks;i++){ /* print each block */
		output << "block_" << i << ": ";
		this->blocks[i]->print(output);
		output << std::endl;
	}
	output.pop();
	
	LOG_FUNC_END	
}

/* print matrix */
template<class VectorBase, class MatrixBase>
void BlockGraphMatrix<VectorBase,MatrixBase>::printcontent(ConsoleOutput &output) const		
{
	LOG_FUNC_BEGIN
	
	output << "Blocks:" << std::endl;

	int i;
	output.push();
	for(i=0;i<this->nmb_blocks;i++){ /* print each block */
		output << "block_" << i << ": " << std::endl;
		this->blocks[i]->printcontent(output);
	}
	output.pop();
	output.synchronize();
	
	LOG_FUNC_END
}

/* matrix-vector multiplication */
template<class VectorBase, class MatrixBase>
void BlockGraphMatrix<VectorBase,MatrixBase>::matmult(VectorBase &y, const VectorBase &x) const { 
	LOG_FUNC_BEGIN
	
	// TODO: maybe y is not initialized, who knows

	// TODO: write block matmult
	LOG_FUNC_END	
}

/* get number of blocks */
template<class VectorBase, class MatrixBase>
int BlockGraphMatrix<VectorBase,MatrixBase>::get_nmb_blocks() const { 
	if(DEBUG_MODE >= 100) coutMaster << "(BlockGraphMatrix)FUNCTION: get_nmb_blocks" << std::endl;

	return this->nmb_blocks;
}


/* get number of blocks */
template<class VectorBase, class MatrixBase>
MatrixBase **BlockGraphMatrix<VectorBase,MatrixBase>::get_blocks() const { 
	if(DEBUG_MODE >= 100) coutMaster << "(BlockGraphMatrix)FUNCTION: get_blocks" << std::endl;

	return this->blocks;
}

/* get blocksize */
template<class VectorBase, class MatrixBase>
int BlockGraphMatrix<VectorBase,MatrixBase>::get_blocksize() const { 
	if(DEBUG_MODE >= 100) coutMaster << "(BlockGraphMatrix)FUNCTION: get_blocksize" << std::endl;

	return this->blocksize;
}



} /* end of namespace */


#endif
