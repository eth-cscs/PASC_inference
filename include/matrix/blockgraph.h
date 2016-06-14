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
		int T; /**< dimension of each block */
		int R; /**< number of vertices = number of blocks in row,col */
		double alpha; /**< general matrix multiplicator */
		
	public:
		BlockGraphMatrix(int T, int R, double alpha=1.0);
		~BlockGraphMatrix(); /* destructor - destroy inner matrix */

		void print(ConsoleOutput &output) const; /* print matrix */
		void printcontent(ConsoleOutput &output) const;
		std::string get_name() const;
		
		void matmult(VectorBase &y, const VectorBase &x) const; /* y = A*x */

		int get_R() const;
		int get_T() const;
		double get_alpha() const;

};

template<class VectorBase, class MatrixBase>
std::string BlockGraphMatrix<VectorBase,MatrixBase>::get_name() const {
	return "BlockGraphMatrix";
}


template<class VectorBase, class MatrixBase>
BlockGraphMatrix<VectorBase,MatrixBase>::BlockGraphMatrix(int T, int R, double alpha){
	LOG_FUNC_BEGIN
	
	this->T = T;
	this->R = R;
	this->alpha = alpha;

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
	
	output << "T:     " << T << std::endl;
	output << "R:     " << R << std::endl;
	output << "alpha: " << alpha << std::endl;
	
	LOG_FUNC_END	
}

/* print matrix */
template<class VectorBase, class MatrixBase>
void BlockGraphMatrix<VectorBase,MatrixBase>::printcontent(ConsoleOutput &output) const		
{
	LOG_FUNC_BEGIN
	
	
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

template<class VectorBase, class MatrixBase>
int BlockGraphMatrix<VectorBase,MatrixBase>::get_R() const { 
	return this->R;
}

template<class VectorBase, class MatrixBase>
int BlockGraphMatrix<VectorBase,MatrixBase>::get_T() const { 
	return this->T;
}

template<class VectorBase, class MatrixBase>
double BlockGraphMatrix<VectorBase,MatrixBase>::get_alpha() const { 
	return this->alpha;
}




} /* end of namespace */


#endif
