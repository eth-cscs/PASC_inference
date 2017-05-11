/** @file blockgraphsparse.h
 *  @brief block graph matrix with tridiag blocks used in GRAPHH1FEM model with matrix-vector multiplication implemented as sparse-matrix operation
 *
 *  @author Lukas Pospisil
 */

#ifndef PASC_BLOCKGRAPHSPARSEMATRIX_H
#define	PASC_BLOCKGRAPHSPARSEMATRIX_H

#include "general/common/decomposition.h"
#include "general/algebra/graph/bgmgraph.h"

namespace pascinference {
using namespace common;

namespace algebra {

/** \class BlockGraphSparseMatrix
 *  \brief sparse graph-based matrix for FEM-H1 model regularisation
 *
*/
template<class VectorBase>
class BlockGraphSparseMatrix: public GeneralMatrix<VectorBase> {
	private:
		Decomposition<VectorBase> *decomposition;

		double alpha; /**< general matrix multiplicator */

		class ExternalContent;
		friend class ExternalContent;
		ExternalContent *externalcontent;			/**< for manipulation with external-specific stuff */

		GeneralVector<VectorBase> *coeffs; /**< vector of coefficient for each block */

	public:
		BlockGraphSparseMatrix(Decomposition<VectorBase> &decomposition, double alpha=1.0, GeneralVector<VectorBase> *new_coeffs=NULL);
		~BlockGraphSparseMatrix(); /* destructor - destroy inner matrix */

		void print(ConsoleOutput &output) const; /* print matrix */
		void print(ConsoleOutput &output_global, ConsoleOutput &output_local) const;
		
		void printcontent(ConsoleOutput &output) const;
		std::string get_name() const;
		
		void matmult(VectorBase &y, const VectorBase &x) const; /* y = A*x */

		int get_R() const;
		int get_Rlocal() const;
		int get_K() const;
		int get_T() const;
		int get_Tlocal() const;
		double get_coeff() const;
		void set_coeff(double coeff);

//TODO: this should be NOT here: (rather return ExternalContent)
#ifdef USE_PETSC
		Mat get_petscmatrix() const;
#endif

};

/* ---------------- IMPLEMENTATION -------------------- */

template<class VectorBase>
std::string BlockGraphSparseMatrix<VectorBase>::get_name() const {
	return "BlockGraphSparseMatrix";
}


template<class VectorBase>
BlockGraphSparseMatrix<VectorBase>::BlockGraphSparseMatrix(Decomposition<VectorBase> &new_decomposition, double alpha, GeneralVector<VectorBase> *new_coeffs){
	LOG_FUNC_BEGIN

	LOG_FUNC_END
}	


template<class VectorBase>
BlockGraphSparseMatrix<VectorBase>::~BlockGraphSparseMatrix(){
	LOG_FUNC_BEGIN
	
	LOG_FUNC_END	
}

/* print matrix */
template<class VectorBase>
void BlockGraphSparseMatrix<VectorBase>::print(ConsoleOutput &output) const		
{
	LOG_FUNC_BEGIN

	output << this->get_name() << std::endl;
	output << " - T:     " << get_T() << std::endl;
	output << " - R:     " << get_R() << std::endl;
	output << " - K:     " << get_K() << std::endl;
	output << " - size:  " << get_T()*get_R()*get_K() << std::endl;
	output << " - alpha: " << alpha << std::endl;

	if(coeffs){
		output << " - coeffs: " << *coeffs << std::endl;
	}
	
	LOG_FUNC_END	
}

template<class VectorBase>
void BlockGraphSparseMatrix<VectorBase>::print(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN

	output_global << this->get_name() << std::endl;
	output_global << " - T:     " << get_T() << std::endl;
	output_global.push();
		output_local << " - Tlocal:  " << get_Tlocal() << " (" << decomposition->get_Tbegin() << "," << decomposition->get_Tend() << ")" << std::endl;
		output_local.synchronize();
	output_global.pop();
	
	output_global << " - R:     " << get_R() << std::endl;
	output_global.push();
		output_local << " - Rlocal:  " << get_Rlocal() << " (" << decomposition->get_Rbegin() << "," << decomposition->get_Rend() << ")" << std::endl;
		output_local.synchronize();
	output_global.pop();

	output_global << " - K:     " << get_K() << std::endl;
	output_global << " - size:  " << get_T()*get_R()*get_K() << std::endl;
	output_global.push();
		output_local << " - sizelocal:  " << get_Tlocal()*get_Rlocal()*get_K() << std::endl;
		output_local.synchronize();
	output_global.pop();

	output_global << " - alpha: " << alpha << std::endl;

	if(coeffs){
		output_local << " - coeffs: " << *coeffs << std::endl;
		output_local.synchronize();
	}

	LOG_FUNC_END	
}

/* print matrix */
template<class VectorBase>
void BlockGraphSparseMatrix<VectorBase>::printcontent(ConsoleOutput &output) const		
{
	LOG_FUNC_BEGIN
	
	LOG_FUNC_END
}

/* matrix-vector multiplication */
template<class VectorBase>
void BlockGraphSparseMatrix<VectorBase>::matmult(VectorBase &y, const VectorBase &x) const { 
	LOG_FUNC_BEGIN

	LOG_FUNC_END
}

template<class VectorBase>
int BlockGraphSparseMatrix<VectorBase>::get_R() const { 
	return decomposition->get_R();
}

template<class VectorBase>
int BlockGraphSparseMatrix<VectorBase>::get_K() const { 
	return decomposition->get_K();
}

template<class VectorBase>
int BlockGraphSparseMatrix<VectorBase>::get_T() const { 
	return decomposition->get_T();
}

template<class VectorBase>
int BlockGraphSparseMatrix<VectorBase>::get_Tlocal() const { 
	return decomposition->get_Tlocal();
}

template<class VectorBase>
int BlockGraphSparseMatrix<VectorBase>::get_Rlocal() const { 
	return decomposition->get_Rlocal();
}

template<class VectorBase>
double BlockGraphSparseMatrix<VectorBase>::get_coeff() const {
	return this->alpha;
}

template<class VectorBase>
void BlockGraphSparseMatrix<VectorBase>::set_coeff(double coeff) {
	this->alpha = coeff;
}

//TODO: this should be somewhere else
#ifdef USE_PETSC
template<class VectorBase>
Mat BlockGraphSparseMatrix<VectorBase>::get_petscmatrix() const {
	return this->A_petsc;
}
#endif

}
} /* end of namespace */


#endif
