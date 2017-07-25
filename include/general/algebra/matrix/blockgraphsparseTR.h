/** @file blockgraphsparseTR.h
 *  @brief block graph matrix for time-space regularization
 *
 *  @author Lukas Pospisil
 */

#ifndef PASC_BLOCKGRAPHSPARSETRMATRIX_H
#define	PASC_BLOCKGRAPHSPARSETRMATRIX_H

#include "general/common/decomposition.h"
#include "general/algebra/graph/bgmgraph.h"

namespace pascinference {
using namespace common;

namespace algebra {

/** \class BlockGraphSparseTRMatrix
 *  \brief sparse graph-based matrix for FEM-H1 model regularization
 *
*/
template<class VectorBase>
class BlockGraphSparseTRMatrix: public GeneralMatrix<VectorBase> {
	public:
		class ExternalContent;	/**< class which includes external content, has to be public because of get_externalcontent() return type */

	private:
		Decomposition<VectorBase> *decomposition;

		double alpha; /**< general matrix multiplicator */
		double sigma; /**< time-space coefficient */

		friend class ExternalContent;
		ExternalContent *externalcontent;			/**< for manipulation with external-specific stuff */

	public:
		BlockGraphSparseTRMatrix(Decomposition<VectorBase> &decomposition, double sigma=0.5, double alpha=1.0);
		~BlockGraphSparseTRMatrix(); /* destructor - destroy inner matrix */

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

		double get_sigma() const;

		ExternalContent *get_externalcontent() const;

};

/* ---------------- IMPLEMENTATION -------------------- */

template<class VectorBase>
std::string BlockGraphSparseTRMatrix<VectorBase>::get_name() const {
	return "BlockGraphSparseTRMatrix";
}


template<class VectorBase>
BlockGraphSparseTRMatrix<VectorBase>::BlockGraphSparseTRMatrix(Decomposition<VectorBase> &new_decomposition, double sigma, double alpha){
	LOG_FUNC_BEGIN

	LOG_FUNC_END
}	


template<class VectorBase>
BlockGraphSparseTRMatrix<VectorBase>::~BlockGraphSparseTRMatrix(){
	LOG_FUNC_BEGIN
	
	LOG_FUNC_END	
}

/* print matrix */
template<class VectorBase>
void BlockGraphSparseTRMatrix<VectorBase>::print(ConsoleOutput &output) const		
{
	LOG_FUNC_BEGIN

	output << this->get_name() << std::endl;
	output << " - T     : " << get_T() << std::endl;
	output << " - R     : " << get_R() << std::endl;
	output << " - K     : " << get_K() << std::endl;
	output << " - size  : " << get_T()*get_R()*get_K() << std::endl;
	output << " - sigma : " << sigma << std::endl;
	output << " - alpha : " << alpha << std::endl;
	
	LOG_FUNC_END	
}

template<class VectorBase>
void BlockGraphSparseTRMatrix<VectorBase>::print(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
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

	output_global << " - sigma: " << sigma << std::endl;
	output_global << " - alpha: " << alpha << std::endl;

	LOG_FUNC_END	
}

/* print matrix */
template<class VectorBase>
void BlockGraphSparseTRMatrix<VectorBase>::printcontent(ConsoleOutput &output) const		
{	LOG_FUNC_BEGIN
	
	LOG_FUNC_END
}

/* matrix-vector multiplication */
template<class VectorBase>
void BlockGraphSparseTRMatrix<VectorBase>::matmult(VectorBase &y, const VectorBase &x) const { 
	LOG_FUNC_BEGIN

	LOG_FUNC_END
}

template<class VectorBase>
int BlockGraphSparseTRMatrix<VectorBase>::get_R() const { 
	return decomposition->get_R();
}

template<class VectorBase>
int BlockGraphSparseTRMatrix<VectorBase>::get_K() const { 
	return decomposition->get_K();
}

template<class VectorBase>
int BlockGraphSparseTRMatrix<VectorBase>::get_T() const { 
	return decomposition->get_T();
}

template<class VectorBase>
int BlockGraphSparseTRMatrix<VectorBase>::get_Tlocal() const { 
	return decomposition->get_Tlocal();
}

template<class VectorBase>
int BlockGraphSparseTRMatrix<VectorBase>::get_Rlocal() const { 
	return decomposition->get_Rlocal();
}

template<class VectorBase>
double BlockGraphSparseTRMatrix<VectorBase>::get_coeff() const {
	return this->alpha;
}

template<class VectorBase>
void BlockGraphSparseTRMatrix<VectorBase>::set_coeff(double coeff) {
	this->alpha = coeff;
}

template<class VectorBase>
double BlockGraphSparseTRMatrix<VectorBase>::get_sigma() const {
	return this->sigma;
}

}
} /* end of namespace */


#endif
