#include "algebra/vector/generalvector.h"

namespace pascinference {
namespace algebra {

template<class VectorBase>
GeneralVector<VectorBase> & GeneralVector<VectorBase>::operator=(GeneralMatrixRHS<VectorBase> rhs){
	rhs.matmult(*this);
	return *this;
}

template<class VectorBase>
void GeneralVector<VectorBase>::set_random(){
}

template<class VectorBase>
void GeneralVector<VectorBase>::set_random2(){
}
			
template<class VectorBase>
std::string GeneralVector<VectorBase>::get_name(){
	return "GeneralVector";
}

			
}
}
