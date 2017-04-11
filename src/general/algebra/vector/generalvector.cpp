#include "algebra/vector/generalvector.h"

namespace pascinference {
namespace algebra {

template<class VectorBase>
GeneralVector<VectorBase> & GeneralVector<VectorBase>::operator=(GeneralMatrixRHS<VectorBase> rhs){
	rhs.matmult(*this);
	return *this;
}
			
template<class VectorBase>
GeneralVector<VectorBase>::GeneralVector(): VectorBase() {
}

template<class VectorBase>
template<class ArgType>
GeneralVector<VectorBase>::GeneralVector(ArgType arg): VectorBase(arg) {
}

template<class VectorBase>
template<class ArgType1, class ArgType2>
GeneralVector<VectorBase>::GeneralVector(ArgType1 arg1, ArgType2 arg2): VectorBase(arg1,arg2) {
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
