#ifndef PETSCVECTOR_WRAPPERMUL_IMPL_H
#define	PETSCVECTOR_WRAPPERMUL_IMPL_H


namespace petscvector {

/* PetscVectorWrapperSub constructor with given IS = create subvector */
PetscVectorWrapperMul::PetscVectorWrapperMul(Vec new_inner_vector1, Vec new_inner_vector2){
	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << "(WrapperMul)CONSTRUCTOR: WrapperMul(inner_vec1, inner_vec2)" << std::endl;

	inner_vector1 = new_inner_vector1; 
	inner_vector2 = new_inner_vector2; 
}

/* PetscVectorWrapperSub destructor */
PetscVectorWrapperMul::~PetscVectorWrapperMul(){

}

/* set all values of the subvector, this function is called from overloaded operator */
void PetscVectorWrapperMul::mul(Vec result){
	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << "(WrapperSub)FUNCTION: mul(Vec result)" << std::endl;

	// TODO: control if vectors were allocated
	TRY( VecPointwiseMult(result, inner_vector1, inner_vector2) );


}

Vec PetscVectorWrapperMul::get_vector1() const {
	return inner_vector1;
}

Vec PetscVectorWrapperMul::get_vector2() const {
	return inner_vector2;
}


} /* end of petscvector namespace */

#endif
