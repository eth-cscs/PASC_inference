#ifndef PETSCVECTOR_WRAPPERSUB_IMPL_H
#define	PETSCVECTOR_WRAPPERSUB_IMPL_H


namespace petscvector {

/* PetscVectorWrapperSub constructor with given IS = create subvector */
PetscVectorWrapperSub::PetscVectorWrapperSub(Vec new_inner_vector, IS new_subvector_is, bool new_free_is){
	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << "(WrapperSub)CONSTRUCTOR: WrapperSub(inner_vec, IS)" << std::endl;

	/* free index set during destruction ? */
	free_is = new_free_is;

	/* store old inner vector - will be used in the destructor to return subvector */
	inner_vector = new_inner_vector; 

	/* copy IS */
	subvector_is = new_subvector_is;

	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << " - get subvector from original vector" << std::endl;

	/* get subvector, restore it in destructor */
	TRY( VecGetSubVector(inner_vector, subvector_is, &subvector) );
	
}

/* PetscVectorWrapperSub destructor */
PetscVectorWrapperSub::~PetscVectorWrapperSub(){
	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << "(WrapperSub)DESTRUCTOR: ~WrapperSub" << std::endl;

	/* if this was a subvector, then restore values */
	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << " - restore subvector" << std::endl;
	TRY( VecRestoreSubVector(inner_vector, subvector_is, &subvector) );

	/* if it is necessary to free IS, then free it */
	if(free_is){
		if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << " - destroy IS" << std::endl;
		TRY( ISDestroy(&subvector_is) );
	}

}

/* set all values of the subvector, this function is called from overloaded operator */
void PetscVectorWrapperSub::set(double new_value){
	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << "(WrapperSub)FUNCTION: set(double)" << std::endl;

	// TODO: control if subvector was allocated

	TRY( VecSet(this->subvector,new_value) );

	valuesUpdate();
}

/* return vector from this node */
Vec PetscVectorWrapperSub::get_subvector(){
	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << "(WrapperSub)FUNCTION: get_subvector()" << std::endl;	
	
	return this->subvector;
}

/* get single value with given id of the vector (works only with local id), really slow */
double PetscVectorWrapperSub::get(int i)
{
	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << "(WrapperSub)FUNCTION: get(int)" << std::endl;

	PetscInt ni = 1;
	PetscInt ix[1];
	PetscScalar y[1];
			
	ix[0] = i;

	TRY( VecGetValues(subvector,ni,ix,y) );
			
	return y[0];
}

/* after update a variable, it is necessary to call asseble begin & end */
void PetscVectorWrapperSub::valuesUpdate() const {
	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << "(WrapperSub)FUNCTION: valuesUpdate()" << std::endl;

	TRY( VecAssemblyBegin(subvector) );
	TRY( VecAssemblyEnd(subvector) );
}

/* subvector = alpha*subvector */
void PetscVectorWrapperSub::scale(PetscScalar alpha) const{
	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << "(PetscVector)FUNCTION: scale(double)" << std::endl;

	//TODO: control subvector

	TRY( VecScale(subvector, alpha) );
	valuesUpdate(); // TODO: has to be called?
}

/* stream insertion << operator */
std::ostream &operator<<(std::ostream &output, const PetscVectorWrapperSub &wrapper)		
{
	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << "(WrapperSub)OPERATOR: <<" << std::endl;

	PetscScalar *arr_vector;
	PetscInt i,local_size;
	const PetscInt *indices;

	output << "[";
	TRY( VecGetLocalSize(wrapper.subvector,&local_size) );

	TRY( VecGetArray(wrapper.subvector,&arr_vector) );
	TRY( ISGetIndices(wrapper.subvector_is,&indices) );

	for (i=0; i<local_size; i++){
		output << "{" << indices[i] << "}=" << arr_vector[i];
		if(i < local_size-1) output << ", ";
	}

	TRY( ISRestoreIndices(wrapper.subvector_is,&indices) );
	TRY( VecRestoreArray(wrapper.subvector,&arr_vector) );
	output << "]";
			
	return output;
}


/* subvec = scalar_value <=> subvec(all) = scalar_value, assignment operator */
PetscVectorWrapperSub &PetscVectorWrapperSub::operator=(double scalar_value){
	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << "(WrapperSub)OPERATOR: (subvec = double)" << std::endl;

	this->set(scalar_value);
	return *this;	
}

/* subvec1 = subvec2, assignment operator (set subvector) */
PetscVectorWrapperSub &PetscVectorWrapperSub::operator=(PetscVectorWrapperSub subvec2){
	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << "(WrapperSub)OPERATOR: (subvec = subvec)" << std::endl;

	/* vec1 is not initialized yet */
	if (!subvector){
		//TODO: give error
	}

	/* else copy the values of inner vectors */
	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << " - copy values" << std::endl;		
	
	VecCopy(subvec2.get_subvector(),subvector);
	this->valuesUpdate(); // TODO: has to be called?
	
	return *this;	
}


/* subvec1 = vec2, assignment operator (set vector) */
PetscVectorWrapperSub &PetscVectorWrapperSub::operator=(const PetscVector &vec2){
	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << "(WrapperSub)OPERATOR: (subvec = vec)" << std::endl;

	/* vec1 is not initialized yet */
	if (!subvector){
		//TODO: give error
	}

	/* else copy the values of inner vectors */
	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << " - copy values" << std::endl;		
	
	VecCopy(vec2.get_vector(),subvector); // TODO: I dont know how to do without this
	this->valuesUpdate(); // TODO: has to be called?
	
	return *this;	
}



/* vec1 = linear_combination_node, perform simple linear combination */
//PetscVectorWrapperSub &PetscVectorWrapperSub::operator=(const PetscVectorWrapperCombNode combnode){
//	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << "(WrapperSub)OPERATOR: (subvec = combnode)" << std::endl;

	/* vec1 is not initialized yet */
//	if (!subvector){
		//TODO: give error
//	}

	/* else copy the vector values and then scale */
//	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << " - copy values" << std::endl;		

//	TRY( VecCopy(combnode.get_vector(),subvector));
	
//    this->scale(combnode.get_coeff());

//	return *this;	
//}

/* vec1 = linear_combination, perform full linear combination */
PetscVectorWrapperSub &PetscVectorWrapperSub::operator=(PetscVectorWrapperComb comb){
	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << "(WrapperSub)OPERATOR: (subvec = comb)" << std::endl;

	/* vec1 is not initialized yet */
	if (!subvector){
		//TODO: give error
	}

	/* vec = comb */
	comb.compute(subvector,0.0);

	return *this;	
}

/* subvec *= alpha */
void operator*=(const PetscVectorWrapperSub &subvec1, double alpha)
{
	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << "(PetscVector)OPERATOR: vec *= double" << std::endl;
	
	subvec1.scale(alpha);
}

/* vec1 += comb */
void operator+=(const PetscVectorWrapperSub &subvec, PetscVectorWrapperComb comb)
{
	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << "(PetscVector)OPERATOR: vec += comb" << std::endl;
	
	comb.compute(subvec.subvector,1.0);
}

/* subvec1 -= comb */
void operator-=(const PetscVectorWrapperSub &subvec1, PetscVectorWrapperComb comb)
{
	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << "(WrapperSub)OPERATOR: subvec -= comb" << std::endl;
	
	subvec1 += (-1.0)*comb;
}

/* vec1 = vec1./subvec2 */
void operator/=(const PetscVectorWrapperSub &subvec1, const PetscVectorWrapperSub subvec2)
{
	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << "(WrapperSub)FUNCTION: vec1/vec2" << std::endl;

	TRY(VecPointwiseDivide(subvec1.subvector,subvec1.subvector,subvec2.subvector) );

	subvec1.valuesUpdate(); // TODO: has to be called?

}


/* subvec == scalar */
bool operator==(PetscVectorWrapperSub subvec, double alpha){
	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << "(WrapperSub)OPERATOR: subvec == double" << std::endl;	
	
	bool return_value = false; // TODO: works only with vector of size 1, otherwise compare only first value
	double vector_value = subvec.get(0);
	
	if(vector_value == alpha){
		return_value = true;
	}	
	
	return return_value;
}

/* vec1 == vec2 */
bool operator==(PetscVectorWrapperSub subvec1, PetscVectorWrapperSub subvec2){
	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << "(WrapperSub)OPERATOR: subvec1 == subvec2" << std::endl;

	PetscBool return_value;

	TRY( VecEqual(subvec1.subvector,subvec2.subvector,&return_value) );
	
	return (bool)return_value;
}

/* subvec1 > subvec2 */
bool operator>(PetscVectorWrapperSub subvec1, PetscVectorWrapperSub subvec2){
	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << "(WrapperSub)OPERATOR: subvec1 > subvec2" << std::endl;

	bool return_value = false; // TODO: works only with vector of size 1, otherwise compare only first value
	
	double vec1_value = subvec1.get(0);
	double vec2_value = subvec2.get(0);
	
	if(vec1_value > vec2_value){
		return_value = true;
	}	
	
	return return_value;

}

/* sum = sum(subvec1) */
double sum(const PetscVectorWrapperSub subvec1)
{
	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << "(WrapperSub)FUNCTION: sum(subvec)" << std::endl;

	double sum_value;
	TRY( VecSum(subvec1.subvector,&sum_value) );
	return sum_value;
}

/* dot = dot(subvec1,subvec2) */
double dot(const PetscVectorWrapperSub subvec1, const PetscVectorWrapperSub subvec2)
{
	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << "(WrapperSub)FUNCTION: dot(subvec1,subvec2)" << std::endl;

	double dot_value;
	TRY( VecDot(subvec1.subvector,subvec2.subvector,&dot_value));
	return dot_value;
}

double dot(const PetscVector &x, const PetscVectorWrapperSub y)
{
	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << "(WrapperSub)FUNCTION: dot(vec,subvec)" << std::endl;

	double dot_value;
	TRY( VecDot(x.inner_vector,y.subvector,&dot_value));
	return dot_value;
}

double dot(const PetscVectorWrapperSub x, const PetscVector &y)
{
	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << "(WrapperSub)FUNCTION: dot(subvec,vec)" << std::endl;

	double dot_value;
	TRY( VecDot(y.inner_vector,x.subvector,&dot_value));
	return dot_value;
}


PetscVectorWrapperMul mul(PetscVectorWrapperSub subvec1, PetscVectorWrapperSub subvec2)
{
	return PetscVectorWrapperMul( subvec1.subvector, subvec2.subvector);
}

} /* end of petscvector namespace */

#endif
