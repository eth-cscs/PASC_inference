namespace minlin {

namespace threx { // TODO: maybe choose the different namespace for my own Petsc stuff
 
/* minlin all */
PetscVectorWrapperSub PetscVector::operator()(minlin::detail::all_type add_in) {
	if(DEBUG_MODE >= 100) std::cout << "operator vec(all)" << std::endl;

	IS new_subvector_is; // TODO: this is quite stupid
	ISCreateStride(PETSC_COMM_WORLD, this->size(), 0,1, &new_subvector_is);
	
	return PetscVectorWrapperSub(this->inner_vector, new_subvector_is, true);
} 
 
 
/*! \fn PetscVector
    \brief Default constructor.
    
    No vector space allocated.
*/
PetscVector::PetscVector(){
	if(DEBUG_MODE >= 100) std::cout << "(PetscVector)CONSTRUCTOR: empty" << std::endl;

	inner_vector = NULL;
}


/* PetscVector constructor with global dimension */
PetscVector::PetscVector(int n){
	if(DEBUG_MODE >= 100) std::cout << "(PetscVector)CONSTRUCTOR: PetscVector(int)" << std::endl;

	TRY( VecCreate(PETSC_COMM_WORLD,&inner_vector) );
	TRY( VecSetSizes(inner_vector,PETSC_DECIDE,n) ); // TODO: there should be more options to set the distribution
	TRY( VecSetFromOptions(inner_vector) );

	valuesUpdate();
}

/* PetscVector copy constructor */
PetscVector::PetscVector(const PetscVector &vec1){
	if(DEBUG_MODE >= 100) std::cout << "(PetscVector)CONSTRUCTOR: PetscVector(&vec) ---- DUPLICATE ----" << std::endl;

	/* inner_vec */
//	inner_vector = vec1.inner_vector;
	// TODO: there is duplicate... this function has to be called as less as possible
	TRY( VecDuplicate(vec1.inner_vector, &inner_vector) );
	TRY( VecCopy(vec1.inner_vector, inner_vector) );
	
}

/* PetscVector constructor with inner_vector */
PetscVector::PetscVector(Vec new_inner_vector){
	if(DEBUG_MODE >= 100) std::cout << "(PetscVector)CONSTRUCTOR: PetscVector(inner_vector)" << std::endl;

	inner_vector = new_inner_vector;
}




/* PetscVector destructor */
PetscVector::~PetscVector(){
	if(DEBUG_MODE >= 100) std::cout << "(PetscVector)DESTRUCTOR" << std::endl;

	/* if there was any inner vector, then destroy it */
	if(inner_vector){
		if(DEBUG_MODE >= 100) std::cout << " - destroy inner vector" << std::endl;

		/* if petsc was finalized in the meantime, then the vector was already destroyed */
		if(PETSC_INITIALIZED){
			/* if the vector wasn't destroyed yet and the petsc is still running, then
			 * destroy the vector */
			TRY( VecDestroy(&inner_vector) );
		}
	}

}

/* after update a variable, it is necessary to call asseble begin & end */
void PetscVector::valuesUpdate(){
	if(DEBUG_MODE >= 100) std::cout << "(PetscVector)FUNCTION: valuesUpdate()" << std::endl;

	TRY( VecAssemblyBegin(inner_vector) );
	TRY( VecAssemblyEnd(inner_vector) );
}

/* set all values of the vector, this function is called from overloaded operator */
void PetscVector::set(double new_value){
	if(DEBUG_MODE >= 100) std::cout << "(PetscVector)FUNCTION: set(double)" << std::endl;

	// TODO: control if inner_vector was allocated

	TRY( VecSet(this->inner_vector,new_value) );

	valuesUpdate();
}

/* set one specific value of the vector, this function is called from overloaded operator */
void PetscVector::set(int index, double new_value){
	if(DEBUG_MODE >= 100) std::cout << "(PetscVector)FUNCTION: set(int,double)" << std::endl;

	// TODO: control if inner_vector was allocated

	TRY( VecSetValue(this->inner_vector,index,new_value, INSERT_VALUES) );
	
	valuesUpdate();
}

/* returns inner vector */
Vec PetscVector::get_vector() const { // TODO: temp
	if(DEBUG_MODE >= 100) std::cout << "(PetscVector)FUNCTION: get_vector()" << std::endl;
		
	return inner_vector;
}

/* get size of the vector */
int PetscVector::size(){
	if(DEBUG_MODE >= 100) std::cout << "(PetscVector)FUNCTION: size()" << std::endl;

	int global_size;

	TRY( VecGetSize(inner_vector,&global_size) );

	return global_size;
}

/* get single value with given id of the vector (works only with local id), really slow */
double PetscVector::get(int i)
{
	if(DEBUG_MODE >= 100) std::cout << "(PetscVector)FUNCTION: get(int)" << std::endl;

	PetscInt ni = 1;
	PetscInt ix[1];
	PetscScalar y[1];
			
	ix[0] = i;

	TRY( VecGetValues(inner_vector,ni,ix,y) );
			
	return y[0];
}

void PetscVector::get_array(double **arr){
	if(DEBUG_MODE >= 100) std::cout << "(PetscVector)FUNCTION: get_array(double **)" << std::endl;

	TRY( VecGetArray(inner_vector,arr) );
}

void PetscVector::restore_array(double **arr){
	if(DEBUG_MODE >= 100) std::cout << "(PetscVector)FUNCTION: restore_array(double **)" << std::endl;

	TRY( VecRestoreArray(inner_vector,arr) );
}

void PetscVector::get_ownership(int *low, int *high){
	if(DEBUG_MODE >= 100) std::cout << "(PetscVector)FUNCTION: get_ownership(int*, int*)" << std::endl;

	//TODO: control inner_vector

	TRY( VecGetOwnershipRange(inner_vector, low, high) );
}

/* inner_vector = alpha*inner_vector */
void PetscVector::scale(PetscScalar alpha){
	if(DEBUG_MODE >= 100) std::cout << "(PetscVector)FUNCTION: scale(double)" << std::endl;

	//TODO: control inner_vector

	TRY( VecScale(inner_vector, alpha) );
	valuesUpdate(); // TODO: has to be called?
}


/* stream insertion << operator */
std::ostream &operator<<(std::ostream &output, const PetscVector &vector)		
{
	if(DEBUG_MODE >= 100) std::cout << "(PetscVector)OPERATOR: <<" << std::endl;

	PetscScalar *arr_vector;
	PetscInt i,local_size;

	output << "[";
	TRY( VecGetLocalSize(vector.inner_vector,&local_size) );
	TRY( VecGetArray(vector.inner_vector,&arr_vector) );
	for (i=0; i<local_size; i++){
		output << arr_vector[i];
		if(i < local_size-1) output << ", ";
	}
	TRY( VecRestoreArray(vector.inner_vector,&arr_vector) );
	output << "]";
			
	return output;
}

/* vec1 = vec2, assignment operator (set vector) */
PetscVector &PetscVector::operator=(const PetscVector &vec2){
	if(DEBUG_MODE >= 100) std::cout << "(PetscVector)OPERATOR: (vec = vec)" << std::endl;

	/* check for self-assignment by comparing the address of the implicit object and the parameter */
	/* vec1 = vec1 */
    if (this == &vec2){
		if(DEBUG_MODE >= 100) std::cout << " - self assignment" << std::endl;		
        return *this;
	}

	/* vec1 is not initialized yet */
	if (!inner_vector){
		if(DEBUG_MODE >= 100) std::cout << " - creating new vector" << std::endl;		
		VecDuplicate(vec2.inner_vector,&(this->inner_vector));
		this->valuesUpdate(); // TODO: has to be called?
	}

	/* else copy the values of inner vectors */
	if(DEBUG_MODE >= 100) std::cout << " - copy values" << std::endl;		
	
	VecCopy(vec2.inner_vector,inner_vector);
	this->valuesUpdate(); // TODO: has to be called?
	
	return *this;	
}

/* vec1 = linear_combination_node, perform simple linear combination */
PetscVector &PetscVector::operator=(PetscVectorWrapperCombNode combnode){
	if(DEBUG_MODE >= 100) std::cout << "(PetscVector)OPERATOR: (vec = combnode)" << std::endl;

	/* vec1 = alpha*vec1 => simple scale */
    if (this->inner_vector == combnode.get_vector()){
		if(DEBUG_MODE >= 100) std::cout << " - scale values" << std::endl;		

        this->scale(combnode.get_coeff());
        return *this;
	}

	/* vec1 is not initialized yet */
	if (!inner_vector){
		if(DEBUG_MODE >= 100) std::cout << " - duplicate vector" << std::endl;		

		TRY( VecDuplicate(combnode.get_vector(),&inner_vector) );
	}

	/* else copy the vector values and then scale */
	if(DEBUG_MODE >= 100) std::cout << " - copy values" << std::endl;		

	TRY( VecCopy(combnode.get_vector(),inner_vector));
	
    this->scale(combnode.get_coeff());

	return *this;	
}

/* vec1 = linear_combination, perform full linear combination */
PetscVector &PetscVector::operator=(PetscVectorWrapperComb comb){
	if(DEBUG_MODE >= 100) std::cout << "(PetscVector)OPERATOR: (vec = comb)" << std::endl;

	/* vec1 is not initialized yet */
	if (!inner_vector){
		if(DEBUG_MODE >= 100) std::cout << " - duplicate vector" << std::endl;		

		TRY( VecDuplicate(comb.get_first_vector(),&inner_vector) );
	}

	/* vec1 = 0, we will perform MAXPY (y += lin_comb) */
	if(DEBUG_MODE >= 100) std::cout << " - set vector LHS = 0" << std::endl;		
	TRY( VecSet(inner_vector,0.0) );

	/* vec += comb */
	*this += comb; 

	return *this;	
}

/* vec1 = scalar_value <=> vec1(all) = scalar_value, assignment operator */
PetscVector &PetscVector::operator=(double scalar_value){
	if(DEBUG_MODE >= 100) std::cout << "(PetscVector)OPERATOR: (vec = double)" << std::endl;

	this->set(scalar_value);
	return *this;	
}

/* return subvector to be able to overload vector(index) = new_value */ 
PetscVectorWrapperSub PetscVector::operator()(int index)
{   
	if(DEBUG_MODE >= 100) std::cout << "(PetscVector)OPERATOR: (int) returns WrapperSub" << std::endl;
	
	/* create new indexset */
	if(DEBUG_MODE >= 100) std::cout << " - create index set" << std::endl;
	IS new_subvector_is;
	PetscInt idxs[1];
	idxs[0] = index;
	
	TRY( ISCreateGeneral(PETSC_COMM_WORLD, 1, idxs, PETSC_COPY_VALUES, &new_subvector_is));
	// TODO: when to delete this index set?
	
	return PetscVectorWrapperSub(this->inner_vector, new_subvector_is, true);
}

/* return subvector vector(index_begin:index_end), i.e. components with indexes: [index_begin, index_begin+1, ..., index_end] */ 
PetscVectorWrapperSub PetscVector::operator()(int index_begin, int index_end)
{   
	if(DEBUG_MODE >= 100) std::cout << "(PetscVector)OPERATOR: (int,int) returns WrapperSub" << std::endl;

	/* create new indexset */
	if(DEBUG_MODE >= 100) std::cout << " - create index set" << std::endl;
	IS new_subvector_is;
	TRY( ISCreateStride(PETSC_COMM_WORLD, index_end-index_begin+1, index_begin,1, &new_subvector_is) );
	// TODO: when to delete this index set?
		
	return PetscVectorWrapperSub(this->inner_vector, new_subvector_is, true);
}

/* return subvector based on provided index set */ 
PetscVectorWrapperSub PetscVector::operator()(const IS new_subvector_is) // TODO: are ju sure that this IS will be nt destroyed?
{   
	if(DEBUG_MODE >= 100) std::cout << "(PetscVector)OPERATOR: vec(IS)" << std::endl;
	
	return PetscVectorWrapperSub(inner_vector,new_subvector_is, false);
}



/* vec1 *= alpha */
void operator*=(PetscVector &vec1, double alpha)
{
	if(DEBUG_MODE >= 100) std::cout << "(PetscVector)OPERATOR: vec *= double" << std::endl;
	
	vec1.scale(alpha);
}

/* vec1 += comb */
void operator+=(PetscVector &vec1, PetscVectorWrapperComb comb)
{
	if(DEBUG_MODE >= 100) std::cout << "(PetscVector)OPERATOR: vec += comb" << std::endl;
	
	int list_size = comb.get_listsize();
	PetscScalar alphas[list_size];
	Vec vectors[list_size];

	/* get array with coefficients and vectors */
	comb.get_arrays(alphas,vectors);

	/* vec1 = vec1 + sum (coeff*vector) */
	if(DEBUG_MODE >= 100) std::cout << " - perform MAXPY" << std::endl;
	VecMAXPY(vec1.inner_vector,list_size,alphas,vectors);
	vec1.valuesUpdate();

}

/* vec1 -= comb */
void operator-=(PetscVector &vec1, PetscVectorWrapperComb comb)
{
	if(DEBUG_MODE >= 100) std::cout << "(PetscVector)OPERATOR: vec -= comb" << std::endl;
	
	vec1 += (-1.0)*comb;
}

/* dot = dot(vec1,vec2) */
double dot(const PetscVector vec1, const PetscVector vec2)
{
	if(DEBUG_MODE >= 100) std::cout << "(PetscVector)FUNCTION: dot(vec1,vec2)" << std::endl;

	double dot_value;
	TRY( VecDot(vec1.inner_vector,vec2.inner_vector,&dot_value));
	return dot_value;
}

/* max = max(vec1) */
double max(const PetscVector vec1)
{
	if(DEBUG_MODE >= 100) std::cout << "(PetscVector)FUNCTION: max(vec)" << std::endl;

	double max_value;
	TRY( VecMax(vec1.inner_vector,NULL, &max_value) );
	return max_value;
}

/* sum = sum(vec1) */
double sum(const PetscVector vec1)
{
	if(DEBUG_MODE >= 100) std::cout << "(PetscVector)FUNCTION: sum(vec)" << std::endl;

	double sum_value;
	TRY( VecSum(vec1.inner_vector,&sum_value) );
	return sum_value;
}

/* vec3 = vec1./vec2 */
const PetscVector operator/(PetscVector vec1, const PetscVector vec2)
{
	if(DEBUG_MODE >= 100) std::cout << "(PetscVector)FUNCTION: vec1/vec2" << std::endl;

	TRY(VecPointwiseDivide(vec1.inner_vector,vec1.inner_vector,vec2.inner_vector) );

	vec1.valuesUpdate(); // TODO: has to be called?
	return vec1;
}




} /* end of namespace */

} /* end of MinLin namespace */

