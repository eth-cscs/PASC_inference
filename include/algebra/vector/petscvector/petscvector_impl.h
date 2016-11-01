#ifndef PETSCVECTOR_IMPL_H
#define	PETSCVECTOR_IMPL_H


namespace petscvector {

PetscVector::PetscVector(){
	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << "(PetscVector)CONSTRUCTOR: empty" << std::endl;

	inner_vector = NULL;
}


PetscVector::PetscVector(int n){
	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << "(PetscVector)CONSTRUCTOR: PetscVector(int)" << std::endl;

	TRY( VecCreate(PETSC_COMM_WORLD,&inner_vector) );
	TRY( VecSetSizes(inner_vector,PETSC_DECIDE,n) ); // TODO: there should be more options to set the distribution
	TRY( VecSetFromOptions(inner_vector) );

	valuesUpdate();
}


PetscVector::PetscVector(double *values, int n){
	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << "(PetscVector)CONSTRUCTOR: PetscVector(values, int)" << std::endl;

	TRY( VecCreateSeqWithArray(PETSC_COMM_SELF, 1, n, values, &inner_vector ) );
	TRY( VecSetFromOptions(inner_vector) );

	valuesUpdate();
}


PetscVector::PetscVector(const PetscVector &vec){
	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << "(PetscVector)CONSTRUCTOR: PetscVector(&vec) ---- DUPLICATE ----" << std::endl;

	/* there is duplicate... this function has to be called as less as possible */
	TRY( VecDuplicate(vec.inner_vector, &inner_vector) );
	TRY( VecCopy(vec.inner_vector, inner_vector) );
	
}


PetscVector::PetscVector(const Vec &new_inner_vector){
	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << "(PetscVector)CONSTRUCTOR: PetscVector(inner_vector)" << std::endl;

	this->inner_vector = new_inner_vector;
}


PetscVector::PetscVector(const PetscVectorWrapperComb &comb){
	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << "(PetscVector)CONSTRUCTOR: PetscVector(comb)" << std::endl;

	inner_vector = NULL;
	*this = comb; /* assemble the linear combination */

}


PetscVector::~PetscVector(){
	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << "(PetscVector)DESTRUCTOR" << std::endl;

	/* if there is any inner vector, then destroy it */
	if(inner_vector){
		if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << " - destroy inner vector" << std::endl;

		/* if petsc was finalized in the meantime, then the vector has been already destroyed */
		if(PETSC_INITIALIZED){
			/* if the vector wasn't destroyed yet and the petsc is still running, then
			 * destroy the vector */
			TRY( VecDestroy(&inner_vector) );
		}
	}

}


void PetscVector::valuesUpdate() const{
	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << "(PetscVector)FUNCTION: valuesUpdate()" << std::endl;

	TRY( VecAssemblyBegin(inner_vector) );
	TRY( VecAssemblyEnd(inner_vector) );
}


void PetscVector::set(double new_value){
	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << "(PetscVector)FUNCTION: set(double)" << std::endl;

	TRY( VecSet(this->inner_vector,new_value) );

	valuesUpdate();
}


void PetscVector::set(int index, double new_value){
	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << "(PetscVector)FUNCTION: set(int,double)" << std::endl;

	TRY( VecSetValue(this->inner_vector,index,new_value, INSERT_VALUES) );
	
	valuesUpdate();
}

void PetscVector::load_local(std::string filename){
	if(!this->inner_vector){
		TRY( VecCreate(PETSC_COMM_SELF,&inner_vector) );
	}

	//TODO: check if file exists

	/* prepare viewer to load from file */
	PetscViewer mviewer;
	TRY( PetscViewerCreate(PETSC_COMM_SELF, &mviewer) );
	TRY( PetscViewerBinaryOpen(PETSC_COMM_SELF ,filename.c_str(), FILE_MODE_READ, &mviewer) );
	
	/* load vector from viewer */
	TRY( VecLoad(this->inner_vector, mviewer) );

	/* destroy the viewer */
	TRY( PetscViewerDestroy(&mviewer) );

	valuesUpdate();
}

void PetscVector::load_global(std::string filename){
	if(!this->inner_vector){
		TRY( VecCreate(PETSC_COMM_WORLD,&inner_vector) );
	}

	//TODO: check if file exists

	/* prepare viewer to load from file */
	PetscViewer mviewer;
	TRY( PetscViewerCreate(PETSC_COMM_WORLD, &mviewer) );
	TRY( PetscViewerBinaryOpen(PETSC_COMM_WORLD ,filename.c_str(), FILE_MODE_READ, &mviewer) );
	
	/* load vector from viewer */
	TRY( VecLoad(this->inner_vector, mviewer) );

	/* destroy the viewer */
	TRY( PetscViewerDestroy(&mviewer) );

	valuesUpdate();
}

void PetscVector::save_binary(std::string filename){
	//TODO: check if vector exists

	/* prepare viewer to save to file */
	PetscViewer mviewer;
	TRY( PetscViewerCreate(PETSC_COMM_WORLD, &mviewer) );
	TRY( PetscViewerBinaryOpen(PETSC_COMM_WORLD ,filename.c_str(), FILE_MODE_WRITE, &mviewer) );
	
	/* load vector from viewer */
	TRY( VecView(this->inner_vector, mviewer) );

	/* destroy the viewer */
	TRY( PetscViewerDestroy(&mviewer) );

	valuesUpdate();
}

void PetscVector::save_ascii(std::string filename){
	//TODO: check if vector exists

	/* prepare viewer to save to file */
	PetscViewer mviewer;
	TRY( PetscViewerCreate(PETSC_COMM_WORLD, &mviewer) );
	TRY( PetscViewerASCIIOpen(PETSC_COMM_WORLD ,filename.c_str(), &mviewer) );
	
	/* load vector from viewer */
	TRY( VecView(this->inner_vector, mviewer) );

	/* destroy the viewer */
	TRY( PetscViewerDestroy(&mviewer) );

	valuesUpdate();
}

Vec PetscVector::get_vector() const { // TODO: temp
	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << "(PetscVector)FUNCTION: get_vector()" << std::endl;
		
	return inner_vector;
}


int PetscVector::size() const{
	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << "(PetscVector)FUNCTION: size()" << std::endl;

	int global_size;

	TRY( VecGetSize(inner_vector,&global_size) );

	return global_size;
}


int PetscVector::local_size() const{
	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << "(PetscVector)FUNCTION: local_size()" << std::endl;

	int local_size;

	TRY( VecGetLocalSize(inner_vector,&local_size) );

	return local_size;
}


double PetscVector::get(int i)
{
	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << "(PetscVector)FUNCTION: get(int)" << std::endl;

	PetscInt ni = 1;
	PetscInt ix[1];
	PetscScalar y[1];
			
	ix[0] = i;

	TRY( VecGetValues(inner_vector,ni,ix,y) );
			
	return y[0];
}


void PetscVector::get_array(double **arr){
	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << "(PetscVector)FUNCTION: get_array(double **)" << std::endl;

	TRY( VecGetArray(inner_vector,arr) );
}


void PetscVector::restore_array(double **arr){
	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << "(PetscVector)FUNCTION: restore_array(double **)" << std::endl;

	TRY( VecRestoreArray(inner_vector,arr) );
}


void PetscVector::get_ownership(int *low, int *high){
	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << "(PetscVector)FUNCTION: get_ownership(int*, int*)" << std::endl;

	//TODO: control inner_vector

	TRY( VecGetOwnershipRange(inner_vector, low, high) );
}


void PetscVector::scale(double alpha){
	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << "(PetscVector)FUNCTION: scale(double)" << std::endl;

	//TODO: control inner_vector

	TRY( VecScale(inner_vector, alpha) );
	valuesUpdate(); // TODO: has to be called?
}


std::ostream &operator<<(std::ostream &output, const PetscVector &vector)		
{
	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << "(PetscVector)OPERATOR: <<" << std::endl;

	PetscScalar *arr_vector;
	PetscInt i,local_size;
	
	// TODO: make more sofisticated for parallel vectors

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
	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << "(PetscVector)OPERATOR: (vec = vec)" << std::endl;

	/* check for self-assignment by comparing the address of the implicit object and the parameter */
	/* vec1 = vec1 */
    if (this == &vec2){
		if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << " - self assignment" << std::endl;		
        return *this;
	}

	/* vec1 is not initialized yet */
	if (!inner_vector){
		if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << " - creating new vector" << std::endl;		
		TRY( VecDuplicate(vec2.inner_vector,&(this->inner_vector)) );
		this->valuesUpdate(); // TODO: has to be called?
	}

	/* else copy the values of inner vectors */
	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << " - copy values" << std::endl;		
	
	TRY( VecCopy(vec2.inner_vector,inner_vector) );
	this->valuesUpdate(); // TODO: has to be called?
	
	return *this;	
}

/* vec1 = linear_combination, perform full linear combination 
 * assemble linear combination and perform maxpy
 * */
PetscVector &PetscVector::operator=(PetscVectorWrapperComb comb){
	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << "(PetscVector)OPERATOR: (vec = comb)" << std::endl;

	/* vec1 is not initialized yet */
	if (!inner_vector){
		if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << " - duplicate vector" << std::endl;		
		TRY( VecDuplicate(comb.get_first_vector(),&inner_vector) );
	}

	/* vec = comb */
	comb.compute(inner_vector,0.0);

	return *this;	
}

/* vec1 = scalar_value <=> vec1(all) = scalar_value, assignment operator */
PetscVector &PetscVector::operator=(double scalar_value){
	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << "(PetscVector)OPERATOR: (vec = double)" << std::endl;

	this->set(scalar_value);
	return *this;	
}

/* vec1 = mul(v1,v2), call vector-vector multiplication on mul 
 * */
PetscVector &PetscVector::operator=(PetscVectorWrapperMul mulinstance){
	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << "(PetscVector)OPERATOR: (vec = mul)" << std::endl;

	/* vec1 is not initialized yet */
	if (!inner_vector){
		if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << " - duplicate vector" << std::endl;		
		TRY( VecDuplicate(mulinstance.get_vector1(),&inner_vector) );
	}

	/* vec = mul */
	mulinstance.mul(inner_vector);

	return *this;	
}

/* return subvector to be able to overload vector(index) = new_value */ 
PetscVectorWrapperSub PetscVector::operator()(int index) const
{   
	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << "(PetscVector)OPERATOR: (int) returns WrapperSub" << std::endl;
	
	/* create new indexset */
	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << " - create index set" << std::endl;
	IS new_subvector_is;
	PetscInt *idxs;
	TRY(PetscMalloc(sizeof(PetscInt),&idxs));
	
	idxs[0] = index;
	
	TRY( ISCreateGeneral(PETSC_COMM_WORLD, 1, idxs, PETSC_OWN_POINTER, &new_subvector_is));
	// TODO: when to delete this index set?
	
	return PetscVectorWrapperSub(this->inner_vector, new_subvector_is, true);
}

/* return subvector vector(index_begin:index_end), i.e. components with indexes: [index_begin, index_begin+1, ..., index_end] */ 
PetscVectorWrapperSub PetscVector::operator()(int index_begin, int index_end) const
{   
	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << "(PetscVector)OPERATOR: (int,int) returns WrapperSub" << std::endl;

	/* create new indexset */
	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << " - create index set" << std::endl;
	IS new_subvector_is;
	TRY( ISCreateStride(PETSC_COMM_WORLD, index_end-index_begin+1, index_begin,1, &new_subvector_is) );
	// TODO: when to delete this index set?
		
	return PetscVectorWrapperSub(this->inner_vector, new_subvector_is, true);
}

/* return subvector based on provided index set */ 
PetscVectorWrapperSub PetscVector::operator()(const IS new_subvector_is) const // TODO: are ju sure that this IS will be nt destroyed?
{   
	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << "(PetscVector)OPERATOR: vec(IS)" << std::endl;
	
	return PetscVectorWrapperSub(inner_vector,new_subvector_is, false);
}

/* define PetscVector(all) */
PetscVectorWrapperSub PetscVector::operator()(petscvector_all_type all_type) const{
	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << "(PetscVector)OPERATOR: vec(all)" << std::endl;

	IS new_subvector_is; // TODO: this is quite stupid, what about returning *this?
	ISCreateStride(PETSC_COMM_WORLD, this->size(), 0,1, &new_subvector_is);
	
	return PetscVectorWrapperSub(this->inner_vector, new_subvector_is, true);
} 



/* vec1 *= alpha */
void operator*=(PetscVector &vec1, double alpha)
{
	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << "(PetscVector)OPERATOR: vec *= double" << std::endl;
	
	vec1.scale(alpha);
}

/* vec1 += comb */
void operator+=(const PetscVector &vec1, PetscVectorWrapperComb comb)
{
	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << "(PetscVector)OPERATOR: vec += comb" << std::endl;
	
	/* vec1.inner_vector should be allocated */
	comb.compute(vec1.inner_vector,1.0);
}

/* vec1 -= comb */
void operator-=(PetscVector &vec1, PetscVectorWrapperComb comb)
{
	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << "(PetscVector)OPERATOR: vec -= comb" << std::endl;
	
	vec1 += (-1.0)*comb;
}

/* dot = dot(vec1,vec2) */
double dot(const PetscVector &vec1, const PetscVector &vec2)
{
	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << "(PetscVector)FUNCTION: dot(vec1,vec2)" << std::endl;

	double dot_value;
	TRY( VecDot(vec1.inner_vector,vec2.inner_vector,&dot_value));
	return dot_value;
}

/* norm = norm_2(vec1) */
double norm(const PetscVector &vec1)
{
	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << "(PetscVector)FUNCTION: norm(vec1)" << std::endl;

	double norm_value;
	TRY( VecNorm(vec1.inner_vector,NORM_2, &norm_value));
	return norm_value;
}

/* max = max(vec1) */
double max(const PetscVector &vec1)
{
	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << "(PetscVector)FUNCTION: max(vec)" << std::endl;

	double max_value;
	TRY( VecMax(vec1.inner_vector,NULL, &max_value) );
	return max_value;
}

/* sum = sum(vec1) */
double sum(const PetscVector &vec1)
{
	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << "(PetscVector)FUNCTION: sum(vec)" << std::endl;

	double sum_value;
	TRY( VecSum(vec1.inner_vector,&sum_value) );
	return sum_value;
}

/* vec3 = vec1./vec2 */
const PetscVector operator/(const PetscVector &vec1, const PetscVector &vec2)
{
	if(DEBUG_MODE_PETSCVECTOR >= 100) std::cout << "(PetscVector)FUNCTION: vec1/vec2" << std::endl;

	TRY(VecPointwiseDivide(vec1.inner_vector,vec1.inner_vector,vec2.inner_vector) );

	vec1.valuesUpdate(); // TODO: has to be called?
	return vec1;
}

PetscVectorWrapperMul mul(const PetscVector &vec1, const PetscVector &vec2)
{
	return PetscVectorWrapperMul( vec1.inner_vector, vec2.inner_vector);
}



} /* end of petscvector namespace */

#endif
