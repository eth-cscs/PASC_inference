#ifndef PETSCVECTOR_H
#define	PETSCVECTOR_H

#include <iostream>
#include <string>
#include <list>

extern int DEBUG_MODE;

namespace minlin {

namespace threx { // TODO: maybe choose the different namespace for my own Petsc stuff
 
class PetscVectorWrapperComb; /* wrapper to allow manipulation with linear combinations of vectors */
class PetscVectorWrapperCombNode; /* one node of previous wrapper */

/* PETSc Vector */
class PetscVector {
		PetscErrorCode ierr; // TODO: I don't know what to do with errors
		Vec inner_vector; /* original Petsc Vector */
		
		/* subvector stuff */
		IS subvector_is; /* is NULL if we work with whole vector */
		Vec inner_vector_orig; /* from which vector we created a subvector? to be able to restore */

	public:

		PetscVector();
		PetscVector(int n);
		PetscVector(const PetscVector &vec1);
		PetscVector(Vec new_inner_vector);
		PetscVector(Vec old_inner_vector, IS new_subvector_is);
		PetscVector(double scalar_value);

		
		~PetscVector();

		void valuesUpdate();
		void scale(PetscScalar alpha);

		Vec get_vector() const; // TODO: temp, direct access to inner vector should be forbidden
		int size();
		double get(int index);
		void get_ownership(int *low, int *high);

		void get_array(double **arr);
		void restore_array(double **arr);

		void set(double new_value);
		void set(int index, double new_value);

		/* assignment */
		PetscVector operator=(const PetscVector vec2);
		PetscVector operator=(PetscVectorWrapperCombNode combnode);	
		PetscVector operator=(PetscVectorWrapperComb comb);	
		PetscVector operator=(double scalar_value);	
		
		/* subvector */
		PetscVector operator()(int index);
		PetscVector operator()(int index_begin,int index_end);
		PetscVector operator()(const IS new_subvector_is);

		friend std::ostream &operator<<(std::ostream &output, const PetscVector &vector);

		friend void operator*=(PetscVector &vec1, double alpha);
		friend void operator+=(PetscVector &vec1, const PetscVectorWrapperComb comb);
		friend void operator-=(PetscVector &vec1, const PetscVectorWrapperComb comb);

		friend bool operator==(PetscVector vec1, double alpha);
		friend bool operator==(PetscVector vec1, PetscVector vec2);

		friend bool operator>(PetscVector vec1, PetscVector vec2);


		friend double dot(const PetscVector vec1, const PetscVector vec2); // TODO: with PetscVectorWrapperComb comb
		friend double max(const PetscVector vec1); // TODO: with PetscVectorWrapperComb comb
		friend double sum(const PetscVector vec1); // TODO: with PetscVectorWrapperComb comb
		friend const PetscVector operator/(PetscVector vec1, PetscVector vec2);
		

		/* define operator PetscVector(all), it is a mixture of petscvector and minlin  */
		PetscVector operator()(minlin::detail::all_type add_in) {
			if(DEBUG_MODE >= 100){
				std::cout << "operator vec(all)" << std::endl;
			}

			/* create new indexset */
			IS new_subvector_is; // TODO: this is quite stupid
			ISCreateStride(PETSC_COMM_WORLD, this->size(), 0,1, &new_subvector_is);
	
			return PetscVector(inner_vector, new_subvector_is);
		}

};

/* wrapper to allow linear combinations of vectors */
class PetscVectorWrapperComb
{
	private:
		std::list<PetscVectorWrapperCombNode> comb_list; /* the list of linear combination nodes */
		
	public:
		PetscVectorWrapperComb();
		PetscVectorWrapperComb(PetscVectorWrapperCombNode comb_node);
		PetscVectorWrapperComb(PetscVector vec);
		PetscVectorWrapperComb(double scalar_value);

		int get_listsize();
		int get_vectorsize();
		Vec get_first_vector();
		void append(PetscVectorWrapperCombNode new_node);
		void merge(PetscVectorWrapperComb comb);
		void get_arrays(PetscScalar *coeffs, Vec *vectors);
		
		friend std::ostream &operator<<(std::ostream &output, PetscVectorWrapperComb wrapper);

		friend const PetscVectorWrapperComb operator*(double alpha, PetscVectorWrapperComb comb);
		friend const PetscVectorWrapperComb operator+(PetscVectorWrapperComb comb1, PetscVectorWrapperComb comb2);
		friend const PetscVectorWrapperComb operator-(PetscVectorWrapperComb comb1, PetscVectorWrapperComb comb2);
		
};

/* wrapper of one node in the linear combinations */
class PetscVectorWrapperCombNode
{
	private:
		Vec inner_vector; /* pointer to vector in linear combination */
		double coeff; /* coefficient in linear combination */
	
	public:
		PetscVectorWrapperCombNode(); // TODO: temp
		PetscVectorWrapperCombNode(double new_coeff, Vec new_vector);
		PetscVectorWrapperCombNode(double scalar_value);
	
		void set_vector(Vec new_vector);
		Vec get_vector();
		int get_size();
		int get_value(int index);
		
		void set_coeff(double new_coeff);
		void scale(double alpha);
		double get_coeff();

		friend std::ostream &operator<<(std::ostream &output, PetscVectorWrapperCombNode wrapper);

};










/* --------------------- PetscVector ----------------------*/

/* PetscVector default constructor */
PetscVector::PetscVector(){
	if(DEBUG_MODE >= 100) std::cout << "empty constructor" << std::endl;

	inner_vector = PETSC_NULL;
	inner_vector_orig = PETSC_NULL; 
	subvector_is = PETSC_NULL;
}


/* PetscVector constructor with global dimension */
PetscVector::PetscVector(int n){
	if(DEBUG_MODE >= 100) std::cout << "constructor PetscVector(int)" << std::endl;

	inner_vector_orig = PETSC_NULL; 
	subvector_is = PETSC_NULL;

	VecCreate(PETSC_COMM_WORLD,&inner_vector);
	VecSetSizes(inner_vector,PETSC_DECIDE,n); // TODO: there should be more options to set the distribution
	VecSetFromOptions(inner_vector);

	valuesUpdate();
}

/* PetscVector copy constructor */
PetscVector::PetscVector(const PetscVector &vec1){
	if(DEBUG_MODE >= 100) std::cout << "constructor PetscVector(&vec)" << std::endl;

	/* inner_vec */
	VecDuplicate(vec1.inner_vector,&inner_vector);
	VecCopy(vec1.inner_vector,inner_vector);

	/* orig_vec */
	if(vec1.subvector_is){
		if(DEBUG_MODE >= 100) std::cout << " - copy also subvector_is" << std::endl;
		
		ISDuplicate(vec1.subvector_is,&subvector_is);
		ISCopy(vec1.subvector_is,subvector_is);

		VecDuplicate(vec1.inner_vector_orig,&inner_vector_orig);
		VecCopy(vec1.inner_vector_orig,inner_vector_orig);

	} else {
		inner_vector_orig = PETSC_NULL; 
		subvector_is = PETSC_NULL;
	}
	
/*	inner_vector = vec1.inner_vector;
	is_subvector = vec1.is_subvector;
	subvector_is = vec1.subvector_is;
	inner_vector_orig = vec1.inner_vector_orig;
*/

//	VecCopy(vec1.inner_vector,inner_vector);
}

/* PetscVector constructor with inner_vector */
PetscVector::PetscVector(Vec new_inner_vector){
	if(DEBUG_MODE >= 100) std::cout << "constructor PetscVector(inner_vector)" << std::endl;

	inner_vector = new_inner_vector;

	subvector_is = PETSC_NULL;
	inner_vector_orig = PETSC_NULL;	
}

/* PetscVector constructor with given IS = create subvector */
PetscVector::PetscVector(Vec old_inner_vector, IS new_subvector_is){
	if(DEBUG_MODE >= 100) std::cout << "constructor PetscVector(inner_vec, IS)" << std::endl;

	/* store old inner vector - will be used in the destructor to return subvector */
	inner_vector_orig = old_inner_vector; 

	/* copy IS */
	ISDuplicate(new_subvector_is, &subvector_is);
	ISCopy(new_subvector_is, subvector_is);

	if(DEBUG_MODE >= 100) std::cout << " - get subvector from original vector" << std::endl;

	/* get subvector, restore it in destructor */
	VecGetSubVector(inner_vector_orig, subvector_is, &inner_vector);
	
}



/* PetscVector destructor */
PetscVector::~PetscVector(){

	/* if this was a subvector, then restore values */
	if(subvector_is){

		/* restore subvector */
		if(DEBUG_MODE >= 100) std::cout << "restore subvector" << std::endl;

		VecRestoreSubVector(inner_vector_orig, subvector_is, &inner_vector);
//		inner_vector_orig = PETSC_NULL; 

		if(DEBUG_MODE >= 100) std::cout << "destroy index set" << std::endl;

		ISDestroy(&subvector_is);
		subvector_is = PETSC_NULL;
		inner_vector = PETSC_NULL;
		
	} else {

		/* if there was any inner vector, then destroy it */
		if(inner_vector){
			if(DEBUG_MODE >= 100) std::cout << "destroy inner vector" << std::endl;
		
			VecDestroy(&inner_vector);
		}
	}

}

/* after update a variable, it is necessary to call asseble begin & end */
void PetscVector::valuesUpdate(){
	VecAssemblyBegin(inner_vector);
	VecAssemblyEnd(inner_vector);
}

/* set all values of the vector, this function is called from overloaded operator */
void PetscVector::set(double new_value){
	if(DEBUG_MODE >= 100) std::cout << "void set(double)" << std::endl;

	VecSet(this->inner_vector,new_value);
	valuesUpdate();
}

/* set one specific value of the vector, this function is called from overloaded operator */
void PetscVector::set(int index, double new_value){
	if(DEBUG_MODE >= 100) std::cout << "void set(int,double)" << std::endl;

	VecSetValue(this->inner_vector,index,new_value, INSERT_VALUES);
	valuesUpdate();
}

/* returns inner vector */
Vec PetscVector::get_vector() const { // TODO: temp
	return inner_vector;
}

/* get size of the vector */
int PetscVector::size(){
	int global_size;
	VecGetSize(inner_vector,&global_size);
	return global_size;
}

/* get single value with given id of the vector (works only with local id), really slow */
double PetscVector::get(int i)
{
	PetscInt ni = 1;
	PetscInt ix[1];
	PetscScalar y[1];
			
	ix[0] = i;
	VecGetValues(inner_vector,ni,ix,y);			
			
	return y[0];
}

void PetscVector::get_array(double **arr){
	VecGetArray(inner_vector,arr);
}

void PetscVector::restore_array(double **arr){
	VecRestoreArray(inner_vector,arr);
}

void PetscVector::get_ownership(int *low, int *high){
	VecGetOwnershipRange(inner_vector, low, high);
}

/* inner_vector = alpha*inner_vector */
void PetscVector::scale(PetscScalar alpha){
	VecScale(inner_vector, alpha);
	valuesUpdate(); // TODO: has to be called?
}


/* stream insertion << operator */
std::ostream &operator<<(std::ostream &output, const PetscVector &vector)		
{
	PetscScalar *arr_vector;
	PetscInt i,local_size;
	
	output << "[";
	VecGetLocalSize(vector.inner_vector,&local_size);
	VecGetArray(vector.inner_vector,&arr_vector);
	for (i=0; i<local_size; i++){
		output << arr_vector[i];
		if(i < local_size-1) output << ",";
	}
	VecRestoreArray(vector.inner_vector,&arr_vector);
	output << "]";
			
	return output;
}

/* vec1 = vec2, assignment operator (set vector) */
PetscVector PetscVector::operator=(const PetscVector vec2){
	if(DEBUG_MODE >= 100) std::cout << "operator (vec = vec)" << std::endl;

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
PetscVector PetscVector::operator=(PetscVectorWrapperCombNode combnode){
	if(DEBUG_MODE >= 100) std::cout << "operator (vec = combnode)" << std::endl;

	/* vec1 = alpha*vec1 => simple scale */
    if (this->inner_vector == combnode.get_vector()){
		if(DEBUG_MODE >= 100) std::cout << " - scale values" << std::endl;		

        this->scale(combnode.get_coeff());
        return *this;
	}

	/* vec1 is not initialized yet */
	if (!inner_vector){
		if(DEBUG_MODE >= 100) std::cout << " - duplicate vector" << std::endl;		

		VecDuplicate(combnode.get_vector(),&inner_vector);
	}

	/* else copy the vector values and then scale */
	if(DEBUG_MODE >= 100) std::cout << " - copy values" << std::endl;		

	VecCopy(combnode.get_vector(),inner_vector);
    this->scale(combnode.get_coeff());

	return *this;	
}

/* vec1 = linear_combination, perform full linear combination */
PetscVector PetscVector::operator=(PetscVectorWrapperComb comb){
	if(DEBUG_MODE >= 100){
		std::cout << "operator (vec = comb)" << std::endl;
	}

	/* vec1 is not initialized yet */
	if (!inner_vector){
		if(DEBUG_MODE >= 100) std::cout << " - duplicate vector" << std::endl;		

		VecDuplicate(comb.get_first_vector(),&inner_vector);
	}

	/* vec1 = 0, we will perform MAXPY (y += lin_comb) */
	VecSet(inner_vector,0.0);

	/* vec += comb */
	*this += comb; 

	return *this;	
}

/* vec1 = scalar_value <=> vec1(all) = scalar_value, assignment operator */
PetscVector PetscVector::operator=(double scalar_value){
	if(DEBUG_MODE >= 100){
		std::cout << "operator (vec = scalar)" << std::endl;
	}

	this->set(scalar_value);
	return *this;	
}

/* return subvector to be able to overload vector(index) = new_value */ 
PetscVector PetscVector::operator()(int index)
{   
	if(DEBUG_MODE >= 100) std::cout << "operator vec(int)" << std::endl;
	
	/* create new indexset */
	IS new_subvector_is;
	PetscInt idxs[1];
	idxs[0] = index;
	
	ISCreateGeneral(PETSC_COMM_WORLD, 1, idxs, PETSC_COPY_VALUES, &new_subvector_is);
	
	return PetscVector(this->inner_vector, new_subvector_is);
}

/* return subvector vector(index_begin:index_end), i.e. components with indexes: [index_begin, index_begin+1, ..., index_end] */ 
PetscVector PetscVector::operator()(int index_begin, int index_end)
{   
	if(DEBUG_MODE >= 100){
		std::cout << "operator vec(int,int)" << std::endl;
	}

	/* create new indexset */
	IS new_subvector_is;
	ISCreateStride(PETSC_COMM_WORLD, index_end-index_begin+1, index_begin,1, &new_subvector_is);
	
	return PetscVector(inner_vector, new_subvector_is);
}

/* return subvector based on provided index set */ 
PetscVector PetscVector::operator()(const IS new_subvector_is)
{   
	if(DEBUG_MODE >= 100) std::cout << "operator vec(IS)" << std::endl;
	
	return PetscVector(inner_vector, new_subvector_is);
}



/* vec1 *= alpha */
void operator*=(PetscVector &vec1, double alpha)
{
	if(DEBUG_MODE >= 100) std::cout << "operator vec *= scalar" << std::endl;
	
	vec1.scale(alpha);
}

/* vec1 += comb */
void operator+=(PetscVector &vec1, PetscVectorWrapperComb comb)
{
	if(DEBUG_MODE >= 100){
		std::cout << "operator vec += comb" << std::endl;
	}
	
	int list_size = comb.get_listsize();
	PetscScalar alphas[list_size];
	Vec vectors[list_size];

	/* get array with coefficients and vectors */
	comb.get_arrays(alphas,vectors);


	std::cout << "maxpy_vector_control before:" << std::endl;
	VecView(vec1.inner_vector,PETSC_VIEWER_STDOUT_SELF);

	/* vec1 = vec1 + sum (coeff*vector) */
	VecMAXPY(vec1.inner_vector,list_size,alphas,vectors);
	vec1.valuesUpdate();

	std::cout << "maxpy_vector_control after:" << std::endl;
	VecView(vec1.inner_vector,PETSC_VIEWER_STDOUT_SELF);

	
}

/* vec1 -= comb */
void operator-=(PetscVector &vec1, PetscVectorWrapperComb comb)
{
	if(DEBUG_MODE >= 100){
		std::cout << "operator vec -= comb" << std::endl;
	}
	
	vec1 += (-1.0)*comb;
}

/* vec1 == scalar */
bool operator==(PetscVector vec1, double alpha){
	bool return_value = false; // TODO: works only with vector of size 1, otherwise compare only first value
	double vector_value = vec1.get(0);
	
	if(vector_value == alpha){
		return_value = true;
	}	
	
	return return_value;
}

/* vec1 == vec2 */
bool operator==(PetscVector vec1, PetscVector vec2){
	PetscBool return_value;

	VecEqual(vec1.inner_vector,vec2.inner_vector,&return_value);
	
	return (bool)return_value;
}

/* vec1 > vec2 */
bool operator>(PetscVector vec1, PetscVector vec2){
	bool return_value = false; // TODO: works only with vector of size 1, otherwise compare only first value
	
	double vec1_value = vec1.get(0);
	double vec2_value = vec2.get(0);
	
	if(vec1_value > vec2_value){
		return_value = true;
	}	
	
	return return_value;

}


/* dot = dot(vec1,vec2) */
double dot(const PetscVector vec1, const PetscVector vec2)
{
	double dot_value;
	VecDot(vec1.inner_vector,vec2.inner_vector,&dot_value);
	return dot_value;
}

/* max = max(vec1) */
double max(const PetscVector vec1)
{
	double max_value;
	VecMax(vec1.inner_vector,NULL, &max_value);
	return max_value;
}

/* sum = sum(vec1) */
double sum(const PetscVector vec1)
{
	double sum_value;
	VecSum(vec1.inner_vector,&sum_value);
	return sum_value;
}

/* vec3 = vec1./vec2 */
const PetscVector operator/(PetscVector vec1, PetscVector vec2)
{
	VecPointwiseDivide(vec1.inner_vector,vec1.inner_vector,vec2.inner_vector);
	vec1.valuesUpdate();
	return vec1;
}


/* --------------------- PetscVectorWrapperComb ----------------------*/

/* default constructor */
PetscVectorWrapperComb::PetscVectorWrapperComb(){
}

/* constructor from node */
PetscVectorWrapperComb::PetscVectorWrapperComb(PetscVectorWrapperCombNode comb_node){
	/* append node */
	this->append(comb_node);
}

/* constructor from vec */
PetscVectorWrapperComb::PetscVectorWrapperComb(PetscVector vec){
	/* create node from vector */
	PetscVectorWrapperCombNode comb_node(1.0,vec.get_vector());

	/* append new node to newly created combination */
	this->append(comb_node);
}

/* constructor from scalar_value - create Node from value */
PetscVectorWrapperComb::PetscVectorWrapperComb(double scalar_value){
	/* create node from scalar_value = create vector of size 1 */
	PetscVectorWrapperCombNode comb_node(scalar_value);

	/* append it to newly created combination */
	this->append(comb_node);
}


/* append new node to the list */
void PetscVectorWrapperComb::append(PetscVectorWrapperCombNode new_node){
	comb_list.push_back(new_node);
}

/* append new list to the end of old list (merge without sort), will be called from overloaded operator+ */
void PetscVectorWrapperComb::merge(PetscVectorWrapperComb comb){
	comb_list.insert(comb_list.end(), comb.comb_list.begin(), comb.comb_list.end());
}

/* get length of the list */
int PetscVectorWrapperComb::get_listsize(){
	return comb_list.size();
}

/* get size of the vectors in the list */
int PetscVectorWrapperComb::get_vectorsize(){
	std::list<PetscVectorWrapperCombNode>::iterator list_iter; /* iterator through list */
	PetscInt vector_size;

	/* get first element and obtain a size of the vector */
	list_iter = comb_list.begin();
	vector_size = list_iter->get_size();

	return vector_size;
}

/* get frist vector from the list */
Vec PetscVectorWrapperComb::get_first_vector(){
	std::list<PetscVectorWrapperCombNode>::iterator list_iter; /* iterator through list */
	Vec vector;

	/* get first element and obtain a size of the vector */
	list_iter = comb_list.begin();
	vector = list_iter->get_vector();

	return vector;
}

/* prepare arrays from linear combination list, I assume that arrays are allocated */
void PetscVectorWrapperComb::get_arrays(PetscScalar *coeffs, Vec *vectors){
	std::list<PetscVectorWrapperCombNode>::iterator list_iter; /* iterator through list */
	int list_size = this->get_listsize();
	int j;

	/* set to the begin of the list */
	list_iter = comb_list.begin();

	/* go through the list and fill the vectors */
	for(j=0;j<list_size;j++){
		coeffs[j] = list_iter->get_coeff();
		vectors[j] = list_iter->get_vector();
		
		if(j < list_size-1){
			/* this is not the last element */
			list_iter++;
		}
	}

}


/* print linear combination */ //TODO: yes, this is really slow, but now I don't care.. who cares? it is only for testing purposes
std::ostream &operator<<(std::ostream &output, PetscVectorWrapperComb wrapper)
{
	PetscInt i,j,vector_size,list_size;
	std::list<PetscVectorWrapperCombNode>::iterator list_iter; /* iterator through list */
	
	output << "[";
	
	list_size = wrapper.get_listsize();
	vector_size = wrapper.get_vectorsize();
	
	/* go through components in the vector */
	for(i=0;i<vector_size;i++){
		list_iter = wrapper.comb_list.begin();
		
		/* for each component go throught the list */
		for(j=0;j<list_size;j++){
			/* print coeff, if coeff < 0, then print it in () */
			if(list_iter->get_coeff() < 0.0){
				output << "(" << list_iter->get_coeff() << ")";
			} else {
				output << list_iter->get_coeff(); 
			}
			output << "*" << list_iter->get_value(i);
			if(j < list_size-1){ 
				/* this is not the last node */
				output << "+";
				list_iter++;
			}
		}
		
		if(i < vector_size-1){ 
			/* this is not the last element */
			output << ", ";
		}
	}

	output << "]";

	return output;
}

/* all nodes in linear combination will be scaled */
const PetscVectorWrapperComb operator*(double alpha, PetscVectorWrapperComb comb){
	/* scale nodes in the list */
	std::list<PetscVectorWrapperCombNode>::iterator list_iter; /* iterator through list */
	int list_size = comb.get_listsize();
	int j;

	/* set to the begin of the list */
	list_iter = comb.comb_list.begin();

	/* go through the list and scale nodes */
	for(j=0;j<list_size;j++){
		list_iter->scale(alpha);
		
		if(j < list_size-1){
			/* this is not the last element */
			list_iter++;
		}
	}
	
	return comb;
}

/* new linear combination created by comb + comb */
const PetscVectorWrapperComb operator+(PetscVectorWrapperComb comb1, PetscVectorWrapperComb comb2){
	/* append second linear combination to the first */
	comb1.merge(comb2);
	
	return comb1;
}

/* new linear combination created by comb + comb */
const PetscVectorWrapperComb operator-(PetscVectorWrapperComb comb1, PetscVectorWrapperComb comb2){
	/* append second linear combination to the first */
	comb1.merge(-1.0*comb2);
	
	return comb1;
}




/* --------------------- PetscVectorWrapperCombNode ----------------------*/

/* default constructor */ //TODO: temp, it is not possible to create node without any vector and scalar
PetscVectorWrapperCombNode::PetscVectorWrapperCombNode(){
}

/* constructor from vector and coefficient */
PetscVectorWrapperCombNode::PetscVectorWrapperCombNode(double new_coeff, Vec new_vector){
	set_vector(new_vector);
	set_coeff(new_coeff);
}

/* create vector from scalar value - to be able to vec(idx) = vec2(idx) + scalar_value */
PetscVectorWrapperCombNode::PetscVectorWrapperCombNode(double scalar_value){
	/* create sequential petsc vector of size 1 */
	VecCreateSeq(PETSC_COMM_SELF,1,&(this->inner_vector));
	VecSet(inner_vector,scalar_value);

	VecAssemblyBegin(inner_vector); // TODO: has to be called?
	VecAssemblyEnd(inner_vector);

	set_coeff(1.0);
}

/* set vector to the node */
void PetscVectorWrapperCombNode::set_vector(Vec new_vector){
	this->inner_vector = new_vector;
}

/* return vector from this node */
Vec PetscVectorWrapperCombNode::get_vector(){
	return this->inner_vector;
}

/* set new coefficient to this node */
void PetscVectorWrapperCombNode::set_coeff(double new_coeff){
	this->coeff = new_coeff;
}

/* node is multiplied by scalar, now multiply only the coefficient of linear combination */
void PetscVectorWrapperCombNode::scale(double alpha){
	this->coeff *= alpha;
}

/* get the coefficient from this node */
double PetscVectorWrapperCombNode::get_coeff(){
	return this->coeff;
}

/* get size of the vector */
int PetscVectorWrapperCombNode::get_size(){
	int global_size;
	VecGetSize(this->inner_vector,&global_size);
	return global_size;
}

/* get value from the vector, really slow */
int PetscVectorWrapperCombNode::get_value(int index){
	PetscInt ni = 1;
	PetscInt ix[1];
	PetscScalar y[1];
			
	ix[0] = index;
	VecGetValues(this->inner_vector,ni,ix,y);			
			
	return y[0];
}

/* stream insertion << operator */
std::ostream &operator<<(std::ostream &output, PetscVectorWrapperCombNode wrapper)
{
	PetscScalar *arr_vector;
	PetscInt i,local_size;
	
	output << "[";
	VecGetLocalSize(wrapper.inner_vector,&local_size);
	VecGetArray(wrapper.inner_vector,&arr_vector);
	for (i=0; i<local_size; i++){
		output << wrapper.coeff << "*" << arr_vector[i];
		if(i < local_size-1) output << ", ";
	}
	VecRestoreArray(wrapper.inner_vector,&arr_vector);
	output << "]";
			
	return output;
}


} /* end of namespace */

} /* end of MinLin namespace */

#endif
