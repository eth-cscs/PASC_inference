#ifndef PETSCVECTOR_H
#define	PETSCVECTOR_H

#include <iostream>
#include <string>
#include <list>


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

		PetscVector(int n);
		PetscVector(const PetscVector &vec1);
		PetscVector(Vec new_inner_vector);
		PetscVector(Vec old_inner_vector, IS new_subvector_is);
		PetscVector(double scalar_value);

		
		~PetscVector();

		void valuesUpdate();
		void scale(PetscScalar alpha);

		Vec get_vector() const; // TODO: temp, direct access to inner vector should be forbidden
		int get_size();
		double get(int index);

		void set(PetscScalar new_value);
		void set(int index, PetscScalar new_value);

		PetscVector& operator=(const PetscVector &vec2);
		PetscVector& operator=(PetscVectorWrapperCombNode combnode);	
		PetscVector& operator=(PetscVectorWrapperComb comb);	
		PetscVector& operator=(double scalar_value);	
		
		PetscVector operator()(int index);
		PetscVector operator()(int index_begin,int index_end);

		friend std::ostream &operator<<(std::ostream &output, const PetscVector &vector);

		friend void operator*=(PetscVector vec1, double alpha);
		friend void operator+=(PetscVector vec1, PetscVectorWrapperComb comb);
		friend void operator-=(PetscVector vec1, PetscVectorWrapperComb comb);
		friend double dot(const PetscVector vec1, const PetscVector vec2); // TODO: with PetscVectorWrapperComb comb


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

/* PetscVector constructor with global dimension */
PetscVector::PetscVector(int n){
	inner_vector_orig = PETSC_NULL; 
	subvector_is = PETSC_NULL;

	VecCreate(PETSC_COMM_WORLD,&inner_vector);
	VecSetSizes(inner_vector,PETSC_DECIDE,n); // TODO: there should be more options to set the distribution
	VecSetFromOptions(inner_vector);

	valuesUpdate();
}

/* PetscVector copy constructor */
PetscVector::PetscVector(const PetscVector &vec1){
	inner_vector = vec1.inner_vector;
//	VecCopy(vec1.inner_vector,inner_vector);
}

/* PetscVector constructor with inner_vector */
PetscVector::PetscVector(Vec new_inner_vector){
	inner_vector = new_inner_vector;
}

/* PetscVector constructor with given IS = create subvector */
PetscVector::PetscVector(Vec old_inner_vector, IS new_subvector_is){
	/* store old inner vector - will be used in the destructor to return subvector */
	inner_vector_orig = old_inner_vector; 
	subvector_is = new_subvector_is;
	
	/* get subvector, restore it in destructor */
	VecGetSubVector(inner_vector_orig, subvector_is, &inner_vector);
	
}



/* PetscVector destructor */
PetscVector::~PetscVector(){

	/* if this was a subvector, then restore values */
	if(subvector_is != PETSC_NULL){
		/* restore subvector */
		VecRestoreSubVector(inner_vector_orig, subvector_is, &inner_vector);
	}
	
}

/* after update a variable, it is necessary to call asseble begin & end */
void PetscVector::valuesUpdate(){
	VecAssemblyBegin(inner_vector);
	VecAssemblyEnd(inner_vector);
}

/* set all values of the vector, this function is called from overloaded operator */
void PetscVector::set(PetscScalar new_value){
	VecSet(inner_vector,new_value);
	valuesUpdate();
}

/* set one specific value of the vector, this function is called from overloaded operator */
void PetscVector::set(int index, PetscScalar new_value){
	VecSetValue(inner_vector,index,new_value, INSERT_VALUES);
	valuesUpdate();
}

/* returns inner vector */
Vec PetscVector::get_vector() const { // TODO: temp
	return inner_vector;
}

/* get size of the vector */
int PetscVector::get_size(){
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
PetscVector& PetscVector::operator=(const PetscVector &vec2){
	/* check for self-assignment by comparing the address of the implicit object and the parameter */
	/* vec1 = vec1 */
    if (this == &vec2)
        return *this;

	/* else copy the values of inner vectors */
	VecCopy(vec2.inner_vector,inner_vector);
	return *this;	
}

/* vec1 = linear_combination_node, perform simple linear combination */
PetscVector& PetscVector::operator=(PetscVectorWrapperCombNode combnode){
	/* vec1 = alpha*vec1 => simple scale */
    if (this->inner_vector == combnode.get_vector()){
        this->scale(combnode.get_coeff());
        return *this;
	}

	/* else copy the vector values and then scale */
	VecCopy(combnode.get_vector(),inner_vector);
    this->scale(combnode.get_coeff());

	return *this;	
}

/* vec1 = linear_combination, perform full linear combination */
PetscVector& PetscVector::operator=(PetscVectorWrapperComb comb){

	/* vec1 = 0, we will perform MAXPY (y += lin_comb) */
	VecSet(inner_vector,0.0);

	/* vec += comb */
	*this += comb; 

	return *this;	
}

/* vec1 = scalar_value <=> vec1(all) = scalar_value, assignment operator */
PetscVector& PetscVector::operator=(double scalar_value){
	this->set(scalar_value);
	return *this;	
}

/* return subvector to be able to overload vector(index) = new_value */ 
PetscVector PetscVector::operator()(int index)
{   
	/* create new indexset */
	IS new_subvector_is;
	PetscInt idxs[1];
	idxs[0] = index;
	
	ISCreateGeneral(PETSC_COMM_WORLD, 1, idxs, PETSC_COPY_VALUES, &new_subvector_is);
	
	return PetscVector(inner_vector, new_subvector_is);
}

/* return subvector vector(index_begin:index_end), i.e. components with indexes: [index_begin, index_begin+1, ..., index_end] */ 
PetscVector PetscVector::operator()(int index_begin, int index_end)
{   
	/* create new indexset */
	IS new_subvector_is;
	ISCreateStride(PETSC_COMM_WORLD, index_end-index_begin+1, index_begin,1, &new_subvector_is);
	
	return PetscVector(inner_vector, new_subvector_is);
}

/* vec1 *= alpha */
void operator*=(PetscVector vec1, double alpha)
{
	vec1.scale(alpha);
}

/* vec1 += comb */
void operator+=(PetscVector vec1, PetscVectorWrapperComb comb)
{
	int list_size = comb.get_listsize();
	PetscScalar alphas[list_size];
	Vec vectors[list_size];

	/* get array with coefficients and vectors */
	comb.get_arrays(alphas,vectors);
	
	/* vec1 = vec1 + sum (coeff*vector) */
	VecMAXPY(vec1.inner_vector,list_size,alphas,vectors);
}

/* vec1 -= comb */
void operator-=(PetscVector vec1, PetscVectorWrapperComb comb)
{
	vec1 += (-1.0)*comb;
}

/* dot = dot(vec1,vec2) */
double dot(const PetscVector vec1, const PetscVector vec2)
{
	double dot_value;
	VecDot(vec1.inner_vector,vec2.inner_vector,&dot_value);
	return dot_value;
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
