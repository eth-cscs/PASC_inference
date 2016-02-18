
namespace minlin {

namespace threx { // TODO: maybe choose the different namespace for my own Petsc stuff




/* --------------------- PetscVectorWrapperComb ----------------------*/

/* default constructor */
PetscVectorWrapperComb::PetscVectorWrapperComb(){
	if(DEBUG_MODE >= 100) std::cout << "(WrapperComb)CONSTRUCTOR: empty" << std::endl;
}

/* constructor from node */
PetscVectorWrapperComb::PetscVectorWrapperComb(PetscVectorWrapperCombNode &comb_node){
	if(DEBUG_MODE >= 100) std::cout << "(WrapperComb)CONSTRUCTOR: comb_node" << std::endl;

	/* append node */
	this->append(comb_node);
}

/* constructor from vec */
PetscVectorWrapperComb::PetscVectorWrapperComb(PetscVector &vec){
	if(DEBUG_MODE >= 100) std::cout << "(WrapperComb)CONSTRUCTOR: vec" << std::endl;

	/* create node from vector */
	PetscVectorWrapperCombNode comb_node(1.0,vec.get_vector());

	/* append new node to newly created combination */
	this->append(comb_node);
}

/* constructor from scalar_value - create Node from value */
PetscVectorWrapperComb::PetscVectorWrapperComb(double scalar_value){
	if(DEBUG_MODE >= 100) std::cout << "(WrapperComb)CONSTRUCTOR: double" << std::endl;

	/* create node from scalar_value = create vector of size 1 */
	PetscVectorWrapperCombNode comb_node(scalar_value);

	/* append it to newly created combination */
	this->append(comb_node);
}

/* constructor from subvector */
PetscVectorWrapperComb::PetscVectorWrapperComb(PetscVectorWrapperSub subvector){
	if(DEBUG_MODE >= 100) std::cout << "(WrapperComb)CONSTRUCTOR: from subvector" << std::endl;

	/* create node from vector */
	PetscVectorWrapperCombNode comb_node(1.0,subvector.get_subvector());

	/* append new node to newly created combination */
	this->append(comb_node);

}

/* destructor */
PetscVectorWrapperComb::~PetscVectorWrapperComb(){
	if(DEBUG_MODE >= 100) std::cout << "(WrapperComb)DESTRUCTOR" << std::endl;

	// TODO: destroy list

}

/* append new node to the list */
void PetscVectorWrapperComb::append(PetscVectorWrapperCombNode &new_node){
	if(DEBUG_MODE >= 100) std::cout << "(WrapperComb)FUNCTION: append(WrapperCombNode)" << std::endl;
	
	comb_list.push_back(new_node);
}

/* append new list to the end of old list (merge without sort), will be called from overloaded operator+ */
void PetscVectorWrapperComb::merge(PetscVectorWrapperComb &comb){
	if(DEBUG_MODE >= 100) std::cout << "(WrapperComb)FUNCTION: merge(WrapperComb)" << std::endl;

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



/* print linear combination without instance, f.x << alpha*vec1 + beta*vec2 */
std::ostream &operator<<(std::ostream &output, PetscVectorWrapperComb &wrapper)
{
	if(DEBUG_MODE >= 100) std::cout << "(WrapperComb)OPERATOR: << wrapper" << std::endl;
		
	PetscInt i,j,vector_size,list_size;
	double value;
	
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

			output << "*";

			/* print value, if value < 0, then print it in () */
			value = list_iter->get_value(i);
			if(value < 0.0){
				output << "(" << value << ")";
			} else {
				output << value; 
			}

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

/* new linear combination created by comb - comb */
const PetscVectorWrapperComb operator-(PetscVectorWrapperComb comb1, PetscVectorWrapperComb comb2){
	/* append second linear combination to the first */
	comb2 = -1*comb2;
	comb1.merge(comb2);
	
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
	TRY( VecGetSize(this->inner_vector,&global_size) );
	return global_size;
}

/* get value from the vector, really slow */
int PetscVectorWrapperCombNode::get_value(int index){
	PetscInt ni = 1;
	PetscInt ix[1];
	PetscScalar y[1];
			
	ix[0] = index;
	TRY( VecGetValues(this->inner_vector,ni,ix,y) );	
			
	return y[0];
}

/* stream insertion << operator */
std::ostream &operator<<(std::ostream &output, PetscVectorWrapperCombNode &wrapper)
{
	if(DEBUG_MODE >= 100) std::cout << "(WrapperCombNode)OPERATOR: << wrapper" << std::endl;
	
	PetscScalar *arr_vector;
	PetscInt i,local_size;
	
	output << "[";
	TRY( VecGetLocalSize(wrapper.inner_vector,&local_size) );
	TRY( VecGetArray(wrapper.inner_vector,&arr_vector) );
	for (i=0; i<local_size; i++){
		output << wrapper.coeff << "*" << arr_vector[i];
		if(i < local_size-1) output << ", ";
	}
	TRY( VecRestoreArray(wrapper.inner_vector,&arr_vector) );
	output << "]";
			
	return output;
}



}

}
