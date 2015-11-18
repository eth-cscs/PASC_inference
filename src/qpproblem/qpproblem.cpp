#include "qpproblem.h"

/* constructor */
QPproblem::QPproblem(Data* data, Gamma *gamma, Theta *theta, PetscScalar eps_sqr){
	this->data = data;
	this->gamma = gamma;
	this->theta = theta;
	this->eps_sqr = eps_sqr;	
}

PetscErrorCode QPproblem::init(){
	PetscFunctionBegin;
    PetscFunctionReturn(0);  
}

PetscErrorCode QPproblem::finalize(){
	PetscFunctionBegin;
    PetscFunctionReturn(0);  
}

PetscErrorCode QPproblem::solve(){
	PetscFunctionBegin;
    PetscFunctionReturn(0);  
}

PetscErrorCode QPproblem::get_function_value(PetscScalar*){
	PetscFunctionBegin;
    PetscFunctionReturn(0);  
}

PetscErrorCode QPproblem::print(PetscViewer){
	PetscFunctionBegin;
    PetscFunctionReturn(0);  
}

PetscErrorCode QPproblem::correct(PetscScalar){
	PetscFunctionBegin;
    PetscFunctionReturn(0);  
}
