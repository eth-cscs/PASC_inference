#include "qpsolver.h"

/* constructor */
QPSolver::QPSolver(Data* data, Gamma *gamma, Theta *theta, PetscScalar eps_sqr){
	this->data = data;
	this->gamma = gamma;
	this->theta = theta;
	this->eps_sqr = eps_sqr;	
}

PetscErrorCode QPSolver::init(){
	PetscFunctionBegin;
    PetscFunctionReturn(0);  
}

PetscErrorCode QPSolver::finalize(){
	PetscFunctionBegin;
    PetscFunctionReturn(0);  
}

PetscErrorCode QPSolver::solve(){
	PetscFunctionBegin;
    PetscFunctionReturn(0);  
}

PetscErrorCode QPSolver::get_function_value(PetscScalar*){
	PetscFunctionBegin;
    PetscFunctionReturn(0);  
}

PetscErrorCode QPSolver::print(PetscViewer){
	PetscFunctionBegin;
    PetscFunctionReturn(0);  
}

PetscErrorCode QPSolver::correct(PetscScalar){
	PetscFunctionBegin;
    PetscFunctionReturn(0);  
}
