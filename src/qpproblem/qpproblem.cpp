#include "qpproblem.h"

/* constructor */
QPproblem::QPproblem(Gamma *gamma, Theta *theta, PetscScalar eps_sqr){
	this->gamma = gamma;
	this->theta = theta;
	this->eps_sqr = eps_sqr;	
}

