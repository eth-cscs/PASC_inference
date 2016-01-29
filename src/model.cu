#include "model.h"

/** 
 * init the model
 *
*/
void Model::init(int dim, int T, int K)
{
	this->dim = dim;
	this->T = T;
	this->K = K;
}

/**
 * finalize the model
 * 
*/ 
void Model::finalize()
{

}

Gamma Model::get_gamma()
{
	return this->gamma;
}

Theta Model::get_theta()
{
	return this->theta;
}


