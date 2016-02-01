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
	
	/* prepare timers */
	this->timer_gamma.restart();
	this->timer_theta.restart();
	
	/* initialize gamma */
	Gamma new_gamma;
	this->timer_gamma.start(); /* start timer for initializing gamma */
	 new_gamma.init(this->dim, this->T, this->K);
	 new_gamma.prepare_random();	/* prepare gammas */
	timer_gamma.stop();

	this->gamma = new_gamma;

	if(DEBUG_MODE >= 2) Message_info_time(" - gamma generated in: ",this->timer_gamma.get_value_last());
	if(DEBUG_MODE >= 3) this->gamma.print();

	/* initialize theta */
	Theta new_theta;
	this->timer_theta.start();
 	 new_theta.init(this->dim,this->K);
	this->timer_theta.stop();

	this->theta = new_theta;
	
	if(DEBUG_MODE >= 2) Message_info_time(" - theta prepared in: ",this->timer_theta.get_value_last());
	if(DEBUG_MODE >= 3) theta.print();
	
	
}

/**
 * finalize the model
 * 
*/ 
void Model::finalize()
{
	this->theta.finalize();
	this->gamma.finalize();

}

/**
 * print info about model
 * 
*/ 
void Model::print() 
{
	
	Message_info("-- MODEL ---");
	Message_info_value("- dim: ",this->dim);
	Message_info_value("- T:   ",this->T);
	Message_info_value("- K:   ",this->K);

}

int Model::get_dim(){
	return this->dim;
}

int Model::get_T(){
	return this->T;
}

int Model::get_K(){
	return this->K;
}

void Model::compute_theta(DataVector<Scalar> data_vec)
{
	this->theta.compute(data_vec,this->gamma);
}

void Model::compute_gamma(DataVector<Scalar> data_vec)
{
	this->gamma.compute(data_vec,this->theta);
}

Gamma Model::get_gamma()
{
	return this->gamma;
}

Theta Model::get_theta()
{
	return this->theta;
}


