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
	 new_gamma.init(this->T, this->K);
	 new_gamma.prepare_random();	/* prepare gammas */
	timer_gamma.stop();
	Message_info_time(" - gamma generated in: ",this->timer_gamma.get_value_last());

	this->gamma = new_gamma;

//	if(DEBUG_PRINTDATA){ /* print gamma */
//		this->gamma.print();
//	}

	/* initialize theta */
	Theta new_theta;
	this->timer_theta.start();
 	 new_theta.init(this->dim,this->K);
	this->timer_theta.stop();
	Message_info_time(" - theta prepared in: ",this->timer_theta.get_value_last());

	this->theta = new_theta;
	
//	if(DEBUG_PRINTDATA){ /* print theta */
//		theta.print();
//	}
	
	
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


