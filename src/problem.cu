#include "problem.h"

/** 
 * init the problem 
 *
*/
void Problem::init()
{
	this->timer_total.restart();
	this->timer_total.start();
}

/**
 * finalize the problem
 * 
*/ 
void Problem::finalize()
{
	this->timer_total.stop();
}

void Problem::set_data(Data new_data)
{
	this->data = new_data;
}

void Problem::set_model(Model new_model)
{
	this->model = new_model;
}


Data Problem::get_data()
{
	return this->data;
}

Model Problem::get_model()
{
	return this->model;
}



