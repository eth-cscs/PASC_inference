#include "problem.h"

/** 
 * init the problem 
 *
*/
void Problem::init()
{
	this->timer_total.restart();


	/* here starts the fun */
	this->timer_total.start();
	
	this->it = 0;
}

/**
 * finalize the problem
 * 
*/ 
void Problem::finalize()
{
	this->model.finalize();
	this->data.finalize();

	this->timer_total.stop();
}

/** 
 * fit model to data
 *
*/
void Problem::solve(int max_s_steps, Scalar deltaL_eps)
{
	/* variables */
	double L, L_old, deltaL; /* object function value */

	/* initialize value of object function */
	L = std::numeric_limits<double>::max(); // TODO: the computation of L should be done in the different way
	
	/* main cycle */
	for(this->it=0;this->it < max_s_steps;this->it++){
		if(DEBUG_MODE >= 2) Message_info_value(" - it = ",this->it);

		/* --- COMPUTE Theta --- */
//		this->model.compute_theta(this->data.get_data_vec());
		
		/* --- COMPUTE gamma --- */
//		this->model.compute_gamma(this->data.get_data_vec());

		/* compute stopping criteria */
		L_old = L;
		L = model.get_function_value();
		deltaL = abs(L - L_old);

		/* print info about cost function */
		if(DEBUG_MODE >= 2){
			Message_info_value("  - L_old       = ",L_old);
			Message_info_value("  - L           = ",L);
			Message_info_value("  - |L - L_old| = ",deltaL);
		}	

		/* end the main cycle if the change of function value is sufficient */
		if (deltaL < deltaL_eps){
			break;
		}
		
	}

}


void Problem::set_data(Data &new_data)
{
	this->data = new_data;
}

void Problem::set_model(Model &new_model)
{
	this->model = new_model;
}

void Problem::print()
{
	Message_info(      "-- PROBLEM --");
	Message_info_time( " - time total: ",this->timer_total.get_value_sum());
//	Message_info_time( "  - time generate data: ",timer_data.get_value_sum());
//	Message_info_time( "  - time gamma:         ",timer_gamma.get_value_sum());
//	Message_info_time( "  - time theta:         ",timer_theta.get_value_sum());
//	Message_info_time( "  - time saveVTK:       ",timer_saveVTK.get_value_sum());
//	Message_info_time( "  - time other:         ",timer_all.get_value_sum() - (timer_data.get_value_sum() + timer_gamma.get_value_sum() + timer_theta.get_value_sum() +  timer_saveVTK.get_value_sum()));

	Message_info_value(" - number of outer iterations: ",this->it+1);
//	Message_info_value(" - |L - L_old| = ",deltaL);

}

void Problem::print_timers()
{
	Message(      "\n- final time info:");
	Message_info_time( " - time total: ",this->timer_total.get_value_sum());
	Message_info_value(" - number of outer iterations: ",this->it+1);

	this->model.print_timers();

	Message_info("");
}

void Problem::saveVTK(std::string name_of_file){
	InputOutput::saveVTK(name_of_file, this->data.get_data_vec(), this->model.get_gamma().get_gamma_vec(), this->model.get_dim(), this->model.get_T(), this->model.get_K());
}


Data Problem::get_data()
{
	return this->data;
}

Model Problem::get_model()
{
	return this->model;
}



