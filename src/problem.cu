#include "problem.h"

/** 
 * init the problem 
 *
*/
void Problem::init()
{
	this->timer_total.restart();
	this->timer_gamma.restart();
	this->timer_theta.restart();

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
	QPSolver qpsolver = this->model.get_gamma().get_qpsolver();

	Scalar L, L_old, deltaL; /* object function value */

	/* initialize value of object function */
	L = std::numeric_limits<Scalar>::max(); // TODO: the computation of L should be done in the different way
	
	/* main cycle */
	for(this->it=0;this->it < max_s_steps;this->it++){
		if(DEBUG_MODE >= 2) Message_info_value(" - it = ",this->it);

		/* --- COMPUTE Theta --- */
		this->timer_theta.start();
		 this->model.compute_theta(this->data.get_data_vec());
		this->timer_theta.stop();

		if(DEBUG_MODE >= 2) Message_info_time("  - theta problem solved in: ",this->timer_theta.get_value_last());
		if(DEBUG_MODE >= 3) this->model.get_theta().print(2);
		
		/* --- COMPUTE gamma --- */
		this->timer_gamma.start();
		 this->model.compute_gamma(this->data.get_data_vec());
		this->timer_gamma.stop(); 

		if(DEBUG_MODE >= 2) Message_info_time("  - gamma problem solved in:  ",this->timer_gamma.get_value_last());
		if(DEBUG_MODE >= 3) this->model.get_gamma().print(2);

		/* compute stopping criteria */
		L_old = L;
		L = qpsolver.get_function_value(this->model.get_gamma().gamma_vec);
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

	qpsolver.finalize();


}


void Problem::set_data(Data new_data)
{
	this->data = new_data;
}

void Problem::set_model(Model new_model)
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

	Message_info_value(" - number of outer iterations: ",this->it);
//	Message_info_value(" - |L - L_old| = ",deltaL);
/*	Message_info(" - QPSolver:");
	Message_info_value("  - it =        ", qpsolver.get_it_all());
	Message_info_value("  - hessmult =  ", qpsolver.get_hessmult_all());
	Message_info_time( "  - time =      ", qpsolver.get_time_total());
	Message_info_time( "   - t_project =  ", qpsolver.get_time_projection());
	Message_info_time( "   - t_matmult =  ", qpsolver.get_time_matmult());
	Message_info_time( "   - t_dot =      ", qpsolver.get_time_dot());
	Message_info_time( "   - t_update =   ", qpsolver.get_time_update());
	Message_info_time( "   - t_stepsize = ", qpsolver.get_time_stepsize());
	Message_info_time( "   - t_fs =       ", qpsolver.get_time_fs());
	Message_info_time( "   - t_other =    ", qpsolver.get_time_total() - (qpsolver.get_time_projection() + qpsolver.get_time_matmult() + qpsolver.get_time_dot() + qpsolver.get_time_update() + qpsolver.get_time_stepsize() + qpsolver.get_time_fs()));
	*/
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



