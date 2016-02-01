/*******************************************************************************
PASC INFERENCE library
Lukas Pospisil, Illia Horenko, Patrick Gagliardini, Will Sawyer
USI Lugano, 2016
lukas.pospisil@usi.ch

*******************************************************************************/

#include "common.h"
#include "problem.h"
#include "data.h"
#include "model.h"

/* PROBLEM SETTINGS */
#define DEFAULT_T 10 /* default length of generated time serie */
#define DEFAULT_K 3 /* default number of clusters */


#define ALGORITHM_deltaL_eps 0.0001 /*stopping criteria of outer main loop */
#define ALGORITHM_max_s_steps 100 /* max number of outer steps */
#define ALGORITHM_EPSSQUARE 10.0 /* FEM regularization parameter */


int main( int argc, char *argv[] )
{
	Initialize(argc,argv); // TODO: load parameters of problem from console input

	Timer timer_all; /* global timer for whole application */
	Timer timer_data; /* for generating the problem */
	Timer timer_model; /* for manipulation with model */

	Timer timer_saveVTK; /* for final export to VTK */

	timer_all.restart();
	timer_data.restart();
	timer_model.restart();
//	timer_gamma.restart();
//	timer_theta.restart();
	timer_saveVTK.restart();

	/* say hello */	
	Message("- start program");
	timer_all.start(); /* here starts the timer for whole application */
	
	/* prepare data */
	Data_kmeans data;
	data.init(2,DEFAULT_T);
	timer_data.start(); 
	 data.generate();
	timer_data.stop();

	Message_info_time(" - problem generated in: ",timer_data.get_value_last());
	/* print problem */
	if(DEBUG_MODE >= 3)	data.print();

	/* prepare model */
	Model_kmeans model;
	timer_model.start(); 
	 model.init(2,DEFAULT_T,DEFAULT_K);
	timer_model.stop();

	Message_info_time(" - model prepared in: ",timer_model.get_value_last());
	/* print problem */
	if(DEBUG_MODE >= 3)	model.print();

	/* prepare problem */
	Problem problem;
	problem.init();
	problem.set_data(data); /* set data to problem */
	problem.set_model(model); /* set model to problem */


	Message("- run main cycle:");
	 problem.solve(ALGORITHM_max_s_steps,ALGORITHM_deltaL_eps);
	Message("- main cycle finished");

	if(DEBUG_MODE >= 2)	problem.print();

	/* save the solution to VTK */
	if(EXPORT_SAVEVTK){
		Message("- save solution to VTK");
//		timer_saveVTK.start();
		problem.saveVTK("output/data.vtk");
//		timer_saveVTK.stop();
	}

	/* here ends the application */
	problem.finalize();
	timer_all.stop();

	/* say bye */	
	Message("- end program");

	
	Finalize();
	return 0;
}

