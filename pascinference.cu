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

#include <boost/program_options.hpp>

#define ALGORITHM_deltaL_eps 0.0001 /*stopping criteria of outer main loop */
#define ALGORITHM_max_s_steps 1 /* max number of outer steps */
#define ALGORITHM_EPSSQUARE 10.0 /* default FEM regularization parameter */

int T = 10; /* default length of generated time serie */
int K = 3; /* default number of clusters */
int dim = 2; /* default dimension of the problem */

/* load options from console arguments */
bool load_from_console(int argc, char *argv[]){
	bool return_value = true; /* continue of not? */

	namespace po = boost::program_options;

	/* define command line options */
	po::options_description description("PASC Inference Usage");
	description.add_options()
		("help,h", "Display this help message")
		("version,v", "Display the version number")
		("debug", po::value<int>(), "Debug mode")
		("length,T", po::value<int>(), "Length of time series")
		("clusters,K", po::value<int>(), "number of clusters");	
	
	/* parse command line arguments */	
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(description).run(), vm);
	po::notify(vm);		

	/* what to do with parsed arguments? */	
	if(vm.count("help")){ // TODO: this can be included in global application
		std::cout << description;
		return_value = false;
	}

	if(vm.count("version")){// TODO: this can be included in global application
		std::cout << "not implemented yet" << std::endl;
		return_value = false;
	}

	if(vm.count("debug")){// TODO: this can be included in global application
		DEBUG_MODE = vm["debug"].as<int>(); /* set global variable */
	}

	if(vm.count("length")){
		T = vm["length"].as<int>(); /* set global variable */
	}

	if(vm.count("clusters")){
		K = vm["clusters"].as<int>(); /* set global variable */
	}

	
	return return_value;
}

/* --- MAIN FUNCTION ---- */
int main( int argc, char *argv[] )
{
	/* load parameters from console input */
	if(!load_from_console(argc, argv)){
		return 0;
	}
		
	Initialize(argc, argv); // TODO: load parameters of problem from console input

	Timer timer_program; /* global timer for whole application */
	Timer timer_data; /* for generating the problem */
	Timer timer_model; /* for manipulation with model */

	timer_program.restart();
	timer_data.restart();
	timer_model.restart();

	/* say hello */	
	Message("- start program");
	timer_program.start(); /* here starts the timer for whole application */
	
	/* prepare data */
	Data_kmeans data;
	
	data.init(dim,T); // TODO: make it more funny using input K

	if(DEBUG_MODE >= 1) Message_info(" - generate data");
	timer_data.start(); // TODO: this timer should be elsewhere
	 data.generate();
	timer_data.stop();

	if(DEBUG_MODE >= 3) Message_info_time(" - problem generated in: ",timer_data.get_value_last());
	if(DEBUG_MODE >= 10)	data.print();

	/* prepare model */
	Model_kmeans model;
	timer_model.start(); 
	 model.init(dim,T,K);
	timer_model.stop();

	if(DEBUG_MODE >= 3) Message_info_time(" - model prepared in: ",timer_model.get_value_last());
	if(DEBUG_MODE >= 10)	model.print();

	/* prepare problem */
	Problem problem;
	problem.init();
	problem.set_data(data); /* set data to problem */
	problem.set_model(model); /* set model to problem */

	return 0;


	if(DEBUG_MODE >= 1) Message("- run main cycle:");
	 problem.solve(ALGORITHM_max_s_steps,ALGORITHM_deltaL_eps);
	if(DEBUG_MODE >= 1) Message("- main cycle finished");

	return 0;

	/* save the solution to VTK */
	if(EXPORT_SAVEVTK){
		if(DEBUG_MODE >= 1) Message("- save solution to VTK");
		problem.saveVTK("output/data.vtk");
	}

	/* here ends the application */
	problem.finalize();
	timer_program.stop();

	/* print final info about problem */
	if(DEBUG_MODE >= 3) problem.print_timers();
	if(DEBUG_MODE >= 10) problem.print();

	/* say bye */	
	Message("- end program");
	if(DEBUG_MODE >= 1) Message_info_time("- elapsed time: ",timer_program.get_value_sum());

	Finalize();
	return 0;
}

