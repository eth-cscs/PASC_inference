/** @file test_bgmgraphgrid2D.cpp
 *  @brief test class and methods: BGMGraphGrid2D
 *
 *  This file tests the class the regular 2D mesh. The distance between nodes is always equal to 1.
 * 
 *  @author Lukas Pospisil
 */

#include <iostream>
#include <list>
#include <algorithm>

#include "pascinference.h"

#ifndef USE_METIS
 #error 'This example is for METIS'
#endif

typedef petscvector::PetscVector PetscVector;

using namespace pascinference;

extern int pascinference::DEBUG_MODE;

int main( int argc, char *argv[] )
{
	/* add local program options */
	boost::program_options::options_description opt_problem("PROBLEM EXAMPLE", consoleArg.get_console_nmb_cols());
	opt_problem.add_options()
		("test_width", boost::program_options::value<int>(), "number of grid nodes on x axis [int]")
		("test_height", boost::program_options::value<int>(), "number of grid nodes on y axis [int]")
		("test_out", boost::program_options::value< std::string >(), "part of name of output file with graph [string]")
		("test_generalprocess", boost::program_options::value<bool>(), "use general processing method (otherwise graph processed as squared grid without using given coeff) [bool]")
		("test_coeff", boost::program_options::value<double>(), "threshold of the graph [double]")
		("test_nmbdomains", boost::program_options::value<int>(), "number of domains for decomposition [int]")
		("test_print", boost::program_options::value<bool>(), "print content of graph or not [bool]");
	consoleArg.get_description()->add(opt_problem);

	/* call initialize */
	if(!Initialize(argc, argv)){
		return 0;
	} 

	/* load console arguments */
	std::string out;
	double coeff;
	bool print, generalprocess;
	int nmbdomains, width, height;
	consoleArg.set_option_value("test_width", &width, 15);
	consoleArg.set_option_value("test_height", &height, 12);
	consoleArg.set_option_value("test_out", &out, "test_graph_out");
	consoleArg.set_option_value("test_coeff", &coeff, 1.1);
	consoleArg.set_option_value("test_nmbdomains", &nmbdomains, 3);
	consoleArg.set_option_value("test_print", &print, false);
	consoleArg.set_option_value("test_generalprocess", &generalprocess, false);

	/* print settings */
	coutMaster << " test_width          = " << std::setw(30) << width << " (number of grid nodes on x axis)" << std::endl;
	coutMaster << " test_height         = " << std::setw(30) << height << " (number of grid nodes on y axis)" << std::endl;
	coutMaster << " test_out            = " << std::setw(30) << out << " (part of name of output file with graph)" << std::endl;
	coutMaster << " test_generalprocess = " << std::setw(30) << generalprocess << " (use general processing method (otherwise graph processed as squared grid without using given coeff))" << std::endl;
	coutMaster << " test_coeff          = " << std::setw(30) << coeff << " (threshold of the graph)" << std::endl;
	coutMaster << " test_nmbdomains     = " << std::setw(30) << nmbdomains << " (number of domains for decomposition)" << std::endl;
	coutMaster << " test_print          = " << std::setw(30) << print << " (print content of graph or not)" << std::endl;

	/* we will measure the time of operations */
	Timer mytimer;
	mytimer.restart();
	
	/* create graph */
	mytimer.start();
	 BGMGraphGrid2D graph(width, height);
	mytimer.stop();
	coutMaster << "- time load      : " << mytimer.get_value_last() << " s" << std::endl;
	
	/* process graph with given coefficient */
	mytimer.start();
	 if(generalprocess){
		graph.process(coeff);
	 } else {
		graph.process_grid();
	 }
	mytimer.stop();
	coutMaster << "- time process   : " << mytimer.get_value_last() << " s" << std::endl;

	/* call domain decomposition */
	mytimer.start();
	 graph.decompose(nmbdomains);
	mytimer.stop();
	coutMaster << "- time decompose : " << mytimer.get_value_last() << " s" << std::endl;

	/* print info about graph */
	if(!print){
		graph.print(coutMaster);
	} else {
		graph.print_content(coutMaster);
	}

	/* save graph into VTK */
	std::ostringstream oss_graph_out_filename;
	oss_graph_out_filename << "results/" << out << ".vtk";
	mytimer.start();
	 graph.saveVTK(oss_graph_out_filename.str());
	mytimer.stop();
	coutMaster << "- time save      : " << mytimer.get_value_last() << " s" << std::endl;

	/* say bye */	
	coutMaster << std::endl;
	coutMaster << "please find generated .vtk file with graph: " << oss_graph_out_filename.str() << std::endl << std::endl;

	Finalize();

	return 0;
}


