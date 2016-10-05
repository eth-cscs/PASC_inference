/** @file test_bgmgraphgrid1D.cpp
 *  @brief test class and methods: BGMGraphGrid1D
 *
 *  This file tests the class the regular 1D mesh (i.e. line). The distance between nodes is always equal to 1.
 *  This mesh is decomposed without metis (since it is trivial problem).
 * 
 *  @author Lukas Pospisil
 */

#include <iostream>
#include <list>
#include <algorithm>

#include "pascinference.h"

typedef petscvector::PetscVector PetscVector;

using namespace pascinference;

extern int pascinference::DEBUG_MODE;

int main( int argc, char *argv[] )
{
	/* add local program options */
	boost::program_options::options_description opt_problem("PROBLEM EXAMPLE", consoleArg.get_console_nmb_cols());
	opt_problem.add_options()
		("test_length", boost::program_options::value<int>(), "number of 1D grid nodes (i.e. line length) [int]")
		("test_out", boost::program_options::value< std::string >(), "part of name of output file with graph [string]")
		("test_generalprocess", boost::program_options::value<bool>(), "use general processing method (otherwise graph processed as grid without use of given coeff) [bool]")
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
	int nmbdomains, length;
	consoleArg.set_option_value("test_length", &length, 10);
	consoleArg.set_option_value("test_out", &out, "test_graph_out");
	consoleArg.set_option_value("test_coeff", &coeff, 1.1);
	consoleArg.set_option_value("test_nmbdomains", &nmbdomains, 3);
	consoleArg.set_option_value("test_print", &print, false);
	consoleArg.set_option_value("test_generalprocess", &generalprocess, false);

	/* print settings */
	coutMaster << " test_length         = " << std::setw(30) << length << " (number of 1D grid nodes (i.e. line length))" << std::endl;
	coutMaster << " test_out            = " << std::setw(30) << out << " (part of name of output file with graph)" << std::endl;
	coutMaster << " test_generalprocess = " << std::setw(30) << generalprocess << " (use general processing method (otherwise graph processed as grid without use of given coeff))" << std::endl;
	coutMaster << " test_coeff          = " << std::setw(30) << coeff << " (threshold of the graph)" << std::endl;
	coutMaster << " test_nmbdomains     = " << std::setw(30) << nmbdomains << " (number of domains for decomposition)" << std::endl;
	coutMaster << " test_print          = " << std::setw(30) << print << " (print content of graph or not)" << std::endl;

	/* we will measure the time of operations */
	Timer mytimer;
	mytimer.restart();
	
	/* create graph */
	mytimer.start();
	 BGMGraphGrid1D graph(length);
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


