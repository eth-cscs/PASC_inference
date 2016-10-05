/** @file test_graph.cpp
 *  @brief test class and methods: BGMGraph 
 *
 *  This file tests the class and metis decomposition of the graph.
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
		("test_graph_filename", boost::program_options::value< std::string >(), "name of file with coordinates [string]")
		("test_graph_dim", boost::program_options::value<int>(), "dimension of graph (1,2,3) [int]")
		("test_graph_out", boost::program_options::value< std::string >(), "part of name of output file with graph [string]")
		("test_graph_coeff", boost::program_options::value<double>(), "threshold of the graph [double]")
		("test_graph_nmbdomains", boost::program_options::value<int>(), "number of domains for decomposition [int]")
		("test_graph_print", boost::program_options::value<bool>(), "print content of graph or not [bool]");
	consoleArg.get_description()->add(opt_problem);

	/* call initialize */
	if(!Initialize(argc, argv)){
		return 0;
	} 

	/* load console arguments */
	std::string filename, out;
	double coeff;
	bool print;
	int nmbdomains, dim;
	consoleArg.set_option_value("test_graph_filename", &filename, "data/test_graph2D.bin");
	consoleArg.set_option_value("test_graph_dim", &dim, 2);
	consoleArg.set_option_value("test_graph_out", &out, "test_graph_out");
	consoleArg.set_option_value("test_graph_coeff", &coeff, 1.1);
	consoleArg.set_option_value("test_graph_nmbdomains", &nmbdomains, 3);
	consoleArg.set_option_value("test_graph_print", &print, false);

	/* print settings */
	coutMaster << " test_graph_filename       = " << std::setw(30) << filename << " (name of file with coordinates)" << std::endl;
	coutMaster << " test_graph_dim            = " << std::setw(30) << dim << " (dimension of graph (1,2,3))" << std::endl;
	coutMaster << " test_graph_out            = " << std::setw(30) << out << " (part of name of output file with graph)" << std::endl;
	coutMaster << " test_graph_coeff          = " << std::setw(30) << coeff << " (threshold of the graph)" << std::endl;
	coutMaster << " test_graph_nmbdomains     = " << std::setw(30) << nmbdomains << " (number of domains for decomposition)" << std::endl;
	coutMaster << " test_graph_print          = " << std::setw(30) << print << " (print content of graph or not)" << std::endl;

	/* we will measure the time of operations */
	Timer mytimer;
	mytimer.restart();
	
	/* create graph (i.e. load from filename) */
	mytimer.start();
	 BGMGraph graph(filename,dim);
	mytimer.stop();
	coutMaster << "- time load      : " << mytimer.get_value_last() << " s" << std::endl;
	
	/* process graph with given coefficient */
	mytimer.start();
	 graph.process(coeff);
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
	coutMaster << "also look into folder data/ to see other graph examples" << std::endl;
	coutMaster << "and switch between them using for example ./test_bgmgraph --test_graph_filename=\"data/test_graph1D\" --test_graph_dim=1\"" << std::endl;

	Finalize();

	return 0;
}


