/** @file test_metis.cu
 *  @brief test graph decomposition using metis
 *
 *  @author Lukas Pospisil
 */

#include <iostream>
#include <list>
#include <algorithm>

#include "pascinference.h"
#include "common/bgmgraph.h"

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
		("test_graph_out", boost::program_options::value< std::string >(), "part of name of output file with graph [string]")
		("test_graph_coeff", boost::program_options::value<double>(), "threshold of the graph [double]")
		("test_graph_nmb_domains", boost::program_options::value<int>(), "number of domains for decomposition [int]")
		("test_view_graph", boost::program_options::value<bool>(), "print content of graph or not [bool]");
	consoleArg.get_description()->add(opt_problem);

	/* call initialize */
	if(!Initialize(argc, argv)){
		return 0;
	} 

/*	if(GlobalManager.get_size() > 1){
		coutMaster << "This example works only on one processor, sorry.\n";
		return 0;		
	}
*/

	/* load console arguments */
	std::string graph_filename, graph_out;
	double graph_coeff;
	bool view_graph;
	int graph_nmb_domains;

	consoleArg.set_option_value("test_graph_filename", &graph_filename, "data/graph_square_10.bin");
	consoleArg.set_option_value("test_graph_out", &graph_out, "graph_square_10");
	consoleArg.set_option_value("test_graph_coeff", &graph_coeff, 1.1);
	consoleArg.set_option_value("test_graph_nmb_domains", &graph_nmb_domains, 3);
	consoleArg.set_option_value("test_view_graph", &view_graph, false);

	/* print settings */
	coutMaster << " test_graph_filename       = " << std::setw(30) << graph_filename << " (name of file with coordinates)\n";
	coutMaster << " test_graph_out            = " << std::setw(30) << graph_filename << " (part of name of output file with graph)\n";
	coutMaster << " test_graph_coeff          = " << std::setw(30) << graph_coeff << " (threshold of the graph)\n";
	coutMaster << " test_graph_nmb_domains    = " << std::setw(30) << graph_nmb_domains << " (number of domains for decomposition)\n";
	coutMaster << " test_view_graph           = " << std::setw(30) << view_graph << " (print content of graph or not)\n";

	/* start logging */
	std::ostringstream oss_name_of_file_log;
	oss_name_of_file_log << "results/test_metis_log_p" << GlobalManager.get_rank() << ".txt";
	logging.begin(oss_name_of_file_log.str());
		
	/* say hello */
	coutMaster << "- start program\n";

	/* create graph */
	BGMGraph graph(graph_filename);
	
	/* process graph */
	graph.process(graph_coeff);

	/* call domain decomposition */
	graph.decompose(graph_nmb_domains);

	/* print info about graph */
	if(!view_graph){
		graph.print(coutMaster);
	} else {
		graph.print_content(coutMaster);
	}

	/* save graph into VTK */
	std::ostringstream oss_graph_out_filename;
	oss_graph_out_filename << "results/" << graph_out << ".vtk";
	graph.saveVTK(oss_graph_out_filename.str());

	/* say bye */	
	coutMaster << "- end program\n";

	logging.end();
	Finalize();

	return 0;
}


