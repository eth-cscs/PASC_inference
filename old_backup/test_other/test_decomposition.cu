/** @file test_decomposition.cu
 *  @brief test decoposition in time and space
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
		("test_T", boost::program_options::value<int>(), "length of time-series [int]")
		("test_DDT", boost::program_options::value<int>(), "decomposition in time [int]")
		("test_DDR", boost::program_options::value<int>(), "decomposition in space [int]")
		("test_graph_filename", boost::program_options::value< std::string >(), "name of file with coordinates [string]")
		("test_graph_out", boost::program_options::value< std::string >(), "part of name of output file with graph [string]")
		("test_graph_coeff", boost::program_options::value<double>(), "threshold of the graph [double]")
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
	int T, DDT_size, DDR_size;
	int nproc = GlobalManager.get_size();

	consoleArg.set_option_value("test_T", &T, 10);
	consoleArg.set_option_value("test_DDT", &DDT_size, nproc);
	consoleArg.set_option_value("test_DDR", &DDR_size, 1);
	consoleArg.set_option_value("test_graph_filename", &graph_filename, "data/graph_square_10.bin");
	consoleArg.set_option_value("test_graph_out", &graph_out, "graph_square_10");
	consoleArg.set_option_value("test_graph_coeff", &graph_coeff, 1.1);
	consoleArg.set_option_value("test_view_graph", &view_graph, false);

	/* print settings */
	coutMaster << " test_T                    = " << std::setw(30) << graph_filename << " (length of time-series)\n";
	coutMaster << " test_graph_filename       = " << std::setw(30) << graph_filename << " (name of file with coordinates)\n";
	coutMaster << " test_graph_out            = " << std::setw(30) << graph_filename << " (part of name of output file with graph)\n";
	coutMaster << " test_graph_coeff          = " << std::setw(30) << graph_coeff << " (threshold of the graph)\n";
	coutMaster << " test_view_graph           = " << std::setw(30) << view_graph << " (print content of graph or not)\n";
	coutMaster << " test_DDT                  = " << std::setw(30) << DDT_size << " (decomposition in time)\n";
	coutMaster << " test_DDR                  = " << std::setw(30) << DDR_size << " (decomposition in space)\n";

	if(DDT_size*DDR_size != nproc){
		coutMaster << "Sorry, DDT*DDR != nproc\n";
		coutMaster << " DDT = " << DDT_size << "\n";
		coutMaster << " DDR = " << DDR_size << "\n";
		coutMaster << " nproc    = " << nproc << "\n";
		
		return 0;
	}

	/* start logging */
	std::ostringstream oss_name_of_file_log;
	oss_name_of_file_log << "results/test_decomposition_log_p" << GlobalManager.get_rank() << ".txt";
	logging.begin(oss_name_of_file_log.str());
		
	/* say hello */
	coutMaster << "- start program\n";

	/* create graph */
	BGMGraph graph(graph_filename);
	
	/* process graph */
	graph.process(graph_coeff);

	/* print info about graph */
	if(!view_graph){
		graph.print(coutMaster);
	} else {
		graph.print_content(coutMaster);
	}

	/* decomposition */
//	Decomposition decomposition(T, DDT_size); /* pure time */
//	Decomposition decomposition(T, graph, DDR_size); /* pure space */
	Decomposition decomposition(T, graph, 1, 1, DDT_size, DDR_size); /* time and space */

	decomposition.print_content(coutMaster, coutAll);

	/* say bye */	
	coutMaster << "- end program\n";

	logging.end();
	Finalize();

	return 0;
}


