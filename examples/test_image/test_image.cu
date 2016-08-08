/** @file test_image.cu
 *  @brief test the kmeans problem solver on image processing problem
 *
 *  Load image file, graph file and solve the problem.
 *
 *  @author Lukas Pospisil
 */

#include "pascinference.h"
#include "solver/tssolver.h"
#include "data/imagedata.h"
#include "model/graphh1fem.h"

#ifndef USE_PETSCVECTOR
 #error 'This example is for PETSCVECTOR'
#endif
 
using namespace pascinference;

typedef petscvector::PetscVector PetscVector;

extern int pascinference::DEBUG_MODE;

int main( int argc, char *argv[] )
{
	/* add local program options */
	boost::program_options::options_description opt_problem("PROBLEM EXAMPLE", consoleArg.get_console_nmb_cols());
	opt_problem.add_options()
		("test_K", boost::program_options::value<int>(), "number of clusters [int]")
		("test_image_filename", boost::program_options::value< std::string >(), "name of input file with image data (vector in PETSc format) [string]")
		("test_image_out", boost::program_options::value< std::string >(), "name of output file with image data (vector in PETSc format) [string]")
		("test_width", boost::program_options::value<int>(), "width of image [int]")
		("test_height", boost::program_options::value<int>(), "height of image [int]")
		("test_graph_filename", boost::program_options::value< std::string >(), "name of input file with graph data (vector in PETSc format) [string]")
		("test_graph_coeff", boost::program_options::value<double>(), "threshold of the graph [double]")
		("test_epssqr", boost::program_options::value<double>(), "penalty parameter [double]")
		("test_annealing", boost::program_options::value<int>(), "number of annealing steps [int]")
		("test_cutgamma", boost::program_options::value<bool>(), "cut gamma to {0,1} [bool]");
	consoleArg.get_description()->add(opt_problem);

	/* call initialize */
	if(!Initialize(argc, argv)){
		return 0;
	} 

	int K, annealing, width, height; 
	double epssqr;
	double graph_coeff; 
	bool cutgamma;

	std::string image_filename;
	std::string image_out;
	std::string graph_filename;

	consoleArg.set_option_value("test_K", &K, 2);
	consoleArg.set_option_value("test_width", &width, 100);
	consoleArg.set_option_value("test_height", &height, 100);
	consoleArg.set_option_value("test_epssqr", &epssqr, 10);
	consoleArg.set_option_value("test_graph_coeff", &graph_coeff, 1.1);
	consoleArg.set_option_value("test_cutgamma", &cutgamma, false);
	consoleArg.set_option_value("test_annealing", &annealing, 1);
	consoleArg.set_option_value("test_image_filename", &image_filename, "data/image1.bin");
	consoleArg.set_option_value("test_image_out", &image_out, "image1");
	consoleArg.set_option_value("test_graph_filename", &graph_filename, "data/graph1.bin");

	/* set decomposition in space */
	int DDR_size = GlobalManager.get_size();

	coutMaster << "- PROBLEM INFO ----------------------------" << std::endl;
	coutMaster << " DDR_size             = " << std::setw(30) << DDR_size << " (decomposition in space)" << std::endl;
	coutMaster << " test_image_filename  = " << std::setw(30) << image_filename << " (name of input file with image data)" << std::endl;
	coutMaster << " test_image_out       = " << std::setw(30) << image_out << " (part of name of output file)" << std::endl;
	coutMaster << " test_width           = " << width << " (width of image)" << std::endl;
	coutMaster << " test_height          = " << height << " (height of image)" << std::endl;
	coutMaster << " test_K               = " << K << " (number of clusters)" << std::endl;
	coutMaster << " test_graph_filename  = " << std::setw(30) << graph_filename << " (name of input file with graph data)" << std::endl;
	coutMaster << " test_graph_coeff     = " << graph_coeff << " (threshold of the graph)" << std::endl;
	coutMaster << " test_epssqr          = " << epssqr << " (penalty)" << std::endl;
	coutMaster << " test_annealing       = " << annealing << " (number of annealing steps)" << std::endl;
	coutMaster << " test_cutgamma        = " << cutgamma << " (cut gamma to {0,1})" << std::endl;
	coutMaster << "-------------------------------------------" << std::endl << std::endl;

	/* start logging */
	std::ostringstream oss_name_of_file_log;
	oss_name_of_file_log << "results/image_p" << GlobalManager.get_rank() << ".txt";
	logging.begin(oss_name_of_file_log.str());
		
	/* say hello */
	coutMaster << "- start program" << std::endl;

/* 1.) prepare graph of image */
	coutMaster << "--- PREPARING GRAPH ---" << std::endl;
	BGMGraphGrid2D graph(width, height);
	graph.process_grid();

//	BGMGraph mygraph(graph_filename);
//	mygraph.process(graph_coeff);

	graph.print(coutMaster);
	
/* 2.) prepare decomposition */
	coutMaster << "--- COMPUTING DECOMPOSITION ---" << std::endl;
	Decomposition decomposition(1, graph, K, 1, DDR_size);
	decomposition.print(coutMaster);

	graph.saveVTK("results/graph.vtk");

/* 3.) prepare time-series data */
	coutMaster << "--- PREPARING DATA ---" << std::endl;
	ImageData<PetscVector> mydata(decomposition, image_filename, width, height);
	mydata.print(coutMaster);

/* 4.) prepare model */
	coutMaster << "--- PREPARING MODEL ---" << std::endl;
	GraphH1FEMModel<PetscVector> mymodel(mydata, epssqr);
	mymodel.print(coutMaster,coutAll);

/* 5.) prepare time-series solver */
	coutMaster << "--- PREPARING SOLVER ---" << std::endl;
	TSSolver<PetscVector> mysolver(mydata, annealing);
	mysolver.print(coutMaster,coutAll);

/* 6.) solve the problem */
	coutMaster << "--- SOLVING THE PROBLEM ---" << std::endl;
	mysolver.solve();

	/* cut gamma */
	if(cutgamma){
		mydata.cut_gamma();
	}

	/* print solution */
	coutMaster << "--- THETA SOLUTION ---" << std::endl;
	mydata.print_thetavector(coutMaster);

/* 7.) save results */
	coutMaster << "--- SAVING OUTPUT ---" << std::endl;
	mydata.saveImage(image_out);

	/* print timers */
	coutMaster << "--- TIMERS INFO ---" << std::endl;
	mysolver.printtimer(coutMaster);

	/* print short info */
	coutMaster << "--- FINAL SOLVER INFO ---" << std::endl;
	mysolver.printstatus(coutMaster);

	/* say bye */	
	coutMaster << "- end program" << std::endl;

	logging.end();
	Finalize();


	return 0;
}

