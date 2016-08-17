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

#include <vector>

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
		("test_cutgamma", boost::program_options::value<bool>(), "cut gamma to set {0;1} [bool]")
		("test_scaledata", boost::program_options::value<bool>(), "scale to interval {-1,1} [bool]")
		("test_cutdata", boost::program_options::value<bool>(), "cut data to interval {0,1} [bool]")
		("test_shiftdata", boost::program_options::value<bool>(), "shift data with -0.5 [bool]")
		("test_shiftdata_coeff", boost::program_options::value<double>(), "coeficient of data shift [double]")
		("test_printstats", boost::program_options::value<bool>(), "print basic statistics of data [bool]")
		("test_Theta", boost::program_options::value<std::vector<double> >()->multitoken(), "given solution Theta [K*int]")
		("test_shortinfo", boost::program_options::value<bool>(), "save shortinfo file after computation [bool]")
		("test_shortinfo_header", boost::program_options::value< std::string >(), "additional header in shortinfo [string]")
		("test_shortinfo_values", boost::program_options::value< std::string >(), "additional values in shortinfo [string]")
		("test_shortinfo_filename", boost::program_options::value< std::string >(), "name of shortinfo file [string]");

	consoleArg.get_description()->add(opt_problem);

	/* call initialize */
	if(!Initialize(argc, argv)){
		return 0;
	} 

	int K, annealing, width, height; 
	double epssqr;
	double graph_coeff, shiftdata_coeff; 
	bool cutgamma, scaledata, cutdata, printstats, shiftdata, shortinfo;

	std::string image_filename;
	std::string image_out;
	std::string graph_filename;
	std::string shortinfo_filename;
	std::string shortinfo_header;
	std::string shortinfo_values;

	consoleArg.set_option_value("test_K", &K, 2);
	consoleArg.set_option_value("test_width", &width, 250);
	consoleArg.set_option_value("test_height", &height, 150);
	consoleArg.set_option_value("test_epssqr", &epssqr, 0.00001);
	consoleArg.set_option_value("test_graph_coeff", &graph_coeff, 1.1);
	consoleArg.set_option_value("test_cutgamma", &cutgamma, false);
	consoleArg.set_option_value("test_scaledata", &scaledata, false);
	consoleArg.set_option_value("test_cutdata", &cutdata, true);
	consoleArg.set_option_value("test_shiftdata", &shiftdata, false);
	consoleArg.set_option_value("test_shiftdata_coeff", &shiftdata_coeff, -0.5);
	consoleArg.set_option_value("test_printstats", &printstats, false);
	consoleArg.set_option_value("test_annealing", &annealing, 1);
	consoleArg.set_option_value("test_image_filename", &image_filename, "data/usi_text/usi_250_150_02.bin");
	consoleArg.set_option_value("test_image_out", &image_out, "usi_test_250_150_02");
	consoleArg.set_option_value("test_shortinfo", &shortinfo, true);
	consoleArg.set_option_value("test_shortinfo_header", &shortinfo_header, "");
	consoleArg.set_option_value("test_shortinfo_values", &shortinfo_values, "");
	consoleArg.set_option_value("test_shortinfo_filename", &shortinfo_filename, "shortinfo/usi_250_150_02.txt");
	
	/* use general graph or 2d grid? */
	bool generalgraph;
	if(consoleArg.set_option_value("test_graph_filename", &graph_filename)){
		generalgraph = true;
	} else {
		generalgraph = false;
	}

	/* maybe theta is given in console parameters */
	bool given_Theta;
	std::vector<double> Theta_list;
	double Theta_solution[K];
	if(consoleArg.set_option_value("test_Theta", &Theta_list)){
		given_Theta = true;
		
		/* control number of provided Theta */
		if(Theta_list.size() != K){
			coutMaster << "number of provided Theta solutions is different then number of clusters!" << std::endl;
			return 0;
		}

		/* store solution in array */
		for(int k=0;k < K;k++){
			Theta_solution[k] = Theta_list[k];
		}
	} else {
		given_Theta = false;
	}	

	/* set decomposition in space */
	int DDR_size = GlobalManager.get_size();

	coutMaster << "- PROBLEM INFO ----------------------------" << std::endl;
	coutMaster << " DDR_size                = " << std::setw(30) << DDR_size << " (decomposition in space)" << std::endl;
	coutMaster << " test_image_filename     = " << std::setw(30) << image_filename << " (name of input file with image data)" << std::endl;
	coutMaster << " test_image_out          = " << std::setw(30) << image_out << " (part of name of output file)" << std::endl;
	coutMaster << " test_width              = " << std::setw(30) << width << " (width of image)" << std::endl;
	coutMaster << " test_height             = " << std::setw(30) << height << " (height of image)" << std::endl;
	coutMaster << " test_K                  = " << std::setw(30) << K << " (number of clusters)" << std::endl;
	if(given_Theta){
		coutMaster << " test_Theta              = " << std::setw(30) << "[";
		for(int k=0;k<K;k++){
			coutMaster << Theta_solution[k];
			if(k<K-1){
				coutMaster << ",";
			}
		}
		coutMaster << "]" << std::endl;
	}
	coutMaster << " test_graph_filename     = " << std::setw(30) << graph_filename << " (name of input file with graph data)" << std::endl;
	coutMaster << " test_graph_coeff        = " << std::setw(30) << graph_coeff << " (threshold of the graph)" << std::endl;
	coutMaster << " test_epssqr             = " << std::setw(30) << epssqr << " (penalty)" << std::endl;
	coutMaster << " test_annealing          = " << std::setw(30) << annealing << " (number of annealing steps)" << std::endl;
	coutMaster << " test_cutgamma           = " << std::setw(30) << cutgamma << " (cut gamma to {0;1})" << std::endl;
	coutMaster << " test_cutdata            = " << std::setw(30) << cutdata << " (cut data to {0,1})" << std::endl;
	coutMaster << " test_scaledata          = " << std::setw(30) << scaledata << " (scale data to {-1,1})" << std::endl;
	coutMaster << " test_shiftdata          = " << std::setw(30) << shiftdata << " (shift data by coeficient)" << std::endl;
	coutMaster << " test_shiftdata_coeff    = " << std::setw(30) << shiftdata_coeff << " (shifting coeficient)" << std::endl;
	coutMaster << " test_shortinfo          = " << std::setw(30) << shortinfo << " (save shortinfo file after computation)" << std::endl;
	coutMaster << " test_shortinfo_header   = " << std::setw(30) << shortinfo_header << " (additional header in shortinfo)" << std::endl;
	coutMaster << " test_shortinfo_values   = " << std::setw(30) << shortinfo_values << " (additional values in shortinfo)" << std::endl;
	coutMaster << " test_shortinfo_filename = " << std::setw(30) << shortinfo_filename << " (name of shortinfo file)" << std::endl;
	coutMaster << "-------------------------------------------" << std::endl << std::endl;

	/* start logging */
	std::ostringstream oss_name_of_file_log;
	oss_name_of_file_log << "results/image_p" << GlobalManager.get_rank() << ".txt";
	logging.begin(oss_name_of_file_log.str());
		
	/* say hello */
	coutMaster << "- start program" << std::endl;

/* 1.) prepare graph of image */
	coutMaster << "--- PREPARING GRAPH ---" << std::endl;

	BGMGraph *graph;
	if(generalgraph){
		graph = new BGMGraph(graph_filename);
		graph->process(graph_coeff);
	} else {
		graph = new BGMGraphGrid2D(width, height);
		((BGMGraphGrid2D*)graph)->process_grid();
	}
	graph->print(coutMaster);
	
/* 2.) prepare decomposition */
	coutMaster << "--- COMPUTING DECOMPOSITION ---" << std::endl;
	Decomposition decomposition(1, *graph, K, 1, DDR_size);
	decomposition.print(coutMaster);

	graph->saveVTK("results/graph.vtk");

/* 3.) prepare time-series data */
	coutMaster << "--- PREPARING DATA ---" << std::endl;
	ImageData<PetscVector> mydata(decomposition, image_filename, width, height);
	mydata.print(coutMaster);

	/* print statistics */
	if(printstats){
		mydata.printstats(coutMaster);
	}

	/* cut data */
	if(cutdata){
		mydata.cutdata(0,1);
	}

	/* shift data */
	if(shiftdata){
		mydata.shiftdata(shiftdata_coeff);
	}

	/* scale data */
	if(scaledata){
		mydata.scaledata(-1,1,0,1);
	}

/* 4.) prepare model */
	coutMaster << "--- PREPARING MODEL ---" << std::endl;
	GraphH1FEMModel<PetscVector> mymodel(mydata, epssqr);
	mymodel.print(coutMaster,coutAll);

/* 5.) prepare time-series solver */
	coutMaster << "--- PREPARING SOLVER ---" << std::endl;
	TSSolver<PetscVector> mysolver(mydata, annealing);
	mysolver.print(coutMaster,coutAll);

	/* set solution if obtained from console */
	if(given_Theta){
		mysolver.set_solution_theta(Theta_solution);
	}
	
/* 6.) solve the problem */
	coutMaster << "--- SOLVING THE PROBLEM ---" << std::endl;
	mysolver.solve();

	/* cut gamma */
	if(cutgamma){
		mydata.cutgamma();
	}

	/* unscale data */
	if(scaledata){
		mydata.scaledata(0,1,-1,1);
	}

	/* unshift data */
	if(shiftdata){
		mydata.shiftdata(-shiftdata_coeff);
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

	/* write short output */
	if(shortinfo){
		std::ostringstream oss_short_output_values;
		std::ostringstream oss_short_output_header;
		
		/* add provided strings from console parameters */
		oss_short_output_header << shortinfo_header;
		oss_short_output_values << shortinfo_values;
		
		mysolver.printshort(oss_short_output_header, oss_short_output_values);

		std::ofstream myfile;
		myfile.open(shortinfo_filename.c_str());
		
		/* master writes the file with short info (used in batch script for quick computation) */
		if( GlobalManager.get_rank() == 0){
			myfile << oss_short_output_header.str() << std::endl;
			myfile << oss_short_output_values.str() << std::endl;
		}
		TRY( PetscBarrier(NULL) );
		
		myfile.close();
	}

	/* say bye */	
	coutMaster << "- end program" << std::endl;

	logging.end();
	Finalize();


	return 0;
}

