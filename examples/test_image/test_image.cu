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
		("test_graph_save", boost::program_options::value<bool>(), "save VTK with graph or not [bool]")
		("test_epssqr", boost::program_options::value<std::vector<double> >()->multitoken(), "penalty parameters [double]")
		("test_annealing", boost::program_options::value<int>(), "number of annealing steps [int]")
		("test_cutgamma", boost::program_options::value<bool>(), "cut gamma to set {0;1} [bool]")
		("test_scaledata", boost::program_options::value<bool>(), "scale to interval {-1,1} [bool]")
		("test_cutdata", boost::program_options::value<bool>(), "cut data to interval {0,1} [bool]")
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

	std::vector<double> epssqr_list;
	if(consoleArg.set_option_value("test_epssqr", &epssqr_list)){
		/* sort list */
		std::sort(epssqr_list.begin(), epssqr_list.end(), std::less<double>());
		
	} else {
		std::cout << "test_epssqr has to be set! Call application with parameter -h to see all parameters\n";
		return 0;
	}

	int K, annealing, width, height; 
	double graph_coeff; 
	bool cutgamma, scaledata, cutdata, printstats, shortinfo_write_or_not, graph_save;

	std::string image_filename;
	std::string image_out;
	std::string graph_filename;
	std::string shortinfo_filename;
	std::string shortinfo_header;
	std::string shortinfo_values;

	consoleArg.set_option_value("test_K", &K, 2);
	consoleArg.set_option_value("test_width", &width, 250);
	consoleArg.set_option_value("test_height", &height, 150);
	consoleArg.set_option_value("test_graph_coeff", &graph_coeff, 1.1);
	consoleArg.set_option_value("test_graph_save", &graph_save, false);
	consoleArg.set_option_value("test_cutgamma", &cutgamma, false);
	consoleArg.set_option_value("test_scaledata", &scaledata, false);
	consoleArg.set_option_value("test_cutdata", &cutdata, true);
	consoleArg.set_option_value("test_printstats", &printstats, false);
	consoleArg.set_option_value("test_annealing", &annealing, 1);
	consoleArg.set_option_value("test_image_filename", &image_filename, "data/usi_text/usi_250_150_02.bin");
	consoleArg.set_option_value("test_image_out", &image_out, "usi_test_250_150_02");
	consoleArg.set_option_value("test_shortinfo", &shortinfo_write_or_not, true);
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
			coutMaster << "number of provided Theta solutions is different then number of clusters!\n";
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

	coutMaster << "- PROBLEM INFO ----------------------------\n";
	coutMaster << " DDR_size                = " << std::setw(30) << DDR_size << " (decomposition in space)\n";
	coutMaster << " test_image_filename     = " << std::setw(30) << image_filename << " (name of input file with image data)\n";
	coutMaster << " test_image_out          = " << std::setw(30) << image_out << " (part of name of output file)\n";
	coutMaster << " test_width              = " << std::setw(30) << width << " (width of image)\n";
	coutMaster << " test_height             = " << std::setw(30) << height << " (height of image)\n";
	coutMaster << " test_K                  = " << std::setw(30) << K << " (number of clusters)\n";
	if(given_Theta){
		coutMaster << " test_Theta              = " << std::setw(30) << print_array(Theta_solution,K) << "\n";
	}
	coutMaster << " test_graph_filename     = " << std::setw(30) << graph_filename << " (name of input file with graph data)\n";
	coutMaster << " test_graph_coeff        = " << std::setw(30) << graph_coeff << " (threshold of the graph)\n";
	coutMaster << " test_graph_save         = " << std::setw(30) << graph_save << " (save VTK with graph or not)\n";
	coutMaster << " test_epssqr             = " << std::setw(30) << print_vector(epssqr_list) << " (penalty)\n";
	coutMaster << " test_annealing          = " << std::setw(30) << annealing << " (number of annealing steps)\n";
	coutMaster << " test_cutgamma           = " << std::setw(30) << cutgamma << " (cut gamma to {0;1})\n";
	coutMaster << " test_cutdata            = " << std::setw(30) << cutdata << " (cut data to {0,1})\n";
	coutMaster << " test_scaledata          = " << std::setw(30) << scaledata << " (scale data to {-1,1})\n";
	coutMaster << " test_shortinfo          = " << std::setw(30) << shortinfo_write_or_not << " (save shortinfo file after computation)\n";
	coutMaster << " test_shortinfo_header   = " << std::setw(30) << shortinfo_header << " (additional header in shortinfo)\n";
	coutMaster << " test_shortinfo_values   = " << std::setw(30) << shortinfo_values << " (additional values in shortinfo)\n";
	coutMaster << " test_shortinfo_filename = " << std::setw(30) << shortinfo_filename << " (name of shortinfo file)\n";
	coutMaster << "-------------------------------------------\n" << "\n";

	/* start logging */
	std::ostringstream oss;
	oss << "log/" << image_out << ".txt";
	logging.begin(oss.str());
	oss.str("");

	/* start shortinfo output */
	if(shortinfo_write_or_not){
		shortinfo.begin(shortinfo_filename);
	}
	std::ostringstream oss_short_output_values;
	std::ostringstream oss_short_output_header;
		
	/* say hello */
	coutMaster << "- start program\n";

/* 1.) prepare graph of image */
	coutMaster << "--- PREPARING GRAPH ---\n";

	BGMGraph *graph;
	if(generalgraph){
		/* if this is general graph, then load it from file and process it with given threshold */
		graph = new BGMGraph(graph_filename);
		graph->process(graph_coeff);
	} else {
		/* this is not general graph, it is "only" 2D grid */
		graph = new BGMGraphGrid2D(width, height);
		((BGMGraphGrid2D*)graph)->process_grid();
	}

	/* print basic info about graph */
	graph->print(coutMaster);
	
/* 2.) prepare decomposition */
	coutMaster << "--- COMPUTING DECOMPOSITION ---\n";

	/* prepare decomposition based on graph, in this case T=1 and DDT_size=1 */
	Decomposition decomposition(1, *graph, K, 1, DDR_size);

	/* print info about decomposition */
	decomposition.print(coutMaster);

	if(graph_save){
		/* save decoposed graph to see if space (graph) decomposition is working */
		oss << "results/" << image_out << "_graph.vtk";
		graph->saveVTK(oss.str());
		oss.str("");
	}

/* 3.) prepare time-series data */
	coutMaster << "--- PREPARING DATA ---\n";
	
	/* load data from file and store it subject to decomposition */
	ImageData<PetscVector> mydata(decomposition, image_filename, width, height);
	
	/* print information about loaded data */
	mydata.print(coutMaster);

	/* print statistics */
	if(printstats) mydata.printstats(coutMaster);

	/* cut data */
	if(cutdata) mydata.cutdata(0,1);

	/* scale data */
	if(scaledata) mydata.scaledata(-1,1,0,1);

/* 4.) prepare model */
	coutMaster << "--- PREPARING MODEL ---\n";

	/* prepare model on the top of given data */
	GraphH1FEMModel<PetscVector> mymodel(mydata, epssqr_list[0]);

	/* print info about model */
	mymodel.print(coutMaster,coutAll);

/* 5.) prepare time-series solver */
	coutMaster << "--- PREPARING SOLVER ---\n";

	/* prepare time-series solver */
	TSSolver<PetscVector> mysolver(mydata, annealing);

	/* print info about solver */
	mysolver.print(coutMaster,coutAll);

	/* set solution if obtained from console */
	if(given_Theta)	mysolver.set_solution_theta(Theta_solution);
	
/* 6.) solve the problem with initial epssqr */
	coutMaster << "--- SOLVING THE PROBLEM with epssqr = " << epssqr_list[0] << " ---\n";
	mysolver.solve();

	/* cut gamma */
	if(cutgamma) mydata.cutgamma();

	/* unscale data before save */
	if(scaledata) mydata.scaledata(0,1,-1,1);

	coutMaster << "--- SAVING OUTPUT ---\n";
	oss << image_out << "_depth0" << "_epssqr" << epssqr_list[0];
	mydata.saveImage(oss.str(),true);
	oss.str("");

	/* write short output */
	if(shortinfo_write_or_not){
		/* add provided strings from console parameters */
		oss_short_output_header << shortinfo_header;
		oss_short_output_values << shortinfo_values;

		/* add info about the problem */
		oss_short_output_header << "width,height,K,depth,epssqr,";
		oss_short_output_values << width << "," << height << "," << K << ",0,0.0,"; 

		/* append Theta solution */
		for(int k=0; k<K; k++) oss_short_output_header << "Theta" << k << ",";
		oss_short_output_values << mydata.print_thetavector(); 

		/* print info from solver */
		mysolver.printshort(oss_short_output_header, oss_short_output_values);

		/* append end of line */
		oss_short_output_header << "\n";
		oss_short_output_values << "\n";

		/* write to shortinfo file */
		shortinfo.write(oss_short_output_header.str());
		shortinfo.write(oss_short_output_values.str());
			
		/* clear streams for next writing */
		oss_short_output_header.str("");
		oss_short_output_values.str("");
	
	}


/* 7.) solve the problems with other epssqr */
	for(int depth = 1; depth < epssqr_list.size();depth++){
		/* set new epssqr */
		mymodel.set_epssqr(epssqr_list[depth]);

		/* decrease the number of annealing steps in TSSolver to 1 */
//		mysolver.set_annealing(1);

		/* scale data before computation */
		if(scaledata) mydata.scaledata(-1,1,0,1);

		coutMaster << "--- SOLVING THE PROBLEM with epssqr = " << epssqr_list[depth] << " ---\n";
		mysolver.solve();

		/* cut gamma */
		if(cutgamma) mydata.cutgamma();

		/* unscale data before export */
		if(scaledata) mydata.scaledata(0,1,-1,1);

		coutMaster << "--- SAVING OUTPUT ---\n";
		oss << image_out << "_depth" << depth << "_epssqr" << epssqr_list[depth];
		mydata.saveImage(oss.str(),false);
		oss.str("");
		
		/* write short output */
		if(shortinfo_write_or_not){
			/* add provided strings from console parameters */
			oss_short_output_values << shortinfo_values << width << "," << height << "," << K << "," << depth << "," << epssqr_list[depth] << ",";

			/* append Theta solution */
			oss_short_output_values << mydata.print_thetavector(); 

			/* append data from solver */
			mysolver.printshort(oss_short_output_header, oss_short_output_values);

			/* append end of line */
			oss_short_output_values << "\n";

			/* write data */
			shortinfo.write(oss_short_output_values.str());

			/* clear streams for next time */
			oss_short_output_values.str("");
		}
	}

	/* print solution */
	coutMaster << "--- THETA SOLUTION ---\n";
	mydata.print_thetavector(coutMaster);

	/* print timers */
	coutMaster << "--- TIMERS INFO ---\n";
	mysolver.printtimer(coutMaster);

	/* print short info */
	coutMaster << "--- FINAL SOLVER INFO ---\n";
	mysolver.printstatus(coutMaster);

	/* say bye */	
	coutMaster << "- end program\n";

	logging.end();
	Finalize();

	return 0;
}

