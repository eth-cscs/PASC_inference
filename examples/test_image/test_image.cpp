/** @file test_image.cu
 *  @brief test the kmeans problem solver on image processing problem
 *
 *  Load image file, graph file and solve the problem.
 *
 *  @author Lukas Pospisil
 */

#include "pascinference.h"
#include <vector>

#ifndef USE_PETSC
 #error 'This example is for PETSC'
#endif
 
using namespace pascinference;

int main( int argc, char *argv[] )
{
	/* add local program options */
	boost::program_options::options_description opt_problem("PROBLEM EXAMPLE", consoleArg.get_console_nmb_cols());
	opt_problem.add_options()
		("test_K", boost::program_options::value<int>(), "number of clusters [int]")
		("test_fem_type", boost::program_options::value<int>(), "type of used FEM to reduce problem [3=FEM2D_SUM/4=FEM2D_HAT]")
		("test_fem_reduce", boost::program_options::value<double>(), "parameter of the reduction of FEM nodes [int,-1=false]")
		("test_image_filename", boost::program_options::value< std::string >(), "name of input file with image data (vector in PETSc format) [string]")
		("test_image_out", boost::program_options::value< std::string >(), "name of output file with image data (vector in PETSc format) [string]")
		("test_width", boost::program_options::value<int>(), "width of image [int]")
		("test_height", boost::program_options::value<int>(), "height of image [int]")
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
	if(!Initialize<PetscVector>(argc, argv)){
		return 0;
	} 

	std::vector<double> epssqr_list;
	if(consoleArg.set_option_value("test_epssqr", &epssqr_list)){
		/* sort list */
		std::sort(epssqr_list.begin(), epssqr_list.end(), std::less<double>());
		
	} else {
		std::cout << "test_epssqr has to be set! Call application with parameter -h to see all parameters" << std::endl;
		return 0;
	}

	int K, annealing, width, height, fem_type; 
	double fem_reduce;
	bool cutgamma, scaledata, cutdata, printstats, shortinfo_write_or_not, graph_save;

	std::string image_filename;
	std::string image_out;
	std::string shortinfo_filename;
	std::string shortinfo_header;
	std::string shortinfo_values;

	consoleArg.set_option_value("test_K", &K, 2);
	consoleArg.set_option_value("test_fem_type", &fem_type, 3);
	consoleArg.set_option_value("test_fem_reduce", &fem_reduce, 1.0);
	consoleArg.set_option_value("test_width", &width, 250);
	consoleArg.set_option_value("test_height", &height, 150);
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
	coutMaster << " DDR_size                = " << std::setw(50) << DDR_size << " (decomposition in space)" << std::endl;
	coutMaster << " test_image_filename     = " << std::setw(50) << image_filename << " (name of input file with image data)" << std::endl;
	coutMaster << " test_image_out          = " << std::setw(50) << image_out << " (part of name of output file)" << std::endl;
	coutMaster << " test_width              = " << std::setw(50) << width << " (width of image)" << std::endl;
	coutMaster << " test_height             = " << std::setw(50) << height << " (height of image)" << std::endl;
	coutMaster << " test_K                  = " << std::setw(50) << K << " (number of clusters)" << std::endl;
	if(given_Theta){
		coutMaster << " test_Theta              = " << std::setw(50) << print_array(Theta_solution,K) << std::endl;
	}
	coutMaster << " test_fem_type           = " << std::setw(50) << fem_type << " (type of used FEM to reduce problem [3=FEM2D_SUM/4=FEM2D_HAT])" << std::endl;
	coutMaster << " test_fem_reduce         = " << std::setw(50) << fem_reduce << " (parameter of the reduction of FEM node)" << std::endl;
	coutMaster << " test_graph_save         = " << std::setw(50) << graph_save << " (save VTK with graph or not)" << std::endl;
	coutMaster << " test_epssqr             = " << std::setw(50) << print_vector(epssqr_list) << " (penalty)" << std::endl;
	coutMaster << " test_annealing          = " << std::setw(50) << annealing << " (number of annealing steps)" << std::endl;
	coutMaster << " test_cutgamma           = " << std::setw(50) << cutgamma << " (cut gamma to {0;1})" << std::endl;
	coutMaster << " test_cutdata            = " << std::setw(50) << cutdata << " (cut data to {0,1})" << std::endl;
	coutMaster << " test_scaledata          = " << std::setw(50) << scaledata << " (scale data to {-1,1})" << std::endl;
	coutMaster << " test_shortinfo          = " << std::setw(50) << shortinfo_write_or_not << " (save shortinfo file after computation)" << std::endl;
	coutMaster << " test_shortinfo_header   = " << std::setw(50) << shortinfo_header << " (additional header in shortinfo)" << std::endl;
	coutMaster << " test_shortinfo_values   = " << std::setw(50) << shortinfo_values << " (additional values in shortinfo)" << std::endl;
	coutMaster << " test_shortinfo_filename = " << std::setw(50) << shortinfo_filename << " (name of shortinfo file)" << std::endl;
	coutMaster << "-------------------------------------------" << std::endl;
	coutMaster << std::endl;

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
	coutMaster << "- start program" << std::endl;

/* 1.) prepare graph of image */
	coutMaster << "--- PREPARING GRAPH ---" << std::endl;

	BGMGraphGrid2D<PetscVector> *graph;
	/* this is not general graph, it is "only" 2D grid */
	graph = new BGMGraphGrid2D<PetscVector>(width, height);
	graph->process_grid();

	/* print basic info about graph */
	graph->print(coutMaster);
	
/* 2.) prepare decomposition */
	coutMaster << "--- COMPUTING DECOMPOSITION ---" << std::endl;

	/* prepare decomposition based on graph, in this case T=1 and DDT_size=1 */
	Decomposition<PetscVector> decomposition(1, *graph, K, 1, DDR_size);

	/* print info about decomposition */
	decomposition.print(coutMaster);

	if(graph_save){
		/* save decoposed graph to see if space (graph) decomposition is working */
		oss << "results/" << image_out << "_graph.vtk";
		graph->saveVTK(oss.str());
		oss.str("");
	}

/* 3.) prepare time-series data */
	coutMaster << "--- PREPARING DATA ---" << std::endl;
	
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
	coutMaster << "--- PREPARING MODEL ---" << std::endl;

	/* prepare FEM reduction */
	Fem<PetscVector> *fem;
	if(fem_type == 3){
		fem = new Fem2D<PetscVector>(fem_reduce);
	}
	if(fem_type == 4){
//		fem = new Fem2DHat<PetscVector>(fem_reduce);
	}

	/* prepare model on the top of given data */
//	GraphH1FEMModel<PetscVector> mymodel(mydata, epssqr_list[0]);
	GraphH1FEMModel<PetscVector> mymodel(mydata, epssqr_list[0], fem);

	/* print info about model */
	mymodel.print(coutMaster,coutAll);

/* 5.) prepare time-series solver */
	coutMaster << "--- PREPARING SOLVER ---" << std::endl;

	/* prepare time-series solver */
	TSSolver<PetscVector> mysolver(mydata, annealing);

	/* print info about solver */
	mysolver.print(coutMaster,coutAll);

	/* set solution if obtained from console */
	if(given_Theta)	mysolver.set_solution_theta(Theta_solution);
	
/* 6.) solve the problem with initial epssqr */
	coutMaster << "--- SOLVING THE PROBLEM with epssqr = " << epssqr_list[0] << " ---" << std::endl;
	mysolver.solve();

	/* cut gamma */
	if(cutgamma) mydata.cutgamma();

	/* unscale data before save */
	if(scaledata) mydata.scaledata(0,1,-1,1);

	coutMaster << "--- SAVING OUTPUT ---" << std::endl;
	oss << image_out << "_epssqr" << epssqr_list[0];
	mydata.saveImage(oss.str(),true);
	oss.str("");

	/* write short output */
	if(shortinfo_write_or_not){
		/* add provided strings from console parameters */
		oss_short_output_header << shortinfo_header;
		oss_short_output_values << shortinfo_values;

		/* add info about the problem */
		oss_short_output_header << "width,height,K,depth,epssqr,";
		oss_short_output_values << width << "," << height << "," << K << ",0," << epssqr_list[0] << ","; 

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

		coutMaster << "--- SOLVING THE PROBLEM with epssqr = " << epssqr_list[depth] << " ---" << std::endl;
		mysolver.solve();

		/* cut gamma */
		if(cutgamma) mydata.cutgamma();

		/* unscale data before export */
		if(scaledata) mydata.scaledata(0,1,-1,1);

		coutMaster << "--- SAVING OUTPUT ---" << std::endl;
		oss << image_out << "_epssqr" << epssqr_list[depth];
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
	coutMaster << "--- THETA SOLUTION ---" << std::endl;
	mydata.print_thetavector(coutMaster);

	/* print timers */
	coutMaster << "--- TIMERS INFO ---" << std::endl;
	mysolver.printtimer(coutMaster);

	/* print short info */
	coutMaster << "--- FINAL SOLVER INFO ---" << std::endl;
	mysolver.printstatus(coutMaster);

	/* say bye */	
	coutMaster << "- end program" << std::endl;

	logging.end();
	Finalize<PetscVector>();

	return 0;
}

