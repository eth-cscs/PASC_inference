/** @file test_image.cpp
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

#define DEFAULT_EPSSQR 1
#define DEFAULT_K 2
#define DEFAULT_FEM_TYPE 2
#define DEFAULT_FEM_REDUCE 1.0
#define DEFAULT_XDIM 1
#define DEFAULT_WIDTH 250
#define DEFAULT_HEIGHT 150
#define DEFAULT_GRAPH_SAVE false
#define DEFAULT_CUTGAMMA false
#define DEFAULT_SCALEDATA false
#define DEFAULT_CUTDATA true
#define DEFAULT_PRINTSTATS false
#define DEFAULT_PRINTINFO false
#define DEFAULT_ANNEALING 1
#define DEFAULT_IMAGE_IN "data/test_image/usi_text/usi_250_150_02.bin"
#define DEFAULT_IMAGE_OUT "usi_test_250_150_02"
#define DEFAULT_IMAGE_SOLUTION "data/test_image/usi_text/usi_250_150_solution.bin"
#define DEFAULT_SHORTINFO true
#define DEFAULT_SHORTINFO_FILENAME "shortinfo/usi_250_150_02.txt"
#define DEFAULT_SAVEALL false
#define DEFAULT_SAVERESULT true

using namespace pascinference;

int main( int argc, char *argv[] )
{
	/* add local program options */
	boost::program_options::options_description opt_problem("PROBLEM EXAMPLE", consoleArg.get_console_nmb_cols());
	opt_problem.add_options()
		("test_K", boost::program_options::value<int>(), "number of clusters [int]")
		("test_fem_type", boost::program_options::value<int>(), "type of used FEM to reduce problem [3=FEM2D_SUM/4=FEM2D_HAT]")
		("test_fem_reduce", boost::program_options::value<double>(), "parameter of the reduction of FEM nodes [int,-1=false]")
		("test_filename_in", boost::program_options::value< std::string >(), "name of input file with image data (vector in PETSc format) [string]")
		("test_filename_out", boost::program_options::value< std::string >(), "name of output file with image data (vector in PETSc format) [string]")
		("test_filename_solution", boost::program_options::value< std::string >(), "name of input file with original image data without noise (vector in PETSc format) [string]")
		("test_width", boost::program_options::value<int>(), "width of image [int]")
		("test_height", boost::program_options::value<int>(), "height of image [int]")
		("test_xdim", boost::program_options::value<int>(), "number of values in every pixel [1=greyscale, 3=rgb]")
		("test_graph_save", boost::program_options::value<bool>(), "save VTK with graph or not [bool]")
		("test_epssqr", boost::program_options::value<std::vector<double> >()->multitoken(), "penalty parameters [double]")
		("test_annealing", boost::program_options::value<int>(), "number of annealing steps [int]")
		("test_cutgamma", boost::program_options::value<bool>(), "cut gamma to set {0;1} [bool]")
		("test_scaledata", boost::program_options::value<bool>(), "scale to interval {-1,1} [bool]")
		("test_cutdata", boost::program_options::value<bool>(), "cut data to interval {0,1} [bool]")
		("test_printstats", boost::program_options::value<bool>(), "print basic statistics of data [bool]")
		("test_printinfo", boost::program_options::value<bool>(), "print informations about created objects [bool]")
		("test_saveall", boost::program_options::value<bool>(), "save results for all epssqr, not only for the best one [bool]")
		("test_saveresult", boost::program_options::value<bool>(), "save the solution [bool]")
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
		/* list is not given, add some value */
		epssqr_list.push_back(DEFAULT_EPSSQR);
	}

	int K, annealing, width, height, xdim, fem_type; 
	double fem_reduce;
	bool cutgamma, scaledata, cutdata, printstats, printinfo, shortinfo_write_or_not, graph_save, saveall, saveresult;

	std::string filename_in;
	std::string filename_out;
	std::string filename_solution;
	std::string shortinfo_filename;
	std::string shortinfo_header;
	std::string shortinfo_values;

	consoleArg.set_option_value("test_K", &K, DEFAULT_K);
	consoleArg.set_option_value("test_fem_type", &fem_type, DEFAULT_FEM_TYPE);
	consoleArg.set_option_value("test_fem_reduce", &fem_reduce, DEFAULT_FEM_REDUCE);
	consoleArg.set_option_value("test_width", &width, DEFAULT_WIDTH);
	consoleArg.set_option_value("test_height", &height, DEFAULT_HEIGHT);
	consoleArg.set_option_value("test_xdim", &xdim, DEFAULT_XDIM);
	consoleArg.set_option_value("test_graph_save", &graph_save, DEFAULT_GRAPH_SAVE);
	consoleArg.set_option_value("test_cutgamma", &cutgamma, DEFAULT_CUTGAMMA);
	consoleArg.set_option_value("test_scaledata", &scaledata, DEFAULT_SCALEDATA);
	consoleArg.set_option_value("test_cutdata", &cutdata, DEFAULT_CUTDATA);
	consoleArg.set_option_value("test_printstats", &printstats, DEFAULT_PRINTSTATS);
	consoleArg.set_option_value("test_printinfo", &printinfo, DEFAULT_PRINTINFO);
	consoleArg.set_option_value("test_annealing", &annealing, DEFAULT_ANNEALING);
	consoleArg.set_option_value("test_filename_in", &filename_in, DEFAULT_IMAGE_IN);
	consoleArg.set_option_value("test_filename_out", &filename_out, DEFAULT_IMAGE_OUT);
	consoleArg.set_option_value("test_saveall", &saveall, DEFAULT_SAVEALL);
	consoleArg.set_option_value("test_saveresult", &saveresult, DEFAULT_SAVERESULT);
	consoleArg.set_option_value("test_shortinfo", &shortinfo_write_or_not, DEFAULT_SHORTINFO);
	consoleArg.set_option_value("test_shortinfo_header", &shortinfo_header, "");
	consoleArg.set_option_value("test_shortinfo_values", &shortinfo_values, "");
	consoleArg.set_option_value("test_shortinfo_filename", &shortinfo_filename, DEFAULT_SHORTINFO_FILENAME);

	/* maybe solution is given */
	bool given_solution;
	if(!consoleArg.set_option_value("test_filename_solution", &filename_solution)){
		given_solution=false;

		/* maybe we run program with default values */
		if(filename_in == DEFAULT_IMAGE_IN){
			filename_solution = DEFAULT_IMAGE_SOLUTION;
			given_solution=true;
		}
	} else {
		given_solution=true;
	}
	
	/* maybe theta is given in console parameters */
	bool given_Theta;
	std::vector<double> Theta_list;
	double Theta_solution[K*xdim];
	if(consoleArg.set_option_value("test_Theta", &Theta_list)){
		given_Theta = true;
		
		/* control number of provided Theta */
		if(Theta_list.size() != K*xdim){
			coutMaster << "number of provided Theta solutions is different then number of clusters (times xdim)!" << std::endl;
			return 0;
		}

		/* store solution in array */
		for(int k=0;k < K*xdim;k++){
			Theta_solution[k] = Theta_list[k];
		}
	} else {
		given_Theta = false;
	}	

	/* set decomposition in space */
	int DDR_size = GlobalManager.get_size();

	coutMaster << "- PROBLEM INFO ----------------------------" << std::endl;
#ifdef USE_CUDA
	coutMaster << " computing on GPU" << std::endl;
#else
	coutMaster << " computing on CPU" << std::endl;
#endif
	coutMaster << " DDR_size                    = " << std::setw(50) << DDR_size << " (decomposition in space)" << std::endl;
	coutMaster << " test_filename_in            = " << std::setw(50) << filename_in << " (name of input file with image data)" << std::endl;
	coutMaster << " test_filename_out           = " << std::setw(50) << filename_out << " (part of name of output file)" << std::endl;
	if(given_solution){
		coutMaster << " test_filename_solution      = " << std::setw(50) << filename_solution << " (name of input file with original image data without noise)" << std::endl;
	} else {
		coutMaster << " test_filename_solution      = " << std::setw(50) << "NO" << " (name of input file with original image data without noise)" << std::endl;
	}
	coutMaster << " test_width                  = " << std::setw(50) << width << " (width of image)" << std::endl;
	coutMaster << " test_height                 = " << std::setw(50) << height << " (height of image)" << std::endl;
	coutMaster << " test_xdim                   = " << std::setw(50) << xdim << " (number of values in every pixel [1=greyscale, 3=rgb])" << std::endl;
	coutMaster << " test_K                      = " << std::setw(50) << K << " (number of clusters)" << std::endl;
	if(given_Theta){
		coutMaster << " test_Theta                  = " << std::setw(50) << print_array(Theta_solution,K) << " (given solution Theta)" << std::endl;
	} else {
		coutMaster << " test_Theta                  = " << std::setw(50) << "NO" << " (given solution Theta)" << std::endl;
	}
	coutMaster << " test_fem_type               = " << std::setw(50) << fem_type << " (type of used FEM to reduce problem [3=FEM2D_SUM/4=FEM2D_HAT])" << std::endl;
	coutMaster << " test_fem_reduce             = " << std::setw(50) << fem_reduce << " (parameter of the reduction of FEM node)" << std::endl;
	coutMaster << " test_graph_save             = " << std::setw(50) << print_bool(graph_save) << " (save VTK with graph or not)" << std::endl;
	coutMaster << " test_saveresult             = " << std::setw(50) << print_bool(saveresult) << " (save reconstructed image)" << std::endl;
	coutMaster << " test_saveall                = " << std::setw(50) << print_bool(saveall) << " (save results for all epssqr, not only for the best one)" << std::endl;
	coutMaster << " test_epssqr                 = " << std::setw(50) << print_vector(epssqr_list) << " (penalty)" << std::endl;
	coutMaster << " test_annealing              = " << std::setw(50) << annealing << " (number of annealing steps)" << std::endl;
	coutMaster << " test_cutgamma               = " << std::setw(50) << print_bool(cutgamma) << " (cut gamma to {0;1})" << std::endl;
	coutMaster << " test_cutdata                = " << std::setw(50) << print_bool(cutdata) << " (cut data to {0,1})" << std::endl;
	coutMaster << " test_scaledata              = " << std::setw(50) << print_bool(scaledata) << " (scale data to {-1,1})" << std::endl;
	coutMaster << " test_printstats             = " << std::setw(50) << print_bool(printstats) << " (print basic statistics of data)" << std::endl;
	coutMaster << " test_printinfo              = " << std::setw(50) << print_bool(printinfo) << " (print informations about created objects)" << std::endl;
	coutMaster << " test_shortinfo              = " << std::setw(50) << print_bool(shortinfo_write_or_not) << " (save shortinfo file after computation)" << std::endl;
	coutMaster << " test_shortinfo_header       = " << std::setw(50) << shortinfo_header << " (additional header in shortinfo)" << std::endl;
	coutMaster << " test_shortinfo_values       = " << std::setw(50) << shortinfo_values << " (additional values in shortinfo)" << std::endl;
	coutMaster << " test_shortinfo_filename     = " << std::setw(50) << shortinfo_filename << " (name of shortinfo file)" << std::endl;
	coutMaster << "-------------------------------------------" << std::endl;
	coutMaster << std::endl;

	/* start logging */
	std::ostringstream oss;
	oss << "log/" << filename_out << ".txt";
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

	/* this is not general graph, it is "only" 2D grid */
	BGMGraphGrid2D<PetscVector> *graph;
	graph = new BGMGraphGrid2D<PetscVector>(width, height);
	graph->process_grid();

	/* print basic info about graph */
	if(printinfo) graph->print(coutMaster);
	
/* 2.) prepare decomposition */
	coutMaster << "--- COMPUTING DECOMPOSITION ---" << std::endl;

	/* prepare decomposition based on graph, in this case T=1 and DDT_size=1 */
	Decomposition<PetscVector> decomposition(1, *graph, K, xdim, DDR_size);

	/* print info about decomposition */
	if(printinfo) decomposition.print(coutMaster);

	if(graph_save){
		/* save decoposed graph to see if space (graph) decomposition is working */
		oss << "results/" << filename_out << "_graph.vtk";
		graph->saveVTK(oss.str());
		oss.str("");
	}

	//TODO: temp
	graph->print_content(coutMaster);


/* 3.) prepare time-series data */
	coutMaster << "--- PREPARING DATA ---" << std::endl;
	
	/* load data from file and store it subject to decomposition */
	ImageData<PetscVector> mydata(decomposition, filename_in, width, height);
	
	/* print information about loaded data */
	if(printinfo) mydata.print(coutMaster);

	/* print statistics */
	if(printstats) mydata.printstats(coutMaster);

/* 4.) prepare and load solution */
	Vec solution_Vec;
	Vec solution_Vec_preload;
	GeneralVector<PetscVector> solution(solution_Vec);
	if(given_solution){
		TRYCXX( VecDuplicate(mydata.get_datavector()->get_vector(),&solution_Vec) );
		TRYCXX( VecDuplicate(mydata.get_datavector()->get_vector(),&solution_Vec_preload) );

		solution.load_global(filename_solution);
		decomposition.permute_bTR_to_dTRb(solution.get_vector(), solution_Vec_preload, decomposition.get_xdim(), false);

		TRYCXX( VecCopy(solution_Vec_preload, solution.get_vector()));
		TRYCXX( VecDestroy(&solution_Vec_preload) );
	}

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
	GraphH1FEMModel<PetscVector> mymodel(mydata, epssqr_list[0], fem);

	/* print info about model */
	if(printinfo) mymodel.print(coutMaster,coutAll);

/* 5.) prepare time-series solver */
	coutMaster << "--- PREPARING SOLVER ---" << std::endl;

	/* prepare time-series solver */
	TSSolver<PetscVector> mysolver(mydata, annealing);

	/* print info about solver */
	if(printinfo) mysolver.print(coutMaster,coutAll);

	/* set solution if obtained from console */
	if(given_Theta)	mysolver.set_solution_theta(Theta_solution);
	
/* 6.) solve the problem with epssqrs and remember best solution */
	double epssqr;
	double epssqr_best = -1;
	double abserr; /* actual error */
	double abserr_best = std::numeric_limits<double>::max(); /* the error of best solution */
	double L; /* actual value of objective function */
	double L_best = std::numeric_limits<double>::max(); /* best value of L from all possible epssqr */

	Vec gammavector_best_Vec; /* here we store solution with best abserr value */
	TRYCXX( VecDuplicate(mydata.get_gammavector()->get_vector(),&gammavector_best_Vec) );

	Vec thetavector_best_Vec; /* here we store solution with best abserr value */
	TRYCXX( VecDuplicate(mydata.get_thetavector()->get_vector(),&thetavector_best_Vec) );

	double tessst;

/* 7.) solve the problems with other epssqr */
	for(int depth = 0; depth < epssqr_list.size();depth++){
		epssqr = epssqr_list[depth];
		coutMaster << "--- SOLVING THE PROBLEM with epssqr = " << epssqr << " ---" << std::endl;

		/* set new epssqr */
		mymodel.set_epssqr(epssqr_list[depth]);

		/* cut data */
		if(cutdata) mydata.cutdata(0,1);

		/* scale data */
		if(scaledata){
			mydata.scaledata(-1,1,0,1);
		}

		/* !!! solve the problem */
		mysolver.solve();

		/* cut gamma */
		if(cutgamma) mydata.cutgamma();

		/* unscale data before save */
		if(scaledata){
			mydata.scaledata(0,1,-1,1);
		}

		std::streamsize ss = coutMaster.precision();
		coutMaster << std::setprecision(17);
		/* compute absolute error of computed solution */
		if(given_solution){
			abserr = mydata.compute_abserr_reconstructed(solution);
		} else {
			abserr = -1.0;
		}
		coutMaster << " - abserr = " << abserr << std::endl;

		/* compute value of objective function */
		L = mysolver.get_L();
		coutMaster << " - L = " << L << std::endl;
		coutMaster << std::setprecision(ss);

		/* store obtained solution */
		if(saveall && saveresult){
			coutMaster << "--- SAVING OUTPUT ---" << std::endl;
			oss << filename_out << "_epssqr" << epssqr;
			mydata.saveImage(oss.str(),false);
			oss.str("");
		}
	
		/* store short info */
		if(shortinfo_write_or_not){
			/* add provided strings from console parameters and info about the problem */
			if(depth==0) oss_short_output_header << shortinfo_header << "width,height,K,depth,epssqr,abserr,";
			oss_short_output_values << shortinfo_values << width << "," << height << "," << K << "," << depth << "," << epssqr_list[depth] << "," << abserr << ",";
			
			/* append Theta solution */
			if(depth==0) for(int k=0; k<K; k++) oss_short_output_header << "Theta" << k << ",";
			oss_short_output_values << mydata.print_thetavector(); 

			/* print info from solver */
			mysolver.printshort(oss_short_output_header, oss_short_output_values);

			/* append end of line */
			if(depth==0) oss_short_output_header << "\n";
			oss_short_output_values << "\n";

			/* write to shortinfo file */
			if(depth==0) shortinfo.write(oss_short_output_header.str());
			shortinfo.write(oss_short_output_values.str());
			
			/* clear streams for next writing */
			oss_short_output_header.str("");
			oss_short_output_values.str("");
		}
	
		/* if this solution is better then previous, then store it */
		if(abserr < abserr_best || depth == 0){
			L_best = L;
			abserr_best = abserr;
			epssqr_best = epssqr;
			TRYCXX(VecCopy(mydata.get_gammavector()->get_vector(),gammavector_best_Vec));
			TRYCXX(VecCopy(mydata.get_thetavector()->get_vector(),thetavector_best_Vec));
		}
	
	}

	/* set best computed solution back to data */
	TRYCXX(VecCopy(gammavector_best_Vec,mydata.get_gammavector()->get_vector()));
	TRYCXX(VecCopy(thetavector_best_Vec, mydata.get_thetavector()->get_vector()));

/* 8.) store best solution */
	if(saveresult && !saveall){
		coutMaster << "--- SAVING OUTPUT ---" << std::endl;
		coutMaster << " - with best epssqr = " << epssqr_best << std::endl;
		oss << filename_out;
		mydata.saveImage(oss.str(),false);
		oss.str("");
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

