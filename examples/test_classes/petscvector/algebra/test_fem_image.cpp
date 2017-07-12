/** @file test_fem_image.cpp
 *  @brief test the fem reduction/prolongation on image problems
 *
 *
 *  @author Lukas Pospisil
 */

#include "pascinference.h"
#include <vector>

#ifndef USE_PETSC
 #error 'This example is for PETSC'
#endif

#define DEFAULT_FEM_TYPE 0
#define DEFAULT_FEM_REDUCE 1.0
#define DEFAULT_XDIM 1
#define DEFAULT_WIDTH 250
#define DEFAULT_HEIGHT 150
#define DEFAULT_GRAPH_SAVE false
#define DEFAULT_PRINTINFO false
#define DEFAULT_FILENAME_IN "data/test_image/usi_text/usi_250_150_02.bin"
#define DEFAULT_FILENAME_OUT "usi_test_250_150_02"
#define DEFAULT_IMAGE_SOLUTION "data/test_image/usi_text/usi_250_150_solution.bin"

using namespace pascinference;

int main( int argc, char *argv[] )
{
	/* add local program options */
	boost::program_options::options_description opt_problem("TEST FEM IMAGE", consoleArg.get_console_nmb_cols());
	opt_problem.add_options()
		("test_fem_type", boost::program_options::value<int>(), "type of used FEM to reduce problem [3=FEM2D_SUM/4=FEM2D_HAT]")
		("test_fem_reduce", boost::program_options::value<double>(), "parameter of the reduction of FEM nodes [int,-1=false]")
		("test_filename_in", boost::program_options::value< std::string >(), "name of input file with image data (vector in PETSc format) [string]")
		("test_filename_out", boost::program_options::value< std::string >(), "name of output file with image data (vector in PETSc format) [string]")
		("test_width", boost::program_options::value<int>(), "width of image [int]")
		("test_height", boost::program_options::value<int>(), "height of image [int]")
		("test_xdim", boost::program_options::value<int>(), "number of values in every pixel [1=greyscale, 3=rgb]")
		("test_graph_save", boost::program_options::value<bool>(), "save VTK with graph or not [bool]")
		("test_printinfo", boost::program_options::value<bool>(), "print informations about created objects [bool]");

	consoleArg.get_description()->add(opt_problem);

	/* call initialize */
	if(!Initialize<PetscVector>(argc, argv)){
		return 0;
	}

	int width, height, xdim, fem_type;
	double fem_reduce;
	bool printinfo, graph_save;

	std::string filename_in;
	std::string filename_out;

	consoleArg.set_option_value("test_fem_type", &fem_type, DEFAULT_FEM_TYPE);
	consoleArg.set_option_value("test_fem_reduce", &fem_reduce, DEFAULT_FEM_REDUCE);
	consoleArg.set_option_value("test_width", &width, DEFAULT_WIDTH);
	consoleArg.set_option_value("test_height", &height, DEFAULT_HEIGHT);
	consoleArg.set_option_value("test_xdim", &xdim, DEFAULT_XDIM);
	consoleArg.set_option_value("test_graph_save", &graph_save, DEFAULT_GRAPH_SAVE);
	consoleArg.set_option_value("test_printinfo", &printinfo, DEFAULT_PRINTINFO);
	consoleArg.set_option_value("test_filename_in", &filename_in, DEFAULT_FILENAME_IN);
	consoleArg.set_option_value("test_filename_out", &filename_out, DEFAULT_FILENAME_OUT);

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

	coutMaster << " test_width                  = " << std::setw(50) << width << " (width of image)" << std::endl;
	coutMaster << " test_height                 = " << std::setw(50) << height << " (height of image)" << std::endl;
	coutMaster << " test_xdim                   = " << std::setw(50) << xdim << " (number of values in every pixel [1=greyscale, 3=rgb])" << std::endl;

	coutMaster << " test_fem_type               = " << std::setw(50) << fem_type << " (type of used FEM to reduce problem [0=FEM2D_SUM/1=FEM2D_HAT])" << std::endl;
	coutMaster << " test_fem_reduce             = " << std::setw(50) << fem_reduce << " (parameter of the reduction of FEM node)" << std::endl;
	coutMaster << " test_graph_save             = " << std::setw(50) << print_bool(graph_save) << " (save VTK with graph or not)" << std::endl;
	coutMaster << " test_printinfo              = " << std::setw(50) << print_bool(printinfo) << " (print informations about created objects)" << std::endl;
	coutMaster << "-------------------------------------------" << std::endl;
	coutMaster << std::endl;

	/* start logging */
	std::ostringstream oss;
	oss << "log/" << filename_out << ".txt";
	logging.begin(oss.str());
	oss.str("");

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
	int K=xdim; /* I am not using Gamma in this test at all, here is a trick: use K instead of xdim */
	Decomposition<PetscVector> decomposition_orig(1, *graph, K, xdim, DDR_size);

	/* print info about decomposition */
	if(printinfo) decomposition_orig.print(coutMaster);

	if(graph_save){
		/* save decoposed graph to see if space (graph) decomposition is working */
		oss << "results/" << filename_out << "_graph_orig.vtk";
		graph->saveVTK(oss.str());
		oss.str("");
	}

/* 3.) prepare time-series data */
	coutMaster << "--- PREPARING DATA ---" << std::endl;

	/* load data from file and store it subject to decomposition */
	ImageData<PetscVector> mydata(decomposition_orig, filename_in, width, height);

	/* print information about loaded data */
	if(printinfo) mydata.print(coutMaster);

	/* prepare FEM reduction */
	Fem<PetscVector> *fem;
	if(fem_type == 0){
		fem = new Fem2DSum<PetscVector>(fem_reduce);
	}
	if(fem_type == 1){
//		fem = new Fem2DHat<PetscVector>(fem_reduce);
	}

	/* set original decomposition to fem */
    fem->set_decomposition_original(&decomposition_orig);

    /* compute reduced decomposition */
    fem->compute_decomposition_reduced();

    /* get width and height of reduced grid */
    int width_reduced, height_reduced;
    if(fem_type == 0){
        width_reduced = ((Fem2DSum<PetscVector>*)fem)->get_grid_reduced()->get_width();
        height_reduced = ((Fem2DSum<PetscVector>*)fem)->get_grid_reduced()->get_height();
    }

	if(printinfo) fem->print(coutMaster, coutAll);

    /* save reduced graph */
	if(graph_save){
		/* save decoposed reduced graph to see if space (graph) decomposition is working */
		oss << "results/" << filename_out << "_graph_reduced.vtk";
		fem->get_decomposition_reduced()->get_graph()->saveVTK(oss.str());
		oss.str("");

        /* save bounding boxes */
        if(fem_type == 0){
            ((Fem2DSum<PetscVector>*)fem)->saveVTK_bounding_box("results/bb1.vtk","results/bb2.vtk");
        }
	}

    /* prepare stuff for reduced data */
	ImageData<PetscVector> mydata_reduced(*(fem->get_decomposition_reduced()), width_reduced, height_reduced);

    /* reduce gamma, however gamma is data now (that's the reason why K=xdim) */
    fem->reduce_gamma(mydata.get_datavector(), mydata_reduced.get_datavector());

    /* save reduced image */
    mydata_reduced.saveImage_datavector(filename_out);

	/* say bye */
	coutMaster << "- end program" << std::endl;

	logging.end();
	Finalize<PetscVector>();

	return 0;
}

