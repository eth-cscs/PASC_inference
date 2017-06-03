/** @file util_print_vec.cpp
 *  @brief print the content of given PETSc vector
 *
 *  @author Lukas Pospisil
 */

#include "pascinference.h"
#include <vector>

#ifndef USE_PETSC
 #error 'This example is for PETSC'
#endif
 
#ifndef USE_DLIB
 #error 'This example is for DLIB'
#endif

#define ENABLE_ASSERTS
#define DLIB_JPEG_SUPPORT

#include <dlib/image_io.h>
#include <dlib/image_transforms.h>

using namespace dlib;

using namespace pascinference;

int main( int argc, char *argv[] )
{
	/* add local program options */
	boost::program_options::options_description opt_problem("UTIL_DLIB_IMAGE_TO_VEC", consoleArg.get_console_nmb_cols());
	opt_problem.add_options()
		("in_filename", boost::program_options::value< std::string >(), "input image [string]")
		("out_filename", boost::program_options::value< std::string >(), "output vector in PETSc format [string]")
		("grayscale", boost::program_options::value< std::string >(), "create greyscale image or not [bool]")
		("new_width", boost::program_options::value< int >(), "width of new image [bool]")
		("noise", boost::program_options::value< double >(), "added noise [bool]");
	consoleArg.get_description()->add(opt_problem);

	/* call initialize */
	if(!Initialize<PetscVector>(argc, argv)){
		return 0;
	} 

	std::string in_filename;
	std::string out_filename;
	bool greyscale;
	int new_width;
	double noise;

	if(!consoleArg.set_option_value("in_filename", &in_filename)){
		std::cout << "in_filename has to be set! Call application with parameter -h to see all parameters\n";
		return 0;
	}

	consoleArg.set_option_value("out_filename", &out_filename, "results/image_output.bin");
	consoleArg.set_option_value("greyscale", &greyscale, true);
	consoleArg.set_option_value("new_width", &new_width, -1);
	consoleArg.set_option_value("noise", &noise, 0.0);

	coutMaster << "- UTIL INFO ----------------------------" << std::endl;
	coutMaster << " in_filename            = " << std::setw(30) << in_filename << " (input image)" << std::endl;
	coutMaster << " out_filename           = " << std::setw(30) << out_filename << " (output vector in PETSc format)" << std::endl;
	coutMaster << " greyscale              = " << std::setw(30) << print_bool(greyscale) << " (create greyscale image or not)" << std::endl;
	coutMaster << " new_width              = " << std::setw(30) << new_width << " (width of new image)" << std::endl;
	coutMaster << " noise                  = " << std::setw(30) << noise << " (added noise)" << std::endl;
	coutMaster << "-------------------------------------------" << std::endl;

	array2d<rgb_pixel> *image_dlib;

	/* rescale image */
	if(new_width > 1){
		/* open image */
		array2d<rgb_pixel> image_orig;
		load_image(image_orig, in_filename);
		
		int width_orig = image_orig.nr();
		int height_orig = image_orig.nc();
		
		image_dlib = new array2d<rgb_pixel>(new_width, height_orig*new_width/(double)width_orig);
		resize_image(image_orig, *image_dlib);

		coutMaster << image_dlib->nr() << "," << image_dlib->nc() << "," << std::endl;

		coutMaster << " orig_width             = " << std::setw(30) << width_orig << " px" << std::endl;
		coutMaster << " orig_height            = " << std::setw(30) << height_orig << " px" << std::endl;
	} else {
		image_dlib = new array2d<rgb_pixel>;
		
		/* open image */
		load_image(*image_dlib, in_filename);
	}

	/* print informations about image */
	int width = image_dlib->nr();
	int height = image_dlib->nc();
	int size =  width * height;
	coutMaster << " width                  = " << std::setw(30) << width << " px" << std::endl;
	coutMaster << " height                 = " << std::setw(30) << height << " px" << std::endl;
	coutMaster << " size                   = " << std::setw(30) << size << std::endl;
	coutMaster << "-------------------------------------------" << std::endl;

	/* prepare vector */
	Vec x_Vec;
	TRYCXX( VecCreate(PETSC_COMM_WORLD,&x_Vec) );
	TRYCXX( VecSetSizes(x_Vec,PETSC_DECIDE,size) );
	TRYCXX( VecSetType(x_Vec, VECSEQ) );
	TRYCXX( VecSetFromOptions(x_Vec) );

	coutMaster << " - convert Dlib image to PETSc Vec" << std::endl;
	double *x_arr;
	TRYCXX( VecGetArray(x_Vec, &x_arr) );
	for(int i=0;i<height;i++){
		for(int j=0;j<width;j++){
//			x_arr[i*width + j] = image_dlib(i)(j);
			if(greyscale){
				x_arr[i*width + j] = ((int)(*image_dlib)[i][j].red + (int)(*image_dlib)[i][j].green + (int)(*image_dlib)[i][j].blue)/(256.0*3.0);
			}

		}
	}
	TRYCXX( VecRestoreArray(x_Vec, &x_arr) );

	/* add noise */
	if(noise > 0){
		coutMaster << " - adding noise" << std::endl;
		Vec x_noise_Vec;
		TRYCXX( VecDuplicate(x_Vec, &x_noise_Vec) );

		PetscRandom rctx;
		TRYCXX( PetscRandomCreate(PETSC_COMM_WORLD,&rctx) );
		TRYCXX( PetscRandomSetFromOptions(rctx) );
		TRYCXX( VecSetRandom(x_noise_Vec,rctx) );
		TRYCXX( PetscRandomDestroy(&rctx) );
		
		TRYCXX( VecAXPY(x_Vec, 1.0, x_noise_Vec) );
		TRYCXX( VecDestroy(&x_noise_Vec) );
		
		coutMaster << " - projection to [0,1]" << std::endl;
		TRYCXX( VecGetArray(x_Vec, &x_arr) );
		for(int i=0;i<size;i++){
			if(x_arr[i] < 0.0){
				x_arr[i] = 0.0;
			}
			if(x_arr[i] > 1.0){
				x_arr[i] = 1.0;
			}
		}
		TRYCXX( VecRestoreArray(x_Vec, &x_arr) );
	}

	/* save vector */
	coutMaster << " - saving vector to file" << std::endl;
	PetscViewer mviewer;
	TRYCXX( PetscViewerCreate(PETSC_COMM_WORLD, &mviewer) );
	TRYCXX( PetscViewerBinaryOpen(PETSC_COMM_WORLD, out_filename.c_str(), FILE_MODE_WRITE, &mviewer) );
	TRYCXX( VecView(x_Vec, mviewer) );
	TRYCXX( PetscViewerDestroy(&mviewer) );
	coutMaster << " - new solution vector saved" << std::endl;

//	GeneralVector<PetscVector> in(in_Vec);

//	in.load_local(in_filename);

//	coutMaster << in << std::endl;

	Finalize<PetscVector>();

	return 0;
}

