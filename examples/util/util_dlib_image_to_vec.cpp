/** @file util_dlib_image_to_vec.cpp
 *  @brief transform image to PETSc vector
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

#define DEFAULT_XDIM 1
#define DEFAULT_TYPE 1
#define DEFAULT_WIDTH -1
#define DEFAULT_NOISE 0.0
#define DEFAULT_OUT_FILENAME "results/image_output.bin"

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
		("xdim", boost::program_options::value< int >(), "number of values in every pixel [1=greyscale, 3=rgb]")
		("new_width", boost::program_options::value< int >(), "width of new image [bool]")
		("type", boost::program_options::value< int >(), "type of output vector [0=Rn, 1=nR]")
		("noise", boost::program_options::value< double >(), "added noise [bool]");
	consoleArg.get_description()->add(opt_problem);

	/* call initialize */
	if(!Initialize<PetscVector>(argc, argv)){
		return 0;
	} 

	std::string in_filename;
	std::string out_filename;
	int xdim;
	int type;
	int new_width;
	double noise;

	if(!consoleArg.set_option_value("in_filename", &in_filename)){
		std::cout << "in_filename has to be set! Call application with parameter -h to see all parameters\n";
		return 0;
	}

	consoleArg.set_option_value("out_filename", &out_filename, DEFAULT_OUT_FILENAME);
	consoleArg.set_option_value("xdim", &xdim, DEFAULT_XDIM);
	consoleArg.set_option_value("type", &type, DEFAULT_TYPE);
	consoleArg.set_option_value("new_width", &new_width, DEFAULT_WIDTH);
	consoleArg.set_option_value("noise", &noise, DEFAULT_NOISE);

	coutMaster << "- UTIL INFO ----------------------------" << std::endl;
	coutMaster << " in_filename            = " << std::setw(30) << in_filename << " (input image)" << std::endl;
	coutMaster << " out_filename           = " << std::setw(30) << out_filename << " (output vector in PETSc format)" << std::endl;
	coutMaster << " xdim                   = " << std::setw(30) << xdim << " (number of values in every pixel [1=greyscale, 3=rgb])" << std::endl;
	coutMaster << " type                   = " << std::setw(30) << type << " (type of output vector [0=Rn, 1=nR])" << std::endl;
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
	int nvalues = xdim * size;
	coutMaster << " width                  = " << std::setw(30) << width << " px" << std::endl;
	coutMaster << " height                 = " << std::setw(30) << height << " px" << std::endl;
	coutMaster << " size                   = " << std::setw(30) << size << std::endl;
	coutMaster << " nvalues                = " << std::setw(30) << nvalues << std::endl;
	coutMaster << "-------------------------------------------" << std::endl;

	/* prepare vector */
	Vec x_Vec;
	TRYCXX( VecCreate(PETSC_COMM_WORLD,&x_Vec) );
	TRYCXX( VecSetSizes(x_Vec,PETSC_DECIDE, nvalues) );
	TRYCXX( VecSetType(x_Vec, VECSEQ) );
	TRYCXX( VecSetFromOptions(x_Vec) );

	coutMaster << " - convert Dlib image to PETSc Vec" << std::endl;
	double *x_arr;
	TRYCXX( VecGetArray(x_Vec, &x_arr) );
	for(int i=0;i<height;i++){
		for(int j=0;j<width;j++){
			
			/* type Rn */
			if(type==0){
				/* greyscale */
				if(xdim == 1){
					x_arr[i*width + j] = ((int)(*image_dlib)[i][j].red + (int)(*image_dlib)[i][j].green + (int)(*image_dlib)[i][j].blue)/(255.0*3.0);
				}

				/* rgb */
				if(xdim == 3){
					x_arr[i*width*xdim + j*xdim + 0] = ((int)(*image_dlib)[i][j].red)/255.0;
					x_arr[i*width*xdim + j*xdim + 1] = ((int)(*image_dlib)[i][j].green)/255.0;
					x_arr[i*width*xdim + j*xdim + 2] = ((int)(*image_dlib)[i][j].blue)/255.0;
				}
			}

			/* type nR */
			if(type==1){
				/* greyscale */
				if(xdim == 1){
					x_arr[i*width + j] = ((int)(*image_dlib)[i][j].red + (int)(*image_dlib)[i][j].green + (int)(*image_dlib)[i][j].blue)/(255.0*3.0);
				}

				/* rgb */
				if(xdim == 3){
					x_arr[0*width*height + i*width + j] = ((int)(*image_dlib)[i][j].red)/255.0;
					x_arr[1*width*height + i*width + j] = ((int)(*image_dlib)[i][j].green)/255.0;
					x_arr[2*width*height + i*width + j] = ((int)(*image_dlib)[i][j].blue)/255.0;
				}
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
		
		coutMaster << " - projection to [0,1]" << std::endl;

		double *x_noise_arr;
		double value;

		TRYCXX( VecGetArray(x_Vec, &x_arr) );
		TRYCXX( VecGetArray(x_noise_Vec, &x_noise_arr) );
		for(int i=0;i<nvalues;i++){
			/* add noise */
			x_arr[i] = x_arr[i] + 2*noise*(x_noise_arr[i] - 0.5);
			
			/* projection */
			if(x_arr[i] < 0.0){
				x_arr[i] = 0.0;
			}
			if(x_arr[i] > 1.0){
				x_arr[i] = 1.0;
			}
		}
		TRYCXX( VecRestoreArray(x_noise_Vec, &x_noise_arr) );
		TRYCXX( VecRestoreArray(x_Vec, &x_arr) );

		TRYCXX( VecDestroy(&x_noise_Vec) );

	}

	/* save vector */
	coutMaster << " - saving vector to file" << std::endl;
	PetscViewer mviewer;
	TRYCXX( PetscViewerCreate(PETSC_COMM_WORLD, &mviewer) );
	TRYCXX( PetscViewerBinaryOpen(PETSC_COMM_WORLD, out_filename.c_str(), FILE_MODE_WRITE, &mviewer) );
	TRYCXX( VecView(x_Vec, mviewer) );
	TRYCXX( PetscViewerDestroy(&mviewer) );
	coutMaster << " - new solution vector saved" << std::endl;

	Finalize<PetscVector>();

	return 0;
}
