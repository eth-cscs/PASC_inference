/** @file util_dlib_seqimage_to_vec.cpp
 *  @brief transform sequence of images to PETSc vector
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

#define DEFAULT_T 1
#define DEFAULT_XDIM 1
#define DEFAULT_TYPE 1
#define DEFAULT_WIDTH -1
#define DEFAULT_HEIGHT -1
#define DEFAULT_NOISE 0.0
#define DEFAULT_OUT_FILENAME "results/movie_output.bin"

#define ENABLE_ASSERTS
#define DLIB_JPEG_SUPPORT

#include <dlib/image_io.h>
#include <dlib/image_transforms.h>
#include <random>

using namespace dlib;

using namespace pascinference;

int compute_length_of_zeros(int T){
    int counter = 0;

    if(T <= 1){
        counter = 1;
    } else {
        double Tcounter = T;
        while(Tcounter >= 1){
            Tcounter = Tcounter/10.0;
            counter++;
        }
    }

    return counter;
}

std::string get_name_with_zeros(int t, int length_with_zeros){
    int t_length = compute_length_of_zeros(t);

	std::ostringstream sout;
    for(int i=0; i<length_with_zeros-t_length;i++){
        sout << "0";
    }
    sout << t;

	return sout.str();
}

int main( int argc, char *argv[] )
{
	/* add local program options */
	boost::program_options::options_description opt_problem("UTIL_DLIB_SEQIMAGE_TO_VEC", consoleArg.get_console_nmb_cols());
	opt_problem.add_options()
		("filename_in", boost::program_options::value< std::string >(), "part of input image <image>0t.jpg (t=1,..,T) [string]")
		("filename_out", boost::program_options::value< std::string >(), "output vector in PETSc format [string]")
		("xdim", boost::program_options::value< int >(), "number of values in every pixel [1=greyscale, 3=rgb]")
		("new_width", boost::program_options::value< int >(), "width of new image [int]")
		("new_height", boost::program_options::value< int >(), "height of new image [int]")
		("T", boost::program_options::value<int>(), "number of frames [int]")
		("type", boost::program_options::value< int >(), "type of output vector [0=TRn, 1=TnR, 2=nTR]")
		("noise", boost::program_options::value< double >(), "added noise [double]");
	consoleArg.get_description()->add(opt_problem);

	/* call initialize */
	if(!Initialize<PetscVector>(argc, argv)){
		return 0;
	} 

	std::string filename_in;
	std::string filename_out;
	int T;
	int xdim;
	int type;
	int new_width, new_height;
	double noise;

	if(!consoleArg.set_option_value("filename_in", &filename_in)){
		std::cout << "filename_in has to be set! Call application with parameter -h to see all parameters\n";
		return 0;
	}

	consoleArg.set_option_value("filename_out", &filename_out, DEFAULT_OUT_FILENAME);
	consoleArg.set_option_value("xdim", &xdim, DEFAULT_XDIM);
	consoleArg.set_option_value("T", &T, DEFAULT_T);
	consoleArg.set_option_value("type", &type, DEFAULT_TYPE);
	consoleArg.set_option_value("new_width", &new_width, DEFAULT_WIDTH);
	consoleArg.set_option_value("new_height", &new_height, DEFAULT_HEIGHT);
	consoleArg.set_option_value("noise", &noise, DEFAULT_NOISE);

	coutMaster << "- UTIL INFO ----------------------------" << std::endl;
	coutMaster << " filename_in            = " << std::setw(30) << filename_in << " (input image)" << std::endl;
	coutMaster << " filename_out           = " << std::setw(30) << filename_out << " (output vector in PETSc format)" << std::endl;
	coutMaster << " xdim                   = " << std::setw(30) << xdim << " (number of values in every pixel [1=greyscale, 3=rgb])" << std::endl;
	coutMaster << " type                   = " << std::setw(30) << Decomposition<PetscVector>::get_type_name(type) << " (type of output vector [" << Decomposition<PetscVector>::get_type_list() << "])" << std::endl;
	coutMaster << " new_width              = " << std::setw(30) << new_width << " (width of new image)" << std::endl;
	coutMaster << " new_height             = " << std::setw(30) << new_height << " (height of new image)" << std::endl;
	coutMaster << " T                      = " << std::setw(30) << T << " (number of frames in movie)" << std::endl;
	coutMaster << " noise                  = " << std::setw(30) << noise << " (added noise)" << std::endl;
	coutMaster << "-------------------------------------------" << std::endl;

	array2d<rgb_pixel> *image_dlib;

	/* construct name of first image */
	int length_with_zeros = compute_length_of_zeros(T);
	std::ostringstream oss;
	oss.str("");
    oss << filename_in << get_name_with_zeros(1, length_with_zeros) << ".jpg";	

	/* open first image set size */
	if(new_width > 1 || new_height > 1){
		/* open image */
		array2d<rgb_pixel> image_orig;
		load_image(image_orig, oss.str());
		
		int width_orig = image_orig.nc();
		int height_orig = image_orig.nr();

		int width_scaled = new_width;
		int height_scaled = new_height;
		
		if(width_scaled <= 0){
			width_scaled = width_orig*new_height/(double)height_orig;
		}

		if(height_scaled <= 0){
			height_scaled = height_orig*new_width/(double)width_orig;
		}
		
		/* scale image */
		image_dlib = new array2d<rgb_pixel>(height_scaled, width_scaled);
		resize_image(image_orig, *image_dlib);

		coutMaster << " orig_width             = " << std::setw(30) << width_orig << " px" << std::endl;
		coutMaster << " orig_height            = " << std::setw(30) << height_orig << " px" << std::endl;
	} else {
		image_dlib = new array2d<rgb_pixel>;
		
		/* open image */
		load_image(*image_dlib, oss.str());
	}

	/* print informations about image */
	int width = image_dlib->nc();
	int height = image_dlib->nr();
	int R =  width * height;
	int nvalues = xdim * R * T;
	coutMaster << " width                  = " << std::setw(30) << width << " px" << std::endl;
	coutMaster << " height                 = " << std::setw(30) << height << " px" << std::endl;
	coutMaster << " R                      = " << std::setw(30) << R << std::endl;
	coutMaster << " nvalues                = " << std::setw(30) << nvalues << std::endl;
	coutMaster << "-------------------------------------------" << std::endl;

	/* prepare vector */
	Vec x_Vec;
	TRYCXX( VecCreate(PETSC_COMM_WORLD,&x_Vec) );
	TRYCXX( VecSetSizes(x_Vec,PETSC_DECIDE, nvalues) );
	TRYCXX( VecSetType(x_Vec, VECSEQ) );
	TRYCXX( VecSetFromOptions(x_Vec) );

	coutMaster << " - convert Dlib images to PETSc Vec" << std::endl;
	double *x_arr;
	TRYCXX( VecGetArray(x_Vec, &x_arr) );

	for(int t=0;t<T;t++){
		/* construct name of file */
		oss.str("");
		oss << filename_in << get_name_with_zeros(t+1, length_with_zeros) << ".jpg";	
		
		/* open image */
		array2d<rgb_pixel> image_orig;
		load_image(image_orig, oss.str());
		
		/* scale image */
		resize_image(image_orig, *image_dlib);

		//TODO: rescale loaded image?
		
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
			
				/* 0=TRn, 1=TnR, 2=nTR */
			
				/* type TRn */
				if(type==0){
					/* greyscale */
					if(xdim == 1){
						x_arr[t*R*xdim + i*width + j] = ((int)(*image_dlib)[i][j].red + (int)(*image_dlib)[i][j].green + (int)(*image_dlib)[i][j].blue)/(255.0*3.0);
					}

					/* rgb */
					if(xdim == 3){
						x_arr[t*R*xdim + i*width*xdim + j*xdim + 0] = ((int)(*image_dlib)[i][j].red)/255.0;
						x_arr[t*R*xdim + i*width*xdim + j*xdim + 1] = ((int)(*image_dlib)[i][j].green)/255.0;
						x_arr[t*R*xdim + i*width*xdim + j*xdim + 2] = ((int)(*image_dlib)[i][j].blue)/255.0;
					}
				}

				/* type TnR */
				if(type==1){
					/* greyscale */
					if(xdim == 1){
						x_arr[t*R*xdim + i*width + j] = ((int)(*image_dlib)[i][j].red + (int)(*image_dlib)[i][j].green + (int)(*image_dlib)[i][j].blue)/(255.0*3.0);
					}

					/* rgb */
					if(xdim == 3){
						x_arr[t*R*xdim + 0*width*height + i*width + j] = ((int)(*image_dlib)[i][j].red)/255.0;
						x_arr[t*R*xdim + 1*width*height + i*width + j] = ((int)(*image_dlib)[i][j].green)/255.0;
						x_arr[t*R*xdim + 2*width*height + i*width + j] = ((int)(*image_dlib)[i][j].blue)/255.0;
					}
				}

				/* type nTR */
				if(type==2){
					/* greyscale */
					if(xdim == 1){
						x_arr[t*R + i*width + j] = ((int)(*image_dlib)[i][j].red + (int)(*image_dlib)[i][j].green + (int)(*image_dlib)[i][j].blue)/(255.0*3.0);
					}

					/* rgb */
					if(xdim == 3){
						x_arr[0*T*R + t*R + i*width + j] = ((int)(*image_dlib)[i][j].red)/255.0;
						x_arr[1*T*R + t*R + i*width + j] = ((int)(*image_dlib)[i][j].green)/255.0;
						x_arr[2*T*R + t*R + i*width + j] = ((int)(*image_dlib)[i][j].blue)/255.0;
					}
				}
				
			}
		}
	}
	
	TRYCXX( VecRestoreArray(x_Vec, &x_arr) );

	/* add noise */
	if(noise > 0){
		coutMaster << " - adding noise" << std::endl;

		std::default_random_engine generator;
		std::normal_distribution<double> distribution(0.0,noise);
		
		coutMaster << " - projection to [0,1]" << std::endl;

		double *x_noise_arr;
		double value;

		TRYCXX( VecGetArray(x_Vec, &x_arr) );
		for(int i=0;i<nvalues;i++){

			double number = distribution(generator);

			/* add noise */
			x_arr[i] += number;
			
			/* projection */
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
	TRYCXX( PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename_out.c_str(), FILE_MODE_WRITE, &mviewer) );
	TRYCXX( VecView(x_Vec, mviewer) );
	TRYCXX( PetscViewerDestroy(&mviewer) );
	coutMaster << " - new vector saved" << std::endl;

	Finalize<PetscVector>();

	return 0;
}

