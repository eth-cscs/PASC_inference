/** @file util_dlib_vec_to_image.cpp
 *  @brief transform PETSc vector to image or sequence of images
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
#define DEFAULT_NOISE 0.0
#define DEFAULT_OUT_FILENAME "results/image_output.jpg"
#define JPEG_QUALITY 95

#define ENABLE_ASSERTS
#define DLIB_JPEG_SUPPORT

#include <dlib/image_io.h>

using namespace dlib;

using namespace pascinference;

std::string cut_filename(const std::string &input);
std::string cut_extension(const std::string &input);

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
	boost::program_options::options_description opt_problem("UTIL_DLIB_IMAGE_TO_VEC", consoleArg.get_console_nmb_cols());
	opt_problem.add_options()
		("filename_in1", boost::program_options::value< std::string >(), "input vector in PETSc format [string]")
		("filename_in2", boost::program_options::value< std::string >(), "optional second input vector in PETSc format [string]")
		("filename_in3", boost::program_options::value< std::string >(), "optional third input vector in PETSc format [string]")
		("K3", boost::program_options::value< int >(), "number of clusters for third image representing gamma [int]")
		("K3only", boost::program_options::value< int >(), "the only cluster which is plotted, otherwise all are plotted [int]")
		("filename_out", boost::program_options::value< std::string >(), "output image in jpeg format [string]")
		("width", boost::program_options::value< int >(), "width of image [int]")
		("T", boost::program_options::value<int>(), "number of frames [int]")
		("xdim", boost::program_options::value< int >(), "number of values in every pixel [1=greyscale, 3=rgb]")
		("type", boost::program_options::value< int >(), "type of input vector [0=TRn, 1=TnR, 2=nTR]");

	consoleArg.get_description()->add(opt_problem);

	/* call initialize */
	if(!Initialize<PetscVector>(argc, argv)){
		return 0;
	}

	std::string filename_in1; bool given_first;
	std::string filename_in2; bool given_second;
	std::string filename_in3; bool given_third;
	std::string filename_out;
    int width, T, xdim, type, K3, K3only;

	if(consoleArg.set_option_value("filename_in1", &filename_in1)){
		given_first = true;
	} else {
        given_first = false;
	}

	if(consoleArg.set_option_value("filename_in2", &filename_in2)){
		given_second = true;
	} else {
        given_second = false;
	}

	if(consoleArg.set_option_value("filename_in3", &filename_in3)){
		given_third = true;
		if(!consoleArg.set_option_value("K3", &K3)){
			std::cout << "K3 for last input vector has to be set! Call application with parameter -h to see all parameters\n";
			return 0;
		}
		if(!consoleArg.set_option_value("K3only", &K3only)){
			K3only = -1;
		}

	} else {
        given_third = false;
	}

	if(!consoleArg.set_option_value("width", &width)){
		std::cout << "width of image has to be set! Call application with parameter -h to see all parameters\n";
		return 0;
	}

	consoleArg.set_option_value("filename_out", &filename_out, DEFAULT_OUT_FILENAME);
	consoleArg.set_option_value("xdim", &xdim, DEFAULT_XDIM);
	consoleArg.set_option_value("type", &type, DEFAULT_TYPE);
	consoleArg.set_option_value("T", &T, DEFAULT_T);

	coutMaster << "- UTIL INFO ----------------------------" << std::endl;

	coutMaster << " filename_in1           = ";
    if(given_first){
        coutMaster << std::setw(30) << filename_in1;
    } else {
        coutMaster << std::setw(30) << print_bool(false);
    }
    coutMaster << " (input vector in PETSc format)" << std::endl;

    coutMaster << " filename_in2           = ";
    if(given_second){
        coutMaster << std::setw(30) << filename_in2;
    } else {
        coutMaster << std::setw(30) << print_bool(false);
    }
    coutMaster << " (optional second input vector in PETSc format)" << std::endl;

    coutMaster << " filename_in3           = ";
    if(given_third){
        coutMaster << std::setw(30) << filename_in3 << " (optional third input vector in PETSc format)" << std::endl;
		coutMaster << " K3                     = " << std::setw(30) << K3 << " (number of clusters for third image representing gamma)" << std::endl;
		coutMaster << " K3only                 = " << std::setw(30) << K3only << " (the only cluster which is plotted, otherwise all are plotted)" << std::endl;
    } else {
        coutMaster << std::setw(30) << print_bool(false) << " (optional third input vector in PETSc format)" << std::endl; 
    }
    coutMaster << " filename_out           = " << std::setw(30) << filename_out << " (output image in jpeg format)" << std::endl;
	coutMaster << " xdim                   = " << std::setw(30) << xdim << " (number of values in every pixel [1=greyscale, 3=rgb])" << std::endl;
	coutMaster << " type                   = " << std::setw(30) << Decomposition<PetscVector>::get_type_name(type) << " (type of output vector [" << Decomposition<PetscVector>::get_type_list() << "])" << std::endl;
	coutMaster << " width                  = " << std::setw(30) << width << " (width of image)" << std::endl;
	coutMaster << " T                      = " << std::setw(30) << T << " (number of frames in movie)" << std::endl;
	coutMaster << "-------------------------------------------" << std::endl;

	int how_many_widths = 0;

	/* prepare and load data vector */
	Vec data1_Vec;
	TRYCXX( VecCreate(PETSC_COMM_WORLD,&data1_Vec) );
	GeneralVector<PetscVector> data1(data1_Vec);
    if(given_first){
        data1.load_local(filename_in1);
        how_many_widths++;
    }

    /* prepare second optional vector */
	Vec data2_Vec;
    TRYCXX( VecCreate(PETSC_COMM_WORLD,&data2_Vec) );
    GeneralVector<PetscVector> data2(data2_Vec);
    if(given_second){
        data2.load_local(filename_in2);
		how_many_widths++;
    }

    /* prepare third optional vector */
	Vec data3_Vec;
    TRYCXX( VecCreate(PETSC_COMM_WORLD,&data3_Vec) );
    GeneralVector<PetscVector> data3(data3_Vec);
    if(given_third){
        data3.load_local(filename_in3);
		if(K3only > 0){
			how_many_widths++;
		} else {
			how_many_widths += K3;
		}
    }

	/* get filename path and extension */
	boost::filesystem::path p(filename_out);
    std::ostringstream oss;
    oss << p.extension();
	std::string extension = oss.str();
	oss.str("");
	if(!extension.compare(".jpg")){
		coutMaster << "ERROR: extension of file '" << extension << "' was not recognized." << std::endl;
        return 0;
	}

    /* remove .jpg from the end of filename */
    std::string dirname_out = filename_out.substr(0, filename_out.size()-4);
	boost::filesystem::path dirname(dirname_out);

	/* if it is movie, then create directory */
    int length_with_zeros; /* length of time image with zeros */
    if(T>1){
        boost::filesystem::create_directory(dirname);
        length_with_zeros = compute_length_of_zeros(T);
    }

    /* compute height of image */
	int data_size;
	TRYCXX( VecGetSize(data1.get_vector(), &data_size) );
	int height = (int)((data_size/(double)xdim)/(double)(width*T));

    /* prepare Dlib stuff for exporting */
	array2d<rgb_pixel> image_dlib;
    image_dlib.set_size(height,how_many_widths*width);

	int R = width*height;

	double *x1_arr;
	if(given_first) TRYCXX( VecGetArray(data1.get_vector(), &x1_arr) );

	double *x2_arr;
    if(given_second) TRYCXX( VecGetArray(data2.get_vector(), &x2_arr) );

	double *x3_arr;
    if(given_third) TRYCXX( VecGetArray(data3.get_vector(), &x3_arr) );

	double left_start;

    for(int t=0; t < T; t++){
        for(int i=0;i<height;i++){
            for(int j=0;j<width;j++){

                /* greyscale */
                if(xdim==1){
                    if(type == 0 | type == 1 | type == 2){

                        left_start = 0.0;

                        if(given_first){
							image_dlib[i][left_start+j].red = x1_arr[t*R + i*width + j]*255;
							image_dlib[i][left_start+j].green = x1_arr[t*R + i*width + j]*255;
							image_dlib[i][left_start+j].blue = x1_arr[t*R + i*width + j]*255;
							
							left_start += width;
						}
						
                        if(given_second){
                            image_dlib[i][left_start+j].red = x2_arr[t*R + i*width + j]*255;
                            image_dlib[i][left_start+j].green = x2_arr[t*R + i*width + j]*255;
                            image_dlib[i][left_start+j].blue = x2_arr[t*R + i*width + j]*255;

							left_start += width;
                        }

                        if(given_third){
							if(K3only < 0){
								left_start += width;
								image_dlib[i][left_start+j].red = x3_arr[t*K3*R + K3only*R + i*width + j]*255;
								image_dlib[i][left_start+j].green = x3_arr[t*K3*R + K3only*R + i*width + j]*255;
								image_dlib[i][left_start+j].blue = x3_arr[t*K3*R + K3only*R + i*width + j]*255;
							} else {
								for(int k=0;k<K3;k++){
									left_start += width;
									image_dlib[i][left_start+j].red = x3_arr[t*K3*R + k*R + i*width + j]*255;
									image_dlib[i][left_start+j].green = x3_arr[t*K3*R + k*R + i*width + j]*255;
									image_dlib[i][left_start+j].blue = x3_arr[t*K3*R + k*R + i*width + j]*255;
								}
							}
                        }
                    }
                }

                if(xdim==3){
                    /* TRn */
                    if(type == 0){
                        image_dlib[i][j].red   = x1_arr[t*xdim*width*height + i*width*xdim + j*xdim + 0]*255;
                        image_dlib[i][j].green = x1_arr[t*xdim*width*height + i*width*xdim + j*xdim + 1]*255;
                        image_dlib[i][j].blue  = x1_arr[t*xdim*width*height + i*width*xdim + j*xdim + 2]*255;

                        if(given_second){
                            image_dlib[i][width+j].red   = x2_arr[t*xdim*width*height + i*width*xdim + j*xdim + 0]*255;
                            image_dlib[i][width+j].green = x2_arr[t*xdim*width*height + i*width*xdim + j*xdim + 1]*255;
                            image_dlib[i][width+j].blue  = x2_arr[t*xdim*width*height + i*width*xdim + j*xdim + 2]*255;
                        }

                        if(given_third){
                            image_dlib[i][2*width+j].red   = x3_arr[t*xdim*width*height + i*width*xdim + j*xdim + 0]*255;
                            image_dlib[i][2*width+j].green = x3_arr[t*xdim*width*height + i*width*xdim + j*xdim + 1]*255;
                            image_dlib[i][2*width+j].blue  = x3_arr[t*xdim*width*height + i*width*xdim + j*xdim + 2]*255;
                        }
                    }

                    /* TnR */
                    if(type == 1){
                        image_dlib[i][j].red   = x1_arr[t*xdim*width*height + 0*width*height + i*width + j]*255;
                        image_dlib[i][j].green = x1_arr[t*xdim*width*height + 1*width*height + i*width + j]*255;
                        image_dlib[i][j].blue  = x1_arr[t*xdim*width*height + 2*width*height + i*width + j]*255;

                        if(given_second){
                            image_dlib[i][width+j].red   = x2_arr[t*xdim*width*height + 0*width*height + i*width + j]*255;
                            image_dlib[i][width+j].green = x2_arr[t*xdim*width*height + 1*width*height + i*width + j]*255;
                            image_dlib[i][width+j].blue  = x2_arr[t*xdim*width*height + 2*width*height + i*width + j]*255;
                        }

                        if(given_third){
                            image_dlib[i][2*width+j].red   = x3_arr[t*xdim*width*height + 0*width*height + i*width + j]*255;
                            image_dlib[i][2*width+j].green = x3_arr[t*xdim*width*height + 1*width*height + i*width + j]*255;
                            image_dlib[i][2*width+j].blue  = x3_arr[t*xdim*width*height + 2*width*height + i*width + j]*255;
                        }
                    }

                    /* nTR */
                    if(type == 2){
                        image_dlib[i][j].red   = x1_arr[0*T*width*height + t*width*height + i*width + j]*255;
                        image_dlib[i][j].green = x1_arr[1*T*width*height + t*width*height + i*width + j]*255;
                        image_dlib[i][j].blue  = x1_arr[2*T*width*height + t*width*height + i*width + j]*255;

                        if(given_second){
                            image_dlib[i][width+j].red   = x2_arr[0*T*width*height + t*width*height + i*width + j]*255;
                            image_dlib[i][width+j].green = x2_arr[1*T*width*height + t*width*height + i*width + j]*255;
                            image_dlib[i][width+j].blue  = x2_arr[2*T*width*height + t*width*height + i*width + j]*255;
                        }

                        if(given_third){
                            image_dlib[i][2*width+j].red   = x3_arr[0*T*width*height + t*width*height + i*width + j]*255;
                            image_dlib[i][2*width+j].green = x3_arr[1*T*width*height + t*width*height + i*width + j]*255;
                            image_dlib[i][2*width+j].blue  = x3_arr[2*T*width*height + t*width*height + i*width + j]*255;
                        }

                    }

                }
            }
        }

        if(T == 1){
            /* save image directly */
            save_jpeg(image_dlib, filename_out, JPEG_QUALITY);
        } else {
            /* compose the name of image */
            oss.str("");
            oss << dirname_out << "/" << get_name_with_zeros(t+1, length_with_zeros) << ".jpg";

            save_jpeg(image_dlib, oss.str(), JPEG_QUALITY);
        }
    }
	if(given_first) TRYCXX( VecRestoreArray(data1.get_vector(), &x1_arr) );
    if(given_second) TRYCXX( VecRestoreArray(data2.get_vector(), &x2_arr) );
    if(given_third) TRYCXX( VecRestoreArray(data3.get_vector(), &x3_arr) );

    //TODO: destroy vectors?

	Finalize<PetscVector>();

	return 0;
}

