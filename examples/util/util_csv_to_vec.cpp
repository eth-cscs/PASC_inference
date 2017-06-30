/** @file util_csv_to_vec.cpp
 *  @brief transform csv to PETSc vector
 *
 *  @author Lukas Pospisil
 */

#include "pascinference.h"
#include <vector>
#include <sstream>
#include <string>

#ifndef USE_PETSC
 #error 'This example is for PETSC'
#endif
 
#define DEFAULT_OUT_FILENAME "results/sample_vec.bin"

using namespace pascinference;

void __mytrim(std::string &str,bool trim_left,bool trim_right)
{
    std::string::size_type begin=0;
    std::string::size_type end=str.size()-1;
    while(trim_left && begin<=end && (str[begin]<=0x20 || str[begin]==0x7f))
        ++begin;
    while(trim_right && end>begin && (str[end]<=0x20 || str[end]==0x7f))
      --end;
    str = str.substr(begin, end - begin + 1);
}
void mytrim(std::string &str)  { __mytrim(str,1,1); }


void ReadCSVFile(std::string s_Path, long T, int xdim, double *values){
	std::ifstream dataFile(s_Path.c_str());

    if (!dataFile)
    {
        std::cerr << "Error opening file:" << s_Path << std::endl;
    }
    
    std::string line;
    std::getline(dataFile, line);
    dataFile.close();
   
	std::stringstream iss(line);
    
    int i = 0; int j = 0;
    while(true){
		std::string val;
		double data;
        
        std::getline(iss, val, ',');
        if(!iss.good())
            break;
        
        mytrim(val);
        std::stringstream convertor(val);
        convertor >> data;

        values[j*T+i] = data;
        if (i == T-1)
        {
            i = 0; j++;
        }
        else
            i++;
    }
}

int main( int argc, char *argv[] )
{
	/* add local program options */
	boost::program_options::options_description opt_problem("UTIL_CSV_TO_VEC", consoleArg.get_console_nmb_cols());
	opt_problem.add_options()
		("filename_in", boost::program_options::value< std::string >(), "input CSV file [string]")
		("filename_out", boost::program_options::value< std::string >(), "output vector in PETSc format [string]")
		("T", boost::program_options::value< int >(), "length of time series [int]")
		("xdim", boost::program_options::value< int >(), "dimension of data [int]")
		("print", boost::program_options::value< bool >(), "print vector [bool]")
		("permute", boost::program_options::value< bool >(), "permute data from nT to Tn format [bool]");
	consoleArg.get_description()->add(opt_problem);

	/* call initialize */
	if(!Initialize<PetscVector>(argc, argv)){
		return 0;
	} 

	std::string filename_in;
	std::string filename_out;
	int xdim;
	int T;
	bool print_or_not, permute_or_not;

	if(!consoleArg.set_option_value("filename_in", &filename_in)){
		std::cout << "filename_in has to be set! Call application with parameter -h to see all parameters\n";
		return 0;
	}

	if(!consoleArg.set_option_value("xdim", &xdim)){
		std::cout << "xdim has to be set! Call application with parameter -h to see all parameters\n";
		return 0;
	}

	if(!consoleArg.set_option_value("T", &T)){
		std::cout << "T has to be set! Call application with parameter -h to see all parameters\n";
		return 0;
	}

	consoleArg.set_option_value("filename_out", &filename_out, DEFAULT_OUT_FILENAME);
	consoleArg.set_option_value("print", &print_or_not, false);
	consoleArg.set_option_value("permute", &permute_or_not, false);


	coutMaster << "- UTIL INFO ----------------------------" << std::endl;
	coutMaster << " filename_in            = " << std::setw(30) << filename_in << " (input image)" << std::endl;
	coutMaster << " filename_out           = " << std::setw(30) << filename_out << " (output vector in PETSc format)" << std::endl;
	coutMaster << " xdim                   = " << std::setw(30) << xdim << " (dimension of data)" << std::endl;
	coutMaster << " T                      = " << std::setw(30) << T << " (length of time series)" << std::endl;
	coutMaster << " print                  = " << std::setw(30) << print_bool(print_or_not) << " (print the content of vector or not)" << std::endl;
	coutMaster << " permute                = " << std::setw(30) << print_bool(permute_or_not) << " (permute data from nT to Tn format)" << std::endl;
	coutMaster << "-------------------------------------------" << std::endl;

	/* prepare vector */
	Vec x_Vec;
	TRYCXX( VecCreate(PETSC_COMM_WORLD,&x_Vec) );
	TRYCXX( VecSetSizes(x_Vec,PETSC_DECIDE, T*xdim) );
	TRYCXX( VecSetType(x_Vec, VECSEQ) );
	TRYCXX( VecSetFromOptions(x_Vec) );

	coutMaster << " - convert CSV file to PETSc Vec" << std::endl;
	double *x_arr;
	TRYCXX( VecGetArray(x_Vec, &x_arr) );

	ReadCSVFile(filename_in, T, xdim, x_arr);
	
	if(print_or_not){
		for(int j=0;j<xdim;j++){
			std::cout << "dim = " << j << std::endl;
			std::cout << "[";
			for(int i=0;i<T;i++){
				std::cout << x_arr[j*T+i];
				if(i < T-1) std::cout << ", ";
			}
			std::cout << "]" << std::endl;
		}
	}
	TRYCXX( VecRestoreArray(x_Vec, &x_arr) );

	if(permute_or_not){
		Vec x_orig_Vec;
		TRYCXX( VecDuplicate(x_Vec,&x_orig_Vec) );
		TRYCXX( VecCopy(x_Vec,x_orig_Vec) );
		
		double *x_orig_arr;
		TRYCXX( VecGetArray(x_Vec, &x_arr) );
		TRYCXX( VecGetArray(x_orig_Vec, &x_orig_arr) );
		
		for(int n=0;n<xdim;n++){
			for(int t=0;t<T;t++){
				x_arr[t*xdim+n] = x_orig_arr[n*T+t];
			}
		}

		TRYCXX( VecRestoreArray(x_Vec, &x_arr) );
		TRYCXX( VecRestoreArray(x_orig_Vec, &x_orig_arr) );

		TRYCXX( VecDestroy(&x_orig_Vec) );
	}

	/* save vector */
	coutMaster << " - saving vector to file" << std::endl;
	PetscViewer mviewer;
	TRYCXX( PetscViewerCreate(PETSC_COMM_WORLD, &mviewer) );
	TRYCXX( PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename_out.c_str(), FILE_MODE_WRITE, &mviewer) );
	TRYCXX( VecView(x_Vec, mviewer) );
	TRYCXX( PetscViewerDestroy(&mviewer) );
	coutMaster << " - new solution vector saved" << std::endl;

	Finalize<PetscVector>();

	return 0;
}

