/** @file process_log.cpp
 *  @brief process log file
 *
 *  @author Lukas Pospisil
 */

#include "pascinference.h"

#ifndef USE_PETSC
 #error 'This example is for PETSC'
#endif

#ifndef USE_BOOST
 #error 'This example is for BOOST'
#endif

#include <iostream>
#include <fstream>
#include <string>

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>

using namespace pascinference;

class log_processor {
	private:
		bool is_opened;
		std::ifstream myfile(in_filename); /* the file reader */
		int line_counter = 0; /* line

		void process_line(std::string myline){
			
		}
	
	public:
		log_processor(std::string filename){
			if (myfile.is_open()){ /* if the file can be opened */
				is_opened = true;
			} else {
				is_opened = false;
			}
		}

		log_processor(){
			if(is_opened){
				/* close the file after reading */
				myfile.close(); 
			}
		}
	
		bool is_opened(){
			return is_opened;
		}
	
		void process(std::delimiter){
			line_counter = 0;
			std::string myline; /* one loaded line from the file */
			while ( getline (myfile,myline) ){ /* we are reading one line */

				/* process the line */
			std::vector<std::string> mywords;
			boost::split(mywords, myline, boost::is_any_of(delimiter), boost::token_compress_on);

			std::cout << myline << std::endl;

			for(int i=0;i<mywords.size();i++){
				std::cout << "  " << mywords[i] << std::endl;
			}
			std::cout << std::endl;


			/* increase counter of lines */
			line_counter++;
			if(line_counter > 5){
				break;
			}
		}
		
	} else {
		/* the file cannot be opened, probably it does not exist */
		coutMaster << "ERROR: Unable to open file" << std::endl;
	}

			
	
}

int main( int argc, char *argv[] )
{
	/* add local program options */
	boost::program_options::options_description opt_problem("UTIL_PROCESS_LOG", consoleArg.get_console_nmb_cols());
	opt_problem.add_options()
		("in_filename", boost::program_options::value< std::string >(), "input log file [string]")
		("delimiter", boost::program_options::value< std::string >(), "delimiter [string]");
	consoleArg.get_description()->add(opt_problem);

	/* call initialize */
	if(!Initialize<PetscVector>(argc, argv)){
		return 0;
	} 

	/* read filename from console */
	std::string in_filename;
	if(!consoleArg.set_option_value("in_filename", &in_filename)){
		std::cout << "in_filename has to be set! Call application with parameter -h to see all parameters\n";
		return 0;
	}

	/* read delimiter from console */
	std::string delimiter;
	consoleArg.set_option_value("delimiter", &delimiter, LOG_SEPARATOR);

	/* give some informations about the settings */
	coutMaster << "- PROCESS LOG ----------------------------------------" << std::endl;
	coutMaster << " in_filename            = " << std::setw(30) << in_filename << " (log file)" << std::endl;
	coutMaster << " delimiter              = " << std::setw(30) << delimiter << std::endl;
	coutMaster << " log_or_not_file_line   = " << std::setw(30) << print_bool(logging.get_log_or_not_file_line()) << std::endl;
	coutMaster << " log_or_not_func_call   = " << std::setw(30) << print_bool(logging.get_log_or_not_func_call()) << std::endl;
	coutMaster << " log_or_not_level       = " << std::setw(30) << print_bool(logging.get_log_or_not_level()) << std::endl;
	coutMaster << " log_or_not_memory      = " << std::setw(30) << print_bool(logging.get_log_or_not_memory()) << std::endl;
	coutMaster << "------------------------------------------------------" << std::endl;
	coutMaster << std::endl;




	Finalize<PetscVector>();

	return 0;
}

