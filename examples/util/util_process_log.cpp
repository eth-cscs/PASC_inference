/** @file util_process_log.cpp
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
		class my_function_node {
			private:
				std::string myname;
				int counter;
				double timer;
				double timer_stamp;
			
			public:
				my_function_node(std::string new_myname){
					myname = new_myname;
					counter = 0;
					timer = 0.0;
					timer_stamp = 0.0;
				}
				
				void begin_call(double new_timer_stamp){
					timer_stamp = new_timer_stamp;
					counter++;
				}

				void end_call(double new_timer_stamp){
					timer += timer_stamp - new_timer_stamp;
					timer_stamp = 0.0;
				}
				
				int get_counter() const {
					return counter;
				}

				double get_timer() const {
					return timer;
				}

				std::string get_name() const {
					return myname;
				}
				
				bool is_me(std::string other_name){
					bool return_value = false;
					int aa = std::strcmp(myname.c_str(),other_name.c_str());
					if(aa == 0){
						return_value = true;
					}
					return return_value;
				}
		};	
	
		bool is_opened;
		std::ifstream *myfile;						/**< the file reader */
		int line_counter; 							/**< number of lines in the file */
		std::vector<my_function_node> myfunctions;	/**< vector of founded functions */

		bool is_there(std::vector<std::string> &mywords, std::string word) {
			bool return_value = false;
			for(int i=0;i<mywords.size();i++){
				int aa = std::strcmp(mywords[i].c_str(),word.c_str());
				if(aa == 0){
					return_value = true;
				}
			}
			return return_value;
		}


		bool is_func_begin(std::vector<std::string> &mywords) {
			return is_there(mywords, "FUNC_BEGIN");
		}

		bool is_func_end(std::vector<std::string> &mywords) {
			return is_there(mywords, "FUNC_END");
		}

		double get_time(std::vector<std::string> &mywords) {
			return std::stod(mywords[0]);
		}

		std::string get_name(std::vector<std::string> &mywords) {
			return mywords[mywords.size()-1];
		}

		void process_line(std::string &myline, std::string &delimiter){
			/* here will be stored words from this line */
			std::vector<std::string> mywords;

			/* split the line with delimiter */
			boost::split(mywords, myline, boost::is_any_of(delimiter), boost::token_compress_on);

			/* is this begin of func call? */
			if(is_func_begin(mywords)){
				/* maybe this function is already in founded functions */
				bool founded = false;
				for(int i=0;i<myfunctions.size();i++){
					if(myfunctions[i].is_me(get_name(mywords))){
						myfunctions[i].begin_call(get_time(mywords));
					}
				}
				
				/* no, it is not here, add it */
				if(!founded){
					my_function_node new_function_node(get_name(mywords));
					new_function_node.begin_call(get_time(mywords));
					myfunctions.push_back(new_function_node);
				}	
				
			}

			/* is this end of func call ? */
			if(is_func_end(mywords)){
				for(int i=0;i<myfunctions.size();i++){
					if(myfunctions[i].is_me(get_name(mywords))){
						myfunctions[i].end_call(get_time(mywords));
					}
				}
			}
			
		}
	
	public:
		log_processor(std::string &filename){
			myfile = new std::ifstream(filename);
			line_counter = 0;
			
			if (myfile->is_open()){ /* if the file can be opened */
				is_opened = true;
			} else {
				is_opened = false;
			}
		}

		log_processor(){
			if(is_opened){
				/* close the file after reading */
				myfile->close(); 
			}
		}
	
		bool get_is_opened(){
			return is_opened;
		}
	
		void process(std::string &delimiter){
			line_counter = 0;
			std::string myline; /* one loaded line from the file */
			while ( getline (*myfile,myline) ){ /* we are reading one line */
				/* process the line */
				process_line(myline, delimiter);

				/* increase counter of lines */
				line_counter++;
				if(line_counter > 500){
					break;
				}
			}
		}
			
		int get_line_counter(){
			return line_counter;
		}
		
		void print_results(){
			std::cout << "Number of processed lines = " << line_counter << std::endl;
			std::cout << std::endl;
			std::cout << "NAME" << std::setw(126) <<  "| COUNTERS" << std::setw(9)  <<  "| TIME" << std::setw(12) << std::endl;
			
			for(int i=0;i<myfunctions.size();i++){
				std::cout << std::setw(120) << myfunctions[i].get_name() << "| ";
				std::cout << std::setw(9) << myfunctions[i].get_counter() << "| ";
				std::cout << std::setw(16) << myfunctions[i].get_timer() << std::endl;
			}
			
		}
};

int main( int argc, char *argv[] )
{
	/* add local program options */
	boost::program_options::options_description opt_problem("UTIL_PROCESS_LOG", consoleArg.get_console_nmb_cols());
	opt_problem.add_options()
		("filename_in", boost::program_options::value< std::string >(), "input log file [string]")
		("delimiter", boost::program_options::value< std::string >(), "delimiter [string]");
	consoleArg.get_description()->add(opt_problem);

	/* call initialize */
	if(!Initialize<PetscVector>(argc, argv)){
		return 0;
	} 

	/* read filename from console */
	std::string filename_in;
	if(!consoleArg.set_option_value("filename_in", &filename_in)){
		std::cout << "ERROR: filename_in has to be set! Call application with parameter -h to see all parameters\n";
		return 0;
	}

	/* read delimiter from console */
	std::string delimiter;
	consoleArg.set_option_value("delimiter", &delimiter, LOG_SEPARATOR);

	/* give some informations about the settings */
	coutMaster << "- PROCESS LOG ----------------------------------------" << std::endl;
	coutMaster << " filename_in            = " << std::setw(30) << filename_in << std::endl;
	coutMaster << " delimiter              = " << std::setw(30) << delimiter << std::endl;
	coutMaster << " log_or_not_file_line   = " << std::setw(30) << print_bool(logging.get_log_or_not_file_line()) << std::endl;
	coutMaster << " log_or_not_func_call   = " << std::setw(30) << print_bool(logging.get_log_or_not_func_call()) << std::endl;
	coutMaster << " log_or_not_level       = " << std::setw(30) << print_bool(logging.get_log_or_not_level()) << std::endl;
	coutMaster << " log_or_not_memory      = " << std::setw(30) << print_bool(logging.get_log_or_not_memory()) << std::endl;
	coutMaster << "------------------------------------------------------" << std::endl;
	coutMaster << std::endl;

	/* create log processor */
	log_processor mylog_processor(filename_in);
	
	if(!mylog_processor.get_is_opened()){
		std::cout << "ERROR: file cannot be opened\n";
		return 0;
	}
	
	/* process the file */
	mylog_processor.process(delimiter);
	
	/* print results */
	mylog_processor.print_results();
	
	

	Finalize<PetscVector>();

	return 0;
}

