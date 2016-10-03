/** @file test_kmeans_simplify.cpp
 *  @brief simplify results by taking only a few records from CSV
 *
 *  @author Lukas Pospisil
 */

#include <iostream>
#include <list>
#include <algorithm>

#include "pascinference.h"

using namespace pascinference;

extern int pascinference::DEBUG_MODE;

int main( int argc, char *argv[] )
{
	/* add local program options */
	boost::program_options::options_description opt_problem("PROBLEM EXAMPLE", consoleArg.get_console_nmb_cols());
	opt_problem.add_options()
		("test_input_filename", boost::program_options::value< std::string >(), "name of file with original CSV data [string]")
		("test_output_filename", boost::program_options::value< std::string >(), "name of file with new CSV data [string]")
		("test_T", boost::program_options::value< int >(), "only every T-th record is taken [int]");
	consoleArg.get_description()->add(opt_problem);

	/* call initialize */
	if(!Initialize(argc, argv)){
		return 0;
	} 

	if(GlobalManager.get_size() > 1){
		coutMaster << "The simplificator works only on one processor, sorry.\n";
		return 0;		
	}

	/* get console parameters */
	int T;
	std::string input_filename;
	std::string output_filename;

	consoleArg.set_option_value("test_input_filename", &input_filename, "results/test_kmeans.csv");
	consoleArg.set_option_value("test_output_filename", &output_filename, "results/test_kmeans_simple.csv");
	consoleArg.set_option_value("test_T", &T, 100000);

	/* print info about what we will compute */
	coutMaster << "- PROBLEM INFO --------------------------------------------------\n";
	coutMaster << " T      = " << std::setw(40) << T << "\n";
	coutMaster << " input  = " << std::setw(40) << input_filename << "\n";
	coutMaster << " output = " << std::setw(40) << output_filename << "\n";

	coutMaster << "------------------------------------------------------------------\n";
	
	/* say hello */
	coutMaster << "- start program\n";

	/* input file */
	std::string line;
	std::ifstream input_file(input_filename);


	if (!input_file.is_open()){
		coutMaster << "Unable to open input file.\n"; 	
		return 1;
	}
	
	/* output file */
	std::ofstream output_file;
	output_file.open(output_filename);

	/* header */
	getline (input_file,line);
	output_file << "t," << line << "\n";

	/* now read every line and T-th line write into new file */
	int T_counter = 0;
	int T_global_couter = 0;
	while ( getline (input_file,line) ){
		if(T_counter == T){
			output_file << T_global_couter << "," << line << "\n";
			T_counter = 0;
		}
		T_counter++;
		T_global_couter++;
	}
	
	/* close files */
	input_file.close();
	output_file.close();

	/* say bye */	
	coutMaster << "- end program\n";

	Finalize();

	return 0;
}

