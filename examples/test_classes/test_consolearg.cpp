/** @file test_consolearg.cpp
 *  @brief test class and methods: ConsoleArg 
 *
 *  This example (as well as whole class) is based on boost::program_options.
 * 
 *  @author Lukas Pospisil
 */

#include "pascinference.h"

using namespace pascinference;

int main( int argc, char *argv[] ){

/* --- ADD LOCAL PROGRAM OPTIONS --- */
	/* create new boost category, run program with "--help" to see them */
	/* here "consoleArg" is a global instance of ConsoleArgClass */
	/* get_console_nmb_cols() returns the actual size of the console to provide the list of arguments in as many columns as possible (set maximum possible width of the list) */
	/* otherwise the list is not through whole console */
	boost::program_options::options_description opt_problem("TEST of ConsoleArgClass", consoleArg.get_console_nmb_cols());

	/* add new program options to the list */
	opt_problem.add_options()
		("test_integer_value", boost::program_options::value<int>(), "this is testing integer value [int]")
		("test_double_value", boost::program_options::value<double>(), "this is testing double value [double]")
		("test_bool_value", boost::program_options::value<bool>(), "this is testing bool value [bool]")
		("test_string_value", boost::program_options::value< std::string >(), "this is testing string value [string]")
		("test_double_more", boost::program_options::value<std::vector<double> >()->multitoken(), "this is testing double multivalue [double]");

	/* append new category of console parameters to the consoleArg (which implicitly contains all console parameters of library) */
	consoleArg.get_description()->add(opt_problem);


/* --- INITIALIZE LIBRARY --- */
	/* call initialize to initialize the library (this method also initializes Petsc if it is used) */
	/* this method also load all pre-defined console parameters of the library */
	if(!Initialize(argc, argv)){
		return 0;
	} 

/* --- LOAD PARAMETERS FROM CONSOLE --- */
	/* define variables for load */
	int test_integer_value; /* here will be integer value loaded */
	double test_double_value; /* here will be double value loaded */
	bool test_bool_value; /* here will be bool value loaded */
	std::string test_string_value; /* here will be string value loaded */
	std::vector<double> test_double_more; /* 

	/* load options from console parameters */
	/* the thirth argument represents the default value of the variable (if it is not set in console arguments) */
	consoleArg.set_option_value("test_integer_value", &test_integer_value, 0);
	consoleArg.set_option_value("test_double_value", &test_double_value, -1.3);
	consoleArg.set_option_value("test_bool_value", &test_bool_value, false);
	consoleArg.set_option_value("test_string_value", &test_string_value, "this is only the test");

	/* if there is given any "--test_double_more=..", then load the list */
	bool given_test_double_more = false; /* are there any "--test_double_more" ? */
	if(consoleArg.set_option_value("test_double_more", &test_double_more)){
		given_test_double_more = true;
	}

/* --- PRINT LOADED OPTIONS --- */
	/* coutMaster denotes overloaded std::cout which prints only by master (i.e. process with rank=0) */
	/* std::setw(int) is used to set fixed length of following variable output */
	coutMaster << "- LOADED ARGUMENTS: ----------------------------" << std::endl;
	coutMaster << " test_integer_value         = " << std::setw(30) << test_integer_value << " (this is testing integer value)\n";
	coutMaster << " test_double_value          = " << std::setw(30) << test_double_value << " (this is testing double value)\n";
	coutMaster << " test_bool_value            = " << std::setw(30) << test_bool_value << " (this is testing bool value)\n";
	coutMaster << " test_string_value          = " << std::setw(30) << test_string_value << " (this is testing string value)\n";

	coutMaster << " test_double_more           = " << std::setw(30);
	if(given_test_double_more){
		/* here I will use our function print_vector implemented in "arrayoperation.h" */
		coutMaster << print_vector(test_double_more);
	} else {
		/* if "--test_double_more" was not provided in console parameters, then print following */
		coutMaster << "not provided";
	}
	coutMaster << "\n";


/* --- FINALIZE LIBRARY --- */
	/* say bye */	
	coutMaster << "- end program" << std::endl;

	/* call Finalize() to finalize Petsc if it was used */
	Finalize();

	return 0;
}
