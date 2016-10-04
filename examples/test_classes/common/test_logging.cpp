/** @file test_logging.cpp
 *  @brief test class and methods: LoggingClass and logging
 *
 *  Test the logging into log file.
 * 
 *  @author Lukas Pospisil
 */

#include "pascinference.h"

using namespace pascinference;

/* define testing functions in separated namespace */
namespace test_namespace {

/* call this function to test log of static function */
static void test_functiontobelogged(std::string something) {
	LOG_FUNC_STATIC_BEGIN
	
	/* do something */
	coutMaster << "performing really useless work: ";
	for(int i=0; i < 20; i++){
		for(int j=0; j < 1000; j++) {};
		coutMaster << something;
	}
	coutMaster << ": done" << std::endl;
	
	LOG_FUNC_STATIC_END
}

class Testclass {
	public:
		/* this function has to be provided, since it is called in LOG_IT as algorithm name */
		std::string get_name(){
			return "ReallyCoolAlgorithm";
		}
	
		/* to test return functions */
		int get_value() {
			LOG_FUNC_BEGIN
			
			/* some funny computation */
			int sum = 0;
			for(int i=0;i<10;i++) sum+=i;
			
			LOG_FUNC_END

			/* return statement has to be after the end of log since after return is not able to call anything */
			return sum;
		}


		/* a testing function, for testing member function call */
		void test_functiontobelogged(std::string something) {
			LOG_FUNC_BEGIN

			/* log function value */
			LOG_FX(3.1416)
			LOG_FX2(6.2932,"two")

			/* log directly something */
			LOG_DIRECT("I just log some function values. Just for fun.")
	
			/* do something */
			coutMaster << "performing another really useless work: ";
			for(int i=0; i < 20; i++){
				for(int j=0; j < 1000; j++) {};
				coutMaster << something;
			}
			coutMaster << ": done" << std::endl;
	
			LOG_FUNC_END
		}
		
		/* a testing function, for testing nested calls and number of it */
		void test_outerfunctiontobelogged(){
			LOG_FUNC_BEGIN
			
			/* log number of iterations (suppose that this class is algorithm) */
			LOG_IT(this->get_value());
			
			/* call first function to test level of log */
			this->test_functiontobelogged("=");
			
			LOG_FUNC_END
		}
}; /* end of class declaration */

} /* end of namespace */

int main( int argc, char *argv[] ){

/* --- INITIALIZE LIBRARY --- */
	/* call initialize to initialize the library (this method also initializes Petsc if it is used) */
	if(!Initialize(argc, argv)){
		/* this happen for example when program is called with "--help" parameter */
		return 0;
	} 

/* --- TEST LOGGING --- */
	/* start to log into logfile, name of file based on proc number, each process has its own log file */
	std::ostringstream oss;
	oss << "log/test_log_p" << GlobalManager.get_rank() << ".txt";
	coutAll << "Test to log into file: " << oss.str() << std::endl;
	coutAll.synchronize();

	/* the following method sets the filename of log file and starts to log into this file */
	logging.begin(oss.str());

	/* turn of logging function calls */
	logging.set_log_or_not_func_call(true);
	logging.set_log_or_not_file_line(true);
	logging.set_log_or_not_level(true);
	logging.set_log_or_not_memory(true);

	/* log directly something */
	LOG_DIRECT("this is a message to log file - log it!")

	/* call the static function to be logged */
	test_namespace::test_functiontobelogged("*");

	/* call the member function to be logged */
	test_namespace::Testclass myinstance;
	myinstance.test_functiontobelogged("-");

	/* call the member function which calls another function to test level */
	myinstance.test_outerfunctiontobelogged();

	/* stop to log */
	logging.end();

/* --- FINALIZE LIBRARY --- */
	/* say bye */	
	coutMaster << "- end program, please check generated log files" << std::endl;

	/* call Finalize() to finalize Petsc if it was used */
	Finalize();

	return 0;
}
