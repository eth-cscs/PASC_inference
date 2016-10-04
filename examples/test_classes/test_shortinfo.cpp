/** @file test_shortinfo.cpp
 *  @brief test class and methods: ShortinfoClass and shortinfo
 *
 *  Test saving informations into shortinfo file.
 * 
 *  @author Lukas Pospisil
 */

#include "pascinference.h"

using namespace pascinference;

/* define testing functions in separated namespace */
namespace test_namespace {

class Testclass {
	private:
		int it;
		double fx;
		int step_counter;
		double computing_time;
	
	public:
		/* constructor */
		Testclass(){
			/* this algorithm doesn't do anything */
			it = 134;
			fx = -1e-6;
			step_counter = 56;
			computing_time = 1e-4;
		}
	
		/* append some values to shortinfo streams */
		void test_printshort(std::ostringstream &header, std::ostringstream &values) const {
			LOG_FUNC_STATIC_BEGIN
	
			/* write some values */
			header << "it, ";
			values << this->it << ", ";

			header << "fx, ";
			values << this->fx << ", ";

			header << "step_counter, ";
			values << this->step_counter << ", ";

			header << "t, ";
			values << this->computing_time << ", ";
	
			LOG_FUNC_STATIC_END
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

/* --- TEST SHORTINFO --- */
	/* prepare shortinfo file */
	std::string shortinfo_filename("shortinfo/test_shortinfo.txt");
	shortinfo.begin(shortinfo_filename);

	/* prepare output streams */
	std::ostringstream oss_short_output_values;
	std::ostringstream oss_short_output_header;
	
	/* call the function to append information to shortinfo streams */
	test_namespace::Testclass myinstance;
	myinstance.test_printshort(oss_short_output_header, oss_short_output_values);

	/* append end of line */
	oss_short_output_header << "\n";
	oss_short_output_values << "\n";

	/* write streams into shortinfo file */
	shortinfo.write(oss_short_output_header.str());
	shortinfo.write(oss_short_output_values.str());
			
	/* clear streams for next writing */
	oss_short_output_header.str("");
	oss_short_output_values.str("");	

/* --- FINALIZE LIBRARY --- */
	/* say bye */	
	coutMaster << "- end program, please check generated shortinfo file: " << shortinfo_filename << std::endl;

	/* call Finalize() to finalize Petsc if it was used */
	Finalize();

	return 0;
}
