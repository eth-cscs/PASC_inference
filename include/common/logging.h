#ifndef PASC_COMMON_LOGGING_H
#define	PASC_COMMON_LOGGING_H

#include "common/globalmanager.h"
#include <fstream>

#define LOG_SEPARATOR "|"

#define LOG_FUNC_BEGIN logging.begin_func(typeid(this).name(),__FUNCTION__,__FILE__,__LINE__);
#define LOG_FUNC_END logging.end_func(typeid(this).name(),__FUNCTION__,__FILE__,__LINE__);

#define LOG_FUNC_STATIC_BEGIN logging.begin_func("static",__FUNCTION__,__FILE__,__LINE__);
#define LOG_FUNC_STATIC_END logging.end_func("static",__FUNCTION__,__FILE__,__LINE__);

namespace pascinference {


class LoggingClass {
	private:
		std::string filename;
		std::ofstream myfile;
		bool log_or_not;
		bool log_or_not_file;

		double reference_time;

		double getUnixTime(void){
			struct timespec tv;
			if(clock_gettime(CLOCK_REALTIME, &tv) != 0) return 0;
			return (((double) tv.tv_sec) + (double) (tv.tv_nsec / 1000000000.0));
		}
		
		void openfile(){
			myfile.open(filename.c_str(), std::fstream::in | std::fstream::out | std::fstream::app);
		}
		
		void closefile(){
			myfile.close();
		}

	public:
		LoggingClass(){
			log_or_not = false;
		};

		~LoggingClass(){
	
		};

		void on(){
			log_or_not = true;
		};
	
		void off(){
			log_or_not = false;
		};
		
		bool get_log_or_not(){
			return log_or_not;
		};

		void begin(std::string new_filename){
			filename = new_filename;
			myfile.open(filename.c_str());

			log_or_not = true;
			log_or_not_file = false; // TODO: this could be turned on
			
			reference_time = getUnixTime();

			myfile << getUnixTime()-reference_time << LOG_SEPARATOR;
			if(log_or_not_file){
				myfile << __FILE__ << LOG_SEPARATOR;
				myfile << __LINE__ << LOG_SEPARATOR;
			}
			myfile << "LOG_OPEN" << LOG_SEPARATOR;
			myfile << "filename=" << filename << ",start time=" << reference_time;
			myfile << "\n";
			closefile();			
		};
		
		void end(){
			if(log_or_not){
				openfile();
				myfile << getUnixTime()-reference_time << LOG_SEPARATOR;
				if(log_or_not_file){
					myfile << __FILE__ << LOG_SEPARATOR;
					myfile << __LINE__ << LOG_SEPARATOR;
				}
				myfile << "LOG_CLOSE" << LOG_SEPARATOR;
				myfile << "filename=" << filename;
				myfile << "\n";
				closefile();			
			}
			log_or_not = false;
		};

		void begin_func(std::string name_class,std::string name_function, std::string file, int line){
			if(log_or_not){
				openfile();
				myfile << getUnixTime()-reference_time << LOG_SEPARATOR;
				if(log_or_not_file){
					myfile << file << LOG_SEPARATOR;
					myfile << line << LOG_SEPARATOR;
				}
				myfile << "FUNC_BEGIN" << LOG_SEPARATOR;
				myfile << name_class << "::" << name_function;
				myfile << "\n";
				closefile();			
			}
		};
		
		void end_func(std::string name_class,std::string name_function, std::string file, int line){
			if(log_or_not){
				openfile();
				myfile << getUnixTime()-reference_time << LOG_SEPARATOR;
				if(log_or_not_file){
					myfile << file << LOG_SEPARATOR;
					myfile << line << LOG_SEPARATOR;
				}
				myfile << "FUNC_END" << LOG_SEPARATOR;
				myfile << name_class << "::" << name_function;
				myfile << "\n";
				closefile();
			}
		};
		
	
};

LoggingClass logging;


}


#endif
