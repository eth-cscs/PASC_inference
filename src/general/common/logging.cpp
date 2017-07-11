#include "general/common/logging.h"

namespace pascinference {
namespace common {

double LoggingClass::getUnixTime(void){
//			struct timespec tv;
//			if(clock_gettime(CLOCK_REALTIME, &tv) != 0) return 0;
//			return (((double) tv.tv_sec) + (double) (tv.tv_nsec / 1000000000.0));
//			MPI_Barrier(MPI_COMM_WORLD);
	return MPI_Wtime();
}

void LoggingClass::openfile(){
	myfile.open(filename->c_str(), std::fstream::in | std::fstream::out | std::fstream::app);
	myfile << std::setprecision(17);
}

void LoggingClass::closefile(){
	myfile.close();
}

LoggingClass::LoggingClass(){
	log_or_not = false;
}

LoggingClass::~LoggingClass(){
}

void LoggingClass::begin(std::string new_filename){
	this->filename = new std::string(new_filename);
	myfile.open(filename->c_str());

	consoleArg.set_option_value("log_or_not", &log_or_not, true);
	consoleArg.set_option_value("log_or_not_file_line", &log_or_not_file_line, DEFAULT_LOG_FILE_LINE);
	consoleArg.set_option_value("log_or_not_func_call", &log_or_not_func_call, DEFAULT_LOG_FUNC_CALL);
	consoleArg.set_option_value("log_or_not_level", &log_or_not_level, DEFAULT_LOG_LEVEL);
	consoleArg.set_option_value("log_or_not_memory", &log_or_not_memory, DEFAULT_LOG_MEMORY);

	/* only master is logging */
	if(GlobalManager.get_rank() != 0){
        log_or_not = false;
	}

	level = -1;
	reference_time = getUnixTime();

	myfile << getUnixTime()-reference_time << LOG_SEPARATOR;
	if(log_or_not_level){
		myfile << level << LOG_SEPARATOR;
	}
	if(log_or_not_memory){
		myfile << MemoryCheck::get_virtual() << LOG_SEPARATOR;
	}
	if(log_or_not_file_line){
		myfile << __FILE__ << LOG_SEPARATOR;
		myfile << __LINE__ << LOG_SEPARATOR;
	}
	myfile << "LOG_OPEN" << LOG_SEPARATOR;
	myfile << "filename=" << *filename << ",start time=" << reference_time;
	myfile << std::endl;
	myfile.close();
}

void LoggingClass::end(){
	if(log_or_not){
		openfile();
		myfile << getUnixTime()-reference_time << LOG_SEPARATOR;
		if(log_or_not_level){
			myfile << "-1" << LOG_SEPARATOR;
		}
		if(log_or_not_memory){
			myfile << MemoryCheck::get_virtual() << LOG_SEPARATOR;
		}
		if(log_or_not_file_line){
			myfile << __FILE__ << LOG_SEPARATOR;
			myfile << __LINE__ << LOG_SEPARATOR;
		}
		myfile << "LOG_CLOSE" << LOG_SEPARATOR;
		myfile << "filename=" << *filename;
		myfile << std::endl;
		closefile();
	}
	log_or_not = false;
}

void LoggingClass::begin_func(std::string name_class,std::string name_function, std::string file, int line){
	level++;

	if(log_or_not && log_or_not_func_call){
		openfile();
		myfile << getUnixTime()-reference_time << LOG_SEPARATOR;


		if(log_or_not_level){
			myfile << level << LOG_SEPARATOR;
		}
		if(log_or_not_memory){
			myfile << MemoryCheck::get_virtual() << LOG_SEPARATOR;
		}
		if(log_or_not_file_line){
			myfile << file << LOG_SEPARATOR;
			myfile << line << LOG_SEPARATOR;
		}
		myfile << "FUNC_BEGIN" << LOG_SEPARATOR;
		myfile << name_class << "::" << name_function;
		myfile << std::endl;
		closefile();
	}
}

void LoggingClass::end_func(std::string name_class,std::string name_function, std::string file, int line){
	if(log_or_not && log_or_not_func_call){
		openfile();
		myfile << getUnixTime()-reference_time << LOG_SEPARATOR;
		if(log_or_not_level){
			myfile << level << LOG_SEPARATOR;
		}
		if(log_or_not_memory){
			myfile << MemoryCheck::get_virtual() << LOG_SEPARATOR;
		}
		if(log_or_not_file_line){
			myfile << file << LOG_SEPARATOR;
			myfile << line << LOG_SEPARATOR;
		}
		myfile << "FUNC_END" << LOG_SEPARATOR;
		myfile << name_class << "::" << name_function;
		myfile << std::endl;
		closefile();
	}

	level--;

}

void LoggingClass::it(std::string name_algorithm, std::string file, int line, int it){
	if(log_or_not){
		openfile();
		myfile << getUnixTime()-reference_time << LOG_SEPARATOR;
		if(log_or_not_level){
			myfile << level << LOG_SEPARATOR;
		}
		if(log_or_not_memory){
			myfile << MemoryCheck::get_virtual() << LOG_SEPARATOR;
		}
		if(log_or_not_file_line){
			myfile << file << LOG_SEPARATOR;
			myfile << line << LOG_SEPARATOR;
		}
		myfile << "IT_" << name_algorithm << LOG_SEPARATOR;
		myfile << it;
		myfile << std::endl;
		closefile();
	}
}

void LoggingClass::fx(std::string name_algorithm, std::string file, int line, double fx_value){
	if(log_or_not){
		openfile();
		myfile << getUnixTime()-reference_time << LOG_SEPARATOR;
		if(log_or_not_level){
			myfile << level << LOG_SEPARATOR;
		}
		if(log_or_not_memory){
			myfile << MemoryCheck::get_virtual() << LOG_SEPARATOR;
		}
		if(log_or_not_file_line){
			myfile << file << LOG_SEPARATOR;
			myfile << line << LOG_SEPARATOR;
		}
		myfile << "FX_" << name_algorithm << LOG_SEPARATOR;

		std::streamsize ss = std::cout.precision();
		myfile << std::setprecision(17) << fx_value << std::setprecision(ss);


		myfile << std::endl;
		closefile();
	}
}

void LoggingClass::fx(std::string name_algorithm, std::string file, int line, double fx_value, std::string name_add){
	if(log_or_not){
		openfile();
		myfile << getUnixTime()-reference_time << LOG_SEPARATOR;
		if(log_or_not_level){
			myfile << level << LOG_SEPARATOR;
		}
		if(log_or_not_memory){
			myfile << MemoryCheck::get_virtual() << LOG_SEPARATOR;
		}
		if(log_or_not_file_line){
			myfile << file << LOG_SEPARATOR;
			myfile << line << LOG_SEPARATOR;
		}
		myfile << "FX_" << name_algorithm << "_" << name_add << LOG_SEPARATOR;

		std::streamsize ss = std::cout.precision();
		myfile << std::setprecision(17) << fx_value << std::setprecision(ss);

		myfile << std::endl;
		closefile();
	}
}

template<class AnyType>
void LoggingClass::direct(AnyType my_string, std::string file, int line){
	if(log_or_not){
		openfile();
		myfile << getUnixTime()-reference_time << LOG_SEPARATOR;
		if(log_or_not_level){
			myfile << level << LOG_SEPARATOR;
		}
		if(log_or_not_memory){
			myfile << MemoryCheck::get_virtual() << LOG_SEPARATOR;
		}
		if(log_or_not_file_line){
			myfile << file << LOG_SEPARATOR;
			myfile << line << LOG_SEPARATOR;
		}
		myfile << my_string << std::endl;
		closefile();
	}
}

void LoggingClass::set_log_or_not(bool new_value){
	log_or_not = new_value;
}

bool LoggingClass::get_log_or_not() const {
	return log_or_not;
}

void LoggingClass::set_log_or_not_func_call(bool new_value){
	log_or_not_func_call = new_value;
}

bool LoggingClass::get_log_or_not_func_call() const {
	return log_or_not_func_call;
}

void LoggingClass::set_log_or_not_file_line(bool new_value){
	log_or_not_file_line = new_value;
}

bool LoggingClass::get_log_or_not_file_line() const {
	return log_or_not_file_line;
}

void LoggingClass::set_log_or_not_level(bool new_value){
	log_or_not_level = new_value;
}

bool LoggingClass::get_log_or_not_level() const {
	return log_or_not_level;
}

void LoggingClass::set_log_or_not_memory(bool new_value){
	log_or_not_memory = new_value;
}

bool LoggingClass::get_log_or_not_memory() const {
	return log_or_not_memory;
}

LoggingClass logging;	/**< global instance of logging class */

}
} /* end of namespace */

