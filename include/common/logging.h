#ifndef PASC_COMMON_LOGGING_H
#define	PASC_COMMON_LOGGING_H

#include "common/globalmanager.h"
#include "common/memorycheck.h"
#include <fstream>

#define LOG_SEPARATOR "|"

#define DEFAULT_LOG_FUNC_CALL	false
#define DEFAULT_LOG_LEVEL		true
#define DEFAULT_LOG_MEMORY		true
#define DEFAULT_LOG_FILE_LINE	false

#define LOG_FUNC_BEGIN logging.begin_func(typeid(this).name(),__FUNCTION__,__FILE__,__LINE__);
#define LOG_FUNC_END logging.end_func(typeid(this).name(),__FUNCTION__,__FILE__,__LINE__);

#define LOG_FUNC_STATIC_BEGIN logging.begin_func("static",__FUNCTION__,__FILE__,__LINE__);
#define LOG_FUNC_STATIC_END logging.end_func("static",__FUNCTION__,__FILE__,__LINE__);

#define LOG_IT(it_num) logging.it(this->get_name(),__FILE__,__LINE__,it_num);
#define LOG_FX(fx_value) logging.fx(this->get_name(),__FILE__,__LINE__,fx_value);
#define LOG_FX2(fx_value,name_add) logging.fx(this->get_name(),__FILE__,__LINE__,fx_value,name_add);

#define LOG_DIRECT(my_string) logging.direct(my_string, __FILE__,__LINE__);


namespace pascinference {


class LoggingClass {
	private:
		std::string *filename;
		std::ofstream myfile;

		bool log_or_not;
		bool log_or_not_func_call;
		bool log_or_not_file_line;
		bool log_or_not_level;
		bool log_or_not_memory;
		
		int level;
		double reference_time;

		double getUnixTime(void){
			struct timespec tv;
			if(clock_gettime(CLOCK_REALTIME, &tv) != 0) return 0;
			return (((double) tv.tv_sec) + (double) (tv.tv_nsec / 1000000000.0));
		}
		
		void openfile(){
			myfile.open(filename->c_str(), std::fstream::in | std::fstream::out | std::fstream::app);
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

		void begin(std::string new_filename){
			this->filename = new std::string(new_filename);
			myfile.open(filename->c_str());

			log_or_not = true;
			log_or_not_file_line = DEFAULT_LOG_FILE_LINE;
			log_or_not_func_call = DEFAULT_LOG_FUNC_CALL;
			log_or_not_level = DEFAULT_LOG_LEVEL;
			log_or_not_memory = DEFAULT_LOG_MEMORY;
			
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
			myfile << "filename=" << filename << ",start time=" << reference_time;
			myfile << "\n";
			myfile.close();			
		};
		
		void end(){
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
				myfile << "\n";
				closefile();			
			}
			log_or_not = false;
		};

		void begin_func(std::string name_class,std::string name_function, std::string file, int line){
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
				myfile << "\n";
				closefile();			
			}
		};
		
		void end_func(std::string name_class,std::string name_function, std::string file, int line){
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
				myfile << "\n";
				closefile();
			}
			
			level--;
		};
		
		void it(std::string name_algorithm, std::string file, int line, int it){
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
				myfile << "\n";
				closefile();
			}
		};

		void fx(std::string name_algorithm, std::string file, int line, double fx_value){
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
				myfile << fx_value;
				myfile << "\n";
				closefile();
			}
		};

		void fx(std::string name_algorithm, std::string file, int line, double fx_value, std::string name_add){
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
				myfile << fx_value;
				myfile << "\n";
				closefile();
			}
		};

		template<class AnyType>
		void direct(AnyType my_string, std::string file, int line){
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
				myfile << my_string << "\n";
				closefile();
			}
		};

};

LoggingClass logging;


}


#endif
