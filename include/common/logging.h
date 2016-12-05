/** @file logging.h
 *  @brief Log files.
 *
 *  @author Lukas Pospisil
 */

#ifndef PASC_COMMON_LOGGING_H
#define	PASC_COMMON_LOGGING_H

#include "common/globalmanager.h"
#include "common/memorycheck.h"
#include <fstream>

/** the separator used in log file */
#define LOG_SEPARATOR "|"

/** default value of logging function calls */
#define DEFAULT_LOG_FUNC_CALL	false 

/** default value of log also the level of called function */
#define DEFAULT_LOG_LEVEL		true

/** default value of log also the state of the memory */
#define DEFAULT_LOG_MEMORY		false

/** default value of log also the file and line of called log function */
#define DEFAULT_LOG_FILE_LINE	false

/** macro used to log begin of the function */
#define LOG_FUNC_BEGIN logging.begin_func(typeid(this).name(),__FUNCTION__,__FILE__,__LINE__);

/** macro used to log end of the function */
#define LOG_FUNC_END logging.end_func(typeid(this).name(),__FUNCTION__,__FILE__,__LINE__);

/** macro used to log begin of the static function */
#define LOG_FUNC_STATIC_BEGIN logging.begin_func("static",__FUNCTION__,__FILE__,__LINE__);

/** macro used to log end of the static function */
#define LOG_FUNC_STATIC_END logging.end_func("static",__FUNCTION__,__FILE__,__LINE__);

/** macro used to log number of iterations, this macro uses this->get_name() as algorithm name */
#define LOG_IT(it_num) logging.it(this->get_name(),__FILE__,__LINE__,it_num);

/** macro used to log function value, this macro uses this->get_name() as algorithm name */
#define LOG_FX(fx_value) logging.fx(this->get_name(),__FILE__,__LINE__,fx_value);

/** macro used to log function value with name of function, this macro uses this->get_name() as algorithm name */
#define LOG_FX2(fx_value,name_add) logging.fx(this->get_name(),__FILE__,__LINE__,fx_value,name_add);

/** macro used to log directly the message represented by string */
#define LOG_DIRECT(my_string) logging.direct(my_string, __FILE__,__LINE__);


namespace pascinference {
namespace common {

/** \class LoggingClass
 *  \brief Manipulation with log file.
 *
 *  Several macros implemented to call logging methods.
 *  Can be used to store the sequence of function calls and/or the progress of algorithms - the function value, stoppig criteria, number of iterations, etc.
 * 
*/
class LoggingClass {
	private:
		std::string *filename;			/**< name of log file */
		std::ofstream myfile;			/**< log file */

		bool log_or_not;				/**< logging (writting into file) is turned on/off */
		bool log_or_not_func_call;		/**< log LOG_FUNC_(STATIC)_BEGIN/LOG_FUNC_(STATIC)_END */
		bool log_or_not_file_line;		/**< log also the file and line of called log function */
		bool log_or_not_level;			/**< log also the level of called function */
		bool log_or_not_memory;			/**< log also the state of the memory */
		
		int level;						/**< level of called function, each LOG_FUNC_(STATIC)_BEGIN increases this counter, LOG_FUNC_(STATIC)_END decreases this counter */
		double reference_time;			/**< used to store start of unix timer */

		/** @brief get actual time in unix format
		 */
		 double getUnixTime(void){
//			struct timespec tv;
//			if(clock_gettime(CLOCK_REALTIME, &tv) != 0) return 0;
//			return (((double) tv.tv_sec) + (double) (tv.tv_nsec / 1000000000.0));
//			MPI_Barrier(MPI_COMM_WORLD);
			return MPI_Wtime();
		}

		/** @brief open log file to append content
		 */
		void openfile(){
			myfile.open(filename->c_str(), std::fstream::in | std::fstream::out | std::fstream::app);
		}
		
		/** @brief close log file
		 */
		void closefile(){
			myfile.close();
		}

	public:
		/** @brief default constructor
		 * 
		 * Logging is turned off until begin() is called.
		 * 
		 */
		LoggingClass(){
			log_or_not = false;
		}

		/** @brief default destructor
		 */ 
		~LoggingClass(){
	
		}

		/** @brief begin to log into log file
		 * 
		 * Set initial values, write "LOG_OPEN" into log file.
		 * 
		 * @param new_filename name of log file
		 */
		void begin(std::string new_filename){
			this->filename = new std::string(new_filename);
			myfile.open(filename->c_str());

			consoleArg.set_option_value("log_or_not", &log_or_not, true);
			consoleArg.set_option_value("log_or_not_file_line", &log_or_not_file_line, DEFAULT_LOG_FILE_LINE);
			consoleArg.set_option_value("log_or_not_func_call", &log_or_not_func_call, DEFAULT_LOG_FUNC_CALL);
			consoleArg.set_option_value("log_or_not_level", &log_or_not_level, DEFAULT_LOG_LEVEL);
			consoleArg.set_option_value("log_or_not_memory", &log_or_not_memory, DEFAULT_LOG_MEMORY);

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

		/** @brief stop to log into file
		 * 
		 * Write "LOG_CLOSE" into log file.
		 * 
		 */
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
				myfile << std::endl;
				closefile();			
			}
			log_or_not = false;
		}

		/** @brief log the begin of called function
		 * 
		 * Increase inner level counter.
		 * 
		 * @param name_class name of class from where the function was called
		 * @param name_function name of called function
		 * @param file name of source file from where the function code is located
		 * @param line number of line in source file where the function is located
		 */
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
				myfile << std::endl;
				closefile();			
			}
		}
		
		/** @brief log the end of called function
		 * 
		 * Decrease inner level counter.
		 * 
		 * @param name_class name of class from where the function was called
		 * @param name_function name of called function
		 * @param file name of source file from where the function code is located
		 * @param line number of line in source file where the function is located
		 */
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
				myfile << std::endl;
				closefile();
			}
			
			level--;

		}

		/** @brief log the number of iterations
		 * 
		 * @param name_algorithm name of algorithm which performed provided number of iterations 
		 * @param file name of source file from where the function code is located
		 * @param line number of line in source file where the function is located
		 * @param it number of iterations to be stored in log file
		 */
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
				myfile << std::endl;
				closefile();
			}
		}

		/** @brief log the function value
		 * 
		 * @param name_algorithm name of algorithm which is connected with this value 
		 * @param file name of source file from where the function code is located
		 * @param line number of line in source file where the function is located
		 * @param fx_value the value of function
		 */
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
				
				std::streamsize ss = std::cout.precision();
				myfile << std::setprecision(17) << fx_value << std::setprecision(ss);
				
				
				myfile << std::endl;
				closefile();
			}
		}

		/** @brief log the function value with additional name
		 * 
		 * @param name_algorithm name of algorithm which is connected with this value 
		 * @param file name of source file from where the function code is located
		 * @param line number of line in source file where the function is located
		 * @param fx_value the value of function
		 * @param name_add the additional string to name_algorithm
		 */
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

				std::streamsize ss = std::cout.precision();
				myfile << std::setprecision(17) << fx_value << std::setprecision(ss);
				
				
				myfile << std::endl;
				closefile();
			}
		}

		/** @brief write directly to log file
		 * 
		 * @param my_string anything which could be stored using << operator
		 * @param file name of source file from where this method was called
		 * @param line number of line in source file from where this method was called
		 */
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
				myfile << my_string << std::endl;
				closefile();
			}
		}

		/** @brief set logging (writting into file) on/off
		 *
		 * @param new_value new value of member variable log_or_not
		 */
		void set_log_or_not(bool new_value){
			log_or_not = new_value;
		}				

		/** @brief set logging LOG_FUNC_(STATIC)_BEGIN/LOG_FUNC_(STATIC)_END on/off
		 *
		 * @param new_value new value of member variable log_or_not_func_call
		 */
		void set_log_or_not_func_call(bool new_value){
			log_or_not_func_call = new_value;
		}				

		/** @brief set log also the file and line of called log function on/off
		 *
		 * @param new_value new value of member variable log_or_not_file_line
		 */
		void set_log_or_not_file_line(bool new_value){
			log_or_not_file_line = new_value;
		}				

		/** @brief set log also the level of called function on/off
		 *
		 * @param new_value new value of member variable log_or_not_level
		 */
		void set_log_or_not_level(bool new_value){
			log_or_not_level = new_value;
		}				

		/** @brief set log also the state of the memory on/off
		 *
		 * @param new_value new value of member variable log_or_not_memory
		 */
		void set_log_or_not_memory(bool new_value){
			log_or_not_memory = new_value;
		}				

};

LoggingClass logging;	/**< global instance of logging class */


}
} /* end of namespace */


#endif
