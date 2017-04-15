/** @file consoleinput.h
 *  @brief Manipulation with console arguments and parameters.
 *
 *  Based on boost::program_options.
 * 
 *  @author Lukas Pospisil
 */

#ifndef PASC_COMMON_CONSOLEINPUT_H
#define	PASC_COMMON_CONSOLEINPUT_H

/* do deal with get console size */
#include <sys/ioctl.h>
#include <unistd.h>
#include <iostream>

#include "external/boost/options.h"

namespace pascinference {
namespace common {

/** @class ConsoleArgClass
 *  @brief for manipulation with console arguments
 * 
 *  Based on boost::program_options.
 * 
 */ 
class ConsoleArgClass {
	private:
		boost::program_options::options_description *description; 	/**< here all the console options (including default library parameters) are stored */
		boost::program_options::variables_map *vm;					/**< used for parsing console arguments */

		int console_nmb_cols;	/**< number of columns of console */
	public:
		/** @brief default constructor
		*
		*   Set number of columns of console.
		*   Set and load the default console parameters of library.
		*/
		ConsoleArgClass() {
			/* get terminal size */
			struct winsize w;
			ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
			console_nmb_cols = w.ws_col;
			
			description = new boost::program_options::options_description("PASC Inference Usage", console_nmb_cols);

			/* define command line options */
			description->add_options()
				("help,h", "Display this help message")
				("version,v", "Display the version number");	

			/* add other options from options.h */
			add_options(description, console_nmb_cols);
		}
		
		/** @brief get the number of columns in console
		*
		*   Can be used to print options through whole console.
		* 
		*  @return number of columns of console
		*/
		int get_console_nmb_cols(){
			return console_nmb_cols;
		}

		/** @brief initialize the option reader
		*
		*  Arguments argc and argv are usually taken from arguments of main function.
		* 
		*  @param argc number of console parameters
		*  @param argv console parameters
		*  @return success of initialisation
		*/
		bool init(int argc, char *argv[]){
			/* parse command line arguments */	
			vm = new boost::program_options::variables_map();
			boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(*description).run(), *vm);
			boost::program_options::notify(*vm);

			/* what to do with parsed arguments? */	
			if(vm->count("help")){
				std::cout << *description << std::endl;
				return false;
			}

			if(vm->count("version")){
				std::cout << "not implemented yet" << std::endl;
				return false;
			}

			return true;
		}

		/** @brief get the description content to add more options
		*
		*  @return the inner instance of description
		*/
		boost::program_options::options_description *get_description(){
			return description;
		}

		/** @brief get the value of option
		*
		*  @param name name of option in string format
		*  @param out_value there will be stored output value of option
		*  @return if option is not provided, then false
		*/
		template<class Type>
		bool set_option_value(std::string name, Type *out_value){
			if(vm->count(name)){
				*out_value = (*vm)[name].as<Type>();
				return true;
			} else {
				return false;
			}
		}

		/** @brief get the value of option or implicit value
		*
		*  @param name name of option in string format
		*  @param out_value there will be stored output value of option
		*  @param implicit_value if the option is not set in console parameters, this value will be used
		*/
		template<class Type, class Type2>
		void set_option_value(std::string name, Type *out_value, Type2 implicit_value){
			if(vm->count(name)){
				*out_value = (*vm)[name].as<Type>();
			} else {
				*out_value = implicit_value;
			}
		}		
};

extern ConsoleArgClass consoleArg; /**< global instance to manipulate with console input arguments */


}
} /* end of namespace */

#endif
