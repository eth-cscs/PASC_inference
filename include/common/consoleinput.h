#ifndef PASC_COMMON_CONSOLEINPUT_H
#define	PASC_COMMON_CONSOLEINPUT_H

/* do deal with get console size */
#include <sys/ioctl.h>
#include <unistd.h>

/* load console parameters with boost */
#include <boost/program_options.hpp>

#include "options.h"

namespace pascinference {

class ConsoleArgClass {
	private:
		boost::program_options::options_description *description;
		boost::program_options::variables_map *vm;

	public:
		ConsoleArgClass() {
			/* get terminal size */
			struct winsize w;
			ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);

			description = new boost::program_options::options_description("PASC Inference Usage", w.ws_row);

			/* define command line options */
			description->add_options()
				("help,h", "Display this help message")
				("version,v", "Display the version number");	

			/* add other options from options.h */
			add_options(description);
		}
		
		bool init(int argc, char *argv[]){
			/* parse command line arguments */	
			vm = new boost::program_options::variables_map();
			boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(*description).run(), *vm);
			boost::program_options::notify(*vm);

			/* what to do with parsed arguments? */	
			if(vm->count("help")){
				std::cout << *description;
				return false;
			}

			if(vm->count("version")){
				std::cout << "not implemented yet" << std::endl;
				return false;
			}

			return true;
		}

		boost::program_options::options_description *get_description(){
			return description;
		}
		
		template<class Type>
		bool set_option_value(std::string name, Type *out_value){
			if(vm->count(name)){
				*out_value = (*vm)[name].as<Type>();
				return true;
			} else {
				return false;
			}
		}
		
		template<class Type, class Type2>
		void set_option_value(std::string name, Type *out_value, Type2 implicit_value){
			if(vm->count(name)){
				*out_value = (*vm)[name].as<Type>();
			} else {
				*out_value = implicit_value;
			}
		}		
};

static ConsoleArgClass consoleArg;

//TODO: define aliases
//typedef boost::program_options::value<int>() ConsoleInt;

} /* end of namespace */

#endif