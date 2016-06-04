#ifndef PASC_COMMON_CONSOLEINPUT_H
#define	PASC_COMMON_CONSOLEINPUT_H

#include <boost/program_options.hpp>
#include "options.h"

namespace pascinference {

class ConsoleArgClass {
	private:
		boost::program_options::options_description *description;
		boost::program_options::variables_map *vm;

	public:
		ConsoleArgClass() {
			description = new boost::program_options::options_description("PASC Inference Usage");

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

//	if(vm.count("debug")){// TODO: this can be included in global application
//		DEBUG_MODE = vm["debug"].as<int>(); /* set global variable */
//	}

//	if(vm.count("length")){
//		T = vm["length"].as<int>(); /* set global variable */
//	}

//	if(vm.count("clusters")){
//		K = vm["clusters"].as<int>(); /* set global variable */
//	}
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
		
		template<class Type>
		void set_option_value(std::string name, Type *out_value, Type implicit_value){
			if(vm->count(name)){
				*out_value = (*vm)[name].as<Type>();
			} else {
				*out_value = implicit_value;
			}
		}		
};

ConsoleArgClass consoleArg;

//TODO: define aliases
//typedef boost::program_options::value<int>() ConsoleInt;

} /* end of namespace */

#endif
