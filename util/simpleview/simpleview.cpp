#include <math.h>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>

#include <vector>
#include <cmath>
#include <boost/tuple/tuple.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include "gnuplot-iostream.h"

#define PI 3.1416

#include "common.h"
#include "sample.h"
#include "showgamma.h"

namespace po = boost::program_options;

int main( int argc, char *argv[] )
{
	/* parse input arguments */
	
	po::options_description desc("Allowed options");
	desc.add_options()
		("help", "produce help message")
		("gammafile", po::value<std::string>(), "set gammafile in binary format");

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);    

	if (vm.count("help")) {
		std::cout << desc << std::endl;
		return 1;
	}

	if (vm.count("gammafile")) {
		std::string filename = vm["gammafile"].as<std::string>();
		std::cout << "Gammafile: " << filename << std::endl;
		show_gamma("sample.txt");

	}	


//	generate_sample("sample.txt", 1000, 4);

	return 0;
}

