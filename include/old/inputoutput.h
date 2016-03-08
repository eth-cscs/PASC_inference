#ifndef INPUTOUTPUT_H
#define	INPUTOUTPUT_H

#include "common.h"
#include <fstream>

namespace pascinference {

class InputOutput {
	public:
		static void saveVTK(std::string name_of_file, DataVector data_vec, GammaVector gamma_vec, int dim, int T, int K);
};

} /* end of namespace */


#endif

