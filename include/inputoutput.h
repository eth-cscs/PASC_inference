#ifndef INPUTOUTPUT_H
#define	INPUTOUTPUT_H

#include "common.h"
#include <fstream>

class InputOutput {
	public:
		static void saveVTK(std::string name_of_file, DataVector<Scalar> data_vec, GammaVector<Scalar> gamma_vec, int dim, int T, int K);
};

#endif

