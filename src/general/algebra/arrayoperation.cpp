#include "general/algebra/arrayoperation.h"

namespace pascinference {
namespace algebra {

bool parse_strings_to_doubles(int K, int Km, std::vector<std::string> Theta_list, double *Theta_solution){
	std::string token;
	size_t pos;
	int counter=0;

	for(int k=0;k < K;k++){
		pos = 0;

		while ((pos = Theta_list[k].find(",")) != std::string::npos) {
			token = Theta_list[k].substr(0, pos);

			if(counter >= K*Km) return false;
			Theta_solution[counter] = atof(token.c_str());
			counter++;

			Theta_list[k].erase(0, pos + 1);
		}

		if(counter >= K*Km) return false;
		Theta_solution[counter] = atof(Theta_list[k].c_str());
		counter++;
	}

	if(counter != K*Km){
		return false;
	} else {
		return true;
	}
}

bool parse_strings_to_doubles(int xdim, std::string mu_string, double *mu){
	std::string token;
	size_t pos=0;
	int counter=0;

    while ((pos = mu_string.find(",")) != std::string::npos) {
        token = mu_string.substr(0, pos);

        if(counter >= xdim) return false;

        mu[counter] = atof(token.c_str());
        counter++;

        mu_string.erase(0, pos + 1);
    }

    /* the last one */
    if(counter >= xdim) return false;
    mu[counter] = atof(mu_string.c_str());
    counter++;

    /* final check of the size */
	if(counter != xdim){
		return false;
	} else {
		return true;
	}
}

void myround(double in, double *out){
	union myUnion {
		double dValue;
		uint64_t iValue;
	} myValue;
	myValue.dValue=in;

//	myValue.iValue = myValue.iValue*0.001;
//	myValue.iValue = myValue.iValue*1000;

	*out = myValue.dValue;
}

std::string printbool(bool input){
	std::string out;
	if(input){
		out = "true";
	} else {
		out = "false";
	}
	return out;
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, std::back_inserter(elems));
    return elems;
}

}
} /* end of namespace */
