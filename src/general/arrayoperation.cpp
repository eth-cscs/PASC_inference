#include "algebra/arrayoperation.h"

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


}
} /* end of namespace */
