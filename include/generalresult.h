#ifndef PASCGENERALRESULT_H
#define	PASCGENERALRESULT_H

#include "common.h"

namespace pascinference {

class GeneralResult {
	protected:

	public:
		GeneralResult() {};
		~GeneralResult() {};

		virtual void print(std::ostream &output) const;
		virtual std::string get_name() const;

		friend std::ostream &operator<<(std::ostream &output, const GeneralResult &result); /* cannot be virtual, therefore it call virtual print() */
	
};

/* general print, call virtual print() */
std::ostream &operator<<(std::ostream &output, const GeneralResult &result){
	if(DEBUG_MODE >= 100) std::cout << "(GeneralResult)OPERATOR: <<" << std::endl;
	result.print(output);
	return output;
}

void GeneralResult::print(std::ostream &output) const {
	output << this->get_name() << std::endl;
}

std::string GeneralResult::get_name() const {
	return "GeneralResult";
}




} /* end of namespace */


#endif
