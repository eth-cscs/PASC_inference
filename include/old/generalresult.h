/** @file generalresult.h
 *  @brief class for manipulation with results
 *
 *  Defines the parent class for manipulation with results - store the content and postprocess.
 *  All specific result implementations should be defined as inherited classes from this class.
 *
 *  @author Lukas Pospisil
 */

#ifndef PASC_GENERALRESULT_H
#define	PASC_GENERALRESULT_H

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
