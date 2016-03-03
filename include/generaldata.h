#ifndef PASCGENERALDATA_H
#define	PASCGENERALDATA_H

#include "common.h"

namespace pascinference {

class GeneralData {
	protected:

	public:
		GeneralData() {};
		~GeneralData() {};

		virtual void print(std::ostream &output) const;
		virtual std::string get_name() const;

		friend std::ostream &operator<<(std::ostream &output, const GeneralData &data); /* cannot be virtual, therefore it call virtual print() */
	
};

/* general print, call virtual print() */
std::ostream &operator<<(std::ostream &output, const GeneralData &data){
	if(DEBUG_MODE >= 100) std::cout << "(GeneralData)OPERATOR: <<" << std::endl;
	data.print(output);
	return output;
}

void GeneralData::print(std::ostream &output) const {
	output << this->get_name() << std::endl;
}

std::string GeneralData::get_name() const {
	return "GeneralData";
}





} /* end of namespace */


#endif
