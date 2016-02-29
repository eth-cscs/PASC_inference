#ifndef PASCGENERALDATA_H
#define	PASCGENERALDATA_H

#include "common.h"

namespace pascinference {

class GeneralData {
	protected:

	public:
		GeneralData() {};
		~GeneralData() {};

		virtual void print(std::ostream &output) const {};

		friend std::ostream &operator<<(std::ostream &output, const GeneralData &data); /* cannot be virtual, therefore it call virtual print() */
	
};

/* general print, call virtual print() */
std::ostream &operator<<(std::ostream &output, const GeneralData &data){
	if(DEBUG_MODE >= 100) std::cout << "(GeneralData)OPERATOR: <<" << std::endl;
	data.print(output);
	return output;
}





} /* end of namespace */


#endif
