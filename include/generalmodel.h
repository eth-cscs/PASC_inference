#ifndef PASCGENERALMODEL_H
#define	PASCGENERALMODEL_H

#include "common.h"

namespace pascinference {

class GeneralModel {
	protected:

	public:
		GeneralModel() {};
		~GeneralModel() {};

		virtual void print(std::ostream &output) const {};

		friend std::ostream &operator<<(std::ostream &output, const GeneralModel &data); /* cannot be virtual, therefore it call virtual print() */
	
};

/* general print, call virtual print() */
std::ostream &operator<<(std::ostream &output, const GeneralModel &model){
	if(DEBUG_MODE >= 100) std::cout << "(GeneralData)OPERATOR: <<" << std::endl;
	model.print(output);
	return output;
}





} /* end of namespace */


#endif
