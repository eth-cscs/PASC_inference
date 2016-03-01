#ifndef PASCGENERALSETTING_H
#define	PASCGENERALSETTING_H

#include "common.h"

namespace pascinference {

class GeneralSetting {
	protected:

	public:
		GeneralSetting() {};
		~GeneralSetting() {};

		virtual void print(std::ostream &output) const {};

		friend std::ostream &operator<<(std::ostream &output, const GeneralSetting &setting); 
	
};


/* general print, call virtual print() */
std::ostream &operator<<(std::ostream &output, const GeneralSetting &setting){
	if(DEBUG_MODE >= 100) std::cout << "(GeneralSolverSetting)OPERATOR: <<" << std::endl;
	setting.print(output);
	return output;
}

} /* end of namespace */

#endif
