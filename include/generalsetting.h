/** @file generalsetting.h
 *  @brief general solver settings
 *
 *  Defines the parent class for manipulation with settings of the solver - i.e. number of iterations, precision, etc.
 *  All specific solvers implementations should contain a setting class defined as inherited classes from this class.
 *
 *  @author Lukas Pospisil
 */

#ifndef PASC_GENERALSETTING_H
#define	PASC_GENERALSETTING_H

#include "common.h"

namespace pascinference {

class GeneralSetting {
	protected:

	public:
		GeneralSetting() {};
		~GeneralSetting() {};

		virtual void print(std::ostream &output) const {};
		virtual std::string get_name() const;

		friend std::ostream &operator<<(std::ostream &output, const GeneralSetting &setting); 
	
};


/* general print, call virtual print() */
std::ostream &operator<<(std::ostream &output, const GeneralSetting &setting){
	if(DEBUG_MODE >= 100) coutMaster << "(GeneralSolverSetting)OPERATOR: <<" << std::endl;
	output << setting.get_name();
	return output;
}

std::string GeneralSetting::get_name() const {
	return "GeneralSetting";
}

} /* end of namespace */

#endif
