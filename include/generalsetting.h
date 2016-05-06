/** @file generalsetting.h
 *  @brief general solver settings
 *
 *  @author Lukas Pospisil
 */

#ifndef PASC_GENERALSETTING_H
#define	PASC_GENERALSETTING_H

#include "common.h"

namespace pascinference {

/** @class GeneralSetting
 *  @brief setting of anything
 * 
 *  Could be used to store any settings.
 *  Parent class for manipulation with settings.
 *  All specific settings implementations should be defined as inherited classes from this class.
 * 
 *  Defines the parent class for manipulation with settings of the solver - i.e. number of iterations, precision, etc.
 *  All specific solvers implementations should contain a setting class defined as inherited classes from this class.
 * 
 */ 
class GeneralSetting {
	protected:

	public:
		/** @brief default constructor
		 */
		GeneralSetting() {};

		/** @brief default destructor
		 */
		~GeneralSetting() {};

		/** @brief print settings
		 * 
		 * @param output where to print
		 */
		virtual void print(ConsoleOutput &output) const {};

		/** @brief return the name of settings
		 * 
		 */
		virtual std::string get_name() const;

		friend ConsoleOutput &operator<<(ConsoleOutput &output, const GeneralSetting &setting); 
	
};


/* general print, call virtual print() */
ConsoleOutput &operator<<(ConsoleOutput &output, const GeneralSetting &setting){
	if(DEBUG_MODE >= 100) coutMaster << "(GeneralSolverSetting)OPERATOR: <<" << std::endl;
	output << setting.get_name();
	return output;
}

std::string GeneralSetting::get_name() const {
	return "GeneralSetting";
}

} /* end of namespace */

#endif
