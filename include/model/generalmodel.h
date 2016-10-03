/** @file generalmodel.h
 *  @brief class for manipulation with models
 *
 *  @author Lukas Pospisil
 */

#ifndef PASC_GENERALMODEL_H
#define	PASC_GENERALMODEL_H

#include "pascinference.h"

namespace pascinference {
namespace model {

/** \class GeneralModel
 *  \brief General class for manipulation with models.
 *
 *  Parent class for manipulation with models - additional information for solving the problem.
 *  All specific model implementations should be defined as inherited classes from this class.
 *	
*/
class GeneralModel {
	protected:

	public:
		/** @brief default constructor
		 */ 
		GeneralModel() {};

		/** @brief default destructor
		 */ 
		~GeneralModel() {};

		/** @brief print info about model
		 */ 
		virtual void print(ConsoleOutput &output) const;

		/** @brief get the name of model
		 */ 
		virtual std::string get_name() const;

		friend ConsoleOutput &operator<<(ConsoleOutput &output, const GeneralModel &model); /* cannot be virtual, therefore it call virtual print() */
	
};

/* general print, call virtual print() */
ConsoleOutput &operator<<(ConsoleOutput &output, const GeneralModel &model){
	if(DEBUG_MODE >= 100) coutMaster << "(GeneralData)OPERATOR: <<" << std::endl;
	output << model.get_name();
	return output;
}

void GeneralModel::print(ConsoleOutput &output) const {
	output <<  this->get_name() << std::endl;
	output.synchronize();
}

std::string GeneralModel::get_name() const {
	return "GeneralModel";
}


}
} /* end of namespace */


#endif
