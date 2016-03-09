/** @file generalmodel.h
 *  @brief class for manipulation with models
 *
 *  Header file which defines the parent class for manipulation with models - additional information for solving the problem.
 *  All specific model implementations should be defined as inherited classes from this class.
 *
 *  @author Lukas Pospisil
 */

#ifndef PASC_GENERALMODEL_H
#define	PASC_GENERALMODEL_H

#include "common.h"

namespace pascinference {

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
		GeneralModel() {};
		~GeneralModel() {};

		virtual void print(std::ostream &output) const;
		virtual std::string get_name() const;

		friend std::ostream &operator<<(std::ostream &output, const GeneralModel &model); /* cannot be virtual, therefore it call virtual print() */
	
};

/* general print, call virtual print() */
std::ostream &operator<<(std::ostream &output, const GeneralModel &model){
	if(DEBUG_MODE >= 100) std::cout << "(GeneralData)OPERATOR: <<" << std::endl;
	model.print(output);
	return output;
}

void GeneralModel::print(std::ostream &output) const {
	output << this->get_name() << std::endl;
}

std::string GeneralModel::get_name() const {
	return "GeneralModel";
}




} /* end of namespace */


#endif
