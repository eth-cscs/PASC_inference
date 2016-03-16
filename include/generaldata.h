/** @file generaldata.h
 *  @brief class for manipulation with data
 *
 *  Defines the parent class for manipulation with data (input and output).
 *  All specific data implementations should be defined as inherited classes from this class.
 *
 *  @author Lukas Pospisil
 */

#ifndef PASC_GENERALDATA_H
#define	PASC_GENERALDATA_H

#include "common.h"

namespace pascinference {

/** \class GeneralData
 *  \brief general class for manipulation with data
 *
 *  Parent class for manipulation with data.
 *  All specific data implementations should be defined as inherited classes from this class.
 *	
*/
class GeneralData {
	protected:

	public:
	
		/** @brief default constructor
		 */ 
		GeneralData() {};

		/** @brief default destructor
		 */ 
		~GeneralData() {};

		/** @brief print basic information about data
		 * 
		 * @param output where to print
		 */
		virtual void print(std::ostream &output) const;

		/** @brief print content of data
		 * 
		 * Print all values in inner data structures.
		 * 
		 * @param output where to print
		 */
		virtual void printcontent(std::ostream &output) const;

		/** @brief get the name of data
		 */ 
		virtual std::string get_name() const;

		/** @brief append name to input stream
		 * 
		 */
		friend std::ostream &operator<<(std::ostream &output, const GeneralData &data); /* cannot be virtual, therefore it call virtual print() */
	
};

/* general print, call virtual print() */
std::ostream &operator<<(std::ostream &output, const GeneralData &data){
	if(DEBUG_MODE >= 100) coutMaster << "(GeneralData)OPERATOR: <<" << std::endl;
	data.print(output);
	return output;
}

void GeneralData::print(std::ostream &output) const {
	output <<  this->get_name() << std::endl;
}

void GeneralData::printcontent(std::ostream &output) const {
	output <<  this->get_name() << std::endl;
}

std::string GeneralData::get_name() const {
	return "GeneralData";
}





} /* end of namespace */


#endif
