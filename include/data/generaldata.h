/** @file generaldata.h
 *  @brief general class for manipulation with data
 *
 *  @author Lukas Pospisil
 */

#ifndef PASC_GENERALDATA_H
#define	PASC_GENERALDATA_H

#include "pascinference.h"

namespace pascinference {
namespace data {

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
		virtual void print(ConsoleOutput &output) const;

		virtual void print(ConsoleOutput &output_global, ConsoleOutput &output_local) const;

		/** @brief print content of data
		 * 
		 * Print all values in inner data structures.
		 * 
		 * @param output where to print
		 */
		virtual void printcontent(ConsoleOutput &output) const;

		virtual void printcontent(ConsoleOutput &output_global, ConsoleOutput &output_local) const;

		/** @brief get the name of data
		 */ 
		virtual std::string get_name() const;

		/** @brief append name to input stream
		 * 
		 */
		friend ConsoleOutput &operator<<(ConsoleOutput &output, const GeneralData &data); /* cannot be virtual, therefore it call virtual print() */
	
};

/* general print, call virtual print() */
ConsoleOutput &operator<<(ConsoleOutput &output, const GeneralData &data){
	if(DEBUG_MODE >= 100) coutMaster << "(GeneralData)OPERATOR: <<\n";
	output << data.get_name();
	return output;
}

void GeneralData::print(ConsoleOutput &output) const {
	output <<  this->get_name() << "\n";
}

void GeneralData::print(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	output_global <<  this->get_name() << "\n";
}

void GeneralData::printcontent(ConsoleOutput &output) const {
	output <<  this->get_name() << "\n";
}

void GeneralData::printcontent(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	output_global <<  this->get_name() << "\n";
}

std::string GeneralData::get_name() const {
	return "GeneralData";
}


}
} /* end of namespace */


#endif
