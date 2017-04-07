/** @file shortinfo.h
 *  @brief Manipulation with output file which includes short final information.
 *
 *  @author Lukas Pospisil
 */

#ifndef PASC_COMMON_SHORTINFO_H
#define	PASC_COMMON_SHORTINFO_H

#include <stdio.h>
#include <fstream>
#include <string>

#include "common/globalmanager.h"

namespace pascinference {
namespace common {

/** \class ShortinfoClass
 *  \brief Manipulation with shortinfo file.
 *
*/
class ShortinfoClass {
	private:
		std::string *filename;			/**< name of log file */
		std::ofstream myfile;			/**< shortinfo file */

		bool shortinfo_or_not;			/**< generate shortinfo file or not */

	public:
		/** @brief basic constructor
		 * 
		 *  Defaultly do not generate shortinfo file, wait until begin() with filename is called.
		 * 
		 */
		ShortinfoClass();

		/** @brief destructor
		 */
		~ShortinfoClass();

		/** @brief start to generate new shortinfo file
		 * 
		 * Set new filename and create shortinfo file.
		 * 
		 * @param new_filename the name of shortinfo file
		 */
		void begin(std::string new_filename);

		/** @brief write sting info shortinfo file
		 * 
		 * The master (process with rank=0) writes input string.
		 * 
		 * @param what_to_write the string which will be written into file
		 */
		void write(std::string what_to_write);
		
};

extern ShortinfoClass shortinfo; /**< global instance of class for manipulating with shortinfo */

}
} /* end of namespace */

#endif
