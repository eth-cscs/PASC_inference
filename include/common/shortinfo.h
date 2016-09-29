#ifndef PASC_COMMON_SHORTINFO_H
#define	PASC_COMMON_SHORTINFO_H

#include "common/globalmanager.h"

namespace pascinference {
namespace common {

class ShortinfoClass {
	private:
		std::string *filename;
		std::ofstream myfile;

		bool shortinfo_or_not;

	public:
		ShortinfoClass(){
			shortinfo_or_not = false;
		}

		~ShortinfoClass(){
	
		}

		void begin(std::string new_filename){
			shortinfo_or_not = true;

			this->filename = new std::string(new_filename);

			/* open file, i.e. create it or delete content */
			myfile.open(filename->c_str());
			myfile.close();			
		}

		void write(std::string what_to_write){
			/* master writes the file with short info (used in batch script for quick computation) */
			if( GlobalManager.get_rank() == 0 && shortinfo_or_not){
				myfile.open(filename->c_str(), std::fstream::in | std::fstream::out | std::fstream::app);

				myfile << what_to_write;
		
				myfile.close();
			}

			/* wait for master */
			TRY( PetscBarrier(NULL) );
		}
		
};

ShortinfoClass shortinfo;


}
} /* end of namespace */

#endif
