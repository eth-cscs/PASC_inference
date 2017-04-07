#include "common/shortinfo.h"

namespace pascinference {
namespace common {

ShortinfoClass::ShortinfoClass(){
	shortinfo_or_not = false;
}

ShortinfoClass::~ShortinfoClass(){
}

void ShortinfoClass::begin(std::string new_filename){
	shortinfo_or_not = true;

	this->filename = new std::string(new_filename);

	/* open file, i.e. create it or delete content */
	if( GlobalManager.get_rank() == 0){
		myfile.open(filename->c_str());
		myfile.close();
	}
}

void ShortinfoClass::write(std::string what_to_write){
	/* master writes the file with short info (used in batch script for quick computation) */
	if( GlobalManager.get_rank() == 0 && shortinfo_or_not){
		myfile.open(filename->c_str(), std::fstream::in | std::fstream::out | std::fstream::app);

		myfile << what_to_write;
		
		myfile.close();
	}

	/* wait for master */
	TRYCXX( PetscBarrier(NULL) );
}
		
ShortinfoClass shortinfo; 


}
} /* end of namespace */


