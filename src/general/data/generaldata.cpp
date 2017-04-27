#include "general/data/generaldata.h"

namespace pascinference {
namespace data {

/* general print, call virtual print() */
ConsoleOutput &operator<<(ConsoleOutput &output, const GeneralData &data){
	output << data.get_name();
	return output;
}

void GeneralData::print(ConsoleOutput &output) const {
	output <<  this->get_name() << std::endl;
}

void GeneralData::print(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	output_global <<  this->get_name() << std::endl;
}

void GeneralData::printcontent(ConsoleOutput &output) const {
	output <<  this->get_name() << std::endl;
}

void GeneralData::printcontent(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	output_global <<  this->get_name() << std::endl;
}

std::string GeneralData::get_name() const {
	return "GeneralData";
}


}
} /* end of namespace */

