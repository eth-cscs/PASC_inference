#include "general/model/generalmodel.h"

namespace pascinference {
namespace model {

/* general print, call virtual print() */
ConsoleOutput &operator<<(ConsoleOutput &output, const GeneralModel &model){
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

