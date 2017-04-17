#include "algebra/integration/entropyintegration.h"

namespace pascinference {
namespace algebra {

/* constructor */
EntropyIntegration::EntropyIntegration(int m_new, int Km_new) {
	LOG_FUNC_BEGIN

	/* set given parameters */
	set_m(m_new);
	set_Km(Km_new);

	LOG_FUNC_END
}

/* destructor */
EntropyIntegration::~EntropyIntegration(){
	LOG_FUNC_BEGIN
	
	LOG_FUNC_END
}

/* print info about integration */
void EntropyIntegration::print(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output << this->get_name() << std::endl;
	
	output <<  " - m                 : " << this->m << std::endl;
	output <<  " - Km                : " << this->Km << std::endl;

	output.synchronize();	

	LOG_FUNC_END
}

/* print info about integration */
void EntropyIntegration::print(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN

	output_global <<  this->get_name() << std::endl;
	
	output_global <<  " - m                 : " << this->m << std::endl;
	output_global <<  " - Km                : " << this->Km << std::endl;

	output_global.synchronize();

	LOG_FUNC_END
}

/* get name of the model */
std::string EntropyIntegration::get_name() const {
	return "Entropy-Integration general";
}


int EntropyIntegration::get_Km() const {
	return this->Km;
}

void EntropyIntegration::set_Km(int Km_new) {
	this->Km = Km_new;
}

int EntropyIntegration::get_m() const {
	return this->m;
}

void EntropyIntegration::set_m(int m_new) {
	this->m = m_new;
}


}
} /* end namespace */




