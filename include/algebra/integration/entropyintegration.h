/** @file entropyintegration.h
 *  @brief Computes numerical integrals in entropy problem
 *
 *  This integral is computed locally using only "local" MPI process.
 * 
 *  @author Lukas Pospisil 
 */
 
#ifndef PASC_ENTROPYINTEGRATION_H
#define	PASC_ENTROPYINTEGRATION_H

namespace pascinference {
namespace algebra {

/** \class EntropyIntegration
 *  \brief Main class providing algorithm for numerical computation of integrals in Entropy problem.
 *
 * \f[
 *  \int\limits_{-1}^{1} x^j e^{-\sum\limits_{k=1}^{m} \lambda_k x^k } ~\textrm{dx}, ~~~ j = 0,\dots,Km
 *	\f] 
 * 
*/
class EntropyIntegration {
	protected:
		int m;		/**< length of lambda vector */
		int Km;		/**< number of integrals (the largest power) */

	public:

		EntropyIntegration(int m_new, int Km_new);
		~EntropyIntegration();

		void print(ConsoleOutput &output) const;
		void print(ConsoleOutput &output_global, ConsoleOutput &output_local) const;
		std::string get_name() const;


		int get_m() const;
		void set_m(int m_new);
		int get_Km() const;
		void set_Km(int Km_new);

};

}
} /* end of namespace */


/* ------------- implementation ----------- */
//TODO: move to impls

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


/* print info about model */
void EntropyIntegration::print(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output << this->get_name() << std::endl;
	
	output <<  " - m                 : " << this->m << std::endl;
	output <<  " - Km                : " << this->Km << std::endl;

	output.synchronize();	

	LOG_FUNC_END
}

/* print info about model */
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

#endif


