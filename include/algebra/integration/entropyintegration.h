/** @file entropyintegration.h
 *  @brief Computes numerical integrals in entropy problem
 *
 *  This integral is computed locally using only "local" MPI process.
 * 
 *  @author Lukas Pospisil 
 */
 
#ifndef PASC_ENTROPYINTEGRATION_H
#define	PASC_ENTROPYINTEGRATION_H

#include <string>
#include <iostream>
#include "common/consoleoutput.h"
#include "common/logging.h"

namespace pascinference {
using namespace common;
	
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

#endif



