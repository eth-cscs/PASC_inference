/** @file entropysolverdlib.h
 *  @brief Solver which solves problem with integrals using DLib library
 *
 *  @author Anna Marchenko & Lukas Pospisil 
 */
 
#ifndef PASC_ENTROPYSOLVERDLIB_H
#define	PASC_ENTROPYSOLVERDLIB_H

#include "general/solver/generalsolver.h"
#include "general/data/entropydata.h"

#define ENTROPYSOLVERDLIB_DEFAULT_MAXIT 1000
#define ENTROPYSOLVERDLIB_DEFAULT_EPS 1e-6
#define ENTROPYSOLVERDLIB_DEFAULT_INTEGRATION_EPS 1e-10
#define ENTROPYSOLVERDLIB_DEFAULT_DEBUGMODE 0

namespace pascinference {
namespace solver {

/** \class EntropySolverDlib
 *  \brief Solver which solves problem with integrals using DLib library
 *
 *  Operates on EntropyData.
*/
template<class VectorBase>
class EntropySolverDlib: public GeneralSolver {
	public:
		class ExternalContent;

	protected:
		friend class ExternalContent;
		ExternalContent *externalcontent;			/**< for manipulation with external-specific stuff */

		Timer timer_solve; /**< total solution time */
		Timer timer_compute_moments; /**< time of moment computation */

		EntropyData<VectorBase> *entropydata; /**< data on which the solver operates */

		/* aux vectors */
		GeneralVector<VectorBase> *moments; /**< vector of computed moments */

		/** @brief set settings of algorithm from arguments in console
		* 
		*/
		void set_settings_from_console();

		double integration_eps;		/**< precision of integration in Dlib adaptive simpson rule */

		/* debug */
		int debugmode;				/**< basic debug mode schema [0/1/2] */
		bool debug_print_it;		/**< print simple info about outer iterations */
		bool debug_print_moments;	/**< print moments during iterations */
		
	public:

		EntropySolverDlib();
		EntropySolverDlib(EntropyData<VectorBase> &new_entropydata); 
		~EntropySolverDlib();

		void solve();

		void print(ConsoleOutput &output) const;
		void print(ConsoleOutput &output_global, ConsoleOutput &output_local) const;
		void printstatus(ConsoleOutput &output) const;
		void printstatus(std::ostringstream &output) const;
		void printcontent(ConsoleOutput &output) const;
		void printtimer(ConsoleOutput &output) const;
		std::string get_name() const;

		EntropyData<VectorBase> *get_data() const;

		void compute_moments();
		
		void compute_residuum(GeneralVector<VectorBase> *residuum) const;
		
		ExternalContent *get_externalcontent() const;
		
		int get_xdim() const;
		int get_K() const;
		int get_Km() const;
		int get_number_of_moments() const;

};


}
} /* end of namespace */

/* ------------- implementation ----------- */
namespace pascinference {
namespace solver {
			
template<class VectorBase>
void EntropySolverDlib<VectorBase>::set_settings_from_console() {
	consoleArg.set_option_value("entropysolverdlib_maxit", &this->maxit, ENTROPYSOLVERDLIB_DEFAULT_MAXIT);
	consoleArg.set_option_value("entropysolverdlib_eps", &this->eps, ENTROPYSOLVERDLIB_DEFAULT_EPS);
	consoleArg.set_option_value("entropysolverdlib_integration_eps", &this->integration_eps, ENTROPYSOLVERDLIB_DEFAULT_INTEGRATION_EPS);

	/* set debug mode */
	consoleArg.set_option_value("entropysolverdlib_debugmode", &this->debugmode, ENTROPYSOLVERDLIB_DEFAULT_DEBUGMODE);
	
	if(debugmode == 1){
		debug_print_it = true;
	}

	if(debugmode == 2){
		debug_print_moments = true;
	}

	consoleArg.set_option_value("entropysolverdlib_debug_print_it",		&debug_print_it, 		debug_print_it);
	consoleArg.set_option_value("entropysolverdlib_debug_print_moments",&debug_print_moments, 		debug_print_moments);

}


/* constructor */
template<class VectorBase>
EntropySolverDlib<VectorBase>::EntropySolverDlib(){
	LOG_FUNC_BEGIN
	
	entropydata = NULL;
	
	/* settings */
	this->maxit = 0;
	this->eps = 0;
	this->debugmode = 0;

	/* settings */
	set_settings_from_console();
	
	/* prepare timers */
	this->timer_solve.restart();	
	this->timer_compute_moments.restart();

	LOG_FUNC_END
}

template<class VectorBase>
EntropySolverDlib<VectorBase>::EntropySolverDlib(EntropyData<VectorBase> &new_entropydata){
	LOG_FUNC_BEGIN

	entropydata = &new_entropydata;

	/* settings */
	this->maxit = 0;
	this->eps = 0;
	this->debugmode = 0;

	/* settings */
	set_settings_from_console();

	/* prepare timers */
	this->timer_solve.restart();	
	this->timer_compute_moments.restart();

	//TODO

	LOG_FUNC_END
}

/* destructor */
template<class VectorBase>
EntropySolverDlib<VectorBase>::~EntropySolverDlib(){
	LOG_FUNC_BEGIN

	LOG_FUNC_END
}


/* print info about problem */
template<class VectorBase>
void EntropySolverDlib<VectorBase>::print(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output <<  this->get_name() << std::endl;
	
	/* print settings */
	output <<  " - maxit             : " << this->maxit << std::endl;
	output <<  " - eps               : " << this->eps << std::endl;
	output <<  " - integration_eps   : " << this->integration_eps << std::endl;

	output <<  " - xdim              : " << this->get_xdim() << std::endl;
	output <<  " - K                 : " << this->get_K() << std::endl;
	output <<  " - Km                : " << this->get_Km() << std::endl;
	output <<  " - number_of_moments : " << this->get_number_of_moments() << std::endl;

	output <<  " - debugmode         : " << this->debugmode << std::endl;

	/* print data */
	if(entropydata){
		coutMaster.push();
		entropydata->print(output);
		coutMaster.pop();
	}
	
	LOG_FUNC_END
}

/* print info about problem */
template<class VectorBase>
void EntropySolverDlib<VectorBase>::print(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN

	output_global <<  this->get_name() << std::endl;
	
	/* print settings */
	output_global <<  " - maxit             : " << this->maxit << std::endl;
	output_global <<  " - eps               : " << this->eps << std::endl;
	output_global <<  " - integration_eps   : " << this->integration_eps << std::endl;

	output_global <<  " - xdim              : " << this->get_xdim() << std::endl;
	output_global <<  " - K                 : " << this->get_K() << std::endl;
	output_global <<  " - Km                : " << this->get_Km() << std::endl;
	output_global <<  " - number_of_moments : " << this->get_number_of_moments() << std::endl;

	output_global <<  " - debugmode         : " << this->debugmode << std::endl;

	/* print data */
	if(entropydata){
		output_global <<  " - data            : " << std::endl;
		coutMaster.push();
		entropydata->print(output_global,output_local);
		coutMaster.pop();
	}
	
	LOG_FUNC_END
}

template<class VectorBase>
void EntropySolverDlib<VectorBase>::printstatus(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output <<  this->get_name() << std::endl;
	output <<  " - used memory: " << MemoryCheck::get_virtual() << "%" << std::endl;

	LOG_FUNC_END
}

template<class VectorBase>
void EntropySolverDlib<VectorBase>::printstatus(std::ostringstream &output) const {
	LOG_FUNC_BEGIN

	std::streamsize ss = std::cout.precision();

	output << std::setprecision(17);
	
	//TODO
	
	output << std::setprecision(ss);

	LOG_FUNC_END
}

/* print content of solver */
template<class VectorBase>
void EntropySolverDlib<VectorBase>::printcontent(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output << this->get_name() << std::endl;
	
	/* print content of data */
	if(entropydata){
		output << "- data:" << std::endl;
		coutMaster.push();
		entropydata->printcontent(output);
		coutMaster.pop();
	}
		
	LOG_FUNC_END
}

template<class VectorBase>
void EntropySolverDlib<VectorBase>::printtimer(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output <<  this->get_name() << std::endl;
	output <<  " - timers" << std::endl;
	output <<  "  - t_solve   = " << this->timer_solve.get_value_sum() << std::endl;
	output <<  "  - t_moments = " << this->timer_compute_moments.get_value_sum() << std::endl;

	LOG_FUNC_END
}

template<class VectorBase>
std::string EntropySolverDlib<VectorBase>::get_name() const {
	std::string return_value = "EntropySolverDlib<" + GeneralVector<VectorBase>::get_name() + ">";
	return return_value;
}

template<class VectorBase>
EntropyData<VectorBase> *EntropySolverDlib<VectorBase>::get_data() const {
	return entropydata;
}

template<class VectorBase>
void EntropySolverDlib<VectorBase>::solve() {
	LOG_FUNC_BEGIN

	//TODO
	
	LOG_FUNC_END
}

template<class VectorBase>
void EntropySolverDlib<VectorBase>::compute_moments() {
	LOG_FUNC_BEGIN

	//TODO
	
	LOG_FUNC_END
}

template<class VectorBase>
void EntropySolverDlib<VectorBase>::compute_residuum(GeneralVector<VectorBase> *residuum) const {
	LOG_FUNC_BEGIN

	//TODO
		
	LOG_FUNC_END
}

template<class VectorBase>
int EntropySolverDlib<VectorBase>::get_xdim() const {
	return entropydata->get_xdim();
}

template<class VectorBase>
int EntropySolverDlib<VectorBase>::get_K() const {
	return entropydata->get_K();
}

template<class VectorBase>
int EntropySolverDlib<VectorBase>::get_Km() const {
	return entropydata->get_Km();
}

template<class VectorBase>
int EntropySolverDlib<VectorBase>::get_number_of_moments() const {
	return entropydata->get_number_of_moments();
}

/* define blank external content for general VectorBase */
template<class VectorBase>
class EntropySolverDlib<VectorBase>::ExternalContent {
};


}
} /* end namespace */

#endif
