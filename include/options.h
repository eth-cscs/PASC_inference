#ifndef PASC_OPTIONS_H
#define	PASC_OPTIONS_H

#include <boost/program_options.hpp>

namespace pascinference {

static void add_options(boost::program_options::options_description *description){
	/* TS_SOLVER_GLOBAL */
	description->add_options()
		("tssolver_global_maxit", boost::program_options::value<int>(), "TSSOLVER_GLOBAL - maximum number of iterations")
		("tssolver_global_eps", boost::program_options::value<double>(), "TSSOLVER_GLOBAL - precision")
		("tssolver_global_debug_mode", boost::program_options::value<int>(), "TSSOLVER_GLOBAL - debug mode");

	/* MULTICG */
	description->add_options()
		("multicgsolver_maxit", boost::program_options::value<int>(), "MULTICGSOLVER - maximum number of iterations")
		("multicgsolver_eps", boost::program_options::value<double>(), "MULTICGSOLVER - precision")
		("multicgsolver_debug_mode", boost::program_options::value<int>(), "MULTICGSOLVER - debug mode");

	/* MULTICG_GLOBAL */
	description->add_options()
		("multicgsolver_global_maxit", boost::program_options::value<int>(), "MULTICGSOLVER_GLOBAL - maximum number of iterations")
		("multicgsolver_global_eps", boost::program_options::value<double>(), "MULTICGSOLVER_GLOBAL - precision")
		("multicgsolver_global_debug_mode", boost::program_options::value<int>(), "MULTICGSOLVER_GLOBAL - debug mode");

	/* QPSOLVER */
	description->add_options()
		("qpsolver_maxit", boost::program_options::value<int>(), "QPSOLVER - maximum number of iterations")
		("qpsolver_eps", boost::program_options::value<double>(), "QPSOLVER - precision")
		("qpsolver_debug_mode", boost::program_options::value<int>(), "QPSOLVER - debug mode");

	/* CGQPSOLVER */
	description->add_options()
		("cgqpsolver_maxit", boost::program_options::value<int>(), "CGQPSOLVER - maximum number of iterations")
		("cgqpsolver_eps", boost::program_options::value<double>(), "CGQPSOLVER - precision")
		("cgqpsolver_debug_mode", boost::program_options::value<int>(), "CGQPSOLVER - debug mode");

	/* SPGQPSOLVER */
	description->add_options()
		("spgqpsolver_maxit", boost::program_options::value<int>(), "SPGQPSOLVER - maximum number of iterations")
		("spgqpsolver_eps", boost::program_options::value<double>(), "SPGQPSOLVER - precision")
		("spgqpsolver_debug_mode", boost::program_options::value<int>(), "SPGQPSOLVER - debug mode")
		("spgqpsolver_m", boost::program_options::value<int>(), "SPGQPSOLVER - parameter of generalized Armijo condition")
		("spgqpsolver_gamma", boost::program_options::value<double>(), "SPGQPSOLVER - parameter of generalized Armijo condition")
		("spgqpsolver_sigma1", boost::program_options::value<double>(), "SPGQPSOLVER - parameter of generalized Armijo condition")
		("spgqpsolver_sigma2", boost::program_options::value<double>(), "SPGQPSOLVER - parameter of generalized Armijo condition")
		("spgqpsolver_alphainit", boost::program_options::value<double>(), "SPGQPSOLVER - initial BB step-size")
		("spgqpsolver_stop_normgp", boost::program_options::value<bool>(), "SPGQPSOLVER - stopping criteria based on norm(gp)")
		("spgqpsolver_stop_normgp_normb", boost::program_options::value<bool>(), "SPGQPSOLVER - stopping criteria based on norm(gp) and norm(b)")
		("spgqpsolver_stop_difff", boost::program_options::value<bool>(), "SPGQPSOLVER - stopping criteria based on difference of object function");

}

}

#endif
