/** @file options.h
 *  @brief defines command line options
 *
 *  Includes the definition of all command line options with description.
 *  Call program with --help to see these options.
 *
 *  @author Lukas Pospisil
 */
 
#ifndef PASC_OPTIONS_H
#define	PASC_OPTIONS_H

#include <boost/program_options.hpp>

namespace pascinference {

/** @brief add all options of common library files
*
* @param description boost program options instance
* @param console_nmb_cols number of columns in console
*/
static void add_options(boost::program_options::options_description *description, int console_nmb_cols){

	/* ----- SOLVERS ------ */
	boost::program_options::options_description opt_solvers("--- SOLVERS ---------------------", console_nmb_cols);
	description->add(opt_solvers);

		/* TSSOLVER */
		boost::program_options::options_description opt_tssolver("TSSOLVER", console_nmb_cols);
		opt_tssolver.add_options()
			("tssolver_maxit", boost::program_options::value<int>(), "maximum number of iterations [int]")
			("tssolver_eps", boost::program_options::value<double>(), "precision [double]")
			("tssolver_debug_mode", boost::program_options::value<int>(), "basic debug mode schema [0/1/2/3]")
			("tssolver_debug_print_annealing", boost::program_options::value<bool>(), "print info about annealing steps")
			("tssolver_debug_print_it", boost::program_options::value<bool>(), "print simple info about outer iterations")
			("tssolver_debug_print_theta", boost::program_options::value<bool>(), "print theta solver info")
			("tssolver_debug_print_theta_solution", boost::program_options::value<bool>(), "print solution of theta problem in each iteration")
			("tssolver_debug_print_gamma", boost::program_options::value<bool>(), "print gamma solver info")
			("tssolver_debug_print_gamma_solution", boost::program_options::value<bool>(), "print solution of gamma problem in each iteration");
		opt_solvers->add(opt_tssolver);

		/* MULTICG */
		boost::program_options::options_description opt_multicgsolver("MULTICGSOLVER", console_nmb_cols);
		opt_multicgsolver.add_options()
			("multicgsolver_maxit", boost::program_options::value<int>(), "maximum number of iterations [int]")
			("multicgsolver_eps", boost::program_options::value<double>(), "precision [double]")
			("multicgsolver_debug_mode", boost::program_options::value<int>(), "debug mode [int]");
		opt_solvers->add(opt_multicgsolver);

		/* MULTICG_GLOBAL */
		boost::program_options::options_description opt_multicgsolver("MULTICGSOLVER_GLOBAL", console_nmb_cols);
		opt_multicgsolver.add_options()
			("multicgsolver_maxit", boost::program_options::value<int>(), "maximum number of iterations [int]")
			("multicgsolver_eps", boost::program_options::value<double>(), "precision [double]")
			("multicgsolver_debug_mode", boost::program_options::value<int>(), "debug mode [int]");
		opt_solvers->add(opt_multicgsolver);

		/* QPSOLVER */
		boost::program_options::options_description opt_qpsolver("QPSOLVER", console_nmb_cols);
		opt_qpsolver.add_options()
			("qpsolver_maxit", boost::program_options::value<int>(), "maximum number of iterations [int]")
			("qpsolver_eps", boost::program_options::value<double>(), "precision [double]")
			("qpsolver_debug_mode", boost::program_options::value<int>(), "debug mode [int]");
		opt_solvers->add(opt_qpsolver);

		/* CGQPSOLVER */
		boost::program_options::options_description opt_cgqpsolver("CGQPSOLVER", console_nmb_cols);
		opt_cgqpsolver.add_options()
			("cgqpsolver_maxit", boost::program_options::value<int>(), "maximum number of iterations [int]")
			("cgqpsolver_eps", boost::program_options::value<double>(), "precision [double]")
			("cgqpsolver_debug_mode", boost::program_options::value<int>(), "debug mode [int]");
		opt_solvers->add(opt_cgqpsolver);

		/* SPGQPSOLVER */
		boost::program_options::options_description opt_spgqpsolver("SPGQPSOLVER", console_nmb_cols);
		opt_spgqpsolver.add_options()
			("spgqpsolver_maxit", boost::program_options::value<int>(), "maximum number of iterations [int]")
			("spgqpsolver_eps", boost::program_options::value<double>(), "precision [double]")
			("spgqpsolver_m", boost::program_options::value<int>(), "parameter of generalized Armijo condition [int]")
			("spgqpsolver_gamma", boost::program_options::value<double>(), "parameter of generalized Armijo condition [double]")
			("spgqpsolver_sigma1", boost::program_options::value<double>(), "parameter of generalized Armijo condition [double]")
			("spgqpsolver_sigma2", boost::program_options::value<double>(), "parameter of generalized Armijo condition [double]")
			("spgqpsolver_alphainit", boost::program_options::value<double>(), "initial BB step-size [double]")
			("spgqpsolver_stop_normgp", boost::program_options::value<bool>(), "stopping criteria based on norm(gp) [bool]")
			("spgqpsolver_stop_Anormgp", boost::program_options::value<bool>(), "stopping criteria based on A-norm(gp) [bool]")
			("spgqpsolver_stop_normgp_normb", boost::program_options::value<bool>(), "stopping criteria based on norm(gp) and norm(b) [bool]")
			("spgqpsolver_stop_Anormgp_normb", boost::program_options::value<bool>(), "stopping criteria based on A-norm(gp) and norm(b) [bool]")
			("spgqpsolver_stop_difff", boost::program_options::value<bool>(), "stopping criteria based on difference of object function [bool]")
			("spgqpsolver_debug_mode", boost::program_options::value<int>(), "basic debug mode schema [0/1/2]")
			("spgqpsolver_debug_print_it", boost::program_options::value<bool>(), "print simple info about outer iterations")
			("spgqpsolver_debug_print_vectors", boost::program_options::value<bool>(), "print content of vectors during iterations")
			("spgqpsolver_debug_print_it", boost::program_options::value<bool>(), "print values of computed scalars during iterations");
		opt_solvers->add(opt_spgqpsolver);

	/* ----- MODELS ------ */
	boost::program_options::options_description opt_models("--- MODELS ---------------------", console_nmb_cols);
	description->add(opt_models);

		/* GRAPHH1FEMMODEL */
		boost::program_options::options_description opt_graphh1femmodel("GRAPHH1FEMMODEL", console_nmb_cols);
		opt_graphh1femmodel.add_options()
			("graphh1femmodel_matrix_type", boost::program_options::value<int>(), "type of used matrix [0=FREE/1=SPARSE]"); //TODO: enum?
		opt_models->add(opt_graphh1femmodel);

		/* KMEANSH1FEMMODEL */
		boost::program_options::options_description opt_kmeansh1femmodel("KMEANSH1FEMMODEL", console_nmb_cols);
		opt_kmeansh1femmodel.add_options()
			("kmeansh1femmodel_matrix_type", boost::program_options::value<int>(), "type of used matrix [0=FREE/1=SPARSE]"); //TODO: enum?
		opt_models->add(opt_kmeansh1femmodel);


}

}

#endif
