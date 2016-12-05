/** @file options.h
 *  @brief define command line options
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

	/* ----- LOG ---- */
	boost::program_options::options_description opt_log("#### LOG ########################", console_nmb_cols);
	opt_log.add_options()
		("log_or_not", boost::program_options::value<bool>(), "logging (writting into log file) is turned on/off [bool]")
		("log_or_not_func_call", boost::program_options::value<bool>(), "log LOG_FUNC_(STATIC)_BEGIN/LOG_FUNC_(STATIC)_END [bool]")
		("log_or_not_file_line", boost::program_options::value<bool>(), "log also the file and line of called log function [bool]")
		("log_or_not_level", boost::program_options::value<bool>(), "log also the level of called function [bool]")
		("log_or_not_memory", boost::program_options::value<bool>(), "log also the state of the memory [bool]");
	description->add(opt_log);

	/* ----- SOLVERS ------ */
	boost::program_options::options_description opt_solvers("#### SOLVERS ########################", console_nmb_cols);

		/* TSSOLVER */
		boost::program_options::options_description opt_tssolver("TSSOLVER", console_nmb_cols);
		opt_tssolver.add_options()
			("tssolver_maxit", boost::program_options::value<int>(), "maximum number of iterations [int]")
			("tssolver_eps", boost::program_options::value<double>(), "precision [double]")
			("tssolver_init_permute", boost::program_options::value<bool>(), "permute initial approximation subject to decomposition [bool]")
			("tssolver_debugmode", boost::program_options::value<int>(), "basic debug mode schema [0/1/2/3]")
			("tssolver_debug_print_annealing", boost::program_options::value<bool>(), "print info about annealing steps [bool]")
			("tssolver_debug_print_it", boost::program_options::value<bool>(), "print simple info about outer iterations [bool]")
			("tssolver_debug_print_theta", boost::program_options::value<bool>(), "print theta solver info [bool]")
			("tssolver_debug_print_theta_solution", boost::program_options::value<bool>(), "print solution of theta problem in each iteration [bool]")
			("tssolver_debug_print_gamma", boost::program_options::value<bool>(), "print gamma solver info [bool]")
			("tssolver_debug_print_gamma_solution", boost::program_options::value<bool>(), "print solution of gamma problem in each iteration [bool]")
			("tssolver_dump", boost::program_options::value<bool>(), "dump solver data [bool]");
		opt_solvers.add(opt_tssolver);

		/* MULTICG */
		boost::program_options::options_description opt_multicgsolver("MULTICGSOLVER", console_nmb_cols);
		opt_multicgsolver.add_options()
			("multicgsolver_maxit", boost::program_options::value<int>(), "maximum number of iterations [int]")
			("multicgsolver_eps", boost::program_options::value<double>(), "precision [double]")
			("multicgsolver_debugmode", boost::program_options::value<int>(), "debug mode [int]")
			("multicgsolver_dump", boost::program_options::value<bool>(), "dump solver data [bool]");
		opt_solvers.add(opt_multicgsolver);

		/* QPSOLVER */
		boost::program_options::options_description opt_qpsolver("QPSOLVER", console_nmb_cols);
		opt_qpsolver.add_options()
			("qpsolver_maxit", boost::program_options::value<int>(), "maximum number of iterations [int]")
			("qpsolver_eps", boost::program_options::value<double>(), "precision [double]")
			("qpsolver_debugmode", boost::program_options::value<int>(), "debug mode [int]")
			("qpsolver_dump", boost::program_options::value<bool>(), "dump solver data [bool]");
		opt_solvers.add(opt_qpsolver);

		/* CGQPSOLVER */
		boost::program_options::options_description opt_cgqpsolver("CGQPSOLVER", console_nmb_cols);
		opt_cgqpsolver.add_options()
			("cgqpsolver_maxit", boost::program_options::value<int>(), "maximum number of iterations [int]")
			("cgqpsolver_eps", boost::program_options::value<double>(), "precision [double]")
			("cgqpsolver_debugmode", boost::program_options::value<int>(), "debug mode [int]")
			("cgqpsolver_dump", boost::program_options::value<bool>(), "dump solver data [bool]");
		opt_solvers.add(opt_cgqpsolver);

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
			("spgqpsolver_debugmode", boost::program_options::value<int>(), "basic debug mode schema [0/1/2]")
			("spgqpsolver_debug_print_it", boost::program_options::value<bool>(), "print simple info about outer iterations")
			("spgqpsolver_debug_print_vectors", boost::program_options::value<bool>(), "print content of vectors during iterations")
			("spgqpsolver_debug_print_scalars", boost::program_options::value<bool>(), "print values of computed scalars during iterations")
			("spgqpsolver_dump", boost::program_options::value<bool>(), "dump solver data [bool]");
		opt_solvers.add(opt_spgqpsolver);

		/* PERMONSOLVER */
		boost::program_options::options_description opt_permonsolver("PERMONSOLVER", console_nmb_cols);
		opt_permonsolver.add_options()
			("permonsolver_maxit", boost::program_options::value<int>(), "maximum number of iterations [int]")
			("permonsolver_eps", boost::program_options::value<double>(), "precision [double]")
			("permonsolver_debugmode", boost::program_options::value<int>(), "debug mode [int]")
			("permonsolver_dump", boost::program_options::value<bool>(), "dump solver data [bool]");
		opt_solvers.add(opt_permonsolver);

	description->add(opt_solvers);


	/* ----- MODELS ------ */
	boost::program_options::options_description opt_models("#### MODELS ########################", console_nmb_cols);

		/* GRAPHH1FEMMODEL */
		boost::program_options::options_description opt_graphh1femmodel("GRAPHH1FEMMODEL", console_nmb_cols);
		opt_graphh1femmodel.add_options()
			("graphh1femmodel_gammasolvertype", boost::program_options::value<int>(), "type of used inner QP solver [0=SOLVER_AUTO/1=SOLVER_SPGQP/2=SOLVER_SPGQPCOEFF/3=SOLVER_PERMON]");
		opt_models.add(opt_graphh1femmodel);

		/* KMEANSH1FEMMODEL */
		boost::program_options::options_description opt_kmeansh1femmodel("KMEANSH1FEMMODEL", console_nmb_cols);
		opt_kmeansh1femmodel.add_options()
			("kmeansh1femmodel_matrix_type", boost::program_options::value<int>(), "type of used matrix [0=FREE/1=SPARSE]"); //TODO: enum?
		opt_models.add(opt_kmeansh1femmodel);

	description->add(opt_models);

}

}

#endif
