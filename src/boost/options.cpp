/** @file options.cpp
 *  @brief define command line options using boost/program_options
 *
 *  Implementation of options.h
 *
 *  @author Lukas Pospisil
 */

#include "options.h"

namespace pascinference {

void add_options(boost::program_options::options_description *description, int console_nmb_cols){

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
			("spgqpsolver_monitor", boost::program_options::value<bool>(), "export the descend of stopping criteria into .m file [bool]")
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
			("permonsolver_use_upperbound", boost::program_options::value<bool>(), "use additional upper bound x<=1 [bool]")
			("permonsolver_use_lambdamax", boost::program_options::value<bool>(), "provide the estimation of max eigen value to permon [bool]")
			("permonsolver_dump", boost::program_options::value<bool>(), "dump solver data [bool]");
		opt_solvers.add(opt_permonsolver);

		/* TAOSOLVER */
		boost::program_options::options_description opt_taosolver("TAOSOLVER", console_nmb_cols);
		opt_taosolver.add_options()
			("taosolver_maxit", boost::program_options::value<int>(), "maximum number of iterations [int]")
			("taosolver_eps", boost::program_options::value<double>(), "precision [double]")
			("taosolver_debugmode", boost::program_options::value<int>(), "debug mode [int]")
			("taosolver_use_upperbound", boost::program_options::value<bool>(), "use additional upper bound x<=1 [bool]")
			("taosolver_dump", boost::program_options::value<bool>(), "dump solver data [bool]");
		opt_solvers.add(opt_taosolver);

		/* ENTROPYSOLVERDLIB */
		boost::program_options::options_description opt_entropysolverdlib("ENTROPYSOLVERDLIB", console_nmb_cols);
		opt_entropysolverdlib.add_options()
			("entropysolverdlib_maxit", boost::program_options::value<int>(), "maximum number of Theta iterations [int]")
			("entropysolverdlib_eps", boost::program_options::value<double>(), "precision [double]")
			("entropysolverdlib_integration_eps", boost::program_options::value<double>(), "precision of integration [double]")
			("entropysolverdlib_integration_type", boost::program_options::value<int>(), "integration type [0=Vegas,1=Suave,2=Divonne,3=Cuhre]")
			("entropysolverdlib_debugmode", boost::program_options::value<int>(), "basic debug mode schema [0/1/2]")
			("entropysolverdlib_debug_print_moments", boost::program_options::value<bool>(), "print computed moments [bool]")
			("entropysolverdlib_debug_print_content", boost::program_options::value<bool>(), "print variables during optimization [bool]")
			("entropysolverdlib_debug_print_it", boost::program_options::value<bool>(), "print simple info about iterations [bool]")
			("entropysolverdlib_debug_print_integration", boost::program_options::value<bool>(), "print CUBA integration output [bool]");
		opt_solvers.add(opt_entropysolverdlib);

		/* ENTROPYSOLVERNEWTON */
		boost::program_options::options_description opt_entropysolvernewton("ENTROPYSOLVERNEWTON", console_nmb_cols);
		opt_entropysolvernewton.add_options()
			("entropysolvernewton_maxit", boost::program_options::value<int>(), "maximum number of iterations [int]")
			("entropysolvernewton_maxit_Axb", boost::program_options::value<int>(), "maximum number of inner KSP iterations [int]")
			("entropysolvernewton_eps", boost::program_options::value<double>(), "precision [double]")
			("entropysolvernewton_eps_Axb", boost::program_options::value<double>(), "precision of inner Ax=b solver [double]")
			("entropysolvernewton_newton_coeff", boost::program_options::value<double>(), "step-size coefficient of Newton update [double]")
			("entropysolvernewton_integration_eps", boost::program_options::value<double>(), "precision of integration [double]")
			("entropysolvernewton_integrationtype", boost::program_options::value<int>(), "type of numerical integration [0=INTEGRATION_AUTO/1=INTEGRATION_DLIB/2=INTEGRATION_MC]")
			("entropysolvernewton_monitor", boost::program_options::value<bool>(), "export the descend of stopping criteria into .m file [bool]")
			("entropysolvernewton_debugmode", boost::program_options::value<int>(), "basic debug mode schema [0/1/2]")
			("entropysolvernewton_debug_print_it", boost::program_options::value<bool>(), "print simple info about outer iterations [bool]")
			("entropysolvernewton_debug_print_Axb", boost::program_options::value<bool>(), "print info about inner Ax=b every outer Newton iteration [bool]")
			("entropysolvernewton_debug_print_vectors", boost::program_options::value<bool>(), "print content of vectors during iterations [bool]")
			("entropysolvernewton_debug_print_scalars", boost::program_options::value<bool>(), "print values of computed scalars during iterations [bool]");
		opt_solvers.add(opt_entropysolvernewton);

	description->add(opt_solvers);


	/* ----- MODELS ------ */
	boost::program_options::options_description opt_models("#### MODELS ########################", console_nmb_cols);

		/* GRAPHH1FEMMODEL */
		boost::program_options::options_description opt_graphh1femmodel("GRAPHH1FEMMODEL", console_nmb_cols);
		opt_graphh1femmodel.add_options()
			("graphh1femmodel_scalef", boost::program_options::value<bool>(), "scale function by 1/T [bool]")
			("graphh1femmodel_gammasolvertype", boost::program_options::value<int>(), "type of used inner QP solver [0=SOLVER_AUTO/1=SOLVER_SPGQP/2=SOLVER_SPGQPCOEFF/3=SOLVER_PERMON/4=SOLVER_TAO]");
		opt_models.add(opt_graphh1femmodel);

		/* KMEANSH1FEMMODEL */
		boost::program_options::options_description opt_kmeansh1femmodel("KMEANSH1FEMMODEL", console_nmb_cols);
		opt_kmeansh1femmodel.add_options()
			("kmeansh1femmodel_matrix_type", boost::program_options::value<int>(), "type of used matrix [0=FREE/1=SPARSE]"); //TODO: enum?
		opt_models.add(opt_kmeansh1femmodel);

		/* ENTROPYH1FEMMODEL */
		boost::program_options::options_description opt_entropyh1femmodel("ENTROPYH1FEMMODEL", console_nmb_cols);
		opt_entropyh1femmodel.add_options()
			("entropyh1femmodel_scalef", boost::program_options::value<bool>(), "scale function by 1/T [bool]")
			("entropyh1femmodel_thetasolvertype", boost::program_options::value<int>(), "type of used inner Entropy solver [0=SOLVER_AUTO/1=SOLVER_ENTROPY_DLIB/2=SOLVER_ENTROPY_NEWTON]")			
			("entropyh1femmodel_gammasolvertype", boost::program_options::value<int>(), "type of used inner QP solver [0=SOLVER_AUTO/1=SOLVER_SPGQP/2=SOLVER_SPGQPCOEFF/3=SOLVER_PERMON/4=SOLVER_TAO]");
		opt_models.add(opt_entropyh1femmodel);


	description->add(opt_models);


	/* ----- PETSC ---- */
#ifdef USE_PETSC
	boost::program_options::options_description opt_petsc("#### PETSC ######################", console_nmb_cols);
	opt_petsc.add_options()
		("petsc_options", boost::program_options::value< std::string >(), "all PETSc options [string]");
	description->add(opt_petsc);
#endif

}

}

