#ifndef PASC_OPTIONS_H
#define	PASC_OPTIONS_H

#include <boost/program_options.hpp>

namespace pascinference {

static void add_options(boost::program_options::options_description *description, int console_nmb_cols){
	/* TSSOLVER */
	boost::program_options::options_description opt_tssolver_global("TSSOLVER", console_nmb_cols);
	opt_tssolver_global.add_options()
		("tssolver_maxit", boost::program_options::value<int>(), "maximum number of iterations [int]")
		("tssolver_eps", boost::program_options::value<double>(), "precision [double]")
		("tssolver_debug_mode", boost::program_options::value<int>(), "debug mode [int]");
	description->add(opt_tssolver_global);

	/* MULTICG */
	boost::program_options::options_description opt_multicgsolver("MULTICGSOLVER", console_nmb_cols);
	opt_multicgsolver.add_options()
		("multicgsolver_maxit", boost::program_options::value<int>(), "maximum number of iterations [int]")
		("multicgsolver_eps", boost::program_options::value<double>(), "precision [double]")
		("multicgsolver_debug_mode", boost::program_options::value<int>(), "debug mode [int]");
	description->add(opt_multicgsolver);

	/* MULTICG_GLOBAL */
	boost::program_options::options_description opt_multicgsolver_global("MULTICGSOLVER_GLOBAL", console_nmb_cols);
	opt_multicgsolver_global.add_options()
		("multicgsolver_global_maxit", boost::program_options::value<int>(), "maximum number of iterations [int]")
		("multicgsolver_global_eps", boost::program_options::value<double>(), "precision [double]")
		("multicgsolver_global_debug_mode", boost::program_options::value<int>(), "debug mode [int]");
	description->add(opt_multicgsolver_global);

	/* QPSOLVER */
	boost::program_options::options_description opt_qpsolver("QPSOLVER", console_nmb_cols);
	opt_qpsolver.add_options()
		("qpsolver_maxit", boost::program_options::value<int>(), "maximum number of iterations [int]")
		("qpsolver_eps", boost::program_options::value<double>(), "precision [double]")
		("qpsolver_debug_mode", boost::program_options::value<int>(), "debug mode [int]");
	description->add(opt_qpsolver);

	/* CGQPSOLVER */
	boost::program_options::options_description opt_cgqpsolver("CGQPSOLVER", console_nmb_cols);
	opt_cgqpsolver.add_options()
		("cgqpsolver_maxit", boost::program_options::value<int>(), "maximum number of iterations [int]")
		("cgqpsolver_eps", boost::program_options::value<double>(), "precision [double]")
		("cgqpsolver_debug_mode", boost::program_options::value<int>(), "debug mode [int]");
	description->add(opt_cgqpsolver);

	/* SPGQPSOLVER */
	boost::program_options::options_description opt_spgqpsolver("SPGQPSOLVER", console_nmb_cols);
	opt_spgqpsolver.add_options()
		("spgqpsolver_maxit", boost::program_options::value<int>(), "maximum number of iterations [int]")
		("spgqpsolver_eps", boost::program_options::value<double>(), "precision [double]")
		("spgqpsolver_debug_mode", boost::program_options::value<int>(), "debug mode [int]")
		("spgqpsolver_m", boost::program_options::value<int>(), "parameter of generalized Armijo condition [int]")
		("spgqpsolver_gamma", boost::program_options::value<double>(), "parameter of generalized Armijo condition [double]")
		("spgqpsolver_sigma1", boost::program_options::value<double>(), "parameter of generalized Armijo condition [double]")
		("spgqpsolver_sigma2", boost::program_options::value<double>(), "parameter of generalized Armijo condition [double]")
		("spgqpsolver_alphainit", boost::program_options::value<double>(), "initial BB step-size [double]")
		("spgqpsolver_stop_normgp", boost::program_options::value<bool>(), "stopping criteria based on norm(gp) [bool]")
		("spgqpsolver_stop_Anormgp", boost::program_options::value<bool>(), "stopping criteria based on A-norm(gp) [bool]")
		("spgqpsolver_stop_normgp_normb", boost::program_options::value<bool>(), "stopping criteria based on norm(gp) and norm(b) [bool]")
		("spgqpsolver_stop_Anormgp_normb", boost::program_options::value<bool>(), "stopping criteria based on A-norm(gp) and norm(b) [bool]")
		("spgqpsolver_stop_difff", boost::program_options::value<bool>(), "stopping criteria based on difference of object function [bool]");
	description->add(opt_spgqpsolver);

	/* GRAPHH1FEMMODEL */
	boost::program_options::options_description opt_graphh1femmodel("GRAPHH1FEMMODEL", console_nmb_cols);
	opt_graphh1femmodel.add_options()
		("graphh1femmodel_matrix_type", boost::program_options::value<int>(), "type of used matrix [0=FREE/1=SPARSE]"); //TODO: enum?
	description->add(opt_graphh1femmodel);

	/* KMEANSH1FEMMODEL */
	boost::program_options::options_description opt_kmeansh1femmodel("KMEANSH1FEMMODEL", console_nmb_cols);
	opt_kmeansh1femmodel.add_options()
		("kmeansh1femmodel_matrix_type", boost::program_options::value<int>(), "type of used matrix [0=FREE/1=SPARSE]"); //TODO: enum?
	description->add(opt_kmeansh1femmodel);


}

}

#endif
