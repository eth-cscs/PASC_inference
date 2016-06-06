#ifndef PASC_OPTIONS_H
#define	PASC_OPTIONS_H

#include <boost/program_options.hpp>

namespace pascinference {

static void add_options(boost::program_options::options_description *description, int console_nmb_cols){
	/* TSSOLVER_GLOBAL */
	boost::program_options::options_description opt_tssolver_global("TSSOLVER_GLOBAL", console_nmb_cols);
	opt_tssolver_global.add_options()
		("tssolver_global_maxit", boost::program_options::value<int>(), "maximum number of iterations")
		("tssolver_global_eps", boost::program_options::value<double>(), "precision")
		("tssolver_global_debug_mode", boost::program_options::value<int>(), "debug mode");
	description->add(opt_tssolver_global);

	/* MULTICG */
	boost::program_options::options_description opt_multicgsolver("MULTICGSOLVER", console_nmb_cols);
	opt_multicgsolver.add_options()
		("multicgsolver_maxit", boost::program_options::value<int>(), "maximum number of iterations")
		("multicgsolver_eps", boost::program_options::value<double>(), "precision")
		("multicgsolver_debug_mode", boost::program_options::value<int>(), "debug mode");
	description->add(opt_multicgsolver);

	/* MULTICG_GLOBAL */
	boost::program_options::options_description opt_multicgsolver_global("MULTICGSOLVER_GLOBAL", console_nmb_cols);
	opt_multicgsolver_global.add_options()
		("multicgsolver_global_maxit", boost::program_options::value<int>(), "maximum number of iterations")
		("multicgsolver_global_eps", boost::program_options::value<double>(), "precision")
		("multicgsolver_global_debug_mode", boost::program_options::value<int>(), "debug mode");
	description->add(opt_multicgsolver_global);

	/* QPSOLVER */
	boost::program_options::options_description opt_qpsolver("QPSOLVER", console_nmb_cols);
	opt_qpsolver.add_options()
		("qpsolver_maxit", boost::program_options::value<int>(), "maximum number of iterations")
		("qpsolver_eps", boost::program_options::value<double>(), "precision")
		("qpsolver_debug_mode", boost::program_options::value<int>(), "debug mode");
	description->add(opt_qpsolver);

	/* CGQPSOLVER */
	boost::program_options::options_description opt_cgqpsolver("CGQPSOLVER", console_nmb_cols);
	opt_cgqpsolver.add_options()
		("cgqpsolver_maxit", boost::program_options::value<int>(), "maximum number of iterations")
		("cgqpsolver_eps", boost::program_options::value<double>(), "precision")
		("cgqpsolver_debug_mode", boost::program_options::value<int>(), "debug mode");
	description->add(opt_cgqpsolver);

	/* SPGQPSOLVER */
	boost::program_options::options_description opt_spgqpsolver("SPGQPSOLVER", console_nmb_cols);
	opt_spgqpsolver.add_options()
		("spgqpsolver_maxit", boost::program_options::value<int>(), "maximum number of iterations")
		("spgqpsolver_eps", boost::program_options::value<double>(), "precision")
		("spgqpsolver_debug_mode", boost::program_options::value<int>(), "debug mode")
		("spgqpsolver_m", boost::program_options::value<int>(), "parameter of generalized Armijo condition")
		("spgqpsolver_gamma", boost::program_options::value<double>(), "parameter of generalized Armijo condition")
		("spgqpsolver_sigma1", boost::program_options::value<double>(), "parameter of generalized Armijo condition")
		("spgqpsolver_sigma2", boost::program_options::value<double>(), "parameter of generalized Armijo condition")
		("spgqpsolver_alphainit", boost::program_options::value<double>(), "initial BB step-size")
		("spgqpsolver_stop_normgp", boost::program_options::value<bool>(), "stopping criteria based on norm(gp)")
		("spgqpsolver_stop_normgp_normb", boost::program_options::value<bool>(), "stopping criteria based on norm(gp) and norm(b)")
		("spgqpsolver_stop_difff", boost::program_options::value<bool>(), "stopping criteria based on difference of object function");
	description->add(opt_spgqpsolver);

	/* VARXH1FEMMODEL */
	boost::program_options::options_description opt_varxh1femmodel("VARXH1FEMMODEL", console_nmb_cols);
	opt_varxh1femmodel.add_options()
		("varxh1femmodel_t_scatter", boost::program_options::value<int>(), "scatter size in assembling Gamma problem, how long time-series to scatter to all processors");
	description->add(opt_varxh1femmodel);


}

}

#endif
