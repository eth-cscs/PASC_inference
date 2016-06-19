/** @file tssolver_global.cu
 *  @brief test the varx global problem solver
 *
 *  Generate random varX problem and solve it using VarX global solver.
 *
 *  @author Lukas Pospisil
 */

#include "pascinference.h"
#include "solver/tssolver_global.h"
#include "data/tsdata_global.h"
#include "model/varxh1fem_global.h"

#ifndef USE_PETSCVECTOR
 #error 'This example is for PETSCVECTOR'
#endif
 
using namespace pascinference;

extern int pascinference::DEBUG_MODE;

int get_cluster_id_solution(int t, int T){
	int id_cluster = 0;
	double coeff = (double)T/3.0;
	if(t >= coeff && t < 2*coeff){
		id_cluster = 1;
	}
	if(t >= 2*coeff){
		id_cluster = 2;
	}

	return id_cluster;
}

int main( int argc, char *argv[] )
{
	/* add local program options */
	boost::program_options::options_description opt_problem("PROBLEM EXAMPLE", consoleArg.get_console_nmb_cols());
	opt_problem.add_options()
		("test_T", boost::program_options::value<int>(), "dimension of the problem")
		("test_K", boost::program_options::value<int>(), "number of clusters")
		("test_xmem", boost::program_options::value<int>(), "VarX parameter")
		("test_epssqr", boost::program_options::value<double>(), "penalty parameter");
	consoleArg.get_description()->add(opt_problem);

	/* call initialize */
	if(!Initialize(argc, argv)){
		return 0;
	} 

	int T, K, xmem; 
	double epssqr; 

	consoleArg.set_option_value("test_T", &T, 1000);
	consoleArg.set_option_value("test_K", &K, 3);
	consoleArg.set_option_value("test_xmem", &xmem, 2);
	consoleArg.set_option_value("test_epssqr", &epssqr, 10);
	
	/* start logging */
	std::ostringstream oss_name_of_file_log;
	oss_name_of_file_log << "results/varx_log_p" << GlobalManager.get_rank() << ".txt";
	logging.begin(oss_name_of_file_log.str());
		
	/* say hello */
	coutMaster << "- start program" << std::endl;

	/* parameters of the model */
	int xdim = 4; /* data dimension */

	/* solution - for generating the problem */
	int K_solution = 3;
	int xmem_solution = 2;
	double coeff = 100;//22564.5;
	
	double theta_solution[108] = {
		 0.0212/coeff,		 1.996973851371985, -0.000487555607403,  0.006182144617033,  0.016898325978254,		-0.997151392835171,  0.000428873596151, -0.006196363191544, -0.016953822805627,		/* K=1,n=1:  mu,A1,A2 */
		-0.0013/coeff,		-0.048771990907163,  1.976630233947632, -0.000172081430470,  0.009037772624055,		 0.049412051999825, -0.977231426020431, -0.000394989636220, -0.008493648930036,		/* K=1,n=2 */
		-0.0436/coeff,		-0.032574927768410, -0.008102611355981,  1.981548648130319, -0.029603366150073,		 0.033237925367826,  0.007963832173957, -0.981922239657675,  0.030099449853664,		/* K=1,n=3 */
		-0.0046/coeff,		 0.004790295629866,  0.001803589031281, -0.001290811055540,  1.992076550632660,		-0.004801590615226, -0.001728325124523,  0.001351998435471, -0.992133547352421,		/* K=1,n=4 */
		 0.0203/coeff,		 1.983547747060682, -0.003037379635604, -0.004431543622847, -0.021803280661815,		-0.984051793112738,  0.003171188444225,  0.004413610739813,  0.021811875222132,		/* K=2,n=1 */
		 0.0358/coeff,		 0.006569191979003,  1.991217087719501,  0.002636081266400,  0.038394708634095,		-0.006384387463868, -0.991557720861810, -0.002841142676061, -0.038556467916247,		/* K=2,n=2 */
		 0.0002/coeff,		-0.009379767734615, -0.015395627701527,  1.980417945671902, -0.041380583760087,		 0.009133059641365,  0.015608835060407, -0.980434238028737,  0.041389794617889,		/* K=2,n=3 */
		-0.0122/coeff,		-0.002906069136457, -0.008809207027806, -0.003884180425821,  1.966082655360234,		 0.002912153314880,  0.008935883850479,  0.003948339785004, -0.966081662311867,		/* K=2,n=4 */
		-0.0163/coeff,		 2.003149401792113, -0.014544606819796, -0.006667554360491,  0.018315146532669,		-1.003661473277484,  0.014877726915554,  0.007011849774012, -0.018321327225863,		/* K=3,n=1 */
		-0.0167/coeff,		 0.015538784060403,  1.975224192754936, -0.013711873751317, -0.010178213477031,		-0.015723160492140, -0.975077491864874,  0.014184789463673,  0.009954113017863,		/* K=3,n=2 */
		 0.0349/coeff,		-0.007494559444909,  0.007469536737359,  1.997085728291710,  0.016291022281772,		 0.007510037938227, -0.007649215150052, -0.997577936599111, -0.016138518043365,		/* K=3,n=3 */
		-0.0041/coeff,		 0.008782704798477, -0.006718795784519, -0.006974446738886,  1.963829240798505,		-0.008423618747776,  0.006554328766455,  0.007167678893911, -0.964105839991197		/* K=3,n=4 */			
	};

	double xstart_solution[8] = {
		50.0/coeff,	50.0503/coeff,	/* n=1 */
		70.0/coeff,	69.8995/coeff,	/* n=2 */
		70.0/coeff,	70.0001/coeff,	/* n=3 */
		90.0/coeff,	89.8995/coeff 	/* n=4 */
	};

	double noise_covariance[4] = {
		1, 1, 1, 2
	};

	/* prepare model */
	coutMaster << "--- PREPARING MODEL ---" << std::endl;
	VarxH1FEMModel_Global mymodel(T, xdim, K, xmem, epssqr);

	/* prepare time-series data */
	coutMaster << "--- PREPARING DATA ---" << std::endl;
	TSData_Global mydata(mymodel);

	/* generate some values to data */
	coutMaster << "--- GENERATING DATA ---" << std::endl;
	mymodel.generate_data(K_solution, xmem_solution, theta_solution, xstart_solution, &get_cluster_id_solution, &mydata, false);
	mymodel.generate_data_add_noise(&mydata, noise_covariance);

	/* prepare time-series solver */
	coutMaster << "--- PREPARING SOLVER ---" << std::endl;
	TSSolver_Global mysolver(mydata);

	mysolver.debug_mode = 2;
//	mysolver.print(coutMaster,coutAll);

	/* solve the problem */
	coutMaster << "--- SOLVING THE PROBLEM ---" << std::endl;
//	mymodel.set_solution_gamma(K_solution, xmem_solution, &get_cluster_id_solution, mydata.get_gammavector());
	mysolver.solve();

	/* print timers */
//	mysolver.printtimer(coutAll);
//	coutAll.synchronize();

	mydata.cut_gamma();

	/* save results into CSV file */
	coutMaster << "--- SAVING CSV ---" << std::endl;
	mydata.saveCSV("results/varx");

	/* say bye */	
	coutMaster << "- end program" << std::endl;

	logging.end();
	Finalize();

	return 0;
}
