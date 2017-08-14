#include "external/petscvector/algebra/integration/entropyintegrationcudavegas.h"

namespace pascinference {
namespace algebra {

template<> EntropyIntegrationCudaVegas<PetscVector>::EntropyIntegrationCudaVegas(EntropyData<PetscVector> *entropydata, double new_eps) : EntropyIntegration<PetscVector>(entropydata, new_eps) {
	LOG_FUNC_BEGIN

	/* load parameters from console */
	set_settings_from_console();

	/* prepare external content */
	#ifdef USE_CUDA
		externalcontent = new ExternalContent();
	
		externalcontent->ncall = 10000; 
		externalcontent->itmx = 10;
		externalcontent->acc = 0.0001; // accuracy
		externalcontent->nBlockSize = 256; // CUDA size of the block
		externalcontent->ndim = entropydata->get_xdim(); 

	#endif
	
	LOG_FUNC_END
}

template<> EntropyIntegrationCudaVegas<PetscVector>::~EntropyIntegrationCudaVegas() {
	LOG_FUNC_BEGIN
	
	/* destroy external content */
	#ifdef USE_CUDA
		free(externalcontent);
	#endif
	
	LOG_FUNC_END
}

template<>
void EntropyIntegrationCudaVegas<PetscVector>::compute(double *integrals_arr, double *lambda_arr, int Km_max) {
	LOG_FUNC_BEGIN

	#ifdef USE_CUDA
		/* call appropriate cuda kernel for integral computation */
		double avgi = 0.;
		double sd = 0.;
		double chi2a = 0.;

		timer.start();
		externalcontent->cuda_gVegas(avgi, sd, chi2a);
		timer.stop();

		double timeTotal = timer.get_value_sum();
		double timeVegasCall = externalcontent->timerVegasCall.get_value_sum();
		double timeVegasMove = externalcontent->timerVegasMove.get_value_sum();
		double timeVegasFill = externalcontent->timerVegasFill.get_value_sum();
		double timeVegasRefine = externalcontent->timerVegasRefine.get_value_sum();

		coutMaster.clear();
		coutMaster << std::setw(10)<<std::setprecision(6)<<std::endl;
		coutMaster << "#==================================" << std::endl;
		coutMaster << "# No. of Thread Block Size  : " << externalcontent->nBlockSize << std::endl;
		coutMaster << "#==================================" << std::endl;
		coutMaster << "# No. of dimensions         : " << externalcontent->ndim << std::endl;
		coutMaster << "# No. of func calls / iter  : " << externalcontent->ncall << std::endl;
		coutMaster << "# No. of max. iterations    : " << externalcontent->itmx << std::endl;
		coutMaster << "# Desired accuracy          : " << externalcontent->acc << std::endl;
		coutMaster << "#==================================" << std::endl;
		coutMaster << std::scientific;
		coutMaster << std::left << std::setfill(' ');
		coutMaster << "# Result                    : "
            << std::setw(12) << std::setprecision(5) << avgi << " +- "
            << std::setw(12) << std::setprecision(5) << sd <<" ( "
            << std::setw(7) << std::setprecision(4)
            << std::fixed << 100.*sd/avgi << "%)" << std::endl;
		coutMaster << std::fixed;
		coutMaster << "# Chisquare                 : " << std::setprecision(4) << chi2a << std::endl;
		coutMaster << "#==================================" << std::endl;
		coutMaster << std::right;
		coutMaster << "# Total Execution Time(sec) : " << std::setw(10) << std::setprecision(4) << timeTotal << std::endl;
		coutMaster << "#==================================" << std::endl;
		coutMaster << "# Time for func calls (sec) : "
            << std::setw(10) << std::setprecision(4) << timeVegasCall
            << " ( " << std::setw(5) << std::setprecision(2)
            << 100.*timeVegasCall/timeTotal << "%)" << std::endl;
		coutMaster << "# Time for data transf (sec): "
            << std::setw(10) << std::setprecision(4) << timeVegasMove
            << " ( " << std::setw(5) << std::setprecision(2)
            << 100.*timeVegasMove/timeTotal << "%)" << std::endl;
		coutMaster << "# Time for data fill (sec)  : "
            << std::setw(10) << std::setprecision(4) << timeVegasFill
            << " ( " << std::setw(5) << std::setprecision(2)
            << 100.*timeVegasFill/timeTotal << "%)" << std::endl;
		coutMaster << "# Time for grid refine (sec): "
            << std::setw(10) << std::setprecision(4) << timeVegasRefine
            << " ( "<< std::setw(5) << std::setprecision(2)
            << 100.*timeVegasRefine/timeTotal << "%)" << std::endl;
		coutMaster << "#==================================" << std::endl;

	#else
		//TODO: throw error - cannot use this method without CUDA
	#endif
	
	LOG_FUNC_END
}



}
}

