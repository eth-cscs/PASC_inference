#ifndef PASC_COMMON_POWERCHECK_H
#define	PASC_COMMON_POWERCHECK_H

#ifdef USE_CRAYPOWER
 #include "common/craypower.h"
#endif

namespace pascinference {
namespace common {

/** @class PowerCheck
 *  @brief get power state (to measure power consumption)
 * 
 *  @todo works only on CRAY system (like Piz Daint)
 * 
 */ 
class PowerCheck {
	public:

		static double get_node_energy() {
			double node_energy = -1.0;
			
			#ifdef USE_CRAYPOWER
				node_energy = Craypower::energy();
			#endif
			
			return node_energy;
		}

		static double get_device_energy() {
			double device_energy = -1.0;
			
			#ifdef USE_CRAYPOWER
				device_energy = Craypower::device_energy();
			#endif
			
			return device_energy;
		}

		static int get_ranks_per_node() {
			int ranks_per_node = 1;

			#ifdef USE_CRAYPOWER
				ranks_per_node = Craypower::ranks_per_node();
			#endif
			
			return ranks_per_node;
		}
		
		static double mpi_sum_reduce(double local_value) {
			double global_sum = 0.;
			MPI_Reduce(&local_value, &global_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			return global_sum;
		}
	
};


}
} /* end of namespace */

#endif
