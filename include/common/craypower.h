#ifndef PASC_COMMON_CRAYPOWER_H
#define	PASC_COMMON_CRAYPOWER_H

namespace pascinference {
namespace common {

/** @class Craypower
 *  @brief get power state on cray machines
 * 
 *  Based on the code provided by Ben Cumming: https://github.com/bcumming/cray-power.git
 * 
 */ 
class Craypower {
	private:
		static double read_pm_file(const std::string &fname) {
			double result = 0.;
			std::ifstream fid(fname.c_str());
			fid >> result;
			return result;
		}

	public:
		static double device_energy(void) {
			return read_pm_file("/sys/cray/pm_counters/accel_energy");
		}

		static double energy() {
			return read_pm_file("/sys/cray/pm_counters/energy");
		}

		static double device_power() {
			return read_pm_file("/sys/cray/pm_counters/accel_power");
		}

		static double power() {
			return read_pm_file("/sys/cray/pm_counters/power");
		}

		static int num_nodes() {
			/* find out the number of nodes using slurm configuration */
			char *ptr = std::getenv("SLURM_JOB_NUM_NODES");
			if(ptr) {
				return atoi(ptr);
			} else {
				return -1;
			}
		}

		static int ranks_per_node() {
			/* returns -1 if unable to determine number of nodes */

			const int maxlen = 512;
			char name[maxlen];

			/* check whether MPI has been initialized */
			int is_initialized = 0;
			MPI_Initialized(&is_initialized);
			if(!is_initialized) {
				std::cerr << "ERROR : MPI not initialized" << std::endl;
				return -1;
			}

			/* get MPI information */
			int rank, size;
			MPI_Comm_size(MPI_COMM_WORLD, &size);
			MPI_Comm_rank(MPI_COMM_WORLD, &rank);

			/* get hostname for this node */
			int result = gethostname(name, maxlen);
			if(result) return -1;

			/* get integer index for this node, by stripping off first 3 characters
			 * on cray systems all compute nodes have hostname set as
			 * nid####### */
			int node = atoi(name+3);

			/* gather list of node identifiers */
			std::vector<int> node_ids(size);
			MPI_Allgather(&node, 1, MPI_INT, &node_ids[0], 1, MPI_INT, MPI_COMM_WORLD);

			/* count the number of mpi ranks that are on the same node as this rank */
			return std::count(node_ids.begin(), node_ids.end(), node);
		}
};

}
}

#endif
