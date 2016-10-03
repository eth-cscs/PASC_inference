#ifndef PASC_COMMON_MEMORYCHECK_H
#define	PASC_COMMON_MEMORYCHECK_H

namespace pascinference {
namespace common {

/** @class MemoryCheck
 *  @brief get memory state
 * 
 *  Several utils for memory management. Could be used to get informations about used memory.
 * 
 */ 
class MemoryCheck {
	public:
		/** @brief get the size of virtual memory in bytes
		 * 
		 * The value of (sysinfo.totalram + sysinfo.totalswap) * sysinfo.mem_unit
		 * 
		 */ 
		static long long get_virtual_all() {
			struct sysinfo memInfo;
			sysinfo (&memInfo);
			long long totalVirtualMem = memInfo.totalram;
			totalVirtualMem += memInfo.totalswap;
			totalVirtualMem *= memInfo.mem_unit;
			
			return totalVirtualMem;
		}

		/** @brief get the size of used virtual memory in bytes
		 * 
		 */ 
		static long long get_virtual_used() {
			struct sysinfo memInfo;
			sysinfo (&memInfo);
		    long long virtualMemUsed = memInfo.totalram - memInfo.freeram;

			virtualMemUsed += memInfo.totalswap - memInfo.freeswap;
			virtualMemUsed *= memInfo.mem_unit;
			
			return virtualMemUsed;
		}

		/** @brief get the percentage usage of virtual memory
		 * 
		 */ 
		static double get_virtual(){
			struct sysinfo memInfo;
			sysinfo (&memInfo);
			
			return 100*((memInfo.totalram - memInfo.freeram) + (memInfo.totalram - memInfo.freeram))/(double)(memInfo.totalram + memInfo.totalram);
		}

		/** @brief get the size of physical memory in bytes
		 * 
		 */ 
		static long long get_physical() {
			struct sysinfo memInfo;
			sysinfo (&memInfo);
			long long totalPhysMem = memInfo.totalram;
			totalPhysMem *= memInfo.mem_unit;
			
			return totalPhysMem;
		}
	
};


}
} /* end of namespace */

#endif
