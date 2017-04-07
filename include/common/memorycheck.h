#ifndef PASC_COMMON_MEMORYCHECK_H
#define	PASC_COMMON_MEMORYCHECK_H

#include <sys/sysinfo.h>

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
		static long long get_virtual_all();

		/** @brief get the size of used virtual memory in bytes
		 * 
		 */ 
		static long long get_virtual_used();

		/** @brief get the percentage usage of virtual memory
		 * 
		 */ 
		static double get_virtual();

		/** @brief get the size of physical memory in bytes
		 * 
		 */ 
		static long long get_physical();
	
};


}
} /* end of namespace */

#endif
