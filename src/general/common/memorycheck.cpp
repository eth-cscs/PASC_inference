#include "common/memorycheck.h"

namespace pascinference {
namespace common {

long long MemoryCheck::get_virtual_all() {
	struct sysinfo memInfo;
	sysinfo (&memInfo);
	long long totalVirtualMem = memInfo.totalram;
	totalVirtualMem += memInfo.totalswap;
	totalVirtualMem *= memInfo.mem_unit;
		
	return totalVirtualMem;
}

long long MemoryCheck::get_virtual_used() {
	struct sysinfo memInfo;
	sysinfo (&memInfo);
    long long virtualMemUsed = memInfo.totalram - memInfo.freeram;

	virtualMemUsed += memInfo.totalswap - memInfo.freeswap;
	virtualMemUsed *= memInfo.mem_unit;
			
	return virtualMemUsed;
}

double MemoryCheck::get_virtual(){
	struct sysinfo memInfo;
	sysinfo (&memInfo);
			
	return 100*((memInfo.totalram - memInfo.freeram) + (memInfo.totalram - memInfo.freeram))/(double)(memInfo.totalram + memInfo.totalram);
}

long long MemoryCheck::get_physical() {
	struct sysinfo memInfo;
	sysinfo (&memInfo);
	long long totalPhysMem = memInfo.totalram;
	totalPhysMem *= memInfo.mem_unit;
			
	return totalPhysMem;
}
	

}
} /* end of namespace */

