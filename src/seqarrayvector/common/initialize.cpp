#include "external/seqarrayvector/common/initialize.h"

namespace pascinference {
namespace common {

template<>
bool Initialize<SeqArrayVector>(int argc, char *argv[]){
	/* console arguments */
	if(!consoleArg.init(argc,argv)){
		return false;
	}

	/* initialize random seed: */
	if(RANDOM_BY_TIME){
		srand(time(NULL));
	} else {
		srand(0);
	}
	
	return true;
}

template<>
void Finalize<SeqArrayVector>(){
}

template<>
void allbarrier<SeqArrayVector>() {
}


}
} /* end of namespace */
