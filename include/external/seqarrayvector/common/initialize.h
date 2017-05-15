#ifndef PASC_SEQARRAYVECTOR_INITIALIZE_H
#define	PASC_SEQARRAYVECTOR_INITIALIZE_H

#include "external/seqarrayvector/algebra/vector/generalvector.h"
#include "general/common/initialize.h"

namespace pascinference {
namespace common {

template<> bool Initialize<SeqArrayVector>(int argc, char *argv[]);
template<> void Finalize<SeqArrayVector>();
template<> void allbarrier<SeqArrayVector>();

}
} /* end of namespace */


#endif
