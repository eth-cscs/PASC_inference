#ifndef PASC_SEQARRAYVECTOR_GENERALVECTOR_H
#define	PASC_SEQARRAYVECTOR_GENERALVECTOR_H

#include "general/algebra/vector/generalvector.h"
#include "external/seqarrayvector/algebra/vector/seqarrayvector.h"

typedef seqarrayvector::SeqArrayVector SeqArrayVector;

namespace pascinference {
namespace algebra {

template<> void GeneralVector<SeqArrayVector>::set_random();
template<> std::string GeneralVector<SeqArrayVector>::get_name();


}
} /* end of namespace */

#endif
