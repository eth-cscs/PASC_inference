#include "external/seqarrayvector/algebra/vector/generalvector.h"

namespace pascinference {
namespace algebra {

template<>
void GeneralVector<SeqArrayVector>::set_random() { 


	double *inner_array = this->get_array();

	for(int i=0;i<this->size();i++){
		inner_array[i] = rand();
	}

}

template<>
std::string GeneralVector<SeqArrayVector>::get_name() {
	return "SeqArrayVector";
}


}
} /* end of namespace */
