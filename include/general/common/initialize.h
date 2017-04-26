/** @file initialize.h
 *  @brief stuff for initialization
 *
 *  @author Lukas Pospisil
 */

#ifndef PASC_INITIALIZE_H
#define	PASC_INITIALIZE_H

#define RANDOM_BY_TIME false  /* if false, then random generator is initialized subject to time, else generated random data are always the same */

#include <string>
#include <vector>
#include "general/common/consoleinput.h"
#include "general/algebra/arrayoperation.h"

namespace pascinference {
using namespace algebra;	
	
namespace common {

/** @brief initialize the library
 * 
 *  Initialize random number generator. 
 *  If we are using Petsc, then Petsc is initialized.
 *
 * @todo process input arguments (use boost library?)
 */
template<class VectorBase> bool Initialize(int, char**){
	//TODO
} 

/** @brief finalize the library
 * 
 *  If we are using Petsc, then Petsc is finalized.
 *
 */
template<class VectorBase> void Finalize(){
	//TODO
}

template<class VectorBase> void allbarrier(){
	//TODO
}


}
} /* end of namespace */


#endif
