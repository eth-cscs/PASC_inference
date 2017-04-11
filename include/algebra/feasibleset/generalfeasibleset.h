/** @file generalfeasibleset.h
 *  @brief class for manipulation with feasible set
 *
 *  @author Lukas Pospisil
 */

#ifndef PASC_GENERALFEASIBLESET_H
#define	PASC_GENERALFEASIBLESET_H

#include "common/logging.h"
#include "common/consoleoutput.h"
#include "algebra/vector/generalvector.h"


namespace pascinference {
using namespace common;
	
namespace algebra {

/** \class GeneralFeasibleSet
 *  \brief general class for manipulation with feasible sets
 *
 *  Parent class for manipulation with feasible sets.
 *  All specific feasible set implementations should be defined as inherited classes from this class.
 *	
*/
template<class VectorBase>
class GeneralFeasibleSet {
	protected:

	public:
		
		/** @brief default constructor
		 */ 
		GeneralFeasibleSet() {};

		/** @brief default destructor
		 */ 
		~GeneralFeasibleSet() {};

		/** @brief print information about feasible set
		 * 
		 * @param output where to print
		 */ 
		virtual void print(ConsoleOutput &output) const;

		/** @brief get the name of feasible set
		 */
		virtual std::string get_name() const;

		/** @brief project given point to feasible set
		 * 
		 * Take input vector and project it into feasible set. The input vector is rewritted with new values.
		 * 
		 * @param x input vector which should be projected
		 * @todo in this case, give an error since the projection is not defined
		 */ 
		virtual void project(GeneralVector<VectorBase> &x) {
			//TODO: give error, the projection is not defined for this type of feasible set
		};

		/** @brief append feasible set name to input stream
		 * 
		 */
		template<class VectorBase2>
		friend ConsoleOutput &operator<<(ConsoleOutput &output, const GeneralFeasibleSet<VectorBase2> &feasibleset); /* cannot be virtual, therefore it call virtual print() */
	
};

/* general print, call virtual print() */
template<class VectorBase2>
ConsoleOutput &operator<<(ConsoleOutput &output, const GeneralFeasibleSet<VectorBase2> &feasibleset){
	output << feasibleset.get_name();
	return output;
}

template<class VectorBase>
void GeneralFeasibleSet<VectorBase>::print(ConsoleOutput &output) const {
	output <<  this->get_name() << std::endl;
}

template<class VectorBase>
std::string GeneralFeasibleSet<VectorBase>::get_name() const {
	return "GeneralFeasibleSet";
}


}
} /* end of namespace */

#endif
