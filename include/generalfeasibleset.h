#ifndef PASCGENERALPROJECTION_H
#define	PASCGENERALPROJECTION_H

#include "common.h"
#include "algebra.h"

namespace pascinference {

/* the template says which vector we project; templated function cannot be virtual */
template<class VectorBase>
class GeneralFeasibleSet {
	protected:

	public:
		
		GeneralFeasibleSet() {};
		~GeneralFeasibleSet() {};

		virtual void print(std::ostream &output) const {};

		virtual void project(GeneralVector<VectorBase> &x) {
			//TODO: give error, the projection is not defined for this type of feasible set
		};

		template<class VectorBase2>
		friend std::ostream &operator<<(std::ostream &output, const GeneralFeasibleSet<VectorBase2> &feasibleset); /* cannot be virtual, therefore it call virtual print() */
	
};

/* general print, call virtual print() */
template<class VectorBase2>
std::ostream &operator<<(std::ostream &output, const GeneralFeasibleSet<VectorBase2> &feasibleset){
	if(DEBUG_MODE >= 100) std::cout << "(GeneralFeasibleSet)OPERATOR: <<" << std::endl;
	feasibleset.print(output);
	return output;
}

} /* end of namespace */

#endif
