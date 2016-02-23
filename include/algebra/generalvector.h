#ifndef GENERALVECTOR_H
#define	GENERALVECTOR_H

namespace pascinference {

	/* define "all" stuff */
/*	enum all_type {
		all
	};
*/
	template<class VectorBase>
	class Vector : public VectorBase {
		
		public:
			/* constructors */
			Vector(): VectorBase() {}
			template<class ArgType> Vector(ArgType arg): VectorBase(arg) {}
			
			/* matrix-vector multiplication with General matrix */
			Vector<VectorBase> &operator=(GeneralMatrixRHS<VectorBase> rhs){
				rhs.matmult(*this);
				return *this;
			}

			/* all stuff */
/*			Vector<VectorBase> operator()(myALL myall) const {
				return *this;
			}
*/
	};

}

#endif
