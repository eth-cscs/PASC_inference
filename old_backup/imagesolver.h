#ifndef PASC_IMAGESOLVER_H
#define	PASC_IMAGESOLVER_H

#ifndef USE_PETSCVECTOR
 #error 'IMAGE is for PETSCVECTOR'
 typedef petscvector::PetscVector PetscVector;
#endif


#include "pascinference.h"

#include "solver/tssolver.h"
#include "data/imagedata.h"
#include "model/tsmodel.h"

#define IMAGESOLVER_DEFAULT_MAXIT 1000
#define IMAGESOLVER_DEFAULT_EPS 0.001
#define IMAGESOLVER_DEFAULT_DEBUG_MODE 0

namespace pascinference {

/* ImageSolver */ 
template<class VectorBase>
class ImageSolver: public TSSolver<VectorBase> {
	protected:
		ImageData<VectorBase> *tsdata; /**< tsdata on which the solver operates */

	public:
		ImageSolver() : TSSolver<VectorBase> () {};
		ImageSolver(ImageData<VectorBase> &new_tsdata, int annealing = 1) : TSSolver<VectorBase> (new_tsdata, annealing) {}; 

		std::string get_name() const;

		ImageData<VectorBase> *get_data() const;

};

} // end of namespace

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {

template<class VectorBase>
std::string ImageSolver<VectorBase>::get_name() const {
	return "ImageSolver<VectorBase>";
}


template<class VectorBase>
ImageData<VectorBase> *ImageSolver<VectorBase>::get_data() const {
	return tsdata;
}

} /* end namespace */

#endif
