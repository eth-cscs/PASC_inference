/** @file tsmodel.h
 *  @brief class for manipulation with time-series models
 *
 *  @author Lukas Pospisil
 */

#ifndef PASC_TSMODEL_H
#define	PASC_TSMODEL_H

#ifndef USE_PETSCVECTOR
 #error 'TSMODEL is for PETSCVECTOR'
#endif

/* for debugging, if >= 100, then print info about ach called function */
extern int DEBUG_MODE;

#include "pascinference.h"
#include "data/tsdata.h"

namespace pascinference {

/* Maybe these classes are not defined yet */ 
namespace data {
	template<class VectorBase>
	class TSData;
}

namespace model {
/** \class TSModel
 *  \brief General class for manipulation with global time-series models.
 *
*/
template<class VectorBase>
class TSModel: public GeneralModel {
	protected:
		Decomposition *decomposition;
		TSData<VectorBase> *tsdata;

		int thetavectorlength_global; /**< global length of thetavector (model parameters) */
		int thetavectorlength_local; /**< local length of thetavector (model parameters) */
		
	public:
		TSModel();
		TSModel(TSData<VectorBase> &tsdata);
		~TSModel();

		virtual void print(ConsoleOutput &output) const;
		virtual std::string get_name() const;

		/* global length */
		virtual int get_thetavectorlength_global();
		virtual int get_thetavectorlength_local();


		virtual GeneralVector<VectorBase> *get_coordinatesVTK() const;
		virtual int get_coordinatesVTK_dim() const;

		/** @brief get value of AIC
		 *  
		 * @param L object function value
		 */ 
		virtual double get_aic(double L) const;
		
		/** @brief alloc memory for gamma solver
		 *  
		 *  Allocate memory for all data of gamma problem.
		 * 
		 */ 
		virtual void initialize_gammasolver(GeneralSolver **gammasolver){
				*gammasolver = NULL;
		};

		/** @brief alloc memory for Theta solver
		 *  
		 *  Allocate memory for all data of Theta problem.
		 * 
		 */ 
		virtual void initialize_thetasolver(GeneralSolver **thetasolver){
				*thetasolver = NULL;
		};

		/** @brief free memory of gamma solver
		 *  
		 *  Deallocate memory for all data of gamma problem.
		 * 
		 */ 
		virtual void finalize_gammasolver(GeneralSolver **gammasolver){
		};

		/** @brief free memory of Theta solver
		 *  
		 *  Deallocate memory for all data of Theta problem.
		 * 
		 */ 
		virtual void finalize_thetasolver(GeneralSolver **thetasolver){
		};

		/** @brief update data values of gamma solver
		 *  
		 *  Update data values of gamma solver, prepare it to the solving.
		 * 
		 */ 
		virtual void update_gammasolver(GeneralSolver *gammasolver){
		};

		/** @brief update data values of Theta solver
		 *  
		 *  Update data values of Theta solver, prepare it to the solving.
		 * 
		 */ 
		virtual void update_thetasolver(GeneralSolver *thetasolver){
		};

		/** @brief get the value of object function L(gamma,Theta,data)
		 *  
		 */ 
		virtual double get_L(GeneralSolver *gammasolver, GeneralSolver *thetasolver){
			return std::numeric_limits<double>::max();
		};

};


}
} /* end of namespace */

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {
namespace model {

/* constructor */
template<class VectorBase>
TSModel<VectorBase>::TSModel(){
	LOG_FUNC_BEGIN

	/* set initial content */
	this->decomposition = NULL;
	this->tsdata = NULL;

	this->thetavectorlength_global = 0;
	this->thetavectorlength_local = 0;

	LOG_FUNC_END
}

template<class VectorBase>
TSModel<VectorBase>::TSModel(TSData<VectorBase> &new_tsdata){
	LOG_FUNC_BEGIN

	/* set initial content */
	this->decomposition = NULL;
	this->tsdata = &new_tsdata;

	this->thetavectorlength_global = 0;
	this->thetavectorlength_local = 0;

	LOG_FUNC_END
}

/* destructor */
template<class VectorBase>
TSModel<VectorBase>::~TSModel(){
	LOG_FUNC_BEGIN
	
	LOG_FUNC_END
}


/* print info about model */
template<class VectorBase>
void TSModel<VectorBase>::print(ConsoleOutput &output) const {
	output <<  this->get_name() << std::endl;
	
	/* give information about presence of the data */
	output <<  " - nproc:  " << GlobalManager.get_size() << std::endl;
	
	output <<  " - T:      " << tsdata->get_T() << std::endl;
	output <<  " - xdim:    " << tsdata->get_xdim() << std::endl;

	output <<  " - K: " << tsdata->get_K() << std::endl;
		
}

/* get name of the model */
template<class VectorBase>
std::string TSModel<VectorBase>::get_name() const {
	return "General Time-Series-Global Model";	
}

template<class VectorBase>
int TSModel<VectorBase>::get_thetavectorlength_global(){
	return this->thetavectorlength_global;
}

template<class VectorBase>
int TSModel<VectorBase>::get_thetavectorlength_local(){
	return this->thetavectorlength_local;
}

template<class VectorBase>
GeneralVector<VectorBase> *TSModel<VectorBase>::get_coordinatesVTK() const{
	return NULL;
}

template<class VectorBase>
int TSModel<VectorBase>::get_coordinatesVTK_dim() const{
	return 0;
}

template<class VectorBase>
double TSModel<VectorBase>::get_aic(double L) const{
	return std::numeric_limits<double>::max();
}


}
} /* end namespace */

#endif
