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
template<class VectorBase>
class TSData;

/** \class TSModel
 *  \brief General class for manipulation with global time-series models.
 *
*/
template<class VectorBase>
class TSModel: public GeneralModel {
	protected:
		int T; /**< length of time-series, including xmem */
		int Tlocal; /**< local part of time-series */
		int Tbegin; /**< range begin of local part of time-series */
		int Tend; /**< range end of local part of time-series */
		
		int xdim; /**< number of components in each time-step */

		int K; /**< number of clusters */

		int datavectorlength_global; /**< global length of datavector (time-series values) */
		int gammavectorlength_global; /**< global length of gammavector (switching functions) */
		int thetavectorlength_global; /**< global length of thetavector (model parameters) */

		int datavectorlength_local; /**< local length of datavector (time-series values) */
		int gammavectorlength_local; /**< local length of gammavector (switching functions) */
		int thetavectorlength_local; /**< local length of thetavector (model parameters) */

	public:
		TSModel();
		~TSModel();

		virtual void print(ConsoleOutput &output) const;
		virtual std::string get_name() const;

		/* global length */
		virtual int get_datavectorlength_global();
		virtual int get_gammavectorlength_global();
		virtual int get_thetavectorlength_global();

		/* local length */
		virtual int get_datavectorlength_local();
		virtual int get_gammavectorlength_local();
		virtual int get_thetavectorlength_local();

		/* get common variables */
		int get_T() const;
		int get_Tlocal() const;
		int get_Tbegin() const;
		int get_Tend() const;
		int get_xdim() const;
		int get_xmem() const;
		virtual GeneralVector<VectorBase> *get_coordinatesVTK() const;
		virtual int get_coordinatesVTK_dim() const;


		int get_K() const;
		
		/** @brief alloc memory for gamma solver
		 *  
		 *  Allocate memory for all data of gamma problem.
		 * 
		 */ 
		virtual void initialize_gammasolver(GeneralSolver **gammasolver, const TSData<VectorBase> *tsdata){
				*gammasolver = NULL;
		};

		/** @brief alloc memory for Theta solver
		 *  
		 *  Allocate memory for all data of Theta problem.
		 * 
		 */ 
		virtual void initialize_thetasolver(GeneralSolver **thetasolver, const TSData<VectorBase> *tsdata){
				*thetasolver = NULL;
		};

		/** @brief free memory of gamma solver
		 *  
		 *  Deallocate memory for all data of gamma problem.
		 * 
		 */ 
		virtual void finalize_gammasolver(GeneralSolver **gammasolver, const TSData<VectorBase> *tsdata){
		};

		/** @brief free memory of Theta solver
		 *  
		 *  Deallocate memory for all data of Theta problem.
		 * 
		 */ 
		virtual void finalize_thetasolver(GeneralSolver **thetasolver, const TSData<VectorBase> *tsdata){
		};

		/** @brief update data values of gamma solver
		 *  
		 *  Update data values of gamma solver, prepare it to the solving.
		 * 
		 */ 
		virtual void update_gammasolver(GeneralSolver *gammasolver, const TSData<VectorBase> *tsdata){
		};

		/** @brief update data values of Theta solver
		 *  
		 *  Update data values of Theta solver, prepare it to the solving.
		 * 
		 */ 
		virtual void update_thetasolver(GeneralSolver *thetasolver, const TSData<VectorBase> *tsdata){
		};

		/** @brief get the value of object function L(gamma,Theta,data)
		 *  
		 */ 
		virtual double get_L(GeneralSolver *gammasolver, GeneralSolver *thetasolver, const TSData<VectorBase> *tsdata){
			return std::numeric_limits<double>::max();
		};

};

} // end of namespace

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {

/* constructor */
template<class VectorBase>
TSModel<VectorBase>::TSModel(){

	/* set initial content */
	this->T = 0;
	this->Tlocal = 0;
	this->Tbegin = 0;
	this->Tend = 0;
		
	this->xdim = 0;

	this->K = 0;

	this->datavectorlength_global = 0;
	this->gammavectorlength_global = 0;
	this->thetavectorlength_global = 0;

	this->datavectorlength_local = 0;
	this->gammavectorlength_local = 0;
	this->thetavectorlength_local = 0;

}

/* destructor */
template<class VectorBase>
TSModel<VectorBase>::~TSModel(){
	
}


/* print info about model */
template<class VectorBase>
void TSModel<VectorBase>::print(ConsoleOutput &output) const {
	output <<  this->get_name() << std::endl;
	
	/* give information about presence of the data */
	output <<  " - nproc:  " << GlobalManager.get_size() << std::endl;
	
	output <<  " - T:      " << this->T << std::endl;
	output <<  " - xdim:    " << this->xdim << std::endl;

	output <<  " - K: " << this->K << std::endl;
		
}

/* get name of the model */
template<class VectorBase>
std::string TSModel<VectorBase>::get_name() const {
	return "General Time-Series-Global Model";	
}

/* ---- GET FUNCTIONS ---- */
template<class VectorBase>
int TSModel<VectorBase>::get_T() const{
	return T; 
}

template<class VectorBase>
int TSModel<VectorBase>::get_Tlocal() const{
	return Tlocal; 
}

template<class VectorBase>
int TSModel<VectorBase>::get_Tbegin() const{
	return Tbegin;
}

template<class VectorBase>
int TSModel<VectorBase>::get_Tend() const{
	return Tend;
}

template<class VectorBase>
int TSModel<VectorBase>::get_xdim() const{
	return xdim; 
}

template<class VectorBase>
int TSModel<VectorBase>::get_K() const{
	return K; 
}

template<class VectorBase>
int TSModel<VectorBase>::get_datavectorlength_global(){
	return this->datavectorlength_global;
}

template<class VectorBase>
int TSModel<VectorBase>::get_gammavectorlength_global(){
	return this->gammavectorlength_global;
}

template<class VectorBase>
int TSModel<VectorBase>::get_thetavectorlength_global(){
	return this->thetavectorlength_global;
}

template<class VectorBase>
int TSModel<VectorBase>::get_datavectorlength_local(){
	return this->datavectorlength_local;
}

template<class VectorBase>
int TSModel<VectorBase>::get_gammavectorlength_local(){
	return this->gammavectorlength_local;
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



} /* end namespace */

#endif
