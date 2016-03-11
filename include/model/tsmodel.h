/** @file tsmodel.h
 *  @brief class for manipulation with time-series models
 *
 *  Header file which defines the parent class for manipulation with time-series models - additional information for solving the time-series problem.
 *  All specific model for time-series problem implementations should be defined as inherited classes from this class.
 *
 *  @author Lukas Pospisil
 */
 
#ifndef PASC_TSMODEL_H
#define	PASC_TSMODEL_H

/* for debugging, if >= 100, then print info about ach called function */
extern int DEBUG_MODE;

#include <iostream>
#include "algebra.h"
#include "generalmodel.h"
#include "generalsolver.h"
#include "generaldata.h"

#include "data/tsdata.h"

namespace pascinference {

/* Maybe these classes are not defined yet */ 
template<class VectorBase> class TSData;

/** \class TSModel
 *  \brief General class for manipulation with time-series models.
 *
 *  Parent class for manipulation with time-series models - additional information for solving the time-series problem by TSSolver.
 *  All specific time-series model implementations should be defined as inherited classes from this class.
 *	
*/
template<class VectorBase>
class TSModel: public GeneralModel {
	protected:
		int T; /**< length of time-series */
		int K; /**< number of clusters */
		int dim; /**< number of components in each time-step */

	public:
		TSModel();
		TSModel(int T, int dim, int K);
		~TSModel();

		virtual void print(std::ostream &output) const;
		virtual std::string get_name() const;

		virtual int get_datavectorlength();
		virtual int get_gammavectorlength();
		virtual int get_thetavectorlength();

		int get_T() const;
		int get_K() const;
		int get_dim() const;
		
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
	if(DEBUG_MODE >= 100) std::cout << "(TSModel)CONSTRUCTOR" << std::endl;

	/* set initial content */
	this->T = 0;
	this->K = 0;
	this->dim = 0;
}


template<class VectorBase>
TSModel<VectorBase>::TSModel(int newT, int newdim, int newK){
	if(DEBUG_MODE >= 100) std::cout << "(TSModel)CONSTRUCTOR" << std::endl;

	/* set initial content */
	this->T = newT;
	this->K = newK;
	this->dim = newdim;
}

/* destructor */
template<class VectorBase>
TSModel<VectorBase>::~TSModel(){
	if(DEBUG_MODE >= 100) std::cout << "(TSModel)DESTRUCTOR" << std::endl;
	
}


/* print info about model */
template<class VectorBase>
void TSModel<VectorBase>::print(std::ostream &output) const {
	output << this->get_name() << std::endl;
	
	/* give information about presence of the data */
	output << "  - T:    " << this->T << std::endl;
	output << "  - K:    " << this->K << std::endl;
	output << "  - dim:  " << this->dim << std::endl;
		
}

/* get name of the model */
template<class VectorBase>
std::string TSModel<VectorBase>::get_name() const {
	return "General Time-Series Model";	
}

/* get the appropriate length of datavector */
template<class VectorBase>
int TSModel<VectorBase>::get_datavectorlength(){
	return 0;
}

/* get the appropriate length of gammavector */
template<class VectorBase>
int TSModel<VectorBase>::get_gammavectorlength(){
	return 0;
}

/* get the appropriate length of thetavector */
template<class VectorBase>
int TSModel<VectorBase>::get_thetavectorlength(){
	return 0; 
}


/* ---- GET FUNCTIONS ---- */
template<class VectorBase>
int TSModel<VectorBase>::get_T() const{
	return T; 
}

template<class VectorBase>
int TSModel<VectorBase>::get_dim() const{
	return dim; 
}

template<class VectorBase>
int TSModel<VectorBase>::get_K() const{
	return K; 
}


} /* end namespace */

#endif
