/** @file tsmodel_global.h
 *  @brief class for manipulation with global time-series models
 *
 *  @author Lukas Pospisil
 */

#ifndef PASC_TSMODEL_GLOBAL_H
#define	PASC_TSMODEL_GLOBAL_H

#ifndef USE_PETSCVECTOR
 #error 'TSMODEL_GLOBAL is for PETSCVECTOR'
#endif

/* for debugging, if >= 100, then print info about ach called function */
extern int DEBUG_MODE;

#include <iostream>
#include "algebra.h"
#include "generalmodel.h"
#include "generalsolver.h"
#include "generaldata.h"

#include "data/tsdata_global.h"

namespace pascinference {

/* Maybe these classes are not defined yet */ 
class TSData_Global;

/** \class TSModel_Global_Global
 *  \brief General class for manipulation with global time-series models.
 *
*/
class TSModel_Global: public GeneralModel {
	protected:
		int T; /**< length of time-series, including xmem */
		int Tlocal; /**< local part of time-series */
		
		int xdim; /**< number of components in each time-step */

		int num; /**< length of K */
		int *K; /**< number of clusters to test */

		int Klocal; // TODO: this should be array

		int datavectorlength_global;
		int gammavectorlength_global;
		int thetavectorlength_global;
		int ulength_global;

		int datavectorlength_local;
		int gammavectorlength_local;
		int thetavectorlength_local;
		int ulength_local;

	public:
		TSModel_Global();
		TSModel_Global(int T, int xdim, int num, int *K);
		~TSModel_Global();

		virtual void print(ConsoleOutput &output) const;
		virtual std::string get_name() const;

		/* global length */
		virtual int get_datavectorlength_global();
		virtual int get_ulength_global();
		virtual int get_gammavectorlength_global();
		virtual int get_thetavectorlength_global();

		/* local length */
		virtual int get_datavectorlength_local();
		virtual int get_ulength_local();
		virtual int get_gammavectorlength_local();
		virtual int get_thetavectorlength_local();

		/* get common variables */
		int get_T() const;
		int get_Tlocal() const;
		int get_xdim() const;

		int get_num() const;
		int* get_K() const;
		int get_Klocal() const;
		
		/** @brief alloc memory for gamma solver
		 *  
		 *  Allocate memory for all data of gamma problem.
		 * 
		 */ 
		virtual void initialize_gammasolver(GeneralSolver **gammasolver, const TSData_Global *tsdata){
				*gammasolver = NULL;
		};

		/** @brief alloc memory for Theta solver
		 *  
		 *  Allocate memory for all data of Theta problem.
		 * 
		 */ 
		virtual void initialize_thetasolver(GeneralSolver **thetasolver, const TSData_Global *tsdata){
				*thetasolver = NULL;
		};

		/** @brief free memory of gamma solver
		 *  
		 *  Deallocate memory for all data of gamma problem.
		 * 
		 */ 
		virtual void finalize_gammasolver(GeneralSolver **gammasolver, const TSData_Global *tsdata){
		};

		/** @brief free memory of Theta solver
		 *  
		 *  Deallocate memory for all data of Theta problem.
		 * 
		 */ 
		virtual void finalize_thetasolver(GeneralSolver **thetasolver, const TSData_Global *tsdata){
		};

		/** @brief update data values of gamma solver
		 *  
		 *  Update data values of gamma solver, prepare it to the solving.
		 * 
		 */ 
		virtual void update_gammasolver(GeneralSolver *gammasolver, const TSData_Global *tsdata){
		};

		/** @brief update data values of Theta solver
		 *  
		 *  Update data values of Theta solver, prepare it to the solving.
		 * 
		 */ 
		virtual void update_thetasolver(GeneralSolver *thetasolver, const TSData_Global *tsdata){
		};

		/** @brief get the value of object function L(gamma,Theta,data)
		 *  
		 */ 
		virtual double get_L(GeneralSolver *gammasolver, GeneralSolver *thetasolver, const TSData_Global *tsdata){
			return std::numeric_limits<double>::max();
		};
};

} // end of namespace

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {

/* constructor */
TSModel_Global::TSModel_Global(){
	if(DEBUG_MODE >= 100) coutMaster << "(TSModel_Global)CONSTRUCTOR" << std::endl;

	/* set initial content */
	this->T = 0;
	this->xdim = 0;

	this->num = 0;
	this->K = NULL;

	this->Klocal = 0;

}


TSModel_Global::TSModel_Global(int newT, int newxdim, int newnum, int *newK){
	if(DEBUG_MODE >= 100) coutMaster << "(TSModel_Global)CONSTRUCTOR" << std::endl;

	/* set initial content */
	this->T = newT;
	this->xdim = newxdim;

	this->num = newnum;
	this->K = newK;
	
	this->Klocal = this->K[GlobalManager.get_rank()];
}

/* destructor */
TSModel_Global::~TSModel_Global(){
	if(DEBUG_MODE >= 100) coutMaster << "(TSModel_Global)DESTRUCTOR" << std::endl;
	
}


/* print info about model */
void TSModel_Global::print(ConsoleOutput &output) const {
	output <<  this->get_name() << std::endl;
	
	/* give information about presence of the data */
	output <<  " - nproc:  " << GlobalManager.get_size() << std::endl;
	
	output <<  " - T:      " << this->T << std::endl;
	output <<  " - xdim:    " << this->xdim << std::endl;

	output <<  " - K:      ";
	print_array(output, this->num, this->K);
	output << std::endl;
	output <<  " - Klocal: " << this->Klocal << std::endl;
		
}

/* get name of the model */
std::string TSModel_Global::get_name() const {
	return "General Time-Series-Global Model";	
}


/* ---- GET FUNCTIONS ---- */
int TSModel_Global::get_T() const{
	return T; 
}

int TSModel_Global::get_Tlocal() const{
	return Tlocal; 
}

int TSModel_Global::get_xdim() const{
	return xdim; 
}

int TSModel_Global::get_Klocal() const{
	return Klocal; 
}

int* TSModel_Global::get_K() const{
	return K; 
}

int TSModel_Global::get_num() const{
	return num; 
}


int TSModel_Global::get_datavectorlength_global(){
	return this->datavectorlength_global;
}

int TSModel_Global::get_gammavectorlength_global(){
	return this->gammavectorlength_global;
}

int TSModel_Global::get_thetavectorlength_global(){
	return this->thetavectorlength_global;
}

int TSModel_Global::get_ulength_global(){
	return this->ulength_global;
}

int TSModel_Global::get_datavectorlength_local(){
	return this->datavectorlength_local;
}

int TSModel_Global::get_gammavectorlength_local(){
	return this->gammavectorlength_local;
}

int TSModel_Global::get_thetavectorlength_local(){
	return this->thetavectorlength_local;
}

int TSModel_Global::get_ulength_local(){
	return this->ulength_local;
}




} /* end namespace */

#endif
