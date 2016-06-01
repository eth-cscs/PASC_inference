/** @file qpdata.h
 *  @brief class for manipulation with QP data
 *
 *  @author Lukas Pospisil
 */

#ifndef PASC_QPDATA_H
#define	PASC_QPDATA_H

/* for debugging, if >= 100, then print info about each called function */
extern int DEBUG_MODE;

#include <iostream>
#include "generaldata.h"
#include "generalfeasibleset.h"
#include "algebra.h"

namespace pascinference {

/** \class QPData
 *  \brief general class for manipulation with QP data
 *
 *  Defines the data for quadratic programming problem, i.e.
 * \f[ \min \frac{1}{2} \langle Ax,x \rangle - \langle b,x \rangle ~\mathrm{s.t.}~ x \in feasible~set \f]
 *	
*/
template<class VectorBase>
class QPData: public GeneralData {
	private:
		GeneralMatrix<VectorBase> *A; /**< Hessian matrix */
		GeneralVector<VectorBase> *b; /**< RHS vector, linear term */
		GeneralVector<VectorBase> *x0; /**< initial approximation */
		GeneralVector<VectorBase> *x; /**< solution */
		GeneralFeasibleSet<VectorBase> *feasibleset; /**< feasible set */

	public:
	
		/** @brief default constructor
		*/ 
		QPData();
		
		/** @brief default destructor
		 */ 
		~QPData();

		/** @brief print basic information about data
		 * 
		 * @param output where to print
		 */
		void print(ConsoleOutput &output) const;

		void print(ConsoleOutput &output_global, ConsoleOutput &output_local) const;
		
		/** @brief print content of data
		 * 
		 * Print all values in inner data structures.
		 * 
		 * @param output where to print
		 */
		void printcontent(ConsoleOutput &output) const;

		void printcontent(ConsoleOutput &output_global, ConsoleOutput &output_local) const;

		
		/** @brief get the name of data
		 */ 		
		std::string get_name() const;

		/** @brief set new Hessian matrix
		 */ 
		void set_A(GeneralMatrix<VectorBase> *A);

		/** @brief get Hessian matrix
		 */ 
		GeneralMatrix<VectorBase> *get_A() const;

		/** @brief set new linear term
		 */ 
		void set_b(GeneralVector<VectorBase> *b);

		/** @brief get linear term
		 */ 
		GeneralVector<VectorBase> *get_b() const;

		/** @brief set new initial approximation
		 */ 
		void set_x0(GeneralVector<VectorBase> *x0);

		/** @brief get initial approximation
		 */ 
		GeneralVector<VectorBase> *get_x0() const;

		/** @brief set new solution vector
		 */ 
		void set_x(GeneralVector<VectorBase> *x);

		/** @brief get solution vector
		 */ 
		GeneralVector<VectorBase> *get_x() const;

		/** @brief set new feasible set
		 */ 
		void set_feasibleset(GeneralFeasibleSet<VectorBase> *x0);

		/** @brief get feasible set
		 */ 
		GeneralFeasibleSet<VectorBase> *get_feasibleset() const;

};

} /* end of namespace */



/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {

/* constructor */
template<class VectorBase>
QPData<VectorBase>::QPData(){
	LOG_FUNC_BEGIN

	/* set initial content */
	this->A = NULL;
	this->b = NULL;
	this->x0 = NULL;
	this->x = NULL;
	this->feasibleset = NULL;

	LOG_FUNC_END
}

/* destructor */
template<class VectorBase>
QPData<VectorBase>::~QPData(){
	LOG_FUNC_BEGIN

	LOG_FUNC_END
}


/* print info about problem */
template<class VectorBase>
void QPData<VectorBase>::print(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output <<  this->get_name() << std::endl;
	
	/* give information about presence of the data */
	output <<  " - A:            ";
	if(this->A){
		output << "YES (" << this->A->get_name() << ")" << std::endl; // TODO: get matrix name 
	} else {
		output << "not set" << std::endl;
	}
	output <<  " - b:            ";
	if(this->b){
		output << "YES (size: " << this->b->size() << ")" << std::endl;
	} else {
		output << "not set" << std::endl;
	}
	output <<  " - x0:           ";
	if(this->x0){
		output << "YES (size: " << this->x0->size() << ")" << std::endl;
	} else {
		output << "not set" << std::endl;
	}
	output <<  " - x:            ";
	if(this->x){
		output << "YES (size: " << this->x->size() << ")" << std::endl;
	} else {
		output << "not set" << std::endl;
	}
	output <<  " - feasible_set: ";
	if(this->feasibleset){
		output << "YES (" << this->feasibleset->get_name() << ")" << std::endl;
	} else {
		output << "not set" << std::endl;
	}
		
	LOG_FUNC_END
}

template<class VectorBase>
void QPData<VectorBase>::print(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN

	output_global <<  this->get_name() << std::endl;
	
	/* give information about presence of the data */
	output_global <<  " - A:            ";
	if(this->A){
		output_global << "YES (" << this->A->get_name() << ")" << std::endl; // TODO: get matrix name 
		output_global.push();
		this->A->print(output_local);
		output_global.pop();
		output_local.synchronize();
	} else {
		output_global << "not set" << std::endl;
	}

	output_global <<  " - b:            ";
	if(this->b){
		output_global << "YES (size: " << this->b->size() << ")" << std::endl;
		output_global.push();
		output_local  <<  "local size: " << this->b->local_size() << std::endl;
		output_global.pop();
		output_local.synchronize();		
	} else {
		output_global << "not set" << std::endl;
	}

	output_global <<  " - x0:           ";
	if(this->x0){
		output_global << "YES (size: " << this->x0->size() << ")" << std::endl;
		output_global.push();
		output_local  <<  "local size: " << this->x0->local_size() << std::endl;
		output_global.pop();
		output_local.synchronize();		
	} else {
		output_global << "not set" << std::endl;
	}
	
	output_global <<  " - x:            ";
	if(this->x){
		output_global << "YES (size: " << this->x->size() << ")" << std::endl;
		output_global.push();
		output_local  <<  "local size: " << this->x->local_size() << std::endl;
		output_global.pop();
		output_local.synchronize();		
	} else {
		output_global << "not set" << std::endl;
	}
	
	output_global <<  " - feasible_set: ";
	if(this->feasibleset){
		output_global << "YES (" << this->feasibleset->get_name() << ")" << std::endl;
		output_global.push();
		this->feasibleset->print(output_local);
		output_global.pop();
		output_local.synchronize();		
	} else {
		output_global << "not set" << std::endl;
	}
		
	output_global.synchronize();

	LOG_FUNC_END
}

template<class VectorBase>
void QPData<VectorBase>::printcontent(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output <<  this->get_name() << std::endl;
	
	/* give information about presence of the data */
	output <<  " - A:            ";
	if(this->A){
		this->A->print(output); 
	} else {
		output << "not set" << std::endl;
	}
	output <<  " - b:            ";
	if(this->b){
		output << *(this->b) << std::endl;
	} else {
		output << "not set" << std::endl;
	}
	output <<  " - x0:           ";
	if(this->x0){
		output << *(this->x0) << std::endl;
	} else {
		output << "not set" << std::endl;
	}
	output <<  " - x:            ";
	if(this->x0){
		output << *(this->x) << std::endl;
	} else {
		output << "not set" << std::endl;
	}
	output <<  " - feasible_set: ";
	if(this->feasibleset){
		output << *(this->feasibleset) << std::endl;
	} else {
		output << "not set" << std::endl;
	}
		
	LOG_FUNC_END
}

template<class VectorBase>
void QPData<VectorBase>::printcontent(ConsoleOutput &output_global,ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN

	output_global <<  this->get_name() << std::endl;
		
	/* give information about presence of the data */
	output_local <<  " - A:            ";
	if(this->A){
		this->A->print(output_local); 
	} else {
		output_local << "not set" << std::endl;
	}
	output_local.synchronize();
	
	output_local <<  " - b:            ";
	if(this->b){
		output_local << *(this->b) << std::endl;
	} else {
		output_local << "not set" << std::endl;
	}
	output_local.synchronize();
	
	output_local <<  " - x0:           ";
	if(this->x0){
		output_local << *(this->x0) << std::endl;
	} else {
		output_local << "not set" << std::endl;
	}
	output_local.synchronize();
	
	output_local <<  " - x:            ";
	if(this->x0){
		output_local << *(this->x) << std::endl;
	} else {
		output_local << "not set" << std::endl;
	}
	output_local.synchronize();
	
	output_local <<  " - feasible_set: ";
	if(this->feasibleset){
		output_local << *(this->feasibleset) << std::endl;
	} else {
		output_local << "not set" << std::endl;
	}
	output_local.synchronize();
		
	LOG_FUNC_END
}

template<class VectorBase>
std::string QPData<VectorBase>::get_name() const {
	return "QPData";
}

/* ----- SET and GET functions --- */
template<class VectorBase>
void QPData<VectorBase>::set_A(GeneralMatrix<VectorBase> *newA){
	this->A = newA;
}

template<class VectorBase>
GeneralMatrix<VectorBase> *QPData<VectorBase>::get_A() const{
	return this->A;
}

template<class VectorBase>
void QPData<VectorBase>::set_b(GeneralVector<VectorBase> *b){
	this->b = b;
}

template<class VectorBase>
GeneralVector<VectorBase> *QPData<VectorBase>::get_b() const{
	return this->b;
}

template<class VectorBase>
void QPData<VectorBase>::set_x0(GeneralVector<VectorBase> *x0){
	this->x0 = x0;
}

template<class VectorBase>
GeneralVector<VectorBase> *QPData<VectorBase>::get_x0() const{
	return this->x0;
}

template<class VectorBase>
void QPData<VectorBase>::set_x(GeneralVector<VectorBase> *x){
	this->x = x;
}

template<class VectorBase>
GeneralVector<VectorBase> *QPData<VectorBase>::get_x() const{
	return this->x;
}

template<class VectorBase>
void QPData<VectorBase>::set_feasibleset(GeneralFeasibleSet<VectorBase> *feasibleset){
	this->feasibleset = feasibleset;
}

template<class VectorBase>
GeneralFeasibleSet<VectorBase> *QPData<VectorBase>::get_feasibleset() const{
	return this->feasibleset;
}


} /* end namespace */

#endif
