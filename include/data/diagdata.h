/** @file diagdata.h
 *  @brief class for manipulation with data for systems of linear equations with diagonal matrix
 *
 *  @author Lukas Pospisil
 */
 
#ifndef PASC_DIAGDATA_H
#define	PASC_DIAGDATA_H

/* for debugging, if >= 100, then print info about ach called function */
extern int DEBUG_MODE;

#include <iostream>
#include "generaldata.h"
#include "algebra.h"

namespace pascinference {

/** class DiagData
 * @brief data of diag(a) x = b 
 * 
 * Class for manipulation with data of systems of linear equations with diagonal matrix.
 *  \f[
 *  \left[
 *   \begin{array}{ccc}
 *    a_1 & & \\
 *    & \ddots & \\
 *    & & a_n
 *   \end{array}
 *  \right]
 *  \left[
 *   \begin{array}{c}
 *    x_1 \\
 *    \vdots \\
 *    x_n
 *   \end{array}
 *  \right] =
 *  \left[
 *   \begin{array}{c}
 *    b_1 \\
 *    \vdots \\
 *    b_n
 *   \end{array}
 *  \right]
 *  \f]
 * 
 */ 
template<class VectorBase>
class DiagData: public GeneralData {
	private:
		/* variables */
		GeneralVector<VectorBase> *a; /**< diagonal of the matrix */
		GeneralVector<VectorBase> *b; /**< right hand-side vector */
		GeneralVector<VectorBase> *x; /**< solution vector */

	public:
	
		/** @brief default constructor
		 */ 
		DiagData();
		
		/** @brief default destructor
		 */ 
		~DiagData();

		/** @brief print information about data
		 * 
		 * @param output where to print
		 */ 
		void print(ConsoleOutput &output) const;

		/** @brief print content of included data
		 * 
		 * @param output where to print
		 */ 
		void printcontent(ConsoleOutput &output) const;

		/** @brief get type of this data
		 */
		std::string get_name() const;

		/** @brief set diagonal of matrix
		 */ 
		void set_a(GeneralVector<VectorBase> *a);

		/** @brief get diagonal of matrix
		 */ 
		GeneralVector<VectorBase> *get_a() const;

		/** @brief set right hand-side vector
		 */ 
		void set_b(GeneralVector<VectorBase> *b);

		/** @brief get right hand-side vector
		 */ 
		GeneralVector<VectorBase> *get_b() const;

		/** @brief set solution vector
		 */ 
		void set_x(GeneralVector<VectorBase> *x);

		/** @brief get solution vector
		 */ 
		GeneralVector<VectorBase> *get_x() const;

};

} // end of namespace

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {

/* constructor */
template<class VectorBase>
DiagData<VectorBase>::DiagData(){
	LOG_FUNC_BEGIN

	/* set initial content */
	this->a = NULL;
	this->b = NULL;
	this->x = NULL;

	LOG_FUNC_END
}

/* destructor */
template<class VectorBase>
DiagData<VectorBase>::~DiagData(){
	LOG_FUNC_BEGIN
	
	LOG_FUNC_END
}


/* print info about problem */
template<class VectorBase>
void DiagData<VectorBase>::print(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output <<  this->get_name() << std::endl;
	
	/* give information about presence of the data */
	output <<  " - a:            ";
	if(this->a){
		output << "YES (size: " << this->a->size() << ")" << std::endl;
	} else {
		output << "not set" << std::endl;
	}
	output <<  " - b:            ";
	if(this->b){
		output << "YES (size: " << this->b->size() << ")" << std::endl;
	} else {
		output << "not set" << std::endl;
	}
	output <<  " - x:            ";
	if(this->x){
		output << "YES (size: " << this->x->size() << ")" << std::endl;
	} else {
		output << "not set" << std::endl;
	}
		
	LOG_FUNC_END
}

template<class VectorBase>
void DiagData<VectorBase>::printcontent(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output <<  this->get_name() << std::endl;
	
	/* give information about presence of the data */
	output <<  " - a:            ";
	if(this->a){
		output << *(this->a) << std::endl;
	} else {
		output << "not set" << std::endl;
	}
	output <<  " - b:            ";
	if(this->b){
		output << *(this->b) << std::endl;
	} else {
		output << "not set" << std::endl;
	}
	output <<  " - x:            ";
	if(this->x){
		output << *(this->x) << std::endl;
	} else {
		output << "not set" << std::endl;
	}

	LOG_FUNC_END
}

template<class VectorBase>
std::string DiagData<VectorBase>::get_name() const {
	return "DiagData";
}

/* ----- SET and GET functions --- */
template<class VectorBase>
void DiagData<VectorBase>::set_a(GeneralVector<VectorBase> *a){
	this->a = a;
}

template<class VectorBase>
GeneralVector<VectorBase> *DiagData<VectorBase>::get_a() const{
	return this->a;
}

template<class VectorBase>
void DiagData<VectorBase>::set_b(GeneralVector<VectorBase> *b){
	this->b = b;
}

template<class VectorBase>
GeneralVector<VectorBase> *DiagData<VectorBase>::get_b() const{
	return this->b;
}

template<class VectorBase>
void DiagData<VectorBase>::set_x(GeneralVector<VectorBase> *x){
	this->x = x;
}

template<class VectorBase>
GeneralVector<VectorBase> *DiagData<VectorBase>::get_x() const{
	return this->x;
}


} /* end namespace */

#endif
