#ifndef PASC_DIAGDATA_H
#define	PASC_DIAGDATA_H

/* for debugging, if >= 100, then print info about ach called function */
extern int DEBUG_MODE;

#include <iostream>
#include "generaldata.h"
#include "algebra.h"

namespace pascinference {

/** class DiagData
 * @brief solve a.*x=b 
 * 
 * solve the system of linear equations with diagonal matrix
 * 
 */ 
template<class VectorBase>
class DiagData: public GeneralData {
	private:
		/* variables */
		GeneralVector<VectorBase> *a;
		GeneralVector<VectorBase> *b;
		GeneralVector<VectorBase> *x; /* solution */

	public:
		DiagData();
		~DiagData();

		void print(std::ostream &output) const;
		std::string get_name() const;

		/* set and get functions */
		void set_a(GeneralVector<VectorBase> *a);
		GeneralVector<VectorBase> *get_a() const;

		void set_b(GeneralVector<VectorBase> *b);
		GeneralVector<VectorBase> *get_b() const;

		void set_x(GeneralVector<VectorBase> *x);
		GeneralVector<VectorBase> *get_x() const;

};

} // end of namespace

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {

/* constructor */
template<class VectorBase>
DiagData<VectorBase>::DiagData(){
	if(DEBUG_MODE >= 100) std::cout << "(DiagData)CONSTRUCTOR" << std::endl;

	/* set initial content */
	this->a = NULL;
	this->b = NULL;
	this->x = NULL;

}

/* destructor */
template<class VectorBase>
DiagData<VectorBase>::~DiagData(){
	if(DEBUG_MODE >= 100) std::cout << "(DiagData)DESTRUCTOR" << std::endl;
	
}


/* print info about problem */
template<class VectorBase>
void DiagData<VectorBase>::print(std::ostream &output) const {
	output << this->get_name() << std::endl;
	
	/* give information about presence of the data */
	output << "  - a:            ";
	if(this->a){
		output << "YES (size: " << this->a->size() << ")" << std::endl;
	} else {
		output << "NO" << std::endl;
	}
	output << "  - b:            ";
	if(this->b){
		output << "YES (size: " << this->b->size() << ")" << std::endl;
	} else {
		output << "NO" << std::endl;
	}
	output << "  - x:            ";
	if(this->x){
		output << "YES (size: " << this->x->size() << ")" << std::endl;
	} else {
		output << "NO" << std::endl;
	}
		
}

template<class VectorBase>
std::string DiagData<VectorBase>::get_name() const {
	return "Diag Data";
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
