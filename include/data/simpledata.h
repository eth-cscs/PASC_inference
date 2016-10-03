/** @file simpledata.h
 *  @brief class for manipulation with data for problems where solution is computed elsewhere
 *
 *  @author Lukas Pospisil
 */
 
#ifndef PASC_SIMPLEDATA_H
#define	PASC_SIMPLEDATA_H

/* for debugging, if >= 100, then print info about ach called function */
extern int DEBUG_MODE;

#include "pascinference.h"

namespace pascinference {
namespace data {

template<class VectorBase>
class SimpleData: public GeneralData {
	private:
		/* variables */
		GeneralVector<VectorBase> *x; /**< solution vector */

	public:
	
		/** @brief default constructor
		 */ 
		SimpleData();
		
		/** @brief default destructor
		 */ 
		~SimpleData();

		/** @brief print information about data
		 * 
		 * @param output where to print
		 */ 
		void print(ConsoleOutput &output) const;

		void print(ConsoleOutput &output_global, ConsoleOutput &output_local) const;

		/** @brief print content of included data
		 * 
		 * @param output where to print
		 */ 
		void printcontent(ConsoleOutput &output) const;

		void printcontent(ConsoleOutput &output_global, ConsoleOutput &output_local) const;

		/** @brief get type of this data
		 */
		std::string get_name() const;

		/** @brief set solution vector
		 */ 
		void set_x(GeneralVector<VectorBase> *x);

		/** @brief get solution vector
		 */ 
		GeneralVector<VectorBase> *get_x() const;

};


}
} /* end of namespace */

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {
namespace data {

/* constructor */
template<class VectorBase>
SimpleData<VectorBase>::SimpleData(){
	LOG_FUNC_BEGIN

	/* set initial content */
	this->x = NULL;

	LOG_FUNC_END
}

/* destructor */
template<class VectorBase>
SimpleData<VectorBase>::~SimpleData(){
	LOG_FUNC_BEGIN
	
	LOG_FUNC_END
}


/* print info about data */
template<class VectorBase>
void SimpleData<VectorBase>::print(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output <<  this->get_name() << "\n";
	
	/* give information about presence of the data */
	output <<  " - x:            ";
	if(this->x){
		output << "YES (size: " << this->x->size() << ")\n";
	} else {
		output << "not set\n";
	}
		
	LOG_FUNC_END
}

template<class VectorBase>
void SimpleData<VectorBase>::print(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN

	output_global <<  this->get_name() << "\n";
	
	/* give information about presence of the data */
	output_global <<  " - x:            ";
	if(this->x){
		output_global << "YES (size: " << this->x->size() << ")\n";
		output_global.push();
		output_local  <<  "local size: " << this->x->local_size() << "\n";
		output_global.pop();
		output_local.synchronize();		
	} else {
		output_global << "not set\n";
	}
		
	LOG_FUNC_END
}

template<class VectorBase>
void SimpleData<VectorBase>::printcontent(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output <<  this->get_name() << "\n";
	
	/* give information about presence of the data */
	output <<  " - x:            ";
	if(this->x){
		output << *(this->x) << "\n";
	} else {
		output << "not set\n";
	}

	LOG_FUNC_END
}

template<class VectorBase>
void SimpleData<VectorBase>::printcontent(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN

	output_global << this->get_name() << "\n";

	// TODO
	LOG_FUNC_END
}

template<class VectorBase>
std::string SimpleData<VectorBase>::get_name() const {
	return "SimpleData";
}

/* ----- SET and GET functions --- */
template<class VectorBase>
void SimpleData<VectorBase>::set_x(GeneralVector<VectorBase> *x){
	this->x = x;
}

template<class VectorBase>
GeneralVector<VectorBase> *SimpleData<VectorBase>::get_x() const{
	return this->x;
}


}
} /* end namespace */

#endif
