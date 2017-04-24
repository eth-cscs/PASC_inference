/** @file entropydata.h
 *  @brief class for manipulation with data for entropy problems
 *
 *  @author Lukas Pospisil
 */
 
#ifndef PASC_ENTROPYDATA_H
#define	PASC_ENTROPYDATA_H

#include <iostream>

namespace pascinference {
namespace data {

template<class VectorBase>
class EntropyData: public GeneralData {
	private:
		/* variables */
		GeneralVector<VectorBase> *lambda; /**< solution vector */
		GeneralVector<VectorBase> *x; /**< vector with data */
		GeneralVector<VectorBase> *gamma; /**< cluster indicator functions */

		Decomposition *decomposition;

		int K; /**< number of clusters */
		int Km; /**< number of moments */
		int T; /**< length of time-series */
		
	public:
	
		/** @brief default constructor
		 */ 
		EntropyData(int T, int K, int Km);
		
		/** @brief default destructor
		 */ 
		~EntropyData();

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

		void set_x(GeneralVector<VectorBase> *x);
		GeneralVector<VectorBase> *get_x() const;

		void set_lambda(GeneralVector<VectorBase> *lambda);
		GeneralVector<VectorBase> *get_lambda() const;

		void set_gamma(GeneralVector<VectorBase> *gamma);
		GeneralVector<VectorBase> *get_gamma() const;

		void set_decomposition(Decomposition *decomposition);
		Decomposition *get_decomposition() const;

		int get_T() const;
		int get_K() const;
		int get_Km() const;

};


}
} /* end of namespace */

/* ------------- implementation ----------- */

namespace pascinference {
namespace data {

/* constructor */
template<class VectorBase>
EntropyData<VectorBase>::EntropyData(int T, int K, int Km){
	LOG_FUNC_BEGIN

	/* set initial content */
	this->lambda = NULL;
	this->x = NULL;
	this->gamma = NULL;

	this->T = T;
	this->K = K;
	this->Km = Km;

	LOG_FUNC_END
}

/* destructor */
template<class VectorBase>
EntropyData<VectorBase>::~EntropyData(){
	LOG_FUNC_BEGIN
	
	LOG_FUNC_END
}


/* print info about data */
template<class VectorBase>
void EntropyData<VectorBase>::print(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output <<  this->get_name() << std::endl;

	output <<  " - T             : " << T << std::endl;
	output <<  " - K             : " << K << std::endl;
	output <<  " - Km            : " << Km << std::endl;
	
	/* give information about presence of the data */
	output <<  " - lambda        : ";
	if(this->lambda){
		output << "YES (size: " << this->lambda->size() << ")" << std::endl;
	} else {
		output << "not set" << std::endl;
	}

	output <<  " - x             : ";
	if(this->x){
		output << "YES (size: " << this->x->size() << ")" << std::endl;
	} else {
		output << "not set" << std::endl;
	}

	output <<  " - gamma         : ";
	if(this->gamma){
		output << "YES (size: " << this->gamma->size() << ")" << std::endl;
	} else {
		output << "not set" << std::endl;
	}
	
	LOG_FUNC_END
}

template<class VectorBase>
void EntropyData<VectorBase>::print(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN

	output_global <<  this->get_name() << std::endl;

	output_global <<  " - T             : " << T << std::endl;
	output_global <<  " - K             : " << K << std::endl;
	output_global <<  " - Km            : " << Km << std::endl;
	
	/* give information about presence of the data */
	output_global <<  " - lambda        : ";
	if(this->lambda){
		output_global << "YES (size: " << this->lambda->size() << ")" << std::endl;
		output_global.push();
		output_local  <<  "local size: " << this->lambda->local_size() << std::endl;
		output_global.pop();
		output_local.synchronize();		
	} else {
		output_global << "not set" << std::endl;
	}

	output_global <<  " - x             : ";
	if(this->x){
		output_global << "YES (size: " << this->x->size() << ")" << std::endl;
		output_global.push();
		output_local  <<  "local size: " << this->x->local_size() << std::endl;
		output_global.pop();
		output_local.synchronize();		
	} else {
		output_global << "not set" << std::endl;
	}

	output_global <<  " - gamma         : ";
	if(this->gamma){
		output_global << "YES (size: " << this->gamma->size() << ")" << std::endl;
		output_global.push();
		output_local  <<  "local size: " << this->gamma->local_size() << std::endl;
		output_global.pop();
		output_local.synchronize();		
	} else {
		output_global << "not set" << std::endl;
	}

		
	LOG_FUNC_END
}

template<class VectorBase>
void EntropyData<VectorBase>::printcontent(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output <<  this->get_name() << std::endl;
	
	/* give information about presence of the data */
	output << " - lambda        : ";
	if(this->lambda){
		output << *(this->lambda) << std::endl;
	} else {
		output << "not set" << std::endl;
	}

	output << " - x             : ";
	if(this->x){
		output << *(this->x) << std::endl;
	} else {
		output << "not set" << std::endl;
	}

	output << " - gamma         : ";
	if(this->gamma){
		output << *(this->gamma) << std::endl;
	} else {
		output << "not set" << std::endl;
	}


	LOG_FUNC_END
}

template<class VectorBase>
void EntropyData<VectorBase>::printcontent(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN

	output_global << this->get_name() << std::endl;

	// TODO
	LOG_FUNC_END
}

template<class VectorBase>
std::string EntropyData<VectorBase>::get_name() const {
	return "EntropyData";
}

/* ----- SET and GET functions --- */
template<class VectorBase>
void EntropyData<VectorBase>::set_lambda(GeneralVector<VectorBase> *lambda){
	this->lambda = lambda;
}

template<class VectorBase>
GeneralVector<VectorBase> *EntropyData<VectorBase>::get_lambda() const{
	return this->lambda;
}

template<class VectorBase>
void EntropyData<VectorBase>::set_x(GeneralVector<VectorBase> *x){
	this->x = x;
}

template<class VectorBase>
GeneralVector<VectorBase> *EntropyData<VectorBase>::get_x() const{
	return this->x;
}

template<class VectorBase>
void EntropyData<VectorBase>::set_gamma(GeneralVector<VectorBase> *gamma){
	this->gamma = gamma;
}

template<class VectorBase>
GeneralVector<VectorBase> *EntropyData<VectorBase>::get_gamma() const{
	return this->gamma;
}

template<class VectorBase>
int EntropyData<VectorBase>::get_T() const {
	return this->T;
}

template<class VectorBase>
int EntropyData<VectorBase>::get_K() const {
	return this->K;
}

template<class VectorBase>
int EntropyData<VectorBase>::get_Km() const {
	return this->Km;
}

template<class VectorBase>
void EntropyData<VectorBase>::set_decomposition(Decomposition *decomposition){
	this->decomposition = decomposition;
}

template<class VectorBase>
Decomposition *EntropyData<VectorBase>::get_decomposition() const {
	return this->decomposition;
}


}
} /* end namespace */

#endif
