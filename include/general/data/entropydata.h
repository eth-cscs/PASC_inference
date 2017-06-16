/** @file entropydata.h
 *  @brief class for manipulation with data for entropy problems
 *
 *  @author Lukas Pospisil
 */
 
#ifndef PASC_ENTROPYDATA_H
#define	PASC_ENTROPYDATA_H

#include <iostream>
#include "general/data/generaldata.h"

namespace pascinference {
namespace data {

template<class VectorBase>
class EntropyData: public GeneralData {
	public:
		class ExternalContent;

	protected:
		friend class ExternalContent;
		ExternalContent *externalcontent;			/**< for manipulation with external-specific stuff */

		/* variables */
		GeneralVector<VectorBase> *lambda; /**< solution vector */
		GeneralVector<VectorBase> *x; /**< vector with data */
		GeneralVector<VectorBase> *gamma; /**< cluster indicator functions */

		Decomposition<VectorBase> *decomposition;

		int K; /**< number of clusters */
		int Km; /**< number of moments */
		int T; /**< length of time-series */
		int xdim; /**< dimension of data */

		GeneralVector<VectorBase> *matrix_D; /**< matrix with power exponents in moments */

		/** @brief fill matrix with power indexes corresponding to moments
		* 
		*/
		void prepare_matrix_D();

		int number_of_moments;		/**< number of moments */
		
	public:
	
		/** @brief default constructor
		 */ 
		EntropyData(int T, int xdim, int K, int Km);
		
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

		void print_matrix_D(ConsoleOutput &output) const;

		/** @brief get type of this data
		 */
		std::string get_name() const;

		void set_x(GeneralVector<VectorBase> *x);
		GeneralVector<VectorBase> *get_x() const;

		void set_lambda(GeneralVector<VectorBase> *lambda);
		GeneralVector<VectorBase> *get_lambda() const;

		void set_gamma(GeneralVector<VectorBase> *gamma);
		GeneralVector<VectorBase> *get_gamma() const;

		void set_decomposition(Decomposition<VectorBase> *decomposition);
		Decomposition<VectorBase> *get_decomposition() const;

		int get_T() const;
		int get_K() const;
		int get_Km() const;
		int get_xdim() const;

		GeneralVector<VectorBase> *get_matrix_D() const;

		int get_number_of_moments() const;
		static int compute_number_of_moments(int xdim, int Km);

		ExternalContent *get_externalcontent() const;
};


}
} /* end of namespace */

/* ------------- implementation ----------- */

namespace pascinference {
namespace data {

/* constructor */
template<class VectorBase>
EntropyData<VectorBase>::EntropyData(int T, int xdim, int K, int Km){
	LOG_FUNC_BEGIN

	/* set initial content */
	this->lambda = NULL;
	this->x = NULL;
	this->gamma = NULL;

	this->T = T;
	this->xdim = xdim;
	this->K = K;
	this->Km = Km;

	this->number_of_moments = compute_number_of_moments(this->xdim, this->Km);

	this->prepare_matrix_D();

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

	output <<  " - T                 : " << T << std::endl;
	output <<  " - xdim              : " << xdim << std::endl;
	output <<  " - K                 : " << K << std::endl;
	output <<  " - Km                : " << Km << std::endl;
	output <<  " - number_of_moments : " << number_of_moments << std::endl;
	
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

	/* print matrix with exponents */
	output <<  " - D             : ";
	output.push();
	this->print_matrix_D(output);
	output.pop();
	
	LOG_FUNC_END
}

template<class VectorBase>
void EntropyData<VectorBase>::print(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN

	output_global <<  this->get_name() << std::endl;

	output_global <<  " - T                 : " << T << std::endl;
	output_global <<  " - xdim              : " << xdim << std::endl;
	output_global <<  " - K                 : " << K << std::endl;
	output_global <<  " - Km                : " << Km << std::endl;
	output_global <<  " - number_of_moments : " << number_of_moments << std::endl;
	
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

	/* print matrix with exponents */
	output_global <<  " - D             : " << std::endl;
	output_global.push();
	this->print_matrix_D(output_global);
	output_global.pop();
		
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

	output <<  " - D            : ";
	this->print_matrix_D(output);

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
int EntropyData<VectorBase>::get_xdim() const {
	return this->xdim;
}

template<class VectorBase>
int EntropyData<VectorBase>::get_number_of_moments() const {
	return this->number_of_moments;
}

template<class VectorBase>
void EntropyData<VectorBase>::set_decomposition(Decomposition<VectorBase> *decomposition){
	this->decomposition = decomposition;
}

template<class VectorBase>
Decomposition<VectorBase> *EntropyData<VectorBase>::get_decomposition() const {
	return this->decomposition;
}

/* number of moments which corresponds to one cluster ! */
template<class VectorBase>
int EntropyData<VectorBase>::compute_number_of_moments(int xdim, int Km){
	/* nmb = 1/xdim! * prod_{i=1}^{xdim} (Km + i) */

    int top = 1;
    int bottom = 1;

	for(int i=1;i<=xdim;i++){
		top *= (Km+i);
		bottom *= i;
	}

	return (int)(top/(double)bottom);
}

template<class VectorBase>
GeneralVector<VectorBase> *EntropyData<VectorBase>::get_matrix_D() const {
	return this->matrix_D;
}

template<class VectorBase>
void EntropyData<VectorBase>::prepare_matrix_D() {
	LOG_FUNC_BEGIN

	//TODO
		
	LOG_FUNC_END
}

template<class VectorBase>
void EntropyData<VectorBase>::print_matrix_D(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	for(int i_moment=0;i_moment < get_number_of_moments(); i_moment++){
		for(int i_xdim=0;i_xdim < get_xdim(); i_xdim++){
			output << (*matrix_D)( i_moment*get_xdim() + i_xdim);
			if(i_xdim < get_xdim()-1) output << ", ";
		}
		output << std::endl;
	}
		
	LOG_FUNC_END
}

/* define blank external content for general VectorBase */
template<class VectorBase>
class EntropyData<VectorBase>::ExternalContent {
};



}
} /* end namespace */

#endif
