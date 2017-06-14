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

		Decomposition<VectorBase> *decomposition;

		int K; /**< number of clusters */
		int Km; /**< number of moments */
		int T; /**< length of time-series */
		int xdim; /**< dimension of data */
		
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

		static int compute_number_of_moments(int xdim, int Km);
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
	output <<  " - xdim          : " << xdim << std::endl;
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
	output_global <<  " - xdim          : " << xdim << std::endl;
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
int EntropyData<VectorBase>::get_xdim() const {
	return this->xdim;
}

template<class VectorBase>
void EntropyData<VectorBase>::set_decomposition(Decomposition<VectorBase> *decomposition){
	this->decomposition = decomposition;
}

template<class VectorBase>
Decomposition<VectorBase> *EntropyData<VectorBase>::get_decomposition() const {
	return this->decomposition;
}

template<class VectorBase>
int EntropyData<VectorBase>::compute_number_of_moments(int d, int k){
    int mp_row = pow(k+1,d); /* upper bound on row number */

    //matrix of powers for calculations on joint moments
    matrix<double> mom_powers; /* = D from THE paper */
    mom_powers.set_size(mp_row, d);
    
    //steps for each dimension
    int step;
    int s;
    int ind = 0;
    //compute powers matrix
    for (int i = d-1; i >= 0; i--){
        step = pow(k+1,ind);
        ind++;
        s = 0;
        for (int j = 0; j< mp_row; j = j + step)
        {
            set_subm(mom_powers,range(j,j+step-1), range(i,i)) = p(s);
            if (s == k)
                s = 0;
            else
                s++;
        }
    }
    
    //remove all rows where the sum of the elements is 0 or > k
    for (long j = 0; j< mp_row; j++)
        if (sum(subm(mom_powers, range(j,j), range(0,d-1))) > k || sum(subm(mom_powers, range(j,j), range(0,d-1))) == 0){
            mom_powers = remove_row(mom_powers, j);
            mp_row--;
            j--;
        }

	return mp_row;
}



}
} /* end namespace */

#endif
