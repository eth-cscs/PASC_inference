/** @file tsdata_global.cu
 *  @brief this is only for PETSC!
 *
 *  @author Lukas Pospisil
 */

#ifndef PASC_TSDATA_H
#define	PASC_TSDATA_H

#include <iostream>
#include "general/data/generaldata.h"
#include "general/common/common.h"

namespace pascinference {

/* Maybe these classes are not defined yet */
namespace model {
	template<class VectorBase>
	class TSModel;
}

using namespace model;

namespace data {

/** \class TSData
 *  \brief General Time-series data.
 *
*/
template<class VectorBase>
class TSData: public GeneralData {
	protected:
		TSModel<VectorBase> *tsmodel; /**< pointer to used time-series model on the data */

		GeneralVector<VectorBase> *datavector; /**< global vector with data of dimension based on model */
		bool destroy_datavector; /**< destroy datavector in destructor? if I am an owner, then TRUE */

		GeneralVector<VectorBase> *gammavector; /**< the characteristic functions of clustered models */
		bool destroy_gammavector;

		GeneralVector<VectorBase> *thetavector; /**< parameters of models */
		bool destroy_thetavector;

		double aic_solution; /**< AIC value in solution */

		Decomposition<VectorBase> *decomposition;

		/* scaling variables */
		double scale_max;
		double scale_min;

	public:
		TSData(Decomposition<VectorBase> &decomposition, GeneralVector<VectorBase> *datavector_new, GeneralVector<VectorBase> *gammavector_new, GeneralVector<VectorBase> *thetavector_new);
		TSData(Decomposition<VectorBase> &decomposition);
		TSData(Decomposition<VectorBase> &decomposition, std::string filename);
		TSData();

		~TSData();

		virtual void print(ConsoleOutput &output) const;
		virtual void print(ConsoleOutput &output_global, ConsoleOutput &output_local) const;

		virtual void printcontent(ConsoleOutput &output) const;
		virtual void printcontent(ConsoleOutput &output_global, ConsoleOutput &output_local) const;
		virtual std::string get_name() const;

		void cutgamma() const;

		/* SET functions */
		void set_model(TSModel<VectorBase> &tsmodel);
		void set_aic(double new_aic);

		/* GET functions */
		int get_T() const;
		int get_Tlocal() const;
		int get_Tbegin() const;
		int get_Tend() const;

		int get_R() const;
		int get_Rlocal() const;
		int get_Rbegin() const;
		int get_Rend() const;

		int get_xdim() const;
		int get_K() const;

		double get_aic() const;

		TSModel<VectorBase> *get_model() const;
		Decomposition<VectorBase> *get_decomposition() const;

		GeneralVector<VectorBase> *get_datavector() const;
		GeneralVector<VectorBase> *get_gammavector() const;
		GeneralVector<VectorBase> *get_thetavector() const;

		virtual void save_datavector(std::string filename, int type = 1) const;
		virtual void save_gammavector(std::string filename) const;
		virtual void save_reconstructed(std::string filename, int type = 1) const;

		void save_thetavector(std::string filename) const;
		void print_thetavector(ConsoleOutput &output) const;
		std::string print_thetavector() const;

		void save_gammavector(std::string filename, int blocksize, int type = 1) const; // impicit type gTbR

		/** @brief print basic statistics about data
		*
		* @param output where to print
		* @param printdetails print details of the data or not
		*/
		void printstats(ConsoleOutput &output, bool printdetails=false) const;

		/** @brief scale data to a,b
		*
		*/
		void scaledata(double a, double b);

		void scaledata(double a, double b, double scale_min, double scale_max);

		/** @brief scale data back to original interval
		*
		*/
		void unscaledata(double a, double b);

		/** @brief cut data to a,b
		*
		* @param a lower threshold
		* @param b upper threshold
		*/
		void cutdata(double a, double b);

		/** @brief shift data by given coeficient
		*
		* @param a shifting value
		*/
		void shiftdata(double a);

		virtual void load_gammavector(std::string filename, int type = 1) const; // impicit type gTbR
		virtual void load_gammavector(VectorBase &gamma0) const;
        virtual double compute_abserr_reconstructed(GeneralVector<VectorBase> &solution) const;

		/** @brief compute nbins of actual gammavector
		 */
		double compute_gammavector_nbins();

        double compute_SNR(double signal_value = -1.0) const;
};


}
} /* end of namespace */

/* ------------- implementation ----------- */
#include "general/model/tsmodel.h"

namespace pascinference {
namespace data {

template<class VectorBase>
TSData<VectorBase>::TSData(){
	LOG_FUNC_BEGIN

	this->decomposition = NULL;
	this->tsmodel = NULL;

	this->datavector = NULL;
	destroy_datavector = false;

	this->gammavector = NULL;
	destroy_gammavector = false;

	this->thetavector = NULL;
	destroy_thetavector = false;

	/* set initial aic */
	this->aic_solution = std::numeric_limits<double>::max();

	LOG_FUNC_END
}

template<class VectorBase>
TSData<VectorBase>::TSData(Decomposition<VectorBase> &new_decomposition, GeneralVector<VectorBase> *datavector_new, GeneralVector<VectorBase> *gammavector_new, GeneralVector<VectorBase> *thetavector_new){
	LOG_FUNC_BEGIN

	this->decomposition = &new_decomposition;
	this->tsmodel = NULL;

	if(datavector_new){
		this->datavector = datavector_new;
	} else {
		this->datavector = NULL;
	}
	destroy_datavector = false;

	if(gammavector_new){
		this->gammavector = gammavector_new;
	} else {
		this->gammavector = NULL;
	}
	destroy_gammavector = false;

	if(thetavector_new){
		this->thetavector = thetavector_new;
	} else {
		this->thetavector = NULL;
	}
	destroy_gammavector = false;

	/* set initial aic */
	this->aic_solution = std::numeric_limits<double>::max();

	LOG_FUNC_END
}


/* no datavector provided - prepare own data vector */
template<class VectorBase>
TSData<VectorBase>::TSData(Decomposition<VectorBase> &new_decomposition){
	LOG_FUNC_BEGIN

	//TODO

	LOG_FUNC_END
}

template<class VectorBase>
TSData<VectorBase>::TSData(Decomposition<VectorBase> &new_decomposition, std::string filename){
	LOG_FUNC_BEGIN

	//TODO: check if file exists
	this->decomposition = &new_decomposition;

	/* gamma and theta vectors are not given */
	this->gammavector = NULL;
	destroy_gammavector = false;

	this->thetavector = NULL;
	destroy_thetavector = false;

	this->tsmodel = NULL;

	LOG_FUNC_END
}

template<class VectorBase>
void TSData<VectorBase>::set_model(TSModel<VectorBase> &tsmodel){
	LOG_FUNC_BEGIN

	/* set initial content */
	this->tsmodel = &tsmodel;

	//TODO

	LOG_FUNC_END
}


/* destructor */
template<class VectorBase>
TSData<VectorBase>::~TSData(){
	LOG_FUNC_BEGIN

	/* if I created a datavector, then I should also be able to destroy it */
	if(this->destroy_datavector){
		free(this->datavector);
	}

	if(this->destroy_gammavector){
		free(this->gammavector);
	}

	if(this->destroy_thetavector){
		free(this->thetavector);
	}

	LOG_FUNC_END
}


/* print info about data */
template<class VectorBase>
void TSData<VectorBase>::print(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output <<  this->get_name() << std::endl;

	/* give information about presence of the data */
	if(this->tsmodel){
		output <<  " - T:           " << get_T() << std::endl;
		output <<  " - xdim:        " << get_xdim() << std::endl;
		output <<  " - K:           " << get_K() << std::endl;
		output <<  " - model:       " << tsmodel->get_name() << std::endl;
	} else {
		output <<  " - model:       NO" << std::endl;
	}

	output <<  " - datavector:  ";
	if(this->datavector){
		output << "YES (size: " << this->datavector->size() << ")" << std::endl;
	} else {
		output << "NO" << std::endl;
	}
	output <<   " - gammavector: ";
	if(this->gammavector){
		output << "YES (size: " << this->gammavector->size() << ")" << std::endl;
	} else {
		output << "NO" << std::endl;
	}
	output <<   " - thetavector: ";
	if(this->thetavector){
		output << "YES (size: " << this->thetavector->size() << ")" << std::endl;
	} else {
		output << "NO" << std::endl;
	}

	output.synchronize();

	LOG_FUNC_END
}

/* print info about data */
template<class VectorBase>
void TSData<VectorBase>::print(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN

	output_global <<  this->get_name() << std::endl;

	/* give information about presence of the data */
	if(this->tsmodel){
		output_global <<  " - T:           " << get_T() << std::endl;
		output_local  <<  "  - Tlocal:     " << get_Tlocal() << std::endl;
		output_local.synchronize();

		output_global <<  " - xdim:        " << get_xdim() << std::endl;
		output_global <<  " - K:           " << get_K() << std::endl;

		output_global <<  " - model:       " << get_name() << std::endl;
	} else {
		output_global <<  " - model:       NO" << std::endl;
	}

	output_global <<  " - datavector:  ";
	if(this->datavector){
		output_global << "YES (size: " << this->datavector->size() << ")" << std::endl;
		output_local  <<  "  - local size: " << this->datavector->local_size() << std::endl;
		output_local.synchronize();
	} else {
		output_global << "NO" << std::endl;
	}

	output_global <<   " - gammavector: ";
	if(this->gammavector){
		output_global << "YES (size: " << this->gammavector->size() << ")" << std::endl;
		output_local  <<  "  - local size: " << this->gammavector->local_size() << std::endl;
		output_local.synchronize();
	} else {
		output_global << "NO" << std::endl;
	}

	output_global <<   " - thetavector: ";
	if(this->thetavector){
		output_global << "YES (size: " << this->thetavector->size() << ")" << std::endl;
		output_local  <<  "  - local size: " << this->thetavector->local_size() << std::endl;
		output_local.synchronize();
	} else {
		output_global << "NO" << std::endl;
	}

	output_global.synchronize();

	LOG_FUNC_END
}

/* print content of all data */
template<class VectorBase>
void TSData<VectorBase>::printcontent(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output <<  this->get_name() << std::endl;

	/* print the content of the data */
	output <<  " - datavector: ";
	if(this->datavector){
		output << *this->datavector << std::endl;
	} else {
		output << "not set" << std::endl;
	}

	output <<  " - gammavector: ";
	if(this->gammavector){
		output << *this->gammavector << std::endl;
	} else {
		output << "not set" << std::endl;
	}

	output <<  " - thetavector: ";
	if(this->thetavector){
		output << *this->thetavector << std::endl;
	} else {
		output << "not set" << std::endl;
	}

	LOG_FUNC_END
}

/* print content of all data */
template<class VectorBase>
void TSData<VectorBase>::printcontent(ConsoleOutput &output_global,ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN

	output_global <<  this->get_name() << std::endl;

	/* print the content of the data */
	output_local <<  " - datavector: ";
	if(this->datavector){
		output_local << *this->datavector << std::endl;
	} else {
		output_local << "not set" << std::endl;
	}
	output_local.synchronize();

	output_local <<  " - gammavector: ";
	if(this->gammavector){
		output_local << *this->gammavector << std::endl;
	} else {
		output_local << "not set" << std::endl;
	}
	output_local.synchronize();

	output_local <<  " - thetavector: ";
	if(this->thetavector){
		output_local << *this->thetavector << std::endl;
	} else {
		output_local << "not set" << std::endl;
	}
	output_local.synchronize();

	output_global.synchronize();

	LOG_FUNC_END
}

template<class VectorBase>
std::string TSData<VectorBase>::get_name() const {
	return "Time-series Data";
}

/* ---------- GET functions --------- */
template<class VectorBase>
int TSData<VectorBase>::get_T() const{
	return decomposition->get_T();
}

template<class VectorBase>
int TSData<VectorBase>::get_Tlocal() const{
	return decomposition->get_Tlocal();
}

template<class VectorBase>
int TSData<VectorBase>::get_Tbegin() const{
	return decomposition->get_Tbegin();
}

template<class VectorBase>
int TSData<VectorBase>::get_Tend() const{
	return decomposition->get_Tend();
}

template<class VectorBase>
int TSData<VectorBase>::get_R() const{
	return decomposition->get_R();
}

template<class VectorBase>
int TSData<VectorBase>::get_Rlocal() const{
	return decomposition->get_Rlocal();
}

template<class VectorBase>
int TSData<VectorBase>::get_Rbegin() const{
	return decomposition->get_Rbegin();
}

template<class VectorBase>
int TSData<VectorBase>::get_Rend() const{
	return decomposition->get_Rend();
}

template<class VectorBase>
int TSData<VectorBase>::get_xdim() const{
	return decomposition->get_xdim();
}

template<class VectorBase>
int TSData<VectorBase>::get_K() const{
	return decomposition->get_K();
}

template<class VectorBase>
TSModel<VectorBase> *TSData<VectorBase>::get_model() const{
	return this->tsmodel;
}

template<class VectorBase>
Decomposition<VectorBase> *TSData<VectorBase>::get_decomposition() const {
	return this->decomposition;
}

template<class VectorBase>
GeneralVector<VectorBase> *TSData<VectorBase>::get_datavector() const{
	return this->datavector;
}

template<class VectorBase>
GeneralVector<VectorBase> *TSData<VectorBase>::get_gammavector() const{
	return this->gammavector;
}

template<class VectorBase>
GeneralVector<VectorBase> *TSData<VectorBase>::get_thetavector() const{
	return this->thetavector;
}

template<class VectorBase>
double TSData<VectorBase>::get_aic() const{
	return this->aic_solution;
}

template<class VectorBase>
void TSData<VectorBase>::set_aic(double new_aic) {
	this->aic_solution = new_aic;
}

template<class VectorBase>
void TSData<VectorBase>::cutgamma() const{
	LOG_FUNC_BEGIN

	//TODO

	LOG_FUNC_END
}

template<class VectorBase>
void TSData<VectorBase>::save_thetavector(std::string filename) const {
	LOG_FUNC_BEGIN

	//TODO

	LOG_FUNC_END
}

template<class VectorBase>
void TSData<VectorBase>::print_thetavector(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	//TODO

	LOG_FUNC_END
}

template<class VectorBase>
std::string TSData<VectorBase>::print_thetavector() const {
	std::ostringstream out;

	//TODO

	return out.str();
}

template<class VectorBase>
void TSData<VectorBase>::save_datavector(std::string filename, int type) const {
	LOG_FUNC_BEGIN

	//TODO

	LOG_FUNC_END
}

template<class VectorBase>
void TSData<VectorBase>::save_gammavector(std::string filename) const {
	LOG_FUNC_BEGIN

	//TODO

	LOG_FUNC_END
}

template<class VectorBase>
void TSData<VectorBase>::save_reconstructed(std::string filename, int type) const {
	LOG_FUNC_BEGIN

	//TODO

	LOG_FUNC_END
}

template<class VectorBase>
void TSData<VectorBase>::printstats(ConsoleOutput &output, bool printdetails) const {
	LOG_FUNC_BEGIN

	//TODO

	LOG_FUNC_END
}

template<class VectorBase>
void TSData<VectorBase>::scaledata(double a, double b){
	LOG_FUNC_BEGIN

	//TODO

	LOG_FUNC_END
}

template<class VectorBase>
void TSData<VectorBase>::unscaledata(double a, double b){
	LOG_FUNC_BEGIN

	//TODO

	LOG_FUNC_END
}

template<class VectorBase>
void TSData<VectorBase>::cutdata(double a, double b){
	LOG_FUNC_BEGIN

	//TODO

	LOG_FUNC_END
}

template<class VectorBase>
void TSData<VectorBase>::shiftdata(double a){
	LOG_FUNC_BEGIN

	//TODO

	LOG_FUNC_END
}

template<class VectorBase>
void TSData<VectorBase>::scaledata(double a, double b, double scale_min, double scale_max){
	LOG_FUNC_BEGIN

	//TODO

	LOG_FUNC_END
}

template<class VectorBase>
void TSData<VectorBase>::load_gammavector(VectorBase &gamma0) const {
	LOG_FUNC_BEGIN

	//TODO

	LOG_FUNC_END
}

template<class VectorBase>
void TSData<VectorBase>::load_gammavector(std::string filename, int type) const {
	LOG_FUNC_BEGIN

	//TODO

	LOG_FUNC_END
}

template<class VectorBase>
double TSData<VectorBase>::compute_gammavector_nbins() {
	LOG_FUNC_BEGIN

	//TODO

	LOG_FUNC_END

	return 0.0;
}

template<class VectorBase>
double TSData<VectorBase>::compute_abserr_reconstructed(GeneralVector<VectorBase> &solution) const {
	LOG_FUNC_BEGIN

	//TODO

	LOG_FUNC_END

	return -1.0;
}

template<class VectorBase>
double TSData<VectorBase>::compute_SNR(double signal_value) const{
	LOG_FUNC_BEGIN

    double theta_value;
    if(signal_value < 0.0){
        double max_theta = max(*thetavector);
        double min_theta = min(*thetavector);
        theta_value = max_theta - min_theta;

        coutMaster << "min theta: " << min_theta << std::endl;
        coutMaster << "max theta: " << max_theta << std::endl;

    } else {
        theta_value = signal_value;
    }

    double max_data = max(*datavector);
    double min_data = min(*datavector);

    coutMaster << "min data: " << min_data << std::endl;
    coutMaster << "max data: " << max_data << std::endl;

    double return_value = theta_value/(max_data - min_data);

	LOG_FUNC_END

    return return_value;
}


}
} /* end namespace */

#endif
