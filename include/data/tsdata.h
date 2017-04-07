/** @file tsdata_global.cu
 *  @brief this is only for PETSC!
 * 
 *  @author Lukas Pospisil
 */

#ifndef PASC_TSDATA_H
#define	PASC_TSDATA_H

#ifndef USE_PETSC
 #error 'TSDATA is for PETSC'
#endif

typedef petscvector::PetscVector PetscVector;

#include <iostream>
#include "common/common.h"
#include "model/tsmodel.h"

namespace pascinference {

/* Maybe these classes are not defined yet */ 
namespace model {
	template<class VectorBase>
	class TSModel;
}

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

		Decomposition *decomposition;

		/* scaling variables */
		double scale_max;
		double scale_min;
		
	public:
		TSData(Decomposition &decomposition, GeneralVector<VectorBase> *datavector_new, GeneralVector<VectorBase> *gammavector_new, GeneralVector<VectorBase> *thetavector_new);
		TSData(Decomposition &decomposition);
		TSData(Decomposition &decomposition, std::string filename);
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
		Decomposition *get_decomposition() const;

		GeneralVector<VectorBase> *get_datavector() const;
		GeneralVector<VectorBase> *get_gammavector() const;
		GeneralVector<VectorBase> *get_thetavector() const;

		void save_datavector(std::string filename) const;

		void save_thetavector(std::string filename) const;
		void print_thetavector(ConsoleOutput &output) const;
		std::string print_thetavector() const;

		void save_gammavector(std::string filename, int blocksize) const;

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

		virtual void load_gammavector(std::string filename) const;
		virtual void load_gammavector(VectorBase &gamma0) const;
		
};


}
} /* end of namespace */

/* ------------- implementation ----------- */
//TODO: move to impls

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

template<>
TSData<PetscVector>::TSData(Decomposition &new_decomposition, GeneralVector<PetscVector> *datavector_new, GeneralVector<PetscVector> *gammavector_new, GeneralVector<PetscVector> *thetavector_new){
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
template<>
TSData<PetscVector>::TSData(Decomposition &new_decomposition){
	LOG_FUNC_BEGIN

	this->decomposition = &new_decomposition;
	
	/* we are ready to prepare datavector */
	Vec data_Vec;
	this->decomposition->createGlobalVec_data(&data_Vec);
	this->datavector = new GeneralVector<PetscVector>(data_Vec);
	destroy_datavector = true;

	/* gamma and theta vectors are not given */
	this->gammavector = NULL;
	destroy_gammavector = false;

	this->thetavector = NULL;
	destroy_thetavector = false;

	/* we don't know anything about model */
	this->tsmodel = NULL;

	/* set initial aic */
	this->aic_solution = std::numeric_limits<double>::max();

	LOG_FUNC_END
}

template<>
TSData<PetscVector>::TSData(Decomposition &new_decomposition, std::string filename){
	LOG_FUNC_BEGIN

	//TODO: check if file exists
	this->decomposition = &new_decomposition;
	
	//TODO: implement loader of vector into decomposition?
	///* new data vector only for loading data */
	//Vec dataPreLoad_Vec;
	//TRYCXX( VecCreate(PETSC_COMM_WORLD,&dataPreLoad_Vec) );

	///* prepare viewer to load from file */
	//PetscViewer mviewer;
	//TRYCXX( PetscViewerCreate(PETSC_COMM_WORLD, &mviewer) );
	//TRYCXX( PetscViewerBinaryOpen(PETSC_COMM_WORLD ,filename.c_str(), FILE_MODE_READ, &mviewer) );
	
	///* load vector from viewer */
	//TRYCXX( VecLoad(dataPreLoad_Vec, mviewer) );

	///* destroy the viewer */
	//TRYCXX( PetscViewerDestroy(&mviewer) );

	///* get T */
	//int vec_size;
	//TRYCXX( VecGetSize(dataPreLoad_Vec,&vec_size) );
	//int T = vec_size/(double)blocksize;

	///* now we know the length of vector, we will load it again on right layout */
	//TRYCXX( VecDestroy(&dataPreLoad_Vec) );

	///* now prepare new layout */
	//Vec layout;
	//TRYCXX( VecCreate(PETSC_COMM_WORLD,&layout) );
	//TRYCXX( VecSetSizes(layout, PETSC_DECIDE, T ));
	//TRYCXX( VecSetFromOptions(layout) );

	//int Tbegin, Tend;
	//TRYCXX( VecGetOwnershipRange(layout,&Tbegin,&Tend) );
	
	///* get Tlocal */
	//int Tlocal;
	//TRYCXX( VecGetLocalSize(layout,&Tlocal) );
	
	///* destroy layout vector - now we know everything what is necessary */
	//TRYCXX( VecDestroy(&layout) );
	
	///* we are ready to prepare real datavector */
	//Vec data_Vec;
	//TRYCXX( VecCreate(PETSC_COMM_WORLD,&data_Vec) );
	//TRYCXX( VecSetSizes(data_Vec, Tlocal*blocksize, T*blocksize ) );
	//TRYCXX( VecSetFromOptions(data_Vec) );	

	//TRYCXX( PetscViewerCreate(PETSC_COMM_WORLD, &mviewer) );
	//TRYCXX( PetscViewerBinaryOpen(PETSC_COMM_WORLD ,filename.c_str(), FILE_MODE_READ, &mviewer) );

	///* load data to vector with right layout */
	//TRYCXX( VecLoad(data_Vec, mviewer) );
	
	//TRYCXX( VecAssemblyBegin(data_Vec) );
	//TRYCXX( VecAssemblyEnd(data_Vec) );

	///* destroy the viewer */
	//TRYCXX( PetscViewerDestroy(&mviewer) );

	///* prepare general vector */
	//this->datavector = new GeneralVector<PetscVector>(data_Vec);
	//destroy_datavector = true;

	/* gamma and theta vectors are not given */
	this->gammavector = NULL;
	destroy_gammavector = false;

	this->thetavector = NULL;
	destroy_thetavector = false;

	this->tsmodel = NULL;
	
	LOG_FUNC_END
}

template<>
void TSData<PetscVector>::set_model(TSModel<PetscVector> &tsmodel){
	LOG_FUNC_BEGIN

	/* set initial content */
	this->tsmodel = &tsmodel;

	/* prepare new vectors based on model */
	if(!this->datavector){
		Vec datavector_Vec;
		decomposition->createGlobalVec_data(&datavector_Vec);
		this->datavector = new GeneralVector<PetscVector>(datavector_Vec);
		this->destroy_datavector = true;
	}

	if(!this->gammavector){
		Vec gammavector_Vec;

		decomposition->createGlobalVec_gamma(&gammavector_Vec);

		this->gammavector = new GeneralVector<PetscVector>(gammavector_Vec);
		this->destroy_gammavector = true;
	}

	/* Theta vector is sequential */
	if(!this->thetavector){
		Vec thetavector_Vec;

		TRYCXX( VecCreate(PETSC_COMM_SELF,&thetavector_Vec) );
		#ifdef USE_CUDA
			TRYCXX(VecSetType(thetavector_Vec, VECSEQCUDA));
		#else
			TRYCXX(VecSetType(thetavector_Vec, VECSEQ));
		#endif
		TRYCXX( VecSetSizes(thetavector_Vec,this->tsmodel->get_thetavectorlength_local(),PETSC_DECIDE) );
		TRYCXX( VecSetFromOptions(thetavector_Vec) );
		
		this->thetavector = new GeneralVector<PetscVector>(thetavector_Vec);
		this->destroy_thetavector = true;
	}

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
template<>
void TSData<PetscVector>::print(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
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
template<>
void TSData<PetscVector>::printcontent(ConsoleOutput &output_global,ConsoleOutput &output_local) const {
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
Decomposition *TSData<VectorBase>::get_decomposition() const {
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

template<>
void TSData<PetscVector>::cutgamma() const{
	LOG_FUNC_BEGIN

	int max_id;
	double max_value;
	
	int K = get_K();
	int gamma_t = decomposition->get_Tlocal()*decomposition->get_Rlocal();
	
	double *gamma_arr;
	TRYCXX( VecGetArray(gammavector->get_vector(),&gamma_arr) );
	
	int t,k;
	for(t = 0; t < gamma_t; t++){
		/* find max value */
		max_id = 0;
		max_value = gamma_arr[t];
		for(k = 1; k < K; k++){
			if(gamma_arr[t*K + k] > max_value){
				max_id = k;
				max_value = gamma_arr[t*K + k];
			}
		}
		
		/* set new values */
		for(k = 0; k < K; k++){
			if(k == max_id){
				gamma_arr[t*K + k] = 1.0;
			} else {
				gamma_arr[t*K + k] = 0.0;
			}
		}


	}

	TRYCXX( VecRestoreArray(gammavector->get_vector(),&gamma_arr) );
	
	LOG_FUNC_END
}

template<>
void TSData<PetscVector>::save_datavector(std::string filename) const {
	LOG_FUNC_BEGIN

	PetscViewer viewer_out;
	TRYCXX( PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename.c_str(),FILE_MODE_WRITE,&viewer_out) );
	TRYCXX( VecView( datavector->get_vector(), viewer_out) );
	TRYCXX( PetscViewerDestroy(&viewer_out) );

	LOG_FUNC_END
}

template<>
void TSData<PetscVector>::save_thetavector(std::string filename) const {
	LOG_FUNC_BEGIN

	LOG_FUNC_END
}

template<>
void TSData<PetscVector>::print_thetavector(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output << "Theta:" << std::endl;

	int theta_size = this->tsmodel->get_thetavectorlength_local();
	double *theta_arr;
	TRYCXX( VecGetArray(thetavector->get_vector(),&theta_arr));
		
	for(int i=0;i<theta_size;i++){
		output << theta_arr[i];
		if(i < theta_size-1){
			output << ", ";
		}
	}
	output << std::endl;
		
	TRYCXX( VecRestoreArray(thetavector->get_vector(),&theta_arr));


	LOG_FUNC_END
}

template<>
std::string TSData<PetscVector>::print_thetavector() const {
	std::ostringstream out;

	int theta_size = this->tsmodel->get_thetavectorlength_local();
	
	double *theta_arr;
	TRYCXX( VecGetArray(thetavector->get_vector(),&theta_arr));
	for(int i=0;i<theta_size;i++){
		out << theta_arr[i] << ",";
	}	

	TRYCXX( VecRestoreArray(thetavector->get_vector(),&theta_arr));

	return out.str();
}


template<>
void TSData<PetscVector>::save_gammavector(std::string filename, int blocksize) const {
	LOG_FUNC_BEGIN

	LOG_FUNC_END
}

template<>
void TSData<PetscVector>::printstats(ConsoleOutput &output, bool printdetails) const {
	LOG_FUNC_BEGIN

	std::streamsize ss = std::cout.precision();
	output << std::setprecision(17);

	output <<  "STATS: " << this->get_name() << std::endl;
	output.push();
		int x_size = this->datavector->size();
		int blocksize = get_R()*get_K();
		output << " - total length:    " << std::setw(25) << x_size << std::endl;
		output << " - nmb of blocks:   " << std::setw(25) << blocksize << std::endl;
		output << " - length of block: " << std::setw(25) << get_T() << std::endl;
		
		/* compute basic statistics: */
		Vec x_Vec = datavector->get_vector();

		double x_sum;
		double x_max;
		double x_min;
		double x_avg;

		TRYCXX( VecSum(x_Vec, &x_sum) );
		TRYCXX( VecMax(x_Vec, NULL, &x_max) );
		TRYCXX( VecMin(x_Vec, NULL, &x_min) );
		x_avg = x_sum/(double)x_size;

		output <<  " - sum:             " << std::setw(25) << x_sum << std::endl;
		output <<  " - max:             " << std::setw(25) << x_max << std::endl;
		output <<  " - min:             " << std::setw(25) << x_min << std::endl;
		output <<  " - avg:             " << std::setw(25) << x_avg << std::endl;

		/* for each dimension compute basic statistics: */
		if(printdetails){
			Vec xk_Vec;
			IS xk_is;
			int xk_size = get_T();
		
			double xk_sum;
			double xk_max;
			double xk_min;
			double xk_avg;
		
			for(int k=0;k<blocksize;k++){
				output << "x_" << k << std::endl;
				output.push();
					TRYCXX( ISCreateStride(PETSC_COMM_WORLD, xk_size, k, blocksize, &xk_is) );
					TRYCXX( VecGetSubVector(x_Vec, xk_is, &xk_Vec) );

					TRYCXX( VecSum(xk_Vec, &xk_sum) );
					TRYCXX( VecMax(xk_Vec, NULL, &xk_max) );
					TRYCXX( VecMin(xk_Vec, NULL, &xk_min) );
					xk_avg = xk_sum/(double)xk_size;

					output <<  " - length: " << std::setw(25) << xk_size << std::endl;
					output <<  " - sum:    " << std::setw(25) << xk_sum << std::endl;
					output <<  " - max:    " << std::setw(25) << xk_max << std::endl;
					output <<  " - min:    " << std::setw(25) << xk_min << std::endl;
					output <<  " - avg:    " << std::setw(25) << xk_avg << std::endl;

					TRYCXX( VecRestoreSubVector(x_Vec, xk_is, &xk_Vec) );
					TRYCXX( ISDestroy(&xk_is) );
				output.pop();
			}
		}
		
	output.pop();
	output << std::setprecision(ss);
		
	LOG_FUNC_END
}

/*
template<>
void TSData<PetscVector>::cutdata(double threshold_down, double threshold_up){
	LOG_FUNC_BEGIN

	int x_size = this->datavector->local_size();
	double *x_arr;
	TRYCXX( VecGetArray(datavector->get_vector(), &x_arr) );	

	double value;
	for(int i=0; i<x_size; i++){
		if(x_arr[i] < threshold_down || x_arr[i] > threshold_up){
			if(x_arr[i] < threshold_down){
				value = threshold_down;
			}
			if(x_arr[i] > threshold_up){
				value = threshold_up;
			}

			if(i>0){
				x_arr[i] = x_arr[i-1]; 
			} else {
				if(i<x_size-1){
					x_arr[i] = x_arr[i+1]; 
				} else {
					x_arr[i] = value;
				}
			}
		}
	}

	TRYCXX( VecRestoreArray(datavector->get_vector(), &x_arr) );	
	
	LOG_FUNC_END
}
*/

template<>
void TSData<PetscVector>::scaledata(double a, double b){
	LOG_FUNC_BEGIN

	Vec x_Vec = datavector->get_vector();

	/* compute max and min for scaling */
	TRYCXX( VecMax(x_Vec, NULL, &scale_max) );
	TRYCXX( VecMin(x_Vec, NULL, &scale_min) );

	/* linear transformation y=k*x + q; */
	double k = (b-a)/((double)(scale_max-scale_min));
	double q = a-k*scale_min;

	/* compute x = (x-scale_min)/(scale_max-scale_min) */
	TRYCXX( VecScale(x_Vec, k) );
	TRYCXX( VecShift(x_Vec, q) );
	
	TRYCXX( VecAssemblyBegin(x_Vec) );
	TRYCXX( VecAssemblyEnd(x_Vec) );

	/* scale also Theta */
	if(thetavector){
		int theta_size = this->tsmodel->get_thetavectorlength_local();
		double *theta_arr;
		TRYCXX( VecGetArray(thetavector->get_vector(),&theta_arr));
		
		for(int i=0;i<theta_size;i++){
			theta_arr[i] = k*theta_arr[i] + q;
		}
		
		TRYCXX( VecRestoreArray(thetavector->get_vector(),&theta_arr));

		TRYCXX( VecAssemblyBegin(thetavector->get_vector()) );
		TRYCXX( VecAssemblyEnd(thetavector->get_vector()) );
	}
	
	LOG_FUNC_END
}

template<>
void TSData<PetscVector>::unscaledata(double a, double b){
	LOG_FUNC_BEGIN

	Vec x_Vec = datavector->get_vector();

	/* linear transformation y=k*x + q; */
	/* inverse 1/k*(y - q) = x; */
	double k = (b-a)/((double)(scale_max-scale_min));
	double q = a-k*scale_min;

	/* compute x = (x-scale_min)/(scale_max-scale_min) */
	TRYCXX( VecShift(x_Vec, -q) );
	TRYCXX( VecScale(x_Vec, 1.0/k) );

	TRYCXX( VecAssemblyBegin(x_Vec) );
	TRYCXX( VecAssemblyEnd(x_Vec) );

	/* scale also computed Theta */
	if(thetavector){
		int theta_size = this->tsmodel->get_thetavectorlength_local();
		double *theta_arr;
		TRYCXX( VecGetArray(thetavector->get_vector(),&theta_arr));
		
		for(int i=0;i<theta_size;i++){
			theta_arr[i] = (theta_arr[i] - q)/k;
		}
		
		TRYCXX( VecRestoreArray(thetavector->get_vector(),&theta_arr));

		TRYCXX( VecAssemblyBegin(thetavector->get_vector()) );
		TRYCXX( VecAssemblyEnd(thetavector->get_vector()) );
	}
	
	LOG_FUNC_END
}

template<>
void TSData<PetscVector>::cutdata(double a, double b){
	LOG_FUNC_BEGIN

	double *data_arr;

	TRYCXX( VecGetArray(datavector->get_vector(),&data_arr));
	for(int i=0;i<datavector->local_size();i++){
		if(data_arr[i] > b){
			data_arr[i] = b;
		} 
		if(data_arr[i] < a){
			data_arr[i] = a;
		} 
	}
	TRYCXX( VecRestoreArray(datavector->get_vector(),&data_arr));
	
	LOG_FUNC_END
}


template<>
void TSData<PetscVector>::shiftdata(double a){
	LOG_FUNC_BEGIN

	Vec x_Vec = datavector->get_vector();

	TRYCXX( VecShift(x_Vec, a) );
	
	TRYCXX( VecAssemblyBegin(x_Vec) );
	TRYCXX( VecAssemblyEnd(x_Vec) );

	/* scale also computed Theta */
	if(thetavector){
		int theta_size = this->tsmodel->get_thetavectorlength_local();
		double *theta_arr;
		TRYCXX( VecGetArray(thetavector->get_vector(),&theta_arr));
		
		for(int i=0;i<theta_size;i++){
			theta_arr[i] = theta_arr[i] + a;
		}
		
		TRYCXX( VecRestoreArray(thetavector->get_vector(),&theta_arr));
		TRYCXX( VecAssemblyBegin(thetavector->get_vector()) );
		TRYCXX( VecAssemblyEnd(thetavector->get_vector()) );
	}
	
	LOG_FUNC_END
}

template<>
void TSData<PetscVector>::scaledata(double a, double b, double scale_min, double scale_max){
	LOG_FUNC_BEGIN

	Vec x_Vec = datavector->get_vector();
	this->scale_max = scale_max;
	this->scale_min = scale_min;

	/* linear transformation y=k*x + q; */
	double k = (b-a)/((double)(scale_max-scale_min));
	double q = a-k*scale_min;

	/* compute x = (x-scale_min)/(scale_max-scale_min) */
	TRYCXX( VecScale(x_Vec, k) );
	TRYCXX( VecShift(x_Vec, q) );
	
	TRYCXX( VecAssemblyBegin(x_Vec) );
	TRYCXX( VecAssemblyEnd(x_Vec) );

	/* scale also Theta */
	if(thetavector){
		int theta_size = this->tsmodel->get_thetavectorlength_local();
		double *theta_arr;
		TRYCXX( VecGetArray(thetavector->get_vector(),&theta_arr));
		
		for(int i=0;i<theta_size;i++){
			theta_arr[i] = k*theta_arr[i] + q;
		}
		
		TRYCXX( VecRestoreArray(thetavector->get_vector(),&theta_arr));

		TRYCXX( VecAssemblyBegin(thetavector->get_vector()) );
		TRYCXX( VecAssemblyEnd(thetavector->get_vector()) );
	}
	
	LOG_FUNC_END
}

template<>
void TSData<PetscVector>::load_gammavector(PetscVector &gamma0) const {
	LOG_FUNC_BEGIN
	
	/* get petsc Vec from provided vector - this vector is stride */
	Vec gamma0_Vec = gamma0.get_vector();
	
	/* variables */
	int K = this->get_K();
	int T = this->get_T();
	int Tlocal = this->get_Tlocal();
	int Tbegin = this->get_Tbegin();
	int xdim = this->get_xdim();

	/* prepare IS with my indexes in provided vector */
	IS gamma_sub_IS;
	
	/* fill the index sets */
	TRYCXX( ISCreateStride(PETSC_COMM_WORLD, Tlocal*K, Tbegin*K, 1, &(gamma_sub_IS)) );

	/* now get subvector with my local values from provided stride vector */
	Vec gamma_sub;
	TRYCXX( VecGetSubVector(gamma0_Vec, gamma_sub_IS, &gamma_sub) );

	/* prepare local vector */
	Vec gamma_local;
	#ifndef USE_CUDA
		TRYCXX( VecCreateSeq(PETSC_COMM_SELF, K*Tlocal, &gamma_local) );
	#else
		TRYCXX( VecCreateSeqCUDA(PETSC_COMM_SELF, K*Tlocal, &gamma_local) );
	#endif

	/* get the vector where I will store my values */
	TRYCXX( VecGetLocalVector(gammavector->get_vector(), gamma_local) );

	/* now copy values from subvector to local vector */
	TRYCXX( VecCopy(gamma_sub, gamma_local) );

	/* restore subvector */
	TRYCXX( VecRestoreLocalVector(gammavector->get_vector(), gamma_local) );
	TRYCXX( VecRestoreSubVector(gamma0_Vec, gamma_sub_IS, &gamma_sub) );

	/* destroy auxiliary index sets */
	TRYCXX( ISDestroy(&gamma_sub_IS) );

	LOG_FUNC_END
}

template<>
void TSData<PetscVector>::load_gammavector(std::string filename) const {
	LOG_FUNC_BEGIN
	
	//TODO: control existence of file

	/* aux vector, we first oad data and then distribute values to procs */
	Vec gamma_preload_Vec;
	TRYCXX( VecCreate(PETSC_COMM_WORLD, &gamma_preload_Vec) );

	#ifdef USE_CUDA
		TRYCXX(VecSetType(gamma_preload_Vec, VECMPICUDA));
	#endif

	/* prepare viewer to load from file */
	PetscViewer mviewer;
	TRYCXX( PetscViewerCreate(PETSC_COMM_WORLD, &mviewer) );
	TRYCXX( PetscViewerBinaryOpen(PETSC_COMM_WORLD ,filename.c_str(), FILE_MODE_READ, &mviewer) );
	
	/* load vector from viewer */
	TRYCXX( VecLoad(gamma_preload_Vec, mviewer) );

	/* destroy the viewer */
	TRYCXX( PetscViewerDestroy(&mviewer) );	

	PetscVector gamma_preload(gamma_preload_Vec);
	this->load_gammavector(gamma_preload);

//	TRYCXX( VecDestroy(&gamma_preload_Vec) );
	
	LOG_FUNC_END
}



}
} /* end namespace */

#endif
