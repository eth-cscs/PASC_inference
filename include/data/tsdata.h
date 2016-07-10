/** @file tsdata_global.cu
 *  @brief this is only for PETSC!
 * 
 *  @author Lukas Pospisil
 */

#ifndef PASC_TSDATA_H
#define	PASC_TSDATA_H

#ifndef USE_PETSCVECTOR
 #error 'TSDATA is for PETSCVECTOR'
#endif

typedef petscvector::PetscVector PetscVector;

/* for debugging, if >= 100, then print info about ach called function */
extern int DEBUG_MODE;

#include <iostream>
#include "common/common.h"
#include "model/tsmodel.h"

namespace pascinference {

/* Maybe these classes are not defined yet */ 
template<class VectorBase>
class TSModel;

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

		// TODO: move variables from model here
		int T;
		int Tlocal;
		int Tbegin;
		int Tend;
		int blocksize;

	public:
		TSData(GeneralVector<VectorBase> *datavector_new, GeneralVector<VectorBase> *gammavector_new, GeneralVector<VectorBase> *thetavector_new, int T);
		TSData(int T, int block_size=1);
		TSData(std::string filename , int block_size=1);
		TSData();

		~TSData();

		virtual void print(ConsoleOutput &output) const;
		virtual void print(ConsoleOutput &output_global, ConsoleOutput &output_local) const;

		virtual void printcontent(ConsoleOutput &output) const;
		virtual void printcontent(ConsoleOutput &output_global, ConsoleOutput &output_local) const;
		virtual std::string get_name() const;

		void cut_gamma() const;

		/* SET functions */
		void set_model(TSModel<VectorBase> &tsmodel);
		void set_aic(double new_aic);

		/* GET functions */
		int get_T() const;
		int get_Tlocal() const;
		int get_Tbegin() const;
		int get_Tend() const;
		int get_xdim() const;
		int get_xmem() const;
		int get_K() const;
		double get_aic() const;
		
		TSModel<VectorBase> *get_model() const;
		GeneralVector<VectorBase> *get_datavector() const;
		GeneralVector<VectorBase> *get_gammavector() const;
		GeneralVector<VectorBase> *get_thetavector() const;

		void save_datavector(std::string filename) const;
		void save_thetavector(std::string filename) const;
		void save_gammavector(std::string filename, int blocksize) const;

};

} // end of namespace

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {

template<class VectorBase>
TSData<VectorBase>::TSData(){
	LOG_FUNC_BEGIN

	this->tsmodel = NULL;

	this->datavector = NULL;
	destroy_datavector = false;

	this->gammavector = NULL;
	destroy_gammavector = false;

	this->thetavector = NULL;
	destroy_thetavector = false;

	this->T = 0;
	this->Tlocal = 0;
	this->Tbegin = 0;
	this->Tend = 0;
	this->blocksize = 0;

	/* set initial aic */
	this->aic_solution = std::numeric_limits<double>::max();

	LOG_FUNC_END
}

template<>
TSData<PetscVector>::TSData(GeneralVector<PetscVector> *datavector_new, GeneralVector<PetscVector> *gammavector_new, GeneralVector<PetscVector> *thetavector_new, int T){
	LOG_FUNC_BEGIN

	this->tsmodel = NULL;
	this->T = T;

	if(datavector_new){
		this->datavector = datavector_new;
	} else {
		this->datavector = NULL;
	}
	destroy_datavector = false;

	if(gammavector_new){
		this->gammavector = gammavector_new;
		
		/* compute distribution of gamma from gamma */
		int global_size, local_size, low, high;
		TRY( VecGetSize(this->gammavector->get_vector(), &global_size) );
		TRY( VecGetLocalSize(this->gammavector->get_vector(), &local_size) );
		TRY( VecGetOwnershipRange(this->gammavector->get_vector(),&low,&high) );
	
		this->blocksize = global_size/(double)T;
		this->Tlocal = local_size/(double)this->blocksize;
		this->Tbegin = low/(double)this->blocksize;
		this->Tend = high/(double)this->blocksize;

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
TSData<PetscVector>::TSData(int T, int block_size){
	LOG_FUNC_BEGIN

	/* prepare new layout */
	Vec layout;
	TRY( VecCreate(PETSC_COMM_WORLD,&layout) );
	TRY( VecSetSizes(layout, PETSC_DECIDE, T ));
	TRY( VecSetFromOptions(layout) );

	int Tbegin, Tend;
	TRY( VecGetOwnershipRange(layout,&Tbegin,&Tend) );
	
	/* get Tlocal */
	int Tlocal;
	TRY( VecGetLocalSize(layout,&Tlocal) );
	
	/* destroy layout vector - now we know everything what is necessary */
	TRY( VecDestroy(&layout) );

	/* we are ready to prepare real datavector */
	Vec data_Vec;
	TRY( VecCreate(PETSC_COMM_WORLD,&data_Vec) );
	TRY( VecSetSizes(data_Vec, Tlocal*block_size, T*block_size ) );
	TRY( VecSetFromOptions(data_Vec) );	

	TRY( VecAssemblyBegin(data_Vec) );
	TRY( VecAssemblyEnd(data_Vec) );

	/* prepare general vector */
	this->datavector = new GeneralVector<PetscVector>(data_Vec);
	destroy_datavector = true;

	/* gamma and theta vectors are not given */
	this->gammavector = NULL;
	destroy_gammavector = false;

	this->thetavector = NULL;
	destroy_thetavector = false;

	/* we don't know anything about model */
	this->tsmodel = NULL;

	/* store provided values */
	this->T = T;
	this->Tlocal = Tlocal;
	this->Tbegin = Tbegin;
	this->Tend = Tend;
	this->blocksize = blocksize;

	/* set initial aic */
	this->aic_solution = std::numeric_limits<double>::max();

	LOG_FUNC_END
}

template<>
TSData<PetscVector>::TSData(std::string filename, int block_size){
	LOG_FUNC_BEGIN

	//TODO: check if file exists

	/* new data vector only for loading data */
	Vec dataPreLoad_Vec;
	TRY( VecCreate(PETSC_COMM_WORLD,&dataPreLoad_Vec) );

	/* prepare viewer to load from file */
	PetscViewer mviewer;
	TRY( PetscViewerCreate(PETSC_COMM_WORLD, &mviewer) );
	TRY( PetscViewerBinaryOpen(PETSC_COMM_WORLD ,filename.c_str(), FILE_MODE_READ, &mviewer) );
	
	/* load vector from viewer */
	TRY( VecLoad(dataPreLoad_Vec, mviewer) );

	/* destroy the viewer */
	TRY( PetscViewerDestroy(&mviewer) );

	/* get T */
	int vec_size;
	TRY( VecGetSize(dataPreLoad_Vec,&vec_size) );
	int T = vec_size/(double)block_size;

	/* now we know the length of vector, we will load it again on right layout */
	TRY( VecDestroy(&dataPreLoad_Vec) );

	/* now prepare new layout */
	Vec layout;
	TRY( VecCreate(PETSC_COMM_WORLD,&layout) );
	TRY( VecSetSizes(layout, PETSC_DECIDE, T ));
	TRY( VecSetFromOptions(layout) );

	int Tbegin, Tend;
	TRY( VecGetOwnershipRange(layout,&Tbegin,&Tend) );
	
	/* get Tlocal */
	int Tlocal;
	TRY( VecGetLocalSize(layout,&Tlocal) );
	
	/* destroy layout vector - now we know everything what is necessary */
	TRY( VecDestroy(&layout) );
	
	/* we are ready to prepare real datavector */
	Vec data_Vec;
	TRY( VecCreate(PETSC_COMM_WORLD,&data_Vec) );
	TRY( VecSetSizes(data_Vec, Tlocal*block_size, T*block_size ) );
	TRY( VecSetFromOptions(data_Vec) );	

	TRY( PetscViewerCreate(PETSC_COMM_WORLD, &mviewer) );
	TRY( PetscViewerBinaryOpen(PETSC_COMM_WORLD ,filename.c_str(), FILE_MODE_READ, &mviewer) );

	/* load data to vector with right layout */
	TRY( VecLoad(data_Vec, mviewer) );
	
	TRY( VecAssemblyBegin(data_Vec) );
	TRY( VecAssemblyEnd(data_Vec) );

	/* destroy the viewer */
	TRY( PetscViewerDestroy(&mviewer) );

	/* prepare general vector */
	this->datavector = new GeneralVector<PetscVector>(data_Vec);
	destroy_datavector = true;

	/* gamma and theta vectors are not given */
	this->gammavector = NULL;
	destroy_gammavector = false;

	this->thetavector = NULL;
	destroy_thetavector = false;

	this->tsmodel = NULL;

	/* store provided values */
	this->T = T;
	this->Tlocal = Tlocal;
	this->Tbegin = Tbegin;
	this->Tend = Tend;
	this->blocksize = blocksize;
	
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
		TRY( VecCreate(PETSC_COMM_WORLD,&datavector_Vec) );
		TRY( VecSetSizes(datavector_Vec,this->tsmodel->get_datavectorlength_local(),this->tsmodel->get_datavectorlength_global()) );
		TRY( VecSetFromOptions(datavector_Vec) );
		this->datavector = new GeneralVector<PetscVector>(datavector_Vec);
		this->destroy_datavector = true;
	}

	if(!this->gammavector){
		Vec gammavector_Vec;

		#ifndef USE_CUDA
			/* classic MPI vector */
			TRY( VecCreate(PETSC_COMM_WORLD,&gammavector_Vec) );
			TRY( VecSetSizes(gammavector_Vec,this->tsmodel->get_gammavectorlength_local(),this->tsmodel->get_gammavectorlength_global()) );
			TRY( VecSetFromOptions(gammavector_Vec) );
		#else
			/* CudaMPI vector */
			TRY( VecCreate(PETSC_COMM_WORLD,&gammavector_Vec) );
			TRY( VecSetType(gammavector_Vec, VECMPICUDA) );
			TRY( VecSetSizes(gammavector_Vec,this->tsmodel->get_gammavectorlength_local(),this->tsmodel->get_gammavectorlength_global()) );
			TRY( VecSetFromOptions(gammavector_Vec) );
		#endif

		this->gammavector = new GeneralVector<PetscVector>(gammavector_Vec);
		this->destroy_gammavector = true;
	}

	if(!this->thetavector){
		Vec thetavector_Vec;

		/* classic MPI vector */
		TRY( VecCreate(PETSC_COMM_WORLD,&thetavector_Vec) );
		TRY( VecSetSizes(thetavector_Vec,this->tsmodel->get_thetavectorlength_local(),this->tsmodel->get_thetavectorlength_global()) );
		TRY( VecSetFromOptions(thetavector_Vec) );
		
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
		output <<  " - T:           " << this->get_T() << std::endl;
		output <<  " - xdim:        " << this->get_xdim() << std::endl;
		output <<  " - K:           " << this->tsmodel->get_K() << std::endl;
		output <<  " - model:       " << this->tsmodel->get_name() << std::endl;
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
		output_global <<  " - T:           " << this->get_T() << std::endl;
		output_local  <<  "  - Tlocal:     " << this->tsmodel->get_Tlocal() << std::endl;
		output_local.synchronize();

		output_global <<  " - xdim:        " << this->get_xdim() << std::endl;
		output_global <<  " - K:           " << this->tsmodel->get_K() << std::endl;

		output_global <<  " - model:       " << this->tsmodel->get_name() << std::endl;
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
	if(this->tsmodel){
		return this->tsmodel->get_T();
	} else {
		return this->T;
	}
}

template<class VectorBase>
int TSData<VectorBase>::get_Tlocal() const{
	if(this->tsmodel){
		return this->tsmodel->get_Tlocal();
	} else {
		return this->Tlocal;
	}
}

template<class VectorBase>
int TSData<VectorBase>::get_Tbegin() const{
	if(this->tsmodel){
		return this->tsmodel->get_Tbegin();
	} else {
		return this->Tbegin;
	}
}

template<class VectorBase>
int TSData<VectorBase>::get_Tend() const{
	if(this->tsmodel){
		return this->tsmodel->get_Tend();
	} else {
		return this->Tend;
	}
}

template<class VectorBase>
int TSData<VectorBase>::get_xdim() const{
	if(this->tsmodel){
		return this->tsmodel->get_xdim();
	} else {
		return 0;
	}
}

template<class VectorBase>
int TSData<VectorBase>::get_K() const{
	if(this->tsmodel){
		return this->tsmodel->get_K();
	} else {
		return 0;
	}
}

template<class VectorBase>
int TSData<VectorBase>::get_xmem() const{
	if(this->tsmodel){
		return this->tsmodel->get_xmem();
	} else {
		return 0;
	}
}

template<class VectorBase>
TSModel<VectorBase> *TSData<VectorBase>::get_model() const{
	return this->tsmodel;
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
void TSData<PetscVector>::cut_gamma() const{
	LOG_FUNC_BEGIN

	int max_id;
	double max_value;
	
	int K = get_K();
	int gamma_t = this->tsmodel->get_gammavectorlength_local()/(double)K;
	
	double *gamma_arr;
	TRY( VecGetArray(gammavector->get_vector(),&gamma_arr) );
	
	int t,k;
	for(t = 0; t < gamma_t; t++){
		/* find max value */
		max_id = 0;
		max_value = gamma_arr[t];
		for(k = 1; k < K; k++){
			if(gamma_arr[k*gamma_t + t] > max_value){
				max_id = k;
				max_value = gamma_arr[k*gamma_t + t];
			}
		}
		
		/* set new values */
		for(k = 0; k < K; k++){
			if(k == max_id){
				gamma_arr[k*gamma_t + t] = 1.0;
			} else {
				gamma_arr[k*gamma_t + t] = 0.0;
			}
		}


	}

	TRY( VecRestoreArray(gammavector->get_vector(),&gamma_arr) );
	
	LOG_FUNC_END
}

template<>
void TSData<PetscVector>::save_datavector(std::string filename) const {
	LOG_FUNC_BEGIN

	PetscViewer viewer_out;
	TRY( PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename.c_str(),FILE_MODE_WRITE,&viewer_out) );
	TRY( VecView( datavector->get_vector(), viewer_out) );
	TRY( PetscViewerDestroy(&viewer_out) );

	LOG_FUNC_END
}

template<>
void TSData<PetscVector>::save_thetavector(std::string filename) const {
	
}

template<>
void TSData<PetscVector>::save_gammavector(std::string filename, int blocksize) const {
	
}


} /* end namespace */

#endif
