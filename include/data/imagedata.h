/** @file imagedata.cu
 *  @brief this is only for PETSC!
 * 
 *  @author Lukas Pospisil
 */

#ifndef PASC_IMAGEDATA_H
#define	PASC_IMAGEDATA_H

#ifndef USE_PETSCVECTOR
 #error 'IMAGEDATA is for PETSCVECTOR'
#endif

typedef petscvector::PetscVector PetscVector;

/* for debugging, if >= 100, then print info about ach called function */
extern int DEBUG_MODE;

#include <iostream>
#include "common/common.h"
#include "common/bgmgraph.h"
#include "matrix/blockgraphfree.h"
#include "model/tsmodel.h"
#include "data/tsdata.h"

namespace pascinference {

template<class VectorBase>
class ImageData: public TSData<VectorBase> {
	protected:
		int width;
		int height;
		int R;
	public:
		ImageData(std::string filename_data, int width, int height, int T=1);
		~ImageData();

		virtual void print(ConsoleOutput &output) const;
		virtual void print(ConsoleOutput &output_global, ConsoleOutput &output_local) const;

		virtual void printcontent(ConsoleOutput &output) const;
		virtual void printcontent(ConsoleOutput &output_global, ConsoleOutput &output_local) const;
		virtual std::string get_name() const;

		void saveImage(std::string filename) const;

		int get_R() const;

};

/* for simplier manipulation with graph of image */
class BGMGraphGrid2D: public BGMGraph {
	protected:
		int width;
		int height;
	public:
		BGMGraphGrid2D(int width, int height);
		BGMGraphGrid2D(std::string filename, int dim=2) : BGMGraph(filename, dim) {};
		BGMGraphGrid2D(const double *coordinates_array, int n, int dim) : BGMGraph(coordinates_array, n, dim) {};

		~BGMGraphGrid2D();
		
		virtual void process_grid();
};

class BGMGraphGrid1D: public BGMGraph {
	protected:
		int width;
	public:
		BGMGraphGrid1D(int width);
		BGMGraphGrid1D(std::string filename, int dim=2) : BGMGraph(filename, dim) {};
		BGMGraphGrid1D(const double *coordinates_array, int n, int dim) : BGMGraph(coordinates_array, n, dim) {};

		~BGMGraphGrid1D();
		
		virtual void process_grid();
};

} // end of namespace

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {

/* from filename */
template<class VectorBase>
ImageData<VectorBase>::ImageData(std::string filename_data, int width, int height, int T){
	LOG_FUNC_BEGIN

	this->width = width;
	this->height = height;
	this->R = (width*height)/(double)T;

	/* prepare new layout */
	Vec layout;
	TRY( VecCreate(PETSC_COMM_WORLD,&layout) );
	TRY( VecSetSizes(layout, PETSC_DECIDE, T ));
	TRY( VecSetFromOptions(layout) );

	/* get Tlocal */
	int Tlocal;
	TRY( VecGetLocalSize(layout,&Tlocal) );
	
	/* destroy layout vector - now we know everything what is necessary */
	TRY( VecDestroy(&layout) );
	
	/* ------ PREPARE DATAVECTOR ------ */
	Vec datavector_Vec;
	TRY( VecCreate(PETSC_COMM_WORLD,&datavector_Vec) );
	TRY( VecSetSizes(datavector_Vec, Tlocal*R, T*R ) );
	TRY( VecSetFromOptions(datavector_Vec) );

	this->datavector = new GeneralVector<PetscVector>(datavector_Vec);
	this->destroy_datavector = true;

	/* load image from file */
	this->datavector->load_global(filename_data);

	this->destroy_gammavector = false;
	this->destroy_thetavector = false;

	/* compute vector lengths */
	this->T = T;

	LOG_FUNC_END
}

/* destructor */
template<class VectorBase>
ImageData<VectorBase>::~ImageData(){
	LOG_FUNC_BEGIN
	
	
	LOG_FUNC_END
}


/* print info about data */
template<class VectorBase>
void ImageData<VectorBase>::print(ConsoleOutput &output) const {
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
	output <<  " - R:           " << this->get_R() << std::endl;
	
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
void ImageData<VectorBase>::print(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
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
	output_global <<  " - R:           " << this->get_R() << std::endl;
	
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
void ImageData<VectorBase>::printcontent(ConsoleOutput &output) const {
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
void ImageData<VectorBase>::printcontent(ConsoleOutput &output_global,ConsoleOutput &output_local) const {
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
std::string ImageData<VectorBase>::get_name() const {
	return "Image Time-series Data";
}

template<class VectorBase>
int ImageData<VectorBase>::get_R() const{
	return this->R;
}

template<>
void ImageData<PetscVector>::saveImage(std::string filename) const{
	Timer timer_saveImage; 
	timer_saveImage.restart();
	timer_saveImage.start();

	std::ostringstream oss_name_of_file;

	/* save datavector - just for fun; to see if it was loaded in a right way */
	oss_name_of_file << "results/" << filename << "_original.bin";
	datavector->save_binary(oss_name_of_file.str());
	oss_name_of_file.str("");

	/* save gammas and compute recovered image */
	Vec gamma_Vec = gammavector->get_vector();
	Vec gammak_Vec;
	GeneralVector<PetscVector> *gammak;
	IS gammak_is;

	Vec data_recovered_Vec;
	TRY( VecDuplicate(datavector->get_vector(), &data_recovered_Vec) );
	TRY( VecSet(data_recovered_Vec,0.0));
	GeneralVector<PetscVector> data_recovered(data_recovered_Vec);

	double *theta_arr;
	TRY( VecGetArray(thetavector->get_vector(),&theta_arr) );

	int K = get_K();
	int R = get_R();
	int Tlocal = get_Tlocal();
	int Tbegin = get_Tbegin();
	int k;

	for(k=0;k<K;k++){ 
		/* get gammak */
		TRY( ISCreateStride(PETSC_COMM_WORLD, R*Tlocal, Tbegin*K*R + k, K, &gammak_is) );
		TRY( VecGetSubVector(gamma_Vec, gammak_is, &gammak_Vec) );

		gammak = new GeneralVector<PetscVector>(gammak_Vec);

		/* save gammak */
		oss_name_of_file << "results/" << filename << "_gamma" << k << ".bin";
		gammak->save_binary(oss_name_of_file.str());
		oss_name_of_file.str("");

		oss_name_of_file << "results/" << filename << "_gamma" << k << ".txt";
		gammak->save_ascii(oss_name_of_file.str());
		oss_name_of_file.str("");


		/* add to recovered image */
		TRY( VecAXPY(data_recovered_Vec, theta_arr[k], gammak_Vec) );

		free(gammak);
	
		TRY( VecRestoreSubVector(gamma_Vec, gammak_is, &gammak_Vec) );
		TRY( ISDestroy(&gammak_is) );
	}	

	TRY( VecRestoreArray(thetavector->get_vector(),&theta_arr) );

	/* save recovered data */
	oss_name_of_file << "results/" << filename << "_recovered.bin";
	data_recovered.save_binary(oss_name_of_file.str());
	oss_name_of_file.str("");

	timer_saveImage.stop();
	coutAll <<  " - problem saved in: " << timer_saveImage.get_value_sum() << std::endl;
	coutAll.synchronize();
}


/* --------------- GraphImage implementation -------------- */
BGMGraphGrid2D::BGMGraphGrid2D(int width, int height) : BGMGraph(){
	this->width = width;
	this->height = height;

	this->dim = 2;
	this->n = width*height;
	
	/* fill coordinates */
	Vec coordinates_Vec;
	TRY( VecCreateSeq(PETSC_COMM_SELF, this->n*this->dim, &coordinates_Vec) );
	
	double *coordinates_arr;
	TRY( VecGetArray(coordinates_Vec, &coordinates_arr) );
	for(int j=0;j<height;j++){
		for(int i=0;i<width;i++){
			coordinates_arr[j*width + i] = i;
			coordinates_arr[j*width + i + this->n] = j;
		}
	}
	TRY( VecRestoreArray(coordinates_Vec, &coordinates_arr) );
	
	this->coordinates = new GeneralVector<PetscVector>(coordinates_Vec);

	this->threshold = -1;
	processed = false;
}

BGMGraphGrid2D::~BGMGraphGrid2D(){
	
}

void BGMGraphGrid2D::process_grid(){
	this->threshold = 1.1;
	this->m = height*(width-1) + width*(height-1);
	this->m_max = 4;

	/* prepare array for number of neighbors */
	neighbor_nmbs = (int*)malloc(n*sizeof(int));
	neighbor_ids = (int**)malloc(n*sizeof(int*));

//	#pragma omp parallel for
	for(int j=0;j<height;j++){
		for(int i=0;i<width;i++){
			int idx = j*width+i;

			/* compute number of neighbors */
			int nmb = 0;
			if(i>0){
				nmb+=1;				
			}
			if(i<width-1){
				nmb+=1;				
			}
			if(j>0){
				nmb+=1;				
			}
			if(j<height-1){
				nmb+=1;				
			}
			neighbor_nmbs[idx] = nmb;
			neighbor_ids[idx] = (int*)malloc(neighbor_nmbs[idx]*sizeof(int));
			
			/* fill neighbors */
			nmb = 0;
			if(i>0){ /* left */
				neighbor_ids[idx][nmb] = idx-1;
				nmb++;
			}
			if(i<width-1){ /* right */
				neighbor_ids[idx][nmb] = idx+1;
				nmb++;
			}
			if(j>0){ /* down */
				neighbor_ids[idx][nmb] = idx-width;
				nmb++;
			}
			if(j<height-1){ /* up */
				neighbor_ids[idx][nmb] = idx+width;
				nmb++;
			}

		}
	}

	#ifdef USE_GPU
		/* copy data to gpu */
		gpuErrchk( cudaMalloc((void **)&neighbor_nmbs_gpu, n*sizeof(int)) );	
		gpuErrchk( cudaMemcpy( neighbor_nmbs_gpu, neighbor_nmbs, n*sizeof(int), cudaMemcpyHostToDevice) );
		
		gpuErrchk( cudaMalloc((void **)&neighbor_ids_gpu, n*sizeof(int)) );	
		for(int i=0;i<n;i++){
			gpuErrchk( cudaMalloc((void **)&(neighbor_ids_gpu[i]), neighbor_nmbs[i]*sizeof(int)) );
			gpuErrchk( cudaMemcpy( neighbor_ids_gpu[i], neighbor_ids[i], n*sizeof(int), cudaMemcpyHostToDevice) );
		}

		gpuErrchk( cudaDeviceSynchronize() );
	#endif
	
	processed = true;
}

BGMGraphGrid1D::BGMGraphGrid1D(int width) : BGMGraph(){
	this->width = width;

	this->dim = 2;
	this->n = width;
	
	/* fill coordinates */
	Vec coordinates_Vec;
	TRY( VecCreateSeq(PETSC_COMM_SELF, this->n*this->dim, &coordinates_Vec) );
	
	double *coordinates_arr;
	TRY( VecGetArray(coordinates_Vec, &coordinates_arr) );
	for(int i=0;i<width;i++){
		coordinates_arr[i] = i;
		coordinates_arr[i + this->n] = 0;
	}
	TRY( VecRestoreArray(coordinates_Vec, &coordinates_arr) );
	
	this->coordinates = new GeneralVector<PetscVector>(coordinates_Vec);

	this->threshold = -1;
	processed = false;
}

BGMGraphGrid1D::~BGMGraphGrid1D(){
	
}

void BGMGraphGrid1D::process_grid(){
	this->threshold = 1.1;
	this->m = width-1;
	this->m_max = 2;

	/* prepare array for number of neighbors */
	neighbor_nmbs = (int*)malloc(n*sizeof(int));
	neighbor_ids = (int**)malloc(n*sizeof(int*));

//	#pragma omp parallel for
	for(int i=0;i<width;i++){
		int idx = i;

		/* compute number of neighbors */
		int nmb = 0;
		if(i>0){
			nmb+=1;				
		}
		if(i<width-1){
			nmb+=1;				
		}
		neighbor_nmbs[idx] = nmb;
		neighbor_ids[idx] = (int*)malloc(neighbor_nmbs[idx]*sizeof(int));
			
		/* fill neighbors */
		nmb = 0;
		if(i>0){ /* left */
			neighbor_ids[idx][nmb] = idx-1;
			nmb++;
		}
		if(i<width-1){ /* right */
			neighbor_ids[idx][nmb] = idx+1;
			nmb++;
		}
	}

	#ifdef USE_GPU
		/* copy data to gpu */
		gpuErrchk( cudaMalloc((void **)&neighbor_nmbs_gpu, n*sizeof(int)) );	
		gpuErrchk( cudaMemcpy( neighbor_nmbs_gpu, neighbor_nmbs, n*sizeof(int), cudaMemcpyHostToDevice) );
		
		gpuErrchk( cudaMalloc((void **)&neighbor_ids_gpu, n*sizeof(int)) );	
		for(int i=0;i<n;i++){
			gpuErrchk( cudaMalloc((void **)&(neighbor_ids_gpu[i]), neighbor_nmbs[i]*sizeof(int)) );
			gpuErrchk( cudaMemcpy( neighbor_ids_gpu[i], neighbor_ids[i], n*sizeof(int), cudaMemcpyHostToDevice) );
		}

		gpuErrchk( cudaDeviceSynchronize() );
	#endif
	
	processed = true;
}


} /* end namespace */

#endif
