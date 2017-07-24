#include "external/petscvector/data/tsdata.h"

namespace pascinference {
namespace data {

template<>
TSData<PetscVector>::TSData(Decomposition<PetscVector> &new_decomposition, GeneralVector<PetscVector> *datavector_new, GeneralVector<PetscVector> *gammavector_new, GeneralVector<PetscVector> *thetavector_new){
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
TSData<PetscVector>::TSData(Decomposition<PetscVector> &new_decomposition){
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
TSData<PetscVector>::TSData(Decomposition<PetscVector> &new_decomposition, std::string filename){
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

template<>
void TSData<PetscVector>::cutgamma() const {
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
void TSData<PetscVector>::save_thetavector(std::string filename) const {
	LOG_FUNC_BEGIN

	//TODO

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
void TSData<PetscVector>::save_datavector(std::string filename, int type) const {
	LOG_FUNC_BEGIN

	Timer timer_save;
	timer_save.restart();
	timer_save.start();

	std::ostringstream oss_name_of_file;

	/* prepare vectors to save as a permutation to original layout */
	Vec datasave_Vec;
	this->decomposition->createGlobalVec_data(&datasave_Vec);
	GeneralVector<PetscVector> datasave(datasave_Vec);

	/* save datavector - just for fun; to see if it was loaded in a right way */
	oss_name_of_file << "results/" << filename << "_datavector.bin";
    this->decomposition->permute_to_pdTRb(datasave_Vec, datavector->get_vector(), decomposition->get_xdim(), type, true);

	datasave.save_binary(oss_name_of_file.str());
	oss_name_of_file.str("");

	timer_save.stop();
	coutAll <<  " - datavector saved in: " << timer_save.get_value_sum() << std::endl;
	coutAll.synchronize();

	LOG_FUNC_END
}

template<>
void TSData<PetscVector>::save_gammavector(std::string filename) const {
	LOG_FUNC_BEGIN

	Timer timer_save;
	timer_save.restart();
	timer_save.start();

	std::ostringstream oss_name_of_file;

	oss_name_of_file << "results/" << filename << "_gamma.bin";

	Vec gammasave_Vec;
    TRYCXX( VecDuplicate(gammavector->get_vector(), &gammasave_Vec) );
	this->decomposition->permute_gTbR_to_pdTRb(gammasave_Vec, gammavector->get_vector(), decomposition->get_K(), true);
	GeneralVector<PetscVector> gammasave(gammasave_Vec);
	gammasave.save_binary(oss_name_of_file.str());

	oss_name_of_file.str("");

	timer_save.stop();
	coutAll <<  " - gammavector saved in: " << timer_save.get_value_sum() << std::endl;
	coutAll.synchronize();

	LOG_FUNC_END
}

template<>
void TSData<PetscVector>::save_reconstructed(std::string filename, int type) const {
	LOG_FUNC_BEGIN

	Timer timer_save;
	timer_save.restart();
	timer_save.start();

	std::ostringstream oss_name_of_file;

	/* compute recovered image */
	Vec gammak_Vec;
	IS gammak_is;
	Vec datan_recovered_Vec;
	IS datan_is;

	Vec data_recovered_Vec;
	TRYCXX( VecDuplicate(datavector->get_vector(), &data_recovered_Vec) );
	TRYCXX( VecSet(data_recovered_Vec,0.0));

	double *theta_arr;
	TRYCXX( VecGetArray(thetavector->get_vector(),&theta_arr) );

	int K = this->get_K();
	int xdim = this->get_xdim();

	for(int k=0;k<K;k++){
		/* get gammak */
		this->decomposition->createIS_gammaK(&gammak_is, k);
		TRYCXX( VecGetSubVector(gammavector->get_vector(), gammak_is, &gammak_Vec) );

		/* add to recovered image */
		for(int n=0;n<xdim;n++){
			/* get datan */
			this->decomposition->createIS_datan(&datan_is, n);
			TRYCXX( VecGetSubVector(data_recovered_Vec, datan_is, &datan_recovered_Vec) );

			/* add to recovered image */
			TRYCXX( VecAXPY(datan_recovered_Vec, theta_arr[k*xdim + n], gammak_Vec) );

			/* restore data */
			TRYCXX( VecRestoreSubVector(data_recovered_Vec, datan_is, &datan_recovered_Vec) );

			TRYCXX( ISDestroy(&datan_is) );
		}

		TRYCXX( VecRestoreSubVector(gammavector->get_vector(), gammak_is, &gammak_Vec) );
		TRYCXX( ISDestroy(&gammak_is) );
	}

	TRYCXX( VecRestoreArray(thetavector->get_vector(),&theta_arr) );

	/* permute recovered data */
	/* prepare vectors to save as a permutation to original layout */
	Vec datasave_Vec;
    TRYCXX( VecDuplicate(data_recovered_Vec, &datasave_Vec) );

    this->decomposition->permute_to_pdTRb(datasave_Vec, data_recovered_Vec, decomposition->get_xdim(), type, true);

	/* save recovered data */
	oss_name_of_file << "results/" << filename << "_recovered.bin";
	GeneralVector<PetscVector> datasave(datasave_Vec);
	datasave.save_binary(oss_name_of_file.str());
	oss_name_of_file.str("");

	timer_save.stop();
	coutAll <<  " - reconstructed signal saved in: " << timer_save.get_value_sum() << std::endl;
	coutAll.synchronize();

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
		output << " - total length:    " << std::setw(25) << x_size << std::endl;
		output << " - xdim:            " << std::setw(25) << get_xdim() << std::endl;
		output << " - R:               " << std::setw(25) << get_R() << std::endl;
		output << " - K:               " << std::setw(25) << get_K() << std::endl;
		output << " - T:               " << std::setw(25) << get_T() << std::endl;

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
/*		if(printdetails){
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
*/
	output.pop();
	output << std::setprecision(ss);

	LOG_FUNC_END
}


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

	//TODO: temp
	int test;
	std::cout << "prdel" << std::endl;
	TRYCXX( VecView(gamma0_Vec, PETSC_VIEWER_STDOUT_WORLD) );

	std::cin >> test;


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
void TSData<PetscVector>::load_gammavector(std::string filename, int type) const {
	LOG_FUNC_BEGIN

	gammavector->load_global(filename);

	Vec gammavector_preload_Vec;
	Vec gammavector_Vec = gammavector->get_vector();

	TRYCXX( VecDuplicate(gammavector_Vec,&gammavector_preload_Vec) );

    this->decomposition->permute_to_pdTRb(gammavector->get_vector(), gammavector_preload_Vec, decomposition->get_K(), type, false);

	TRYCXX( VecCopy(gammavector_preload_Vec, gammavector->get_vector()));
	TRYCXX( VecDestroy(&gammavector_preload_Vec) );

	LOG_FUNC_END
}

template<>
double TSData<PetscVector>::compute_gammavector_nbins() {
	LOG_FUNC_BEGIN

	double nbins = 1;
	bool is_changed;

	int K = this->get_K();
	int T = this->get_T();

	double *gamma_arr;
	TRYCXX( VecGetArray(gammavector->get_vector(), &gamma_arr) );

	for(int t=0;t<T-1;t++){
		is_changed = false;
		for(int k=0;k<K;k++){
			if(std::abs(gamma_arr[t*K + k] - gamma_arr[(t+1)*K+k]) > 1e-10 ){
				is_changed = true;
			}
		}

		if(is_changed) nbins++;
	}
	TRYCXX( VecRestoreArray(gammavector->get_vector(), &gamma_arr) );

	LOG_FUNC_END

	return nbins;
}

template<>
double TSData<PetscVector>::compute_abserr_reconstructed(GeneralVector<PetscVector> &solution) const {
	LOG_FUNC_BEGIN

	double abserr;

	Vec data_Vec = this->datavector->get_vector();
	Vec solution_Vec = solution.get_vector();
	Vec gammavector_Vec = this->gammavector->get_vector();

	/* vector of absolute error */
	Vec data_abserr_Vec;
	TRYCXX( VecDuplicate(data_Vec, &data_abserr_Vec) );

	/* compute recovered signal */
	Vec datan_abserr_Vec;
	IS datan_is;
	Vec gammak_Vec;
	IS gammak_is;

	/* abserr = -solution */
	TRYCXX( VecCopy(solution_Vec,data_abserr_Vec));
	TRYCXX( VecScale(data_abserr_Vec,-1.0));

	double *theta_arr;
	TRYCXX( VecGetArray(this->thetavector->get_vector(),&theta_arr) );

	int K = this->get_K();
	int xdim = this->get_xdim();

	for(int k=0;k<K;k++){
		/* get gammak */
		this->decomposition->createIS_gammaK(&gammak_is, k);
		TRYCXX( VecGetSubVector(gammavector_Vec, gammak_is, &gammak_Vec) );

		for(int n=0;n<xdim;n++){
			/* get datan */
			this->decomposition->createIS_datan(&datan_is, n);
			TRYCXX( VecGetSubVector(data_abserr_Vec, datan_is, &datan_abserr_Vec) );

			/* add to recovered image */
			TRYCXX( VecAXPY(datan_abserr_Vec, theta_arr[k*xdim + n], gammak_Vec) );

			/* restore data */
			TRYCXX( VecRestoreSubVector(data_abserr_Vec, datan_is, &datan_abserr_Vec) );
			TRYCXX( ISDestroy(&datan_is) );

		}

		TRYCXX( VecRestoreSubVector(gammavector_Vec, gammak_is, &gammak_Vec) );
		TRYCXX( ISDestroy(&gammak_is) );
	}

	TRYCXX( VecRestoreArray(this->thetavector->get_vector(),&theta_arr) );

	/* compute mean(abs(solution - data_recovered) */
//	TRYCXX( VecAbs(data_abserr_Vec) );
//	TRYCXX( VecSum(data_abserr_Vec, &abserr) );
//	int T = this->get_T();
//	abserr = abserr/(double)T;

	TRYCXX( VecNorm(data_abserr_Vec, NORM_1, &abserr) );

	LOG_FUNC_END

	return abserr;
}

}
} /* end namespace */


