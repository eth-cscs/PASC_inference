#include "external/petscvector/data/edfdata.h"

namespace pascinference {
namespace data {

template<>
void EdfData<PetscVector>::edfRead(std::string filename, int max_record_nmb){
	LOG_FUNC_BEGIN

	/* open file */
	std::ifstream myfile(filename.c_str(), std::ios::in | std::ios::binary);

	myfile.seekg(0, std::ios::beg);

	int i;
	char buffer[100];

	/* ------ HEADER ------ */
	myfile.read(buffer, 8);
	hdr_ver = atoi(buffer);

	myfile.read(buffer, 80);
	hdr_patientID = std::string(buffer);

	myfile.read(buffer, 80);
	hdr_recordID = std::string(buffer);

	myfile.read(buffer, 8);
	hdr_startdate = std::string(buffer);

	myfile.read(buffer, 8);
	hdr_starttime = std::string(buffer);

	myfile.read(buffer, 8);
	hdr_bytes = atoi(buffer);

	myfile.read(buffer, 44);
	/* reserved */
	
	myfile.read(buffer, 8);
	hdr_records = atoi(buffer);

	/* cut the dataset if user provided max number of records */
	if(max_record_nmb > 0 && hdr_records > max_record_nmb){
		hdr_records = max_record_nmb;
	}

	myfile.read(buffer, 8);
	hdr_duration = atoi(buffer);

	myfile.read(buffer, 4);
	hdr_ns = atoi(buffer);

	/* arrays */
	hdr_records_detail = (Record*)malloc(sizeof(Record)*hdr_ns);
	free_hdr_records_detail = true;
	
	for(i=0;i<hdr_ns;i++){
		myfile.read(buffer, 16);
		hdr_records_detail[i].hdr_label = new std::string(buffer);
	}
			
	for(i=0;i<hdr_ns;i++){
		myfile.read(buffer, 80);
		hdr_records_detail[i].hdr_transducer = new std::string(buffer);
	}

	for(i=0;i<hdr_ns;i++){
		myfile.read(buffer, 8);
		hdr_records_detail[i].hdr_units = new std::string(buffer);
	}
	
	for(i=0;i<hdr_ns;i++){
		myfile.read(buffer, 8);
		hdr_records_detail[i].hdr_physicalMin = atof(buffer);
	}
	
	for(i=0;i<hdr_ns;i++){
		myfile.read(buffer, 8);
		hdr_records_detail[i].hdr_physicalMax = atof(buffer);
	}

	for(i=0;i<hdr_ns;i++){
		myfile.read(buffer, 8);
		hdr_records_detail[i].hdr_digitalMin = atof(buffer);
	}

	for(i=0;i<hdr_ns;i++){
		myfile.read(buffer, 8);
		hdr_records_detail[i].hdr_digitalMax = atof(buffer);
	}

	for(i=0;i<hdr_ns;i++){
		myfile.read(buffer, 80);
		hdr_records_detail[i].hdr_prefilter = new std::string(buffer);
	}

	for(i=0;i<hdr_ns;i++){
		myfile.read(buffer, 8);
		hdr_records_detail[i].hdr_samples = atoi(buffer);
	}	

	for(i=0;i<hdr_ns;i++){
		myfile.read(buffer, 32);
		/* reserved */
	}	

	/* ------ PREPARE DATAVECTOR ------ */
	/* compute vector lengths */
	this->Tpreliminary = hdr_records_detail[0].hdr_samples*hdr_records;
	int R = hdr_ns-1;

	/* prepare preliminary datavector and load data */
	Vec datapreload_Vec;
	TRYCXX( VecCreate(PETSC_COMM_WORLD,&datapreload_Vec) );
	#ifdef USE_CUDA
		TRYCXX(VecSetType(datepreload_Vec, VECMPICUDA));
	#endif

	TRYCXX( VecSetSizes(datapreload_Vec,PETSC_DECIDE,Tpreliminary*R) );
	TRYCXX( VecSetFromOptions(datapreload_Vec) );
	this->datavectorpreliminary = new GeneralVector<PetscVector>(datapreload_Vec);

	/* ------ RECORDS ------ */
	int recnum, ii, samplei, index;
	double scalefac, dc;
    
	int16_t value;

    for(recnum = 0; recnum < hdr_records; recnum++){
		for(ii = 0; ii < R; ii++){
			scalefac = (hdr_records_detail[ii].hdr_physicalMax - hdr_records_detail[ii].hdr_physicalMin)/(double)(hdr_records_detail[ii].hdr_digitalMax - hdr_records_detail[ii].hdr_digitalMin);
			dc = hdr_records_detail[ii].hdr_physicalMax - scalefac*hdr_records_detail[ii].hdr_digitalMax;

			for(samplei=0; samplei < hdr_records_detail[ii].hdr_samples; samplei++){
				myfile.read((char *)&value, sizeof(int16_t)); /* read block of memory */
				value =  value * scalefac + dc;
				index = ii*this->Tpreliminary + recnum*hdr_records_detail[ii].hdr_samples + samplei;
				TRYCXX( VecSetValue(datapreload_Vec, index, value, INSERT_VALUES) );
			}
        }
    }

	/* vector is prepared */
	TRYCXX( VecAssemblyBegin(datapreload_Vec) );
	TRYCXX( VecAssemblyEnd(datapreload_Vec) );

	/* close file */
    myfile.close();		

	TRYCXX( PetscBarrier(NULL) );

	LOG_FUNC_END
}

/* set decomposition - from preliminary to real data */
template<>
void EdfData<PetscVector>::set_decomposition(Decomposition<PetscVector> &new_decomposition) {
	LOG_FUNC_BEGIN

	this->decomposition = &new_decomposition;

	/* prepare real datavector */
	Vec data_Vec;
	this->decomposition->createGlobalVec_data(&data_Vec);
	this->datavector = new GeneralVector<PetscVector>(data_Vec);
	this->destroy_datavector = true;
	
	/* permute orig to new using parallel layout */
	Vec datapreload_Vec = datavectorpreliminary->get_vector();
	this->decomposition->permute_TRxdim(datapreload_Vec, data_Vec);
	
	/* destroy preliminary data */
	TRYCXX(VecDestroy(&datapreload_Vec));
	
	LOG_FUNC_END
}


/* from filename */
template<>
EdfData<PetscVector>::EdfData(std::string filename_data, int max_record_nmb){
	LOG_FUNC_BEGIN

	/* read data from input file */
	edfRead(filename_data, max_record_nmb);

	this->destroy_gammavector = false;
	this->destroy_thetavector = false;

	LOG_FUNC_END
}

/* destructor */
template<>
EdfData<PetscVector>::~EdfData(){
	LOG_FUNC_BEGIN
	
	/* if I created a datavector, then I should also be able to destroy it */
	if(this->free_hdr_records_detail){
		free(this->hdr_records_detail);
	}

	LOG_FUNC_END
}

template<>
void EdfData<PetscVector>::saveVTK(std::string filename) const{
	Timer timer_saveVTK; 
	timer_saveVTK.restart();
	timer_saveVTK.start();

	int T = get_T();
	int Tlocal = get_Tlocal();
	int Tbegin = get_Tbegin();
	int Tend = get_Tend();
	const int *Tranges = decomposition->get_DDT_ranges();

	int K = get_K();
	int R = get_R();
	int Rlocal = get_Rlocal();
	int *DDR_affiliation = decomposition->get_DDR_affiliation();
	int DDR_rank = decomposition->get_DDR_rank();
	int DDR_size = decomposition->get_DDR_size();

	int xdim = get_xdim();

	int prank = GlobalManager.get_rank();
	int psize = GlobalManager.get_size();

	/* to manipulate with filename */
	std::ostringstream oss_filename;

	/* to manipulate with file */
	std::ofstream myfile;

	/* master writes the main file */
	if(prank == 0){
		/* create folder */
		oss_filename << "results/" << filename << "_vtk";
		boost::filesystem::path dir(oss_filename.str().c_str());
		boost::filesystem::create_directory(dir);
		oss_filename.str("");
		
		/* write to the name of file */
		oss_filename << "results/" << filename << "_vtk/" << filename << ".pvd";
		myfile.open(oss_filename.str().c_str());
		oss_filename.str("");

		/* write header to file */
		myfile << "<?xml version=\"1.0\"?>" << std::endl;
		myfile << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << std::endl;
		myfile << "<Collection>" << std::endl;
		for(int t=0;t<T;t++){
			for(int r=0;r<DDR_size;r++){
				myfile << " <DataSet timestep=\"" << t << "\" group=\"\" part=\"r\" file=\"edf_" << r << "_" << t <<".vtu\"/>" << std::endl;
			}
		}

		myfile << "</Collection>" << std::endl;
		myfile << "</VTKFile>";
		
		myfile.close();
	}

	TRYCXX( PetscBarrier(NULL));

	/* compute recovered vector */
	Vec gammak_Vec;
	IS gammak_is;
	
	Vec data_recovered_Vec;
	TRYCXX( VecDuplicate(datavector->get_vector(), &data_recovered_Vec) );
	TRYCXX( VecSet(data_recovered_Vec,0.0));
	GeneralVector<PetscVector> data_recovered(data_recovered_Vec);

	double *theta_arr;
	TRYCXX( VecGetArray(thetavector->get_vector(),&theta_arr) );

	for(int k=0;k<K;k++){ 
		/* get gammak */
		this->decomposition->createIS_gammaK(&gammak_is, k);
		TRYCXX( VecGetSubVector(gammavector->get_vector(), gammak_is, &gammak_Vec) );

		/* add to recovered image */
		TRYCXX( VecAXPY(data_recovered_Vec, theta_arr[k], gammak_Vec) );

		TRYCXX( VecRestoreSubVector(gammavector->get_vector(), gammak_is, &gammak_Vec) );
		TRYCXX( ISDestroy(&gammak_is) );
	}	

	double *data_arr;
	TRYCXX( VecGetArray(datavector->get_vector(), &data_arr) );

	double *data_recovered_arr;
	TRYCXX( VecGetArray(data_recovered_Vec, &data_recovered_arr) );

	double *gamma_arr;
	TRYCXX( VecGetArray(gammavector->get_vector(), &gamma_arr) );

	int coordinates_dim = tsmodel->get_coordinatesVTK_dim();
	double *coordinates_arr;
	TRYCXX( VecGetArray(tsmodel->get_coordinatesVTK()->get_vector(), &coordinates_arr) );

	double gamma_max;
	int gamma_maxk;

	/* each processor writes its own portion of data */
	for(int t=Tbegin;t < Tend;t++){
		oss_filename.str("");
		oss_filename << "results/" << filename << "_vtk/edf_" << DDR_rank << "_" << t << ".vtu";

		myfile.open(oss_filename.str().c_str());

		myfile << "<?xml version=\"1.0\"?>" << std::endl;
		myfile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
		myfile << "  <UnstructuredGrid>" << std::endl;
		myfile << "	  <Piece NumberOfPoints=\"" << Rlocal << "\" NumberOfCells=\"0\" >" << std::endl;
		myfile << "      <PointData Scalars=\"scalars\">" << std::endl;

		/* original data */
		myfile << "        <DataArray type=\"Float32\" Name=\"original\" format=\"ascii\">" << std::endl;
		for(int r=0;r<Rlocal;r++){
			myfile << data_arr[(t-Tbegin)*Rlocal+r] << std::endl;
		}
		myfile << "        </DataArray>" << std::endl;

		/* value of gamma */
		for(int k=0;k<K;k++){
			myfile << "        <DataArray type=\"Float32\" Name=\"gamma_" << k << "\" format=\"ascii\">" << std::endl;
			for(int r=0;r<Rlocal;r++){
				myfile << gamma_arr[(t-Tbegin)*Rlocal*K+r*K+k] << std::endl;
			}
			myfile << "        </DataArray>" << std::endl;
		}

		/* cluster affiliation */
		myfile << "        <DataArray type=\"Float32\" Name=\"gamma_max\" format=\"ascii\">" << std::endl;
		for(int r=0;r<Rlocal;r++){
			gamma_max = 0.0;
			gamma_maxk = 0;
			for(int k=0;k<K;k++){
				if(gamma_arr[(t-Tbegin)*Rlocal*K + r*K + k] > gamma_max){
					gamma_max = gamma_arr[(t-Tbegin)*Rlocal*K + r*K + k];
					gamma_maxk = k;
				}
			}
			myfile << gamma_maxk << std::endl;
		}
		myfile << "        </DataArray>" << std::endl;

		/* original data */
		myfile << "        <DataArray type=\"Float32\" Name=\"recovered\" format=\"ascii\">" << std::endl;
		for(int r=0;r<Rlocal;r++){
			myfile << data_recovered_arr[(t-Tbegin)*Rlocal+r] << std::endl;
		}
		myfile << "        </DataArray>" << std::endl;
		myfile << "      </PointData>" << std::endl;

		myfile << "      <CellData Scalars=\"scalars\">" << std::endl;
		myfile << "      </CellData>" << std::endl;
		myfile << "      <Points>" << std::endl;
		myfile << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;

		for(int r=0;r<R;r++){
			if(DDR_rank == DDR_affiliation[r]){
				/* 1D */
				if(coordinates_dim == 1){
					myfile << coordinates_arr[r] << " 0 0" << std::endl;
				}

				/* 2D */
				if(coordinates_dim == 2){
					myfile << coordinates_arr[r] << " " << coordinates_arr[r+R] << " 0" << std::endl;
				}

				/* 3D */
				if(coordinates_dim == 3){
					myfile << coordinates_arr[r] << coordinates_arr[r+R] << " " << coordinates_arr[r+2*R] << std::endl;
				}
			}
		}
		
		myfile << "        </DataArray>" << std::endl;
		myfile << "      </Points>" << std::endl;
		myfile << "      <Cells>" << std::endl;
		myfile << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" << std::endl;
		myfile << "        </DataArray>" << std::endl;
		myfile << "		<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << std::endl;
		myfile << "        </DataArray>" << std::endl;
		myfile << "		<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << std::endl;
		myfile << "        </DataArray>" << std::endl;
		myfile << "      </Cells>" << std::endl;
		myfile << "    </Piece>" << std::endl;
		myfile << "  </UnstructuredGrid>" << std::endl;
		myfile << "</VTKFile>" << std::endl;

		myfile.close();
	}

	TRYCXX( VecRestoreArray(gammavector->get_vector(), &gamma_arr) );
	TRYCXX( VecRestoreArray(thetavector->get_vector(), &theta_arr) );
	TRYCXX( VecRestoreArray(datavector->get_vector(), &data_arr) );
	TRYCXX( VecRestoreArray(data_recovered_Vec, &data_arr) );
	TRYCXX( VecRestoreArray(tsmodel->get_coordinatesVTK()->get_vector(), &coordinates_arr) );

	timer_saveVTK.stop();
	coutAll <<  " - problem saved to VTK in: " << timer_saveVTK.get_value_sum() << std::endl;
	coutAll.synchronize();
}

template<>
void EdfData<PetscVector>::saveVector(std::string filename, bool save_original) const{
	Timer timer_saveVector; 
	timer_saveVector.restart();
	timer_saveVector.start();

	std::ostringstream oss_name_of_file;

	/* prepare vectors to save as a permutation to original layout */
	Vec datasave_Vec;
	this->decomposition->createGlobalVec_data(&datasave_Vec);
	GeneralVector<PetscVector> datasave(datasave_Vec);

	Vec gammasave_Vec;
	this->decomposition->createGlobalVec_gamma(&gammasave_Vec);
	GeneralVector<PetscVector> gammasave(gammasave_Vec);

	/* save datavector - just for fun; to see if it was loaded in a right way */
	if(save_original){
		oss_name_of_file << "results/" << filename << "_original.bin";
		this->decomposition->permute_TRxdim(datasave_Vec, datavector->get_vector(), true);
		datasave.save_binary(oss_name_of_file.str());
		oss_name_of_file.str("");
	}

	/* save gamma */
	oss_name_of_file << "results/" << filename << "_gamma.bin";
	this->decomposition->permute_TRK(gammasave_Vec, gammavector->get_vector(), true);
	gammasave.save_binary(oss_name_of_file.str());
	oss_name_of_file.str("");

	/* compute recovered signal */
	Vec gammak_Vec;
	IS gammak_is;

	Vec data_recovered_Vec;
	TRYCXX( VecDuplicate(datavector->get_vector(), &data_recovered_Vec) );
	TRYCXX( VecSet(data_recovered_Vec,0.0));
	GeneralVector<PetscVector> data_recovered(data_recovered_Vec);

	double *theta_arr;
	TRYCXX( VecGetArray(thetavector->get_vector(),&theta_arr) );

	int K = this->get_K();

	for(int k=0;k<K;k++){ 
		/* get gammak */
		this->decomposition->createIS_gammaK(&gammak_is, k);
		TRYCXX( VecGetSubVector(gammavector->get_vector(), gammak_is, &gammak_Vec) );

		/* add to recovered image */
		TRYCXX( VecAXPY(data_recovered_Vec, theta_arr[k], gammak_Vec) );

		TRYCXX( VecRestoreSubVector(gammavector->get_vector(), gammak_is, &gammak_Vec) );
		TRYCXX( ISDestroy(&gammak_is) );
	}	

	TRYCXX( VecRestoreArray(thetavector->get_vector(),&theta_arr) );

	/* save recovered data */
	oss_name_of_file << "results/" << filename << "_recovered.bin";
	
	/* but at first, permute recovered data, datasave can be used */
	this->decomposition->permute_TRxdim(datasave_Vec, data_recovered_Vec, true);
	datasave.save_binary(oss_name_of_file.str());
	oss_name_of_file.str("");

	/* destroy vectors with original layout */
//	TRYCXX( VecDestroy(&datasave_Vec) );

	timer_saveVector.stop();
	coutAll <<  " - problem saved in: " << timer_saveVector.get_value_sum() << std::endl;
	coutAll.synchronize();
}


}
} /* end namespace */

