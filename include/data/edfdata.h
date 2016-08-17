/** @file edfdata.cu
 *  @brief this is only for PETSC!
 * 
 *  @author Lukas Pospisil
 */

#ifndef PASC_EDFDATA_H
#define	PASC_EDFDATA_H

#ifndef USE_PETSCVECTOR
 #error 'EDFDATA is for PETSCVECTOR'
#endif

typedef petscvector::PetscVector PetscVector;

#include <iostream>
#include "common/common.h"
#include "common/bgmgraph.h"
#include "model/tsmodel.h"
#include "data/tsdata.h"

namespace pascinference {

template<class VectorBase>
class EdfData: public TSData<VectorBase> {
	protected:
		/* informations from EDF file */
		struct Record {
			std::string *hdr_label;
			std::string *hdr_transducer;
			std::string *hdr_units;
			double hdr_physicalMin; // TODO: double or int?
			double hdr_physicalMax;
			double hdr_digitalMin;
			double hdr_digitalMax;
			std::string *hdr_prefilter;
			int hdr_samples;
		};
		int hdr_ver;
		std::string hdr_patientID;
		std::string hdr_recordID;
		std::string hdr_startdate;
		std::string hdr_starttime;
		int hdr_bytes;
		int hdr_records;
		int hdr_duration;
		int hdr_ns;
		Record *hdr_records_detail;
		bool free_hdr_records_detail;

		void edfRead(std::string filename, int max_record_nmb = -1);

		/* preliminary data */
		int Tpreliminary;
		GeneralVector<VectorBase> *datavectorpreliminary;

	public:
		EdfData(std::string filename_data, int max_record_nmb = -1);
		~EdfData();

		virtual void print(ConsoleOutput &output) const;
		virtual void print(ConsoleOutput &output_global, ConsoleOutput &output_local) const;

		virtual void printcontent(ConsoleOutput &output) const;
		virtual void printcontent(ConsoleOutput &output_global, ConsoleOutput &output_local) const;
		virtual std::string get_name() const;

		void saveVTK(std::string filename) const;

		int get_Tpreliminary() const;
		void set_decomposition(Decomposition &decomposition);

};

} // end of namespace

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {

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
	TRY( VecCreate(PETSC_COMM_WORLD,&datapreload_Vec) );
	TRY( VecSetSizes(datapreload_Vec,PETSC_DECIDE,Tpreliminary*R) );
	TRY( VecSetFromOptions(datapreload_Vec) );
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
				TRY( VecSetValue(datapreload_Vec, index, value, INSERT_VALUES) );
			}
        }
    }

	/* vector is prepared */
	TRY( VecAssemblyBegin(datapreload_Vec) );
	TRY( VecAssemblyEnd(datapreload_Vec) );

	/* close file */
    myfile.close();		

	LOG_FUNC_END
}

/* set decomposition - from preliminary to real data */
template<class VectorBase>
void EdfData<VectorBase>::set_decomposition(Decomposition &new_decomposition) {
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
	TRY(VecDestroy(&datapreload_Vec));
	
	LOG_FUNC_END
}


/* from filename */
template<class VectorBase>
EdfData<VectorBase>::EdfData(std::string filename_data, int max_record_nmb){
	LOG_FUNC_BEGIN

	/* read data from input file */
	edfRead(filename_data, max_record_nmb);

	this->destroy_gammavector = false;
	this->destroy_thetavector = false;

	LOG_FUNC_END
}

/* destructor */
template<class VectorBase>
EdfData<VectorBase>::~EdfData(){
	LOG_FUNC_BEGIN
	
	/* if I created a datavector, then I should also be able to destroy it */
	if(this->free_hdr_records_detail){
		free(this->hdr_records_detail);
	}

	LOG_FUNC_END
}


/* print info about data */
template<class VectorBase>
void EdfData<VectorBase>::print(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output <<  this->get_name() << std::endl;
	
	/* give information about presence of the data */
	output <<  " - version of this data format:            " << hdr_ver << std::endl;
	output <<  " - local patient identification:           " << hdr_patientID << std::endl;
	output <<  " - local recording identification:         " << hdr_recordID << std::endl;
	output <<  " - startdate of recording (dd.mm.yy):      " << hdr_startdate << std::endl;
	output <<  " - starttime of recording (hh.mm.ss):      " << hdr_starttime << std::endl;
	output <<  " - number of bytes in header record:       " << hdr_bytes << std::endl;
	output <<  " - number of data records (-1 if unknown): " << hdr_records << std::endl;
	output <<  " - duration of a data record, in seconds:  " << hdr_duration << std::endl;
	output <<  " - number of signals (ns) in data record:  " << hdr_ns << std::endl;
/*
	output <<  " - record details:" << std::endl;
	for(int i=0;i<hdr_ns;i++){
		output <<  "   - id:                 " << i << std::endl;
		output <<  "     label:              " << *hdr_records_detail[i].hdr_label << std::endl;
		output <<  "     transducer type:    " << *hdr_records_detail[i].hdr_transducer << std::endl;
		output <<  "     physical dimension: " << *hdr_records_detail[i].hdr_units << std::endl;
		output <<  "     physical minimum:   " << hdr_records_detail[i].hdr_physicalMin << std::endl;
		output <<  "     physical maximum:   " << hdr_records_detail[i].hdr_physicalMax << std::endl;
		output <<  "     digital minimum:    " << hdr_records_detail[i].hdr_digitalMin << std::endl;
		output <<  "     digital maximum:    " << hdr_records_detail[i].hdr_digitalMax << std::endl;
		output <<  "     prefiltering:       " << *hdr_records_detail[i].hdr_prefilter << std::endl;
		output <<  "     nr of samples:      " << hdr_records_detail[i].hdr_samples << std::endl;
	}
*/
	output <<  "----------------------------------------------------------------" << std::endl;

	output <<  " - Tpreliminary: " << this->get_T() << std::endl;

	if(this->decomposition){
		output <<  " - T           : " << this->get_T() << std::endl;
		output <<  " - xdim        : " << this->get_xdim() << std::endl;
		output <<  " - K           : " << this->get_K() << std::endl;
		output <<  " - R           : " << this->get_R() << std::endl;
	}

	if(this->tsmodel){
		output <<  " - model       : " << this->tsmodel->get_name() << std::endl;
	} else {
		output <<  " - model       : NO" << std::endl;
	}
	
	output <<  " - datavector  : ";
	if(this->datavector){
		output << "YES (size: " << this->datavector->size() << ")" << std::endl;
	} else {
		output << "NO" << std::endl;
	}
	output <<  " - gammavector : ";
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
void EdfData<VectorBase>::print(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN

	output_global <<  this->get_name() << std::endl;
	
	/* give information about presence of the data */
	output_global <<  " - version of this data format:            " << hdr_ver << std::endl;
	output_global <<  " - local patient identification:           " << hdr_patientID << std::endl;
	output_global <<  " - local recording identification:         " << hdr_recordID << std::endl;
	output_global <<  " - startdate of recording (dd.mm.yy):      " << hdr_startdate << std::endl;
	output_global <<  " - starttime of recording (hh.mm.ss):      " << hdr_starttime << std::endl;
	output_global <<  " - number of bytes in header record:       " << hdr_bytes << std::endl;
	output_global <<  " - number of data records (-1 if unknown): " << hdr_records << std::endl;
	output_global <<  " - duration of a data record, in seconds:  " << hdr_duration << std::endl;
	output_global <<  " - number of signals (ns) in data record:  " << hdr_ns << std::endl;
/*
	output_global <<  " - record details:" << std::endl;
	for(int i=0;i<hdr_ns;i++){
		output_global <<  "   - id:                 " << i << std::endl;
		output_global <<  "     label:              " << *hdr_records_detail[i].hdr_label << std::endl;
		output_global <<  "     transducer type:    " << *hdr_records_detail[i].hdr_transducer << std::endl;
		output_global <<  "     physical dimension: " << *hdr_records_detail[i].hdr_units << std::endl;
		output_global <<  "     physical minimum:   " << hdr_records_detail[i].hdr_physicalMin << std::endl;
		output_global <<  "     physical maximum:   " << hdr_records_detail[i].hdr_physicalMax << std::endl;
		output_global <<  "     digital minimum:    " << hdr_records_detail[i].hdr_digitalMin << std::endl;
		output_global <<  "     digital maximum:    " << hdr_records_detail[i].hdr_digitalMax << std::endl;
		output_global <<  "     prefiltering:       " << *hdr_records_detail[i].hdr_prefilter << std::endl;
		output_global <<  "     nr of samples:      " << hdr_records_detail[i].hdr_samples << std::endl;
	}
*/
	output_global <<  "----------------------------------------------------------------" << std::endl;
	
	/* give information about presence of the data */
	output_global <<  " - Tpreliminary: " << this->get_T() << std::endl;

	if(this->decomposition){
		output_global <<  " - T           : " << this->get_T() << std::endl;
		output_local  <<  "  - Tlocal     : " << this->get_Tlocal() << std::endl;
		output_local.synchronize();
		output_global <<  " - xdim        : " << this->get_xdim() << std::endl;
		output_global <<  " - K           : " << this->get_K() << std::endl;
		output_global <<  " - R           : " << this->get_R() << std::endl;
		output_local  <<  "  - Rlocal     : " << this->get_Rlocal() << std::endl;
		output_local.synchronize();
	}

	if(this->tsmodel){
		output_global <<  " - model       : " << this->tsmodel->get_name() << std::endl;
	} else {
		output_global <<  " - model       : NO" << std::endl;
	}
	
	output_global <<  " - datavector  : ";
	if(this->datavector){
		output_global << "YES (size: " << this->datavector->size() << ")" << std::endl;
		output_local  <<  "  - local size : " << this->datavector->local_size() << std::endl;
		output_local.synchronize();
	} else {
		output_global << "NO" << std::endl;
	}
	
	output_global <<  " - gammavector : ";
	if(this->gammavector){
		output_global << "YES (size: " << this->gammavector->size() << ")" << std::endl;
		output_local  <<  "  - local size : " << this->gammavector->local_size() << std::endl;
		output_local.synchronize();
	} else {
		output_global << "NO" << std::endl;
	}
	
	output_global << " - thetavector : ";
	if(this->thetavector){
		output_global << "YES (size: " << this->thetavector->size() << ")" << std::endl;
		output_local  <<  "  - local size : " << this->thetavector->local_size() << std::endl;
		output_local.synchronize();
	} else {
		output_global << "NO" << std::endl;
	}

	output_global.synchronize();

	LOG_FUNC_END
}

/* print content of all data */
template<class VectorBase>
void EdfData<VectorBase>::printcontent(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output << this->get_name() << std::endl;
	
	/* print the content of the data */
	output <<  " - datavector  : ";
	if(this->datavector){
		output << *this->datavector << std::endl;
	} else {
		output << "not set" << std::endl;
	}

	output <<  " - gammavector : ";
	if(this->gammavector){
		output << *this->gammavector << std::endl;
	} else {
		output << "not set" << std::endl;
	}

	output <<  " - thetavector : ";
	if(this->thetavector){
		output << *this->thetavector << std::endl;
	} else {
		output << "not set" << std::endl;
	}

	LOG_FUNC_END
}

/* print content of all data */
template<class VectorBase>
void EdfData<VectorBase>::printcontent(ConsoleOutput &output_global,ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN

	output_global <<  this->get_name() << std::endl;
	
	/* print the content of the data */
	output_local <<  " - datavector : ";
	if(this->datavector){
		output_local << *this->datavector << std::endl;
	} else {
		output_local << "not set" << std::endl;
	}
	output_local.synchronize();

	output_local <<  " - gammavector : ";
	if(this->gammavector){
		output_local << *this->gammavector << std::endl;
	} else {
		output_local << "not set" << std::endl;
	}
	output_local.synchronize();

	output_local <<  " - thetavector : ";
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
std::string EdfData<VectorBase>::get_name() const {
	return "EDF Time-series Data";
}

template<class VectorBase>
int EdfData<VectorBase>::get_Tpreliminary() const{
	return this->Tpreliminary;
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
		/* write to the name of file */
		oss_filename << "results/" << filename << ".pvd";
		myfile.open(oss_filename.str().c_str());
		oss_filename.str("");

		/* write header to file */
		myfile << "<?xml version=\"1.0\"?>\n";
		myfile << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">\n";
		myfile << "<Collection>\n";
		for(int t=0;t<T;t++){
			for(int r=0;r<DDR_size;r++){
				myfile << " <DataSet timestep=\"" << t << "\" group=\"\" part=\"r\" file=\"" << filename << "_" << r << "_" << t <<".vtu\"/>\n";
			}
		}

		myfile << "</Collection>\n";
		myfile << "</VTKFile>";
		
		myfile.close();
	}

	TRY( PetscBarrier(NULL));

	/* compute recovered vector */
	Vec gammak_Vec;
	IS gammak_is;
	
	Vec data_recovered_Vec;
	TRY( VecDuplicate(datavector->get_vector(), &data_recovered_Vec) );
	TRY( VecSet(data_recovered_Vec,0.0));
	GeneralVector<PetscVector> data_recovered(data_recovered_Vec);

	double *theta_arr;
	TRY( VecGetArray(thetavector->get_vector(),&theta_arr) );

	for(int k=0;k<K;k++){ 
		/* get gammak */
		this->decomposition->createIS_gammaK(&gammak_is, k);
		TRY( VecGetSubVector(gammavector->get_vector(), gammak_is, &gammak_Vec) );

		/* add to recovered image */
		TRY( VecAXPY(data_recovered_Vec, theta_arr[k], gammak_Vec) );

		TRY( VecRestoreSubVector(gammavector->get_vector(), gammak_is, &gammak_Vec) );
		TRY( ISDestroy(&gammak_is) );
	}	

	double *data_arr;
	TRY( VecGetArray(datavector->get_vector(), &data_arr) );

	double *data_recovered_arr;
	TRY( VecGetArray(data_recovered_Vec, &data_recovered_arr) );

	double *gamma_arr;
	TRY( VecGetArray(gammavector->get_vector(), &gamma_arr) );

	int coordinates_dim = tsmodel->get_coordinatesVTK_dim();
	double *coordinates_arr;
	TRY( VecGetArray(tsmodel->get_coordinatesVTK()->get_vector(), &coordinates_arr) );

	double gamma_max;
	int gamma_maxk;

	/* each processor writes its own portion of data */
	for(int t=Tbegin;t < Tend;t++){
		oss_filename.str("");
		oss_filename << "results/" << filename << "_" << DDR_rank << "_" << t << ".vtu";

		myfile.open(oss_filename.str().c_str());

		myfile << "<?xml version=\"1.0\"?>\n";
		myfile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
		myfile << "  <UnstructuredGrid>\n";
		myfile << "	  <Piece NumberOfPoints=\"" << Rlocal << "\" NumberOfCells=\"0\" >\n";
		myfile << "      <PointData Scalars=\"scalars\">\n";

		/* original data */
		myfile << "        <DataArray type=\"Float32\" Name=\"original\" format=\"ascii\">\n";
		for(int r=0;r<Rlocal;r++){
			myfile << data_arr[(t-Tbegin)*Rlocal+r] << "\n";
		}
		myfile << "        </DataArray>\n";

		/* value of gamma */
		for(int k=0;k<K;k++){
			myfile << "        <DataArray type=\"Float32\" Name=\"gamma_" << k << "\" format=\"ascii\">\n";
			for(int r=0;r<Rlocal;r++){
				myfile << gamma_arr[(t-Tbegin)*Rlocal*K+r*K+k] << "\n";
			}
			myfile << "        </DataArray>\n";
		}

		/* cluster affiliation */
		myfile << "        <DataArray type=\"Float32\" Name=\"gamma_max\" format=\"ascii\">\n";
		for(int r=0;r<Rlocal;r++){
			gamma_max = 0.0;
			gamma_maxk = 0;
			for(int k=0;k<K;k++){
				if(gamma_arr[(t-Tbegin)*Rlocal*K + r*K + k] > gamma_max){
					gamma_max = gamma_arr[(t-Tbegin)*Rlocal*K + r*K + k];
					gamma_maxk = k;
				}
			}
			myfile << gamma_maxk << "\n";
		}
		myfile << "        </DataArray>\n";

		/* original data */
		myfile << "        <DataArray type=\"Float32\" Name=\"recovered\" format=\"ascii\">\n";
		for(int r=0;r<Rlocal;r++){
			myfile << data_recovered_arr[(t-Tbegin)*Rlocal+r] << "\n";
		}
		myfile << "        </DataArray>\n";
		myfile << "      </PointData>\n";

		myfile << "      <CellData Scalars=\"scalars\">\n";
		myfile << "      </CellData>\n";
		myfile << "      <Points>\n";
		myfile << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";

		for(int r=0;r<R;r++){
			if(DDR_rank == DDR_affiliation[r]){
				/* 1D */
				if(coordinates_dim == 1){
					myfile << coordinates_arr[r] << " 0 0\n";
				}

				/* 2D */
				if(coordinates_dim == 2){
					myfile << coordinates_arr[r] << " " << coordinates_arr[r+R] << " 0\n";
				}

				/* 3D */
				if(coordinates_dim == 3){
					myfile << coordinates_arr[r] << coordinates_arr[r+R] << " " << coordinates_arr[r+2*R] << "\n";
				}
			}
		}
		
		myfile << "        </DataArray>\n";
		myfile << "      </Points>\n";
		myfile << "      <Cells>\n";
		myfile << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
		myfile << "        </DataArray>\n";
		myfile << "		<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
		myfile << "        </DataArray>\n";
		myfile << "		<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
		myfile << "        </DataArray>\n";
		myfile << "      </Cells>\n";
		myfile << "    </Piece>\n";
		myfile << "  </UnstructuredGrid>\n";
		myfile << "</VTKFile>\n";

		myfile.close();
	}

	TRY( VecRestoreArray(gammavector->get_vector(), &gamma_arr) );
	TRY( VecRestoreArray(thetavector->get_vector(), &theta_arr) );
	TRY( VecRestoreArray(datavector->get_vector(), &data_arr) );
	TRY( VecRestoreArray(data_recovered_Vec, &data_arr) );
	TRY( VecRestoreArray(tsmodel->get_coordinatesVTK()->get_vector(), &coordinates_arr) );

	timer_saveVTK.stop();
	coutAll <<  " - problem saved to VTK in: " << timer_saveVTK.get_value_sum() << std::endl;
	coutAll.synchronize();
}


} /* end namespace */

#endif
