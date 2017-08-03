#include "external/petscvector/algebra/graph/bgmgraphgrid2D.h"

namespace pascinference {
namespace algebra {

template<>
BGMGraphGrid2D<PetscVector>::BGMGraphGrid2D(int width, int height) : BGMGraph<PetscVector>(){
	LOG_FUNC_BEGIN

	this->width = width;
	this->height = height;

	this->dim = 2;
	this->n = width*height;

	/* fill coordinates */
	Vec coordinates_Vec;
	TRYCXX( VecCreateSeq(PETSC_COMM_SELF, this->n*this->dim, &coordinates_Vec) );

	double *coordinates_arr;
	TRYCXX( VecGetArray(coordinates_Vec, &coordinates_arr) );

	for(int idx=0;idx<width*height;idx++){
		int i = idx/(double)width; /* index of row */
		int j = idx - i*width; /* index of column */

		coordinates_arr[idx] = j;
		coordinates_arr[idx + this->n] = i;
	}

	TRYCXX( VecRestoreArray(coordinates_Vec, &coordinates_arr) );

	this->coordinates = new GeneralVector<PetscVector>(coordinates_Vec);

	this->threshold = -1;
	this->processed = false;

	LOG_FUNC_END
}

template<>
void BGMGraphGrid2D<PetscVector>::process_grid(){
	LOG_FUNC_BEGIN

	this->threshold = 1.1;
	this->m = height*(width-1) + width*(height-1);
	this->m_max = 4;

	/* prepare array for number of neighbors */
	this->neighbor_nmbs = (int*)malloc(this->n*sizeof(int));
	this->neighbor_ids = (int**)malloc(this->n*sizeof(int*));

//	#pragma omp parallel for
	for(int idx=0;idx<width*height;idx++){
		int i = idx/(double)width; /* index of row */
		int j = idx - i*width; /* index of column */

		/* compute number of neighbors */
		int nmb = 0;
		if(j>0){
			nmb+=1;
		}
		if(j<width-1){
			nmb+=1;
		}
		if(i>0){
			nmb+=1;
		}
		if(i<height-1){
			nmb+=1;
		}
		this->neighbor_nmbs[idx] = nmb;
		this->neighbor_ids[idx] = (int*)malloc(this->neighbor_nmbs[idx]*sizeof(int));

		/* fill neighbors */
		nmb = 0;
		if(j>0){ /* left */
			this->neighbor_ids[idx][nmb] = idx-1;
			nmb+=1;
		}
		if(j<width-1){ /* right */
			this->neighbor_ids[idx][nmb] = idx+1;
			nmb+=1;
		}
		if(i>0){ /* down */
			this->neighbor_ids[idx][nmb] = idx-width;
			nmb+=1;
		}
		if(i<height-1){ /* up */
			this->neighbor_ids[idx][nmb] = idx+width;
			nmb+=1;
		}
	}

	#ifdef USE_CUDA
		externalcontent->process_grid_cuda(neighbor_nmbs, neighbor_ids);
	#endif

	this->processed = true;

	LOG_FUNC_END
}

template<>
void BGMGraphGrid2D<PetscVector>::saveVTK_bounding_box(std::string filename, int *bounding_box_local) const {
	LOG_FUNC_BEGIN

	Timer timer_saveVTK;
	timer_saveVTK.restart();
	timer_saveVTK.start();

	/* to manipulate with file */
	std::ofstream myfile;

    /* gather all values of bounding boxes to master */
    Vec bounding_box_Vec;
    TRYCXX( VecCreate(PETSC_COMM_WORLD, &bounding_box_Vec) );
    TRYCXX( VecSetSizes(bounding_box_Vec,4,4*GlobalManager.get_size()) );
    TRYCXX( VecSetFromOptions(bounding_box_Vec) );

    /* copy local values */
    double *bounding_box_arr;
    TRYCXX( VecGetArray(bounding_box_Vec, &bounding_box_arr) );
    for(int i=0;i<4;i++) bounding_box_arr[i] = bounding_box_local[i];
    TRYCXX( VecRestoreArray(bounding_box_Vec, &bounding_box_arr) );

    /* gather all values */
    IS all_idx;
    TRYCXX( ISCreateStride(PETSC_COMM_WORLD, GlobalManager.get_size()*4, 0,1, &all_idx) );

    Vec all_Vec;
    TRYCXX( VecGetSubVector(bounding_box_Vec, all_idx, &all_Vec) );

    double *all_arr;
    TRYCXX(VecGetArray(all_Vec, &all_arr));

    /* only master writes VTK */
    if(GlobalManager.get_rank() == 0){
        std::ofstream myfile;

		myfile.open(filename.c_str());

		/* write header to file */
		myfile << "# vtk DataFile Version 3.1" << std::endl;
		myfile << "PASCInference: Bounding Box" << std::endl;
		myfile << "ASCII" << std::endl;
		myfile << "DATASET POLYDATA" << std::endl;

		/* write points - coordinates */
		myfile << "POINTS " << 4*GlobalManager.get_size() << " FLOAT" << std::endl;
        for(int j=0;j < GlobalManager.get_size(); j++){
			myfile << all_arr[j*4+0] << " " << all_arr[j*4+2] << " 0" << std::endl;
			myfile << all_arr[j*4+1] << " " << all_arr[j*4+2] << " 0" << std::endl;
			myfile << all_arr[j*4+1] << " " << all_arr[j*4+3] << " 0" << std::endl;
			myfile << all_arr[j*4+0] << " " << all_arr[j*4+3] << " 0" << std::endl;
        }

		myfile << "\nPOLYGONS " << GlobalManager.get_size() << " " << 5*GlobalManager.get_size() << std::endl; /* actually, edges are here twice */
        for(int j=0;j < GlobalManager.get_size(); j++){
			myfile << 4 << " " << 4*j+0 << " " << 4*j+1 << " " << 4*j+2 << " " << 4*j+3 << std::endl;
        }

		/* write domain affiliation */
		myfile << "\nCELL_DATA " << GlobalManager.get_size() << std::endl;
		myfile << "SCALARS domain float 1" << std::endl;
		myfile << "LOOKUP_TABLE default" << std::endl;
        for(int j=0;j < GlobalManager.get_size(); j++){
			myfile << j << std::endl;
		}

		myfile.close();
    }

    TRYCXX(VecRestoreArray(all_Vec, &all_arr));

    TRYCXX( VecRestoreSubVector(bounding_box_Vec, all_idx, &all_Vec) );

    TRYCXX( VecDestroy(&bounding_box_Vec));

	timer_saveVTK.stop();
	coutMaster <<  " - bounding box saved to VTK in: " << timer_saveVTK.get_value_sum() << std::endl;


	LOG_FUNC_END
}

template<> BGMGraphGrid2D<PetscVector>::ExternalContent * BGMGraphGrid2D<PetscVector>::get_externalcontent() const {
	return externalcontent;
}

}
} /* end of namespace */
