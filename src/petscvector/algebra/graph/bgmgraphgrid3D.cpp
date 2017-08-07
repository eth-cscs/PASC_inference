#include "external/petscvector/algebra/graph/bgmgraphgrid3D.h"

namespace pascinference {
namespace algebra {

template<>
BGMGraphGrid3D<PetscVector>::BGMGraphGrid3D(int x_size, int y_size, int z_size) : BGMGraph<PetscVector>(){
	LOG_FUNC_BEGIN

	this->x_size = x_size;
	this->y_size = y_size;
	this->z_size = z_size;

	this->dim = 3;
	this->n = x_size*y_size*z_size;

	/* fill coordinates */
	Vec coordinates_Vec;
	TRYCXX( VecCreateSeq(PETSC_COMM_SELF, this->n*this->dim, &coordinates_Vec) );

	double *coordinates_arr;
	TRYCXX( VecGetArray(coordinates_Vec, &coordinates_arr) );

	for(int idx=0;idx<x_size*y_size*z_size;idx++){
		int z = idx/(double)(x_size*y_size); /* index of z */
		int y = (idx-z*x_size*y_size)/(double)(x_size); /* index of z */
		int x = idx- z*x_size*y_size - y*x_size; /* index of x */

		coordinates_arr[idx] = x;
		coordinates_arr[idx + this->n] = y;
		coordinates_arr[idx + 2*this->n] = z;
	}

	TRYCXX( VecRestoreArray(coordinates_Vec, &coordinates_arr) );

	this->coordinates = new GeneralVector<PetscVector>(coordinates_Vec);

	this->threshold = -1;
	this->processed = false;

	LOG_FUNC_END
}

template<>
void BGMGraphGrid3D<PetscVector>::process_grid(){
	LOG_FUNC_BEGIN

	this->threshold = 1.1;
	this->m = (2*x_size - 1)*(2*y_size - 1)*(2*z_size - 1) - x_size*y_size*z_size; /* number of edges */ 
	this->m_max = 6; /* max number of neighbours */

	/* prepare array for number of neighbors */
	this->neighbor_nmbs = (int*)malloc(this->n*sizeof(int));
	this->neighbor_ids = (int**)malloc(this->n*sizeof(int*));

//	#pragma omp parallel for
	for(int idx=0;idx<x_size*y_size*z_size;idx++){
		int z = idx/(double)(x_size*y_size); /* index of z */
		int y = (idx-z*x_size*y_size)/(double)(x_size); /* index of z */
		int x = idx- z*x_size*y_size - y*x_size; /* index of x */

		/* compute number of neighbors */
		int nmb = 0;
		if(x > 0)	nmb+=1;
		if(x < x_size-1) nmb+=1;
		if(y > 0)	nmb+=1;
		if(y < y_size-1) nmb+=1;
		if(z > 0)	nmb+=1;
		if(z < z_size-1) nmb+=1;

		this->neighbor_nmbs[idx] = nmb;
		this->neighbor_ids[idx] = (int*)malloc(this->neighbor_nmbs[idx]*sizeof(int));

		/* fill neighbors */
		nmb = 0;
		if(x > 0){ /* left */
			this->neighbor_ids[idx][nmb] = idx-1;
			nmb+=1;
		}
		if(x < x_size-1){ /* right */
			this->neighbor_ids[idx][nmb] = idx+1;
			nmb+=1;
		}
		if(y > 0){ /* bottom */
			this->neighbor_ids[idx][nmb] = idx-x_size;
			nmb+=1;
		}
		if(y < y_size-1){ /* top */
			this->neighbor_ids[idx][nmb] = idx+x_size;
			nmb+=1;
		}
		if(z > 0){ /* closer */
			this->neighbor_ids[idx][nmb] = idx-x_size*y_size;
			nmb+=1;
		}
		if(z < z_size-1){ /* distant */
			this->neighbor_ids[idx][nmb] = idx+x_size*y_size;
			nmb+=1;
		}
	}

	#ifdef USE_CUDA
//		externalcontent->process_grid_cuda(neighbor_nmbs, neighbor_ids);
	#endif

	this->processed = true;

	LOG_FUNC_END
}

template<>
void BGMGraphGrid3D<PetscVector>::saveVTK_bounding_box(std::string filename, int *bounding_box_local) const {
	LOG_FUNC_BEGIN

	Timer timer_saveVTK;
	timer_saveVTK.restart();
	timer_saveVTK.start();

	/* to manipulate with file */
	std::ofstream myfile;

    /* gather all values of bounding boxes to master */
    Vec bounding_box_Vec;
    TRYCXX( VecCreate(PETSC_COMM_WORLD, &bounding_box_Vec) );
    TRYCXX( VecSetSizes(bounding_box_Vec,6,6*GlobalManager.get_size()) );
    TRYCXX( VecSetFromOptions(bounding_box_Vec) );

    /* copy local values */
    double *bounding_box_arr;
    TRYCXX( VecGetArray(bounding_box_Vec, &bounding_box_arr) );
    for(int i=0;i<6;i++) bounding_box_arr[i] = bounding_box_local[i];
    TRYCXX( VecRestoreArray(bounding_box_Vec, &bounding_box_arr) );

    /* gather all values */
    IS all_idx;
    TRYCXX( ISCreateStride(PETSC_COMM_WORLD, GlobalManager.get_size()*6, 0,1, &all_idx) );

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
		myfile << "POINTS " << 6*GlobalManager.get_size() << " FLOAT" << std::endl;
        for(int j=0;j < GlobalManager.get_size(); j++){
			myfile << all_arr[j*6+0] << " " << all_arr[j*6+2] << " " << all_arr[j*6+4] << std::endl;
			myfile << all_arr[j*6+1] << " " << all_arr[j*6+2] << " " << all_arr[j*6+4] << std::endl;
			myfile << all_arr[j*6+0] << " " << all_arr[j*6+3] << " " << all_arr[j*6+4] << std::endl;
			myfile << all_arr[j*6+1] << " " << all_arr[j*6+3] << " " << all_arr[j*6+4] << std::endl;
			myfile << all_arr[j*6+0] << " " << all_arr[j*6+2] << " " << all_arr[j*6+5] << std::endl;
			myfile << all_arr[j*6+1] << " " << all_arr[j*6+2] << " " << all_arr[j*6+5] << std::endl;
			myfile << all_arr[j*6+0] << " " << all_arr[j*6+3] << " " << all_arr[j*6+5] << std::endl;
			myfile << all_arr[j*6+1] << " " << all_arr[j*6+3] << " " << all_arr[j*6+5] << std::endl;
        }

		myfile << "\nPOLYGONS " << 6*GlobalManager.get_size() << " " << 5*6*GlobalManager.get_size() << std::endl; /* actually, edges are here twice */
        for(int j=0;j < GlobalManager.get_size(); j++){
			myfile << 4 << " " << 8*j+0 << " " << 8*j+4 << " " << 8*j+5 << " " << 8*j+1 << std::endl;
			myfile << 4 << " " << 8*j+2 << " " << 8*j+3 << " " << 8*j+7 << " " << 8*j+6 << std::endl;
			myfile << 4 << " " << 8*j+0 << " " << 8*j+1 << " " << 8*j+3 << " " << 8*j+2 << std::endl;
			myfile << 4 << " " << 8*j+5 << " " << 8*j+4 << " " << 8*j+6 << " " << 8*j+7 << std::endl;
			myfile << 4 << " " << 8*j+1 << " " << 8*j+5 << " " << 8*j+7 << " " << 8*j+3 << std::endl;
			myfile << 4 << " " << 8*j+4 << " " << 8*j+0 << " " << 8*j+2 << " " << 8*j+6 << std::endl;
        }

		/* write domain affiliation */
		myfile << "\nCELL_DATA " << 6*GlobalManager.get_size() << std::endl;
		myfile << "SCALARS domain float 1" << std::endl;
		myfile << "LOOKUP_TABLE default" << std::endl;
        for(int j=0;j < GlobalManager.get_size(); j++){
			myfile << j << std::endl;
			myfile << j << std::endl;
			myfile << j << std::endl;
			myfile << j << std::endl;
			myfile << j << std::endl;
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

template<> BGMGraphGrid3D<PetscVector>::ExternalContent * BGMGraphGrid3D<PetscVector>::get_externalcontent() const {
	return externalcontent;
}

}
} /* end of namespace */
