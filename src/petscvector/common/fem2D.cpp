#include "external/petscvector/common/fem2D.h"

namespace pascinference {
namespace common {

template<>
Fem2D<PetscVector>::Fem2D(Decomposition<PetscVector> *decomposition1, Decomposition<PetscVector> *decomposition2, double fem_reduce) : Fem<PetscVector>(decomposition1, decomposition2, fem_reduce){
	LOG_FUNC_BEGIN

	this->grid1 = (BGMGraphGrid2D<PetscVector>*)(this->decomposition1->get_graph());
	this->grid2 = (BGMGraphGrid2D<PetscVector>*)(this->decomposition2->get_graph());

	#ifdef USE_CUDA
		/* compute optimal kernel calls */
		externalcontent->cuda_occupancy();
		externalcontent->gridSize_reduce = (decomposition2->get_Tlocal() + externalcontent->blockSize_reduce - 1)/ externalcontent->blockSize_reduce;
		externalcontent->gridSize_prolongate = (decomposition2->get_Tlocal() + externalcontent->blockSize_prolongate - 1)/ externalcontent->blockSize_prolongate;
	#endif

	this->diff = 1; /* time */
	this->diff_x = (grid1->get_width()-1)/(double)(grid2->get_width()-1);
	this->diff_y = (grid1->get_height()-1)/(double)(grid2->get_height()-1);

	if(this->is_reduced()){
		this->bounding_box1 = new int[4];
		set_value_array(4, this->bounding_box1, 0); /* initial values */
		this->bounding_box2 = new int[4];
		set_value_array(4, this->bounding_box2, 0); /* initial values */

		compute_overlaps();
	}

	LOG_FUNC_END
}

template<>
void Fem2D<PetscVector>::compute_decomposition_reduced() {
	LOG_FUNC_BEGIN
	
	/* decomposition1 has to be set */
	this->grid1 = (BGMGraphGrid2D<PetscVector>*)(this->decomposition1->get_graph());

	if(this->is_reduced()){
		int T_reduced = 1;
		int width_reduced = ceil(grid1->get_width()*this->fem_reduce);
		int height_reduced = ceil(grid1->get_height()*this->fem_reduce);
		
		this->grid2 = new BGMGraphGrid2D<PetscVector>(width_reduced, height_reduced);
		this->grid2->process_grid();
		
		/* decompose second grid based on the decomposition of the first grid */
		this->grid2->decompose(this->grid1, this->bounding_box1, this->bounding_box2);
		this->grid2->print(coutMaster);
		
		/* compute new decomposition */
		this->decomposition2 = new Decomposition<PetscVector>(T_reduced, 
				*(this->grid2), 
				this->decomposition1->get_K(), 
				this->decomposition1->get_xdim(), 
				this->decomposition1->get_DDT_size(), 
				this->decomposition1->get_DDR_size());

		compute_overlaps();
	} else {
		/* there is not reduction of the data, we can reuse the decomposition */
		this->set_decomposition_reduced(this->decomposition1);
		this->grid2 = this->grid1;
	}

	#ifdef USE_CUDA
		/* compute optimal kernel calls */
		externalcontent->cuda_occupancy();
		externalcontent->gridSize_reduce = (decomposition2->get_Tlocal() + externalcontent->blockSize_reduce - 1)/ externalcontent->blockSize_reduce;
		externalcontent->gridSize_prolongate = (decomposition2->get_Tlocal() + externalcontent->blockSize_prolongate - 1)/ externalcontent->blockSize_prolongate;
	#endif

	this->diff = 1; /* time */
	this->diff_x = (grid1->get_width()-1)/(double)(grid2->get_width()-1);
	this->diff_y = (grid1->get_height()-1)/(double)(grid2->get_height()-1);
			
	LOG_FUNC_END
}

template<>
void Fem2D<PetscVector>::reduce_gamma(GeneralVector<PetscVector> *gamma1, GeneralVector<PetscVector> *gamma2) const {
	LOG_FUNC_BEGIN

	double *gammak1_arr;
	double *gammak2_arr;

	Vec gamma1_Vec = gamma1->get_vector();
	Vec gamma2_Vec = gamma2->get_vector();

	Vec gammak1_Vec;
	Vec gammak2_Vec;

	IS gammak1_is;
	IS gammak2_is;

	/* stuff for getting subvector for local computation */
	IS gammak1_overlap_is;
	Vec gammak1_overlap_Vec;

	int *DD_permutation1 = grid1->get_DD_permutation(); 
	int *DD_invpermutation1 = grid1->get_DD_invpermutation(); 
	int *DD_permutation2 = grid2->get_DD_permutation(); 
	int *DD_invpermutation2 = grid2->get_DD_invpermutation(); 

	int Rbegin1 = this->decomposition1->get_Rbegin();
	int Rbegin2 = this->decomposition2->get_Rbegin();

	int width1 = grid1->get_width();
	int width2 = grid2->get_width();
	int width_overlap1 = bounding_box1[1] - bounding_box1[0] + 1;
	int height_overlap1 = bounding_box1[3] - bounding_box1[2] + 1;

	for(int k=0;k<this->decomposition2->get_K();k++){

		/* get gammak */
		this->decomposition1->createIS_gammaK(&gammak1_is, k);
		this->decomposition2->createIS_gammaK(&gammak2_is, k);

		TRYCXX( VecGetSubVector(gamma1_Vec, gammak1_is, &gammak1_Vec) );
		TRYCXX( VecGetSubVector(gamma2_Vec, gammak2_is, &gammak2_Vec) );

		/* get local necessary part for local computation */
		TRYCXX( ISCreateGeneral(PETSC_COMM_SELF,overlap1_idx_size,overlap1_idx,PETSC_USE_POINTER,&gammak1_overlap_is) );
		TRYCXX( VecGetSubVector(gammak1_Vec, gammak1_overlap_is, &gammak1_overlap_Vec) );

		#ifndef USE_CUDA
			/* sequential version */
			TRYCXX( VecGetArray(gammak1_overlap_Vec,&gammak1_arr) );
			TRYCXX( VecGetArray(gammak2_Vec,&gammak2_arr) );

			//TODO: OpenMP?
			for(int r2=0; r2 < this->decomposition2->get_Rlocal(); r2++){
				int id2 = DD_permutation2[Rbegin2 + r2];
				int id_y2 = floor(id2/(double)width2);
				int id_x2 = id2 - id_y2*width2;

				/* coordinates in overlap */
				double center_x1 = (id_x2)*this->diff_x - bounding_box1[0];
				double left_x1 = (id_x2-1)*this->diff_x - bounding_box1[0];
				double right_x1 = (id_x2+1)*this->diff_x - bounding_box1[0];

				double center_y1 = (id_y2)*this->diff_y - bounding_box1[2];
				double left_y1 = (id_y2-1)*this->diff_y - bounding_box1[2];
				double right_y1 = (id_y2+1)*this->diff_y - bounding_box1[2];
				
				double mysum = 0.0;
				int counter = 0;
				for(int x1 = floor(left_x1); x1 < right_x1; x1++){
					for(int y1 = floor(left_y1); y1 < right_y1; y1++){
						if(x1 >= 0 && x1 < width_overlap1 && y1 >= 0 && y1 < height_overlap1){
							mysum += gammak1_arr[y1*width_overlap1 + x1];
							counter += 1;
						}
					}
				}
				gammak2_arr[r2] = mysum;// /(double)counter;
				
//				coutAll << "r2 = " << r2 << ", counter = " << counter << ", mysum = " << mysum << std::endl;
			}
//			coutAll.synchronize();

			TRYCXX( VecRestoreArray(gammak1_overlap_Vec,&gammak1_arr) );
			TRYCXX( VecRestoreArray(gammak2_Vec,&gammak2_arr) );
			
		#else
			/* cuda version */
			TRYCXX( VecCUDAGetArrayReadWrite(gammak1_overlap_Vec,&gammak1_arr) );
			TRYCXX( VecCUDAGetArrayReadWrite(gammak2_Vec,&gammak2_arr) );

			externalcontent->cuda_femhat_reduce_data(gammak1_arr, gammak2_arr, this->decomposition1->get_T(), this->decomposition2->get_T(), this->decomposition1->get_Tbegin(), this->decomposition2->get_Tbegin(), this->decomposition1->get_Tlocal(), this->decomposition2->get_Tlocal(), left_t1_idx, left_t2_idx, this->diff);

			TRYCXX( VecCUDARestoreArrayReadWrite(gammak1_overlap_Vec,&gammak1_arr) );
			TRYCXX( VecCUDARestoreArrayReadWrite(gammak2_Vec,&gammak2_arr) );			
		#endif

		/* restore local necessary part for local computation */
		TRYCXX( VecRestoreSubVector(gammak1_Vec, gammak1_overlap_is, &gammak1_overlap_Vec) );
		TRYCXX( ISDestroy(&gammak1_overlap_is) );

		TRYCXX( VecRestoreSubVector(gamma1_Vec, gammak1_is, &gammak1_Vec) );
		TRYCXX( VecRestoreSubVector(gamma2_Vec, gammak2_is, &gammak2_Vec) );

		TRYCXX( ISDestroy(&gammak1_is) );
		TRYCXX( ISDestroy(&gammak2_is) );

	}

	LOG_FUNC_END
}

template<>
void Fem2D<PetscVector>::prolongate_gamma(GeneralVector<PetscVector> *gamma2, GeneralVector<PetscVector> *gamma1) const {
	LOG_FUNC_BEGIN

	double *gammak1_arr;
	double *gammak2_arr;

	Vec gamma1_Vec = gamma1->get_vector();
	Vec gamma2_Vec = gamma2->get_vector();

	Vec gammak1_Vec;
	Vec gammak2_Vec;
	
	IS gammak1_is;
	IS gammak2_is;

	/* stuff for getting subvector for local computation */
	IS gammak2_overlap_is;
	Vec gammak2_overlap_Vec;

	int *DD_permutation1 = grid1->get_DD_permutation(); 
	int *DD_invpermutation1 = grid1->get_DD_invpermutation(); 
	int *DD_permutation2 = grid2->get_DD_permutation(); 
	int *DD_invpermutation2 = grid2->get_DD_invpermutation(); 

	int Rbegin1 = this->decomposition1->get_Rbegin();
	int Rbegin2 = this->decomposition2->get_Rbegin();

	int width1 = grid1->get_width();
	int height1 = grid1->get_height();
	int width2 = grid2->get_width();
	int height2 = grid2->get_height();
	int width_overlap2 = bounding_box2[1] - bounding_box2[0] + 1;
	int height_overlap2 = bounding_box2[3] - bounding_box2[2] + 1;

	for(int k=0;k<this->decomposition1->get_K();k++){

		/* get gammak */
		this->decomposition1->createIS_gammaK(&gammak1_is, k);
		this->decomposition2->createIS_gammaK(&gammak2_is, k);

		TRYCXX( VecGetSubVector(gamma1_Vec, gammak1_is, &gammak1_Vec) );
		TRYCXX( VecGetSubVector(gamma2_Vec, gammak2_is, &gammak2_Vec) );

		/* get local necessary part for local computation */
		TRYCXX( ISCreateGeneral(PETSC_COMM_SELF,overlap2_idx_size,overlap2_idx,PETSC_USE_POINTER,&gammak2_overlap_is) );
		TRYCXX( VecGetSubVector(gammak2_Vec, gammak2_overlap_is, &gammak2_overlap_Vec) );

		#ifndef USE_CUDA
			/* sequential version */
			TRYCXX( VecGetArray(gammak1_Vec,&gammak1_arr) );
			TRYCXX( VecGetArray(gammak2_overlap_Vec,&gammak2_arr) );

			//TODO: OpenMP?
			for(int r1=0; r1 < this->decomposition1->get_Rlocal(); r1++){
				int id1 = DD_invpermutation1[Rbegin1 + r1];
				int id_y1 = floor(id1/(double)width1);
				int id_x1 = id1 - id_y1*width1;

				/* coordinates in overlap */
				int center_x2 = floor((id_x1)/this->diff_x) - bounding_box2[0];
				int center_y2 = floor((id_y1)/this->diff_y) - bounding_box2[2];
				
//				gammak1_arr[r1] = GlobalManager.get_rank()/(double)GlobalManager.get_size();
//				gammak1_arr[r1] = id1/((double)(width1*height1));
				gammak1_arr[r1] = gammak2_arr[center_y2*width_overlap2 + center_x2];
				
//				coutAll << "r1 = " << r1 << ", value = " << gammak2_arr[center_y2*width_overlap2 + center_x2] << std::endl;
			}
//			coutAll.synchronize();

			TRYCXX( VecRestoreArray(gammak1_Vec,&gammak1_arr) );
			TRYCXX( VecRestoreArray(gammak2_overlap_Vec,&gammak2_arr) );
			
		#else
			/* cuda version */
			TRYCXX( VecCUDAGetArrayReadWrite(gammak1_overlap_Vec,&gammak1_arr) );
			TRYCXX( VecCUDAGetArrayReadWrite(gammak2_Vec,&gammak2_arr) );

			externalcontent->cuda_femhat_reduce_data(gammak1_arr, gammak2_arr, this->decomposition1->get_T(), this->decomposition2->get_T(), this->decomposition1->get_Tbegin(), this->decomposition2->get_Tbegin(), this->decomposition1->get_Tlocal(), this->decomposition2->get_Tlocal(), left_t1_idx, left_t2_idx, this->diff);

			TRYCXX( VecCUDARestoreArrayReadWrite(gammak1_overlap_Vec,&gammak1_arr) );
			TRYCXX( VecCUDARestoreArrayReadWrite(gammak2_Vec,&gammak2_arr) );			
		#endif

		/* restore local necessary part for local computation */
		TRYCXX( VecRestoreSubVector(gammak2_Vec, gammak2_overlap_is, &gammak2_overlap_Vec) );
		TRYCXX( ISDestroy(&gammak2_overlap_is) );

		TRYCXX( VecRestoreSubVector(gamma1_Vec, gammak1_is, &gammak1_Vec) );
		TRYCXX( VecRestoreSubVector(gamma2_Vec, gammak2_is, &gammak2_Vec) );

		TRYCXX( ISDestroy(&gammak1_is) );
		TRYCXX( ISDestroy(&gammak2_is) );
	}

	TRYCXX( VecRestoreArray(gamma1_Vec,&gammak1_arr) );
	TRYCXX( VecRestoreArray(gamma2_Vec,&gammak2_arr) );

	LOG_FUNC_END
}


}
} /* end of namespace */

