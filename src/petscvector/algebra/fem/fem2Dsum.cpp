#include "external/petscvector/algebra/fem/fem2Dsum.h"

namespace pascinference {
namespace algebra {

template<>
void Fem2DSum<PetscVector>::reduce_gamma(GeneralVector<PetscVector> *gamma1, GeneralVector<PetscVector> *gamma2) const {
	LOG_FUNC_BEGIN

    if(!is_reduced()){
        TRYCXX( VecCopy(gamma1->get_vector(), gamma2->get_vector()));
    } else {
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

        int Rlocal1 = this->decomposition1->get_Rlocal();
        int Rlocal2 = this->decomposition2->get_Rlocal();

        int width1 = grid1->get_width();
        int width2 = grid2->get_width();
        int height1 = grid1->get_height();
        int height2 = grid2->get_height();
        int width_overlap1 = bounding_box1[1] - bounding_box1[0] + 1;
        int height_overlap1 = bounding_box1[3] - bounding_box1[2] + 1;

        for(int k=0;k<this->decomposition2->get_K();k++){

            /* get gammak */
            this->decomposition1->createIS_gammaK(&gammak1_is, k);
            this->decomposition2->createIS_gammaK(&gammak2_is, k);

            TRYCXX( VecGetSubVector(gamma1_Vec, gammak1_is, &gammak1_Vec) );
            TRYCXX( VecGetSubVector(gamma2_Vec, gammak2_is, &gammak2_Vec) );

            /* get local necessary part for local computation */
            TRYCXX( ISCreateGeneral(PETSC_COMM_WORLD,overlap1_idx_size,overlap1_idx,PETSC_USE_POINTER,&gammak1_overlap_is) );
            TRYCXX( VecGetSubVector(gammak1_Vec, gammak1_overlap_is, &gammak1_overlap_Vec) );

            //TODO:temp
//            TRYCXX( ISView(gammak1_overlap_is, PETSC_VIEWER_STDOUT_WORLD));

            /* sequential version */
            TRYCXX( VecGetArray(gammak1_overlap_Vec,&gammak1_arr) );
            TRYCXX( VecGetArray(gammak2_Vec,&gammak2_arr) );

            //TODO: OpenMP? GPU?
            for(int r2=0; r2 <= Rlocal2; r2++){
                /* r2 is in dR format, transfer it to R format */
                int r2R = DD_invpermutation2[Rbegin2+r2];

                /* coordinates in new image */
                int y2 = (int)(r2R/(double)width2);
                int x2 = r2R - y2*width2;

                /* coresponding point in original image */
                double x1 = x2*this->diff_x;
                double y1 = y2*this->diff_y;

                /* go through window and compute something */
                double value = 0.0;
                int nmb = 0;
                for(int xx1 = floor(x1 - 0.5*this->diff_x); xx1 <= ceil(x1 + 0.5*this->diff_x); xx1++){
                    for(int yy1 = floor(y1 - 0.5*this->diff_y); yy1 <= ceil(y1 + 0.5*this->diff_y); yy1++){
                        /* maybe we are out of the image */
                        if(xx1 >= 0 && xx1<width1 && yy1 >= 0 && yy1 < height1){
                            /* compute coordinate in overlap */
                            int r1_overlap = (yy1 - bounding_box1[2])*width_overlap1 + (xx1 - bounding_box1[0]);

                            value += gammak1_arr[r1_overlap];
                            nmb++;
                        }
                    }
                }
                /* write value */
                gammak2_arr[r2] = value/(double)nmb;

            }


            TRYCXX( VecRestoreArray(gammak1_overlap_Vec,&gammak1_arr) );
            TRYCXX( VecRestoreArray(gammak2_Vec,&gammak2_arr) );

            /* restore local necessary part for local computation */
            TRYCXX( VecRestoreSubVector(gammak1_Vec, gammak1_overlap_is, &gammak1_overlap_Vec) );
            TRYCXX( ISDestroy(&gammak1_overlap_is) );

            TRYCXX( VecRestoreSubVector(gamma1_Vec, gammak1_is, &gammak1_Vec) );
            TRYCXX( VecRestoreSubVector(gamma2_Vec, gammak2_is, &gammak2_Vec) );

            TRYCXX( ISDestroy(&gammak1_is) );
            TRYCXX( ISDestroy(&gammak2_is) );

        }
    }

	LOG_FUNC_END
}

template<>
void Fem2DSum<PetscVector>::prolongate_gamma(GeneralVector<PetscVector> *gamma2, GeneralVector<PetscVector> *gamma1) const {
	LOG_FUNC_BEGIN

    if(!is_reduced()){
        TRYCXX( VecCopy(gamma2->get_vector(), gamma1->get_vector()));
    } else {

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

		/* sequential version */
		TRYCXX( VecGetArray(gammak1_Vec,&gammak1_arr) );
		TRYCXX( VecGetArray(gammak2_overlap_Vec,&gammak2_arr) );

		//TODO: OpenMP? GPU?
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

		}

		TRYCXX( VecRestoreArray(gammak1_Vec,&gammak1_arr) );
		TRYCXX( VecRestoreArray(gammak2_overlap_Vec,&gammak2_arr) );

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

    }

	LOG_FUNC_END
}



}
} /* end of namespace */

