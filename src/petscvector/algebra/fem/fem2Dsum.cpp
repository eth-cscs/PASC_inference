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

        int Tbegin = this->decomposition1->get_Tbegin();
        int Tlocal = this->decomposition1->get_Tlocal();

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
            TRYCXX( ISCreateGeneral(PETSC_COMM_WORLD,overlap1_idx_size,overlap1_idx,PETSC_COPY_VALUES,&gammak1_overlap_is) );
            TRYCXX( VecGetSubVector(gammak1_Vec, gammak1_overlap_is, &gammak1_overlap_Vec) );

            /* sequential version */
            TRYCXX( VecGetArray(gammak1_overlap_Vec,&gammak1_arr) );
            TRYCXX( VecGetArray(gammak2_Vec,&gammak2_arr) );

            //TODO: OpenMP? GPU?
            for(int t=0; t < Tlocal; t++){
                for(int r2=0; r2 < Rlocal2; r2++){
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

                    for(int xx1 = floor(x1 - this->diff_x); xx1 <= ceil(x1 + this->diff_x); xx1++){
                        for(int yy1 = floor(y1 - this->diff_y); yy1 <= ceil(y1 + this->diff_y); yy1++){
                            /* compute coordinate in overlap */
                            int r1_overlap = (yy1 - bounding_box1[2])*width_overlap1 + (xx1 - bounding_box1[0]);
//                          int r1_overlap = ((int)y1 - bounding_box1[2])*width_overlap1 + ((int)x1 - bounding_box1[0]);

                            if(xx1 >= 0 & xx1 < width1 & yy1 >= 0 & yy1 <= height1 & r1_overlap >= 0 & r1_overlap < overlap1_idx_size){
                                value += gammak1_arr[t*this->overlap1_idx_size + r1_overlap];
                                nmb++;
                            }

                        }
                    }

                    /* write value */
                    gammak2_arr[t*Rlocal2 + r2] = value/(double)nmb; ////0.0001*r1_overlap;
                }
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
        IS gammak2_overlap_is;
        Vec gammak2_overlap_Vec;

        int *DD_permutation1 = grid1->get_DD_permutation();
        int *DD_invpermutation1 = grid1->get_DD_invpermutation();
        int *DD_permutation2 = grid2->get_DD_permutation();
        int *DD_invpermutation2 = grid2->get_DD_invpermutation();

        int Rbegin1 = this->decomposition1->get_Rbegin();
        int Rbegin2 = this->decomposition2->get_Rbegin();

        int Rlocal1 = this->decomposition1->get_Rlocal();
        int Rlocal2 = this->decomposition2->get_Rlocal();

        int Tlocal1 = this->decomposition1->get_Tlocal();

        int width1 = grid1->get_width();
        int width2 = grid2->get_width();
        int height1 = grid1->get_height();
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
            TRYCXX( ISCreateGeneral(PETSC_COMM_WORLD,overlap2_idx_size,overlap2_idx,PETSC_COPY_VALUES,&gammak2_overlap_is) ); //TODO:: PETSC_USE_POINTER ?
            TRYCXX( VecGetSubVector(gammak2_Vec, gammak2_overlap_is, &gammak2_overlap_Vec) );

            /* sequential version */
            TRYCXX( VecGetArray(gammak1_Vec,&gammak1_arr) );
            TRYCXX( VecGetArray(gammak2_overlap_Vec,&gammak2_arr) );

            //TODO: OpenMP? GPU?
            for(int t=0; t < Tlocal1; t++){
                for(int r1=0; r1 < Rlocal1; r1++){
                    /* r1 is in dR format, transfer it to R format */
                    int r1R = DD_invpermutation1[Rbegin1+r1];

                    /* coordinates in original image */
                    int y1 = (int)(r1R/(double)width1);
                    int x1 = r1R - y1*width1;

                    /* coresponding point in original image */
                    double x2 = x1*(double)(grid2->get_width())/((double)grid1->get_width());//this->diff_x;
                    double y2 = y1*(double)(grid2->get_height())/((double)grid1->get_height());//this->diff_y;

                    /* in this case, the window is only one point */
                    double value = 0.0;
                    int xx2 = floor(x2);
                    int yy2 = floor(y2);

                    /* compute coordinate in overlap */
                    int r2_overlap = (yy2 - bounding_box2[2])*width_overlap2 + (xx2 - bounding_box2[0]);

                    if(xx2 >= 0 & xx2 < width2 & yy2 >= 0 & yy2 <= height2 & r2_overlap >= 0 & r2_overlap < overlap2_idx_size){
						value += gammak2_arr[t*this->overlap2_idx_size + r2_overlap];
                    }

                    /* write value */
                    gammak1_arr[t*Rlocal1 + r1] = value;
                }
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
    }

	LOG_FUNC_END
}



}
} /* end of namespace */

