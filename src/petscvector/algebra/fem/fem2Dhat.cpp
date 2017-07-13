#include "external/petscvector/algebra/fem/fem2Dhat.h"

namespace pascinference {
namespace algebra {

template<>
void Fem2DHat<PetscVector>::reduce_gamma(GeneralVector<PetscVector> *gamma1, GeneralVector<PetscVector> *gamma2) const {
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
            TRYCXX( ISCreateGeneral(PETSC_COMM_WORLD,overlap1_idx_size,overlap1_idx,PETSC_COPY_VALUES,&gammak1_overlap_is) );
            TRYCXX( VecGetSubVector(gammak1_Vec, gammak1_overlap_is, &gammak1_overlap_Vec) );

            //TODO:temp
//            TRYCXX( ISView(gammak1_overlap_is, PETSC_VIEWER_STDOUT_WORLD));

            /* sequential version */
            TRYCXX( VecGetArray(gammak1_overlap_Vec,&gammak1_arr) );
            TRYCXX( VecGetArray(gammak2_Vec,&gammak2_arr) );

            //TODO: OpenMP? GPU?
            for(int r2=0; r2 < Rlocal2; r2++){
                /* r2 is in dR format, transfer it to R format */
                int r2R = DD_invpermutation2[Rbegin2+r2];

                /* coordinates in new image */
                int y2 = (int)(r2R/(double)width2);
                int x2 = r2R - y2*width2;

                /* coresponding point in original image */
                double x1 = x2*this->diff_x;
                double y1 = y2*this->diff_y;

                /* prepare hat vectors */
                /*
                 * P3 ------ P2
                 *  | \ p2 / |
                 *  |  \  /  |
                 *  |p3 P4 p1|
                 *  |  /  \  |
                 *  | / p0 \ |
                 * P0 ------ P1
                 */
                double P0[3], P1[3], P2[3], P3[3], P4[3]; /* x,y,f(x,y) */
                compute_window_values1(gammak1_arr, x1, y1, P0, P1, P2, P3, P4);

                /* go through window and compute something */
                double value = 0.0;
                int nmb = 0;

                for(double xx1 = P0[0]; xx1 <= P1[0]; xx1++){
                    for(double yy1 = P0[1]; yy1 <= P3[1]; yy1++){

                        double new_value;
                        double alpha, beta;

                        double counted = false;

                        /* p0 */
                        compute_plane_interpolation(&alpha,&beta,&new_value, xx1, yy1, P4, P0, P1);
                        if(alpha >= 0 & beta >= 0 & !counted){
                            value += new_value;
                            nmb++;

                            counted = true;
                        }

                        /* p1 */
                        compute_plane_interpolation(&alpha,&beta,&new_value, xx1, yy1, P4, P1, P2);
                        if(alpha >= 0 & beta >= 0 & !counted){
                            value += new_value;
                            nmb++;

                            counted = true;
                        }

                        /* p2 */
                        compute_plane_interpolation(&alpha,&beta,&new_value, xx1, yy1, P4, P2, P3);
                        if(alpha >= 0 & beta >= 0 & !counted){
                            value += new_value;
                            nmb++;

                            counted = true;
                        }

                        /* p3 */
                        compute_plane_interpolation(&alpha,&beta,&new_value, xx1, yy1, P4, P3, P0);
                        if(alpha >= 0 & beta >= 0 & !counted){
                            value += new_value;
                            nmb++;

                            counted = true;
                        }

                    }
                }

                /* write value */
//				gammak2_arr[r2] = 0.25*(P0[2] + P1[2] + P2[2] + P3[2] ); //value/(double)nmb; ////0.0001*r1_overlap;
//				gammak2_arr[r2] = 0.2*(P0[2] + P1[2] + P2[2] + P3[2] + P4[2] );
				gammak2_arr[r2] = value/(double)nmb; ////0.0001*r1_overlap;

//                coutMaster << gammak2_arr[r2] << "," << value << ", " << nmb << std::endl;

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
void Fem2DHat<PetscVector>::prolongate_gamma(GeneralVector<PetscVector> *gamma2, GeneralVector<PetscVector> *gamma1) const {
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
            for(int r1=0; r1 < Rlocal1; r1++){
                /* r1 is in dR format, transfer it to R format */
                int r1R = DD_invpermutation1[Rbegin1+r1];

                /* coordinates in original image */
                int y1 = (int)(r1R/(double)width1);
                int x1 = r1R - y1*width1;

                /* coresponding point in reduced image */
//                double x2 = x1*(double)(grid2->get_width())/((double)grid1->get_width());
//                double y2 = y1*(double)(grid2->get_height())/((double)grid1->get_height());
                double x2 = x1/this->diff_x;
                double y2 = y1/this->diff_y;

				/* there are 4 points which influencing the value in original image */
                /*
                 * P3 ------ P2
                 *  |        |
                 *  |        |
                 *  |   P    |
                 *  |        |
                 *  |        |
                 * P0 ------ P1
                 */
				double P0[3], P1[3], P2[3], P3[3]; /* x,y,f(x,y) */
                compute_window_values2(gammak2_arr, x2, y2, P0, P1, P2, P3);
                
                double alpha, beta, new_value;
                double value = 0.0;
				int nmb = 0;

				/* P0 */
				compute_plane_interpolation(&alpha,&beta,&new_value, x2, y2, P0, P1, P2);
				if(alpha >= 0 & beta >= 0){
					value += new_value;
					nmb++;
                }
				compute_plane_interpolation(&alpha,&beta,&new_value, x2, y2, P0, P2, P3);
				if(alpha >= 0 & beta >= 0){
					value += new_value;
					nmb++;
                }

				/* P1 */
				compute_plane_interpolation(&alpha,&beta,&new_value, x2, y2, P1, P0, P3);
				if(alpha >= 0 & beta >= 0){
					value += new_value;
					nmb++;
                }
				compute_plane_interpolation(&alpha,&beta,&new_value, x2, y2, P1, P2, P3);
				if(alpha >= 0 & beta >= 0){
					value += new_value;
					nmb++;
                }

				/* P2 */
				compute_plane_interpolation(&alpha,&beta,&new_value, x2, y2, P2, P0, P1);
				if(alpha >= 0 & beta >= 0){
					value += new_value;
					nmb++;
                }
				compute_plane_interpolation(&alpha,&beta,&new_value, x2, y2, P2, P0, P3);
				if(alpha >= 0 & beta >= 0){
					value += new_value;
					nmb++;
                }

				/* P3 */
				compute_plane_interpolation(&alpha,&beta,&new_value, x2, y2, P3, P0, P1);
				if(alpha >= 0 & beta >= 0){
					value += new_value;
					nmb++;
                }
				compute_plane_interpolation(&alpha,&beta,&new_value, x2, y2, P3, P1, P2);
				if(alpha >= 0 & beta >= 0){
					value += new_value;
					nmb++;
                }

                /* write value */
				gammak1_arr[r1] = value/(double)nmb;
//				gammak1_arr[r1] = (P0[2] + P1[2] + P2[2] + P3[2])/4.0;//value/(double)nmb;

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

