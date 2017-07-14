/** @file fem2Dhat.h
 *  @brief class for reduction and prolongation on fem meshes using hat functions in 2D
 *
 *  @author Lukas Pospisil
 */

#ifndef PASC_FEM2DHAT_H
#define	PASC_FEM2DHAT_H

#include "general/algebra/fem/fem2Dsum.h"
#include "general/algebra/graph/bgmgraphgrid2D.h"

namespace pascinference {
namespace algebra {

/** \class Fem2DHat
 *  \brief Reduction/prolongation between FEM meshes using linear interpolation.
 *
*/
template<class VectorBase>
class Fem2DHat : public Fem2DSum<VectorBase> {
	public:
		class ExternalContent;

	protected:
		friend class ExternalContent;
		ExternalContent *externalcontent;			/**< for manipulation with external-specific stuff */

        double get_imagevaluefromoverlap1(double* overlap_values, bool *is_inside, double xR, double yR) const;
        void compute_window_values1(double* overlap_values, double x1, double y1, double *PV, double *P0, double *P1, double *P2, double *P3, double *P4, double *P5, double *P6, double *P7) const;
        void compute_plane_interpolation(double *alpha, double *beta, double* value, double x, double y, double *PV, double *PA, double *PB, double *PC) const;

		double get_imagevaluefromoverlap2(double* overlap_values, bool *is_inside, int xR, int yR) const;
		void compute_window_values2(double* overlap_values, double x2, double y2, double *P0, double *P1, double *P2, double *P3) const;


		virtual double get_overlap_extension() const;
	public:
		/** @brief create FEM mapping between two decompositions
		*/
		Fem2DHat(Decomposition<VectorBase> *decomposition1, Decomposition<VectorBase> *decomposition2, double fem_reduce);

		/** @brief create general FEM mapping
		 *
		 * do not forget to call Fem::set_decomposition_original() and afterwards Fem::compute_decomposition_reduced() to compute decomposition2 internally
		 *
		 */
		Fem2DHat(double fem_reduce = 1.0);

		/** @brief destructor
		*/
		~Fem2DHat();

		virtual std::string get_name() const;

		virtual void reduce_gamma(GeneralVector<VectorBase> *gamma1, GeneralVector<VectorBase> *gamma2) const;
		virtual void prolongate_gamma(GeneralVector<VectorBase> *gamma2, GeneralVector<VectorBase> *gamma1) const;

		ExternalContent *get_externalcontent() const;
};

/* ----------------- Fem implementation ------------- */
template<class VectorBase>
Fem2DHat<VectorBase>::Fem2DHat(Decomposition<VectorBase> *decomposition1, Decomposition<VectorBase> *decomposition2, double fem_reduce) : Fem2DSum<VectorBase>(decomposition1, decomposition2, fem_reduce){
	LOG_FUNC_BEGIN

	LOG_FUNC_END
}

template<class VectorBase>
Fem2DHat<VectorBase>::Fem2DHat(double fem_reduce) : Fem2DSum<VectorBase>(fem_reduce){
	LOG_FUNC_BEGIN

	LOG_FUNC_END
};

template<class VectorBase>
Fem2DHat<VectorBase>::~Fem2DHat(){
	LOG_FUNC_BEGIN

	LOG_FUNC_END
}

template<class VectorBase>
std::string Fem2DHat<VectorBase>::get_name() const {
	return "FEM-2D-HAT";
}

template<class VectorBase>
void Fem2DHat<VectorBase>::reduce_gamma(GeneralVector<VectorBase> *gamma1, GeneralVector<VectorBase> *gamma2) const {
	LOG_FUNC_BEGIN

	//TODO

	LOG_FUNC_END
}

template<class VectorBase>
void Fem2DHat<VectorBase>::prolongate_gamma(GeneralVector<VectorBase> *gamma2, GeneralVector<VectorBase> *gamma1) const {
	LOG_FUNC_BEGIN

	//TODO

	LOG_FUNC_END
}

template<class VectorBase>
double Fem2DHat<VectorBase>::get_imagevaluefromoverlap1(double* overlap_values, bool *is_inside, double xR, double yR) const {
    double return_value = 0.0;

    if(floor(xR) >= 0 & ceil(xR) < this->grid1->get_width() & floor(yR) >= 0 & ceil(yR) < this->grid1->get_height() ){ // & r1_overlap >= 0 & r1_overlap < this->overlap1_idx_size
		/* compute average value of neighbour integer values */
		int r1_overlap1 = (floor(yR) - this->bounding_box1[2])*(this->bounding_box1[1] - this->bounding_box1[0] + 1) + (floor(xR) - this->bounding_box1[0]);
		int r1_overlap2 = (floor(yR) - this->bounding_box1[2])*(this->bounding_box1[1] - this->bounding_box1[0] + 1) + (ceil(xR) - this->bounding_box1[0]);
		int r1_overlap3 = (ceil(yR) - this->bounding_box1[2])*(this->bounding_box1[1] - this->bounding_box1[0] + 1) + (floor(xR) - this->bounding_box1[0]);
		int r1_overlap4 = (ceil(yR) - this->bounding_box1[2])*(this->bounding_box1[1] - this->bounding_box1[0] + 1) + (ceil(xR) - this->bounding_box1[0]);

		return_value = (overlap_values[r1_overlap1] + overlap_values[r1_overlap2] + overlap_values[r1_overlap3] + overlap_values[r1_overlap4])/4.0;

        *is_inside = true;
	} else {
        *is_inside = false;
	}

    return return_value;
}

template<class VectorBase>
double Fem2DHat<VectorBase>::get_imagevaluefromoverlap2(double* overlap_values, bool *is_inside, int xR, int yR) const {
    double return_value = 0.0;

    if(xR >= 0 & xR < this->grid2->get_width() & yR >= 0 & yR < this->grid2->get_height()){ 
		int r2_overlap = (yR - this->bounding_box2[2])*(this->bounding_box2[1] - this->bounding_box2[0] + 1) + (xR - this->bounding_box2[0]);
		return_value = overlap_values[r2_overlap];

        *is_inside = true;
	} else {
        *is_inside = false;
	}

    return return_value;
}

template<class VectorBase>
void Fem2DHat<VectorBase>::compute_plane_interpolation(double *alpha_out, double *beta_out,double* value, double x, double y, double *PV, double *PA, double *PB, double *PC) const {
    /* find coefficients alpha,beta in
        [x,y] = P4 + alpha*(P0-P4) + beta*(P1-P4)
       rewritting in form of system of linear equations
       [a,b; c,d] [alpha; beta] = [e,f]
    */

    double a1 = PA[0] - PV[0];
    double b1 = PB[0] - PV[0];
    double c1 = PA[1] - PV[1];
    double d1 = PB[1] - PV[1];
    double e1 = x - PV[0];
    double f1 = y - PV[1];

    double a2 = PB[0] - PV[0];
    double b2 = PC[0] - PV[0];
    double c2 = PB[1] - PV[1];
    double d2 = PC[1] - PV[1];
    double e2 = x - PV[0];
    double f2 = y - PV[1];

    double invDeterminant1 = 1.0/(b1*c1-a1*d1);
    double invDeterminant2 = 1.0/(b2*c2-a2*d2);

    double alpha1 = (-d1*e1 + b1*f1)*invDeterminant1;
    double beta1 = (c1*e1 - a1*f1)*invDeterminant1;

    double alpha2 = (-d2*e2 + b2*f2)*invDeterminant2;
    double beta2 = (c2*e2 - a2*f2)*invDeterminant2;

	if(alpha1 >= 0 & beta1 >= 0){
		*alpha_out = alpha1;
		*beta_out = beta1;
		*value = PV[2] + alpha1*(PA[2] - PV[2]) + beta1*(PB[2] - PV[2]);
	} else {
		*alpha_out = alpha2;
		*beta_out = beta2;
		*value = PV[2] + alpha2*(PB[2] - PV[2]) + beta2*(PC[2] - PV[2]);
	}
	
}

template<class VectorBase>
void Fem2DHat<VectorBase>::compute_window_values1(double* overlap_values, double x1, double y1, double *PV, double *P0, double *P1, double *P2, double *P3, double *P4, double *P5, double *P6, double *P7) const {
	/*
	* P6 -- P5 -- P4
	*  | \   |  / |
	*  |  \  | /  |
	* P7 -  PV  - P3
	*  |  /  | \  |
	*  | /   |  \ |
	* P0 -- P1 -- P2
	*/

    bool PV_inside;
    bool P0_inside;
    bool P1_inside;
    bool P2_inside;
    bool P3_inside;
    bool P4_inside;
    bool P5_inside;
    bool P6_inside;
    bool P7_inside;

    PV[0] = x1;
    PV[1] = y1;

    P0[0] = x1 - this->diff_x;
    P0[1] = y1 - this->diff_y;

    P1[0] = PV[0];
    P1[1] = P0[1];

    P2[0] = x1 + this->diff_x;
    P2[1] = P1[1];

    P3[0] = P2[0];
    P3[1] = PV[1];

    P4[0] = P3[0];
    P4[1] = y1 + this->diff_y;

    P5[0] = PV[0];
    P5[1] = P4[1];

    P6[0] = P0[0];
    P6[1] = P5[1];

    P7[0] = P0[0];
    P7[1] = PV[1];

    PV[2] = get_imagevaluefromoverlap1(overlap_values,&PV_inside,PV[0],PV[1]);
    P0[2] = get_imagevaluefromoverlap1(overlap_values,&P0_inside,P0[0],P0[1]);
    P1[2] = get_imagevaluefromoverlap1(overlap_values,&P1_inside,P1[0],P1[1]);
    P2[2] = get_imagevaluefromoverlap1(overlap_values,&P2_inside,P2[0],P2[1]);
    P3[2] = get_imagevaluefromoverlap1(overlap_values,&P3_inside,P3[0],P3[1]);
    P4[2] = get_imagevaluefromoverlap1(overlap_values,&P4_inside,P4[0],P4[1]);
    P5[2] = get_imagevaluefromoverlap1(overlap_values,&P5_inside,P4[0],P4[1]);
    P6[2] = get_imagevaluefromoverlap1(overlap_values,&P6_inside,P4[0],P4[1]);
    P7[2] = get_imagevaluefromoverlap1(overlap_values,&P7_inside,P4[0],P4[1]);

    if(!PV_inside){
        //TODO: throw error
        coutMaster << "ERROR in FEM2DHAT - center point is not inside the domain" << std::endl;
    }

    /* deal with corners */
    if(!P0_inside) P0[2] = PV[2];
    if(!P1_inside) P1[2] = PV[2];
    if(!P2_inside) P2[2] = PV[2];
    if(!P3_inside) P3[2] = PV[2];
    if(!P4_inside) P4[2] = PV[2];
    if(!P5_inside) P5[2] = PV[2];
    if(!P6_inside) P6[2] = PV[2];
    if(!P7_inside) P7[2] = PV[2];

}

template<class VectorBase>
void Fem2DHat<VectorBase>::compute_window_values2(double* overlap_values, double x2, double y2, double *P0, double *P1, double *P2, double *P3) const {

	/*
	* P3 ------ P2
	*  |        |
	*  |        |
	*  |   P    |
	*  |        |
	*  |        |
	* P0 ------ P1
	*/
                 
    bool P0_inside;
    bool P1_inside;
    bool P2_inside;
    bool P3_inside;

    P0[0] = floor(x2);
    P0[1] = floor(y2);

    P1[0] = P0[0]+1;//ceil(x2);
    P1[1] = P0[1];

    P2[0] = P1[0];
    P2[1] = P0[1]+1;//ceil(y2);

    P3[0] = P0[0];
    P3[1] = P2[1];

    P0[2] = get_imagevaluefromoverlap2(overlap_values,&P0_inside,P0[0],P0[1]);
    P1[2] = get_imagevaluefromoverlap2(overlap_values,&P1_inside,P1[0],P1[1]);
    P2[2] = get_imagevaluefromoverlap2(overlap_values,&P2_inside,P2[0],P2[1]);
    P3[2] = get_imagevaluefromoverlap2(overlap_values,&P3_inside,P3[0],P3[1]);


    if(!P1_inside) P1[2] = P0[2];
    if(!P2_inside) P2[2] = P0[2];
    if(!P3_inside) P3[2] = P0[2];

}


template<class VectorBase>
double Fem2DHat<VectorBase>::get_overlap_extension() const {
    return 1.0;
}


}
} /* end of namespace */

#endif
