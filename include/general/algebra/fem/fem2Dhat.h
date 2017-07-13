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

        double get_imagevaluefromoverlap1(double* overlap_values, bool *is_inside, int xR, int yR) const;
        void compute_window_values1(double* overlap_values, double x1, double y1, double *P0, double *P1, double *P2, double *P3, double *P4) const;
        void compute_plane_interpolation(double *alpha, double *beta, double* value, double x, double y, double *P4, double *P0, double *P1) const;

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
double Fem2DHat<VectorBase>::get_imagevaluefromoverlap1(double* overlap_values, bool *is_inside, int xR, int yR) const {
    double return_value = 0.0;

    int r1_overlap = (yR - this->bounding_box1[2])*(this->bounding_box1[1] - this->bounding_box1[0] + 1) + (xR - this->bounding_box1[0]);
    if(xR >= 0 & xR < this->grid1->get_width() & yR >= 0 & yR <= this->grid1->get_height() & r1_overlap >= 0 & r1_overlap < this->overlap1_idx_size){
		return_value = overlap_values[r1_overlap];
        *is_inside = true;
	} else {
        *is_inside = false;
	}

    return return_value;
}

template<class VectorBase>
void Fem2DHat<VectorBase>::compute_plane_interpolation(double *alpha_out, double *beta_out,double* value, double x, double y, double *P4, double *P0, double *P1) const {
    /* find coefficients alpha,beta in
        [x,y] = P4 + alpha*(P0-P4) + beta*(P1-P4)
       rewritting in form of system of linear equations
       [a,b; c,d] [alpha; beta] = [e,f]
    */

    double a = P0[0] - P4[0];
    double b = P1[0] - P4[0];
    double c = P0[1] - P4[1];
    double d = P1[1] - P4[1];
    double e = x - P4[0];
    double f = y - P4[1];

    double invDeterminant = 1.0/(b*c-a*d);

    double alpha = (d*e + b*f)*invDeterminant;
    double beta = (c*e + a*f)*invDeterminant;

    *alpha_out = alpha;
    *beta_out = beta;
    *value = P4[2] + alpha*(P0[2] - P4[2]) + beta*(P1[2] - P4[2]);
}

template<class VectorBase>
void Fem2DHat<VectorBase>::compute_window_values1(double* overlap_values, double x1, double y1, double *P0, double *P1, double *P2, double *P3, double *P4) const {

    bool P0_inside;
    bool P1_inside;
    bool P2_inside;
    bool P3_inside;
    bool P4_inside;

    P4[0] = x1;
    P4[1] = y1;

    P0[0] = P4[0] - this->diff_x;
    P0[1] = P4[1] - this->diff_y;

    P1[0] = P4[0] + this->diff_x;
    P1[1] = P0[1];

    P2[0] = P1[0];
    P2[1] = P4[1] + this->diff_y;

    P3[0] = P0[0];
    P3[1] = P2[1];

/*
    P4[0] = x1;
    P4[1] = y1;

    P0[0] = P4[0];
    P0[1] = y1 - this->diff_y;

    P1[0] = x1 + this->diff_x;
    P1[1] = P4[1];

    P2[0] = P4[0];
    P2[1] = y1 + this->diff_y;

    P3[0] = x1 - this->diff_x;
    P3[1] = P4[1];
*/
    P0[2] = get_imagevaluefromoverlap1(overlap_values,&P0_inside,(int)P0[0],(int)P0[1]);
    P1[2] = get_imagevaluefromoverlap1(overlap_values,&P1_inside,(int)P1[0],(int)P1[1]);
    P2[2] = get_imagevaluefromoverlap1(overlap_values,&P2_inside,(int)P2[0],(int)P2[1]);
    P3[2] = get_imagevaluefromoverlap1(overlap_values,&P3_inside,(int)P3[0],(int)P3[1]);
    P4[2] = get_imagevaluefromoverlap1(overlap_values,&P4_inside,(int)P4[0],(int)P4[1]);

    if(!P4_inside){
        //TODO: throw error
        coutMaster << "ERROR in FEM2DHAT - center point is not inside the domain" << std::endl;
    }

    /* deal with corners */
    if(!P0_inside) P0[2] = P4[2];
    if(!P1_inside) P1[2] = P4[2];
    if(!P2_inside) P2[2] = P4[2];
    if(!P3_inside) P3[2] = P4[2];

}


}
} /* end of namespace */

#endif
