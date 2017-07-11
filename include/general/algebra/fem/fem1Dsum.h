/** @file fem1Dsum.h
 *  @brief class for reduction and prolongation on fem meshes
 *
 *  @author Lukas Pospisil
 */
 
#ifndef PASC_FEM1DSUM_H
#define	PASC_FEM1DSUM_H

#include "general/algebra/fem/fem.h"

namespace pascinference {
using namespace common;	
	
namespace algebra {

/** \class Fem1DSum
 *  \brief Reduction/prolongation between FEM meshes using averaging.
 *
*/
template<class VectorBase>
class Fem1DSum : public Fem<VectorBase> {
	public:
		class ExternalContent;

	protected:
		friend class ExternalContent;
		ExternalContent *externalcontent;			/**< for manipulation with external-specific stuff */

		double diff;

	public:
		/** @brief create FEM mapping between two decompositions
		*/
		Fem1DSum(Decomposition<VectorBase> *decomposition1, Decomposition<VectorBase> *decomposition2, double fem_reduce);

		/** @brief create general FEM mapping
		 * 
		 * do not forget to call Fem1DSum::set_decomposition_original() and afterwards Fem1DSum::compute_decomposition_reduced() to compute decomposition2 internally
		 * 
		 */
		Fem1DSum(double fem_reduce = 1.0);

		/** @brief destructor
		*/
		~Fem1DSum();

		/** @brief print info about fem
		 * 
		 * @param output where to print
		 */	
		virtual void print(ConsoleOutput &output_global, ConsoleOutput &output_local) const;
		
		virtual void reduce_gamma(GeneralVector<VectorBase> *gamma1, GeneralVector<VectorBase> *gamma2) const;
		virtual void prolongate_gamma(GeneralVector<VectorBase> *gamma2, GeneralVector<VectorBase> *gamma1) const;

		double get_diff() const;
		virtual std::string get_name() const;

		virtual void compute_decomposition_reduced();

		ExternalContent *get_externalcontent() const;		

};


/* ----------------- implementation ------------- */

template<class VectorBase>
Fem1DSum<VectorBase>::Fem1DSum(double fem_reduce) : Fem<VectorBase>(fem_reduce) {
	LOG_FUNC_BEGIN

	this->diff = 0; /* I dont have this information without decompositions */
	
	LOG_FUNC_END
}

template<class VectorBase>
Fem1DSum<VectorBase>::Fem1DSum(Decomposition<VectorBase> *decomposition1, Decomposition<VectorBase> *decomposition2, double fem_reduce) : Fem<VectorBase>(decomposition1, decomposition2, fem_reduce) {
	LOG_FUNC_BEGIN

	diff = (decomposition1->get_T())/(double)(decomposition2->get_T());

	LOG_FUNC_END
}

template<class VectorBase>
Fem1DSum<VectorBase>::~Fem1DSum(){
	LOG_FUNC_BEGIN

	LOG_FUNC_END
}

template<class VectorBase>
std::string Fem1DSum<VectorBase>::get_name() const {
	return "FEM-1D-SUM";
}

template<class VectorBase>
void Fem1DSum<VectorBase>::print(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN

	output_global << this->get_name() << std::endl;
	
	/* information of reduced problem */
	output_global <<  " - is reduced       : " << printbool(this->is_reduced()) << std::endl;
	output_global <<  " - diff             : " << diff << std::endl;
	output_global <<  " - fem_reduce       : " << this->get_fem_reduce() << std::endl;
	output_global <<  " - fem_type         : " << get_name() << std::endl;
	
	if(this->get_decomposition_original() == NULL){
		output_global <<  " - decomposition1   : NO" << std::endl;
	} else {
		output_global <<  " - decomposition1   : YES" << std::endl;
		output_global.push();
		this->get_decomposition_original()->print(output_global);
		output_global.pop();
	}

	if(this->get_decomposition_reduced() == NULL){
		output_global <<  " - decomposition2   : NO" << std::endl;
	} else {
		output_global <<  " - decomposition2   : YES" << std::endl;
		output_global.push();
		this->get_decomposition_reduced()->print(output_global);
		output_global.pop();
	}
	
	output_global.synchronize();	

	LOG_FUNC_END
}

template<class VectorBase>
void Fem1DSum<VectorBase>::reduce_gamma(GeneralVector<VectorBase> *gamma1, GeneralVector<VectorBase> *gamma2) const {
	LOG_FUNC_BEGIN

	LOG_FUNC_END
}

template<class VectorBase>
void Fem1DSum<VectorBase>::prolongate_gamma(GeneralVector<VectorBase> *gamma2, GeneralVector<VectorBase> *gamma1) const {
	LOG_FUNC_BEGIN

	LOG_FUNC_END
}

template<class VectorBase>
double Fem1DSum<VectorBase>::get_diff() const {
	return diff;
}

template<class VectorBase>
void Fem1DSum<VectorBase>::compute_decomposition_reduced() {
	LOG_FUNC_BEGIN
	
	if(this->is_reduced()){
		int T_reduced = ceil(this->decomposition1->get_T()*this->fem_reduce);
		
		/* compute new decomposition */
		this->decomposition2 = new Decomposition<VectorBase>(T_reduced, 
				*(this->decomposition1->get_graph()), 
				this->decomposition1->get_K(), 
				this->decomposition1->get_xdim(), 
				this->decomposition1->get_DDT_size(), 
				this->decomposition1->get_DDR_size());

	} else {
		/* there is not reduction of the data, we can reuse the decomposition */
		this->set_decomposition_reduced(this->decomposition1);
	}

	diff = (this->decomposition1->get_T())/(double)(this->decomposition2->get_T());
	
	LOG_FUNC_END
}


}
} /* end of namespace */


#endif
