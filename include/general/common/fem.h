/** @file fem.h
 *  @brief class for reduction and prolongation on fem meshes
 *
 *  @author Lukas Pospisil
 */
 
#ifndef PASC_FEM_H
#define	PASC_FEM_H

#include "general/common/decomposition.h"
#include "general/common/initialize.h"

namespace pascinference {
using namespace algebra;	
	
namespace common {

/** \class Fem
 *  \brief Reduction/prolongation between FEM meshes using constant functions.
 *
*/
template<class VectorBase>
class Fem {
	public:
		class ExternalContent;

	protected:
		friend class ExternalContent;
		ExternalContent *externalcontent;			/**< for manipulation with external-specific stuff */
	
		Decomposition<VectorBase> *decomposition1; /**< decomposition of the larger problem */
		Decomposition<VectorBase> *decomposition2; /**< decomposition of smaller problem */
		
		double diff;
		
		double fem_reduce;
	public:
		/** @brief create FEM mapping between two decompositions
		*/
		Fem(Decomposition<VectorBase> *decomposition1, Decomposition<VectorBase> *decomposition2, double fem_reduce);

		/** @brief create general FEM mapping
		 * 
		 * do not forget to call Fem::set_decomposition_original() and afterwards Fem::compute_decomposition_reduced() to compute decomposition2 internally
		 * 
		 */
		Fem(double fem_reduce = 1.0);

		/** @brief destructor
		*/
		~Fem();

		/** @brief print info about fem
		 * 
		 * @param output where to print
		 */	
		virtual void print(ConsoleOutput &output_global, ConsoleOutput &output_local) const;
		
		virtual void reduce_gamma(GeneralVector<VectorBase> *gamma1, GeneralVector<VectorBase> *gamma2) const;
		virtual void prolongate_gamma(GeneralVector<VectorBase> *gamma2, GeneralVector<VectorBase> *gamma1) const;

		//TODO: cannot me used in general FEM!
		double get_diff() const;
		double get_fem_reduce() const;
		virtual std::string get_name() const;

		void set_decomposition_original(Decomposition<VectorBase> *decomposition1);
		void set_decomposition_reduced(Decomposition<VectorBase> *decomposition2);

		virtual void compute_decomposition_reduced();

		Decomposition<VectorBase>* get_decomposition_original() const;
		Decomposition<VectorBase>* get_decomposition_reduced() const;

		bool is_reduced() const;
		
		ExternalContent *get_externalcontent() const;		
};


/* ----------------- implementation ------------- */

template<class VectorBase>
Fem<VectorBase>::Fem(double fem_reduce){
	LOG_FUNC_BEGIN

	decomposition1 = NULL;
	decomposition2 = NULL;

	this->fem_reduce = fem_reduce;
	this->diff = 0; /* I dont have this information without decompositions */
	
	LOG_FUNC_END
}

template<class VectorBase>
Fem<VectorBase>::Fem(Decomposition<VectorBase> *decomposition1, Decomposition<VectorBase> *decomposition2, double fem_reduce){
	LOG_FUNC_BEGIN

	this->fem_reduce = fem_reduce;

	this->set_decomposition_original(decomposition1);
	this->set_decomposition_reduced(decomposition2);

	diff = (decomposition1->get_T())/(double)(decomposition2->get_T());

	LOG_FUNC_END
}

template<class VectorBase>
Fem<VectorBase>::~Fem(){
	LOG_FUNC_BEGIN

	LOG_FUNC_END
}

template<class VectorBase>
std::string Fem<VectorBase>::get_name() const {
	return "FEM-SUM";
}

template<class VectorBase>
void Fem<VectorBase>::print(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN

	output_global << this->get_name() << std::endl;
	
	/* information of reduced problem */
	output_global <<  " - is reduced       : " << printbool(is_reduced()) << std::endl;
	output_global <<  " - diff             : " << diff << std::endl;
	output_global <<  " - fem_reduce       : " << fem_reduce << std::endl;
	output_global <<  " - fem_type         : " << get_name() << std::endl;
	
	if(decomposition1 == NULL){
		output_global <<  " - decomposition1   : NO" << std::endl;
	} else {
		output_global <<  " - decomposition1   : YES" << std::endl;
		output_global.push();
		decomposition1->print(output_global);
		output_global.pop();
	}

	if(decomposition2 == NULL){
		output_global <<  " - decomposition2   : NO" << std::endl;
	} else {
		output_global <<  " - decomposition2   : YES" << std::endl;
		output_global.push();
		decomposition2->print(output_global);
		output_global.pop();
	}
	
	output_global.synchronize();	

	LOG_FUNC_END
}

template<class VectorBase>
void Fem<VectorBase>::reduce_gamma(GeneralVector<VectorBase> *gamma1, GeneralVector<VectorBase> *gamma2) const {
	LOG_FUNC_BEGIN

	LOG_FUNC_END
}

template<class VectorBase>
void Fem<VectorBase>::prolongate_gamma(GeneralVector<VectorBase> *gamma2, GeneralVector<VectorBase> *gamma1) const {
	LOG_FUNC_BEGIN

	LOG_FUNC_END
}

template<class VectorBase>
double Fem<VectorBase>::get_diff() const {
	return diff;
}

template<class VectorBase>
double Fem<VectorBase>::get_fem_reduce() const {
	return fem_reduce;
}

template<class VectorBase>
void Fem<VectorBase>::set_decomposition_original(Decomposition<VectorBase> *decomposition1) {
	LOG_FUNC_BEGIN

	this->decomposition1 = decomposition1;

	LOG_FUNC_END
}

template<class VectorBase>
void Fem<VectorBase>::set_decomposition_reduced(Decomposition<VectorBase> *decomposition2) {
	LOG_FUNC_BEGIN

	this->decomposition2 = decomposition2;

	LOG_FUNC_END
}

template<class VectorBase>
Decomposition<VectorBase>* Fem<VectorBase>::get_decomposition_original() const {
	return decomposition1;
}

template<class VectorBase>
Decomposition<VectorBase>* Fem<VectorBase>::get_decomposition_reduced() const {
	return decomposition2;
}

template<class VectorBase>
void Fem<VectorBase>::compute_decomposition_reduced() {
	LOG_FUNC_BEGIN
	
	if(is_reduced()){
		int T_reduced = ceil(decomposition1->get_T()*fem_reduce);
		
		/* compute new decomposition */
		decomposition2 = new Decomposition<VectorBase>(T_reduced, 
				*(decomposition1->get_graph()), 
				decomposition1->get_K(), 
				decomposition1->get_xdim(), 
				decomposition1->get_DDT_size(), 
				decomposition1->get_DDR_size());

	} else {
		/* there is not reduction of the data, we can reuse the decomposition */
		this->set_decomposition_reduced(decomposition1);
	}

	diff = (decomposition1->get_T())/(double)(decomposition2->get_T());
	
	LOG_FUNC_END
}

template<class VectorBase>
bool Fem<VectorBase>::is_reduced() const {
	bool return_value;

	if(fem_reduce < 1.0) {
		return_value = true;
	} else {
		return_value = false;
	}
	
	return return_value;
}


}
} /* end of namespace */


#endif
