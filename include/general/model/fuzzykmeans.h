#ifndef PASC_FUZZYKMEANSMODEL_H
#define PASC_FUZZYKMEANSMODEL_H

#include "general/common/common.h"
#include "general/algebra/fem/fem.h"

/* gamma, theta problem */
#include "general/solver/simplesolver.h"
#include "general/data/simpledata.h"
#include "general/data/tsdata.h"

namespace pascinference {
namespace model {

/** \class FuzzyKmeansModel
 *  \brief time-series model with quadratic penalty time-space regularisation.
 *
*/
template<class VectorBase>
class FuzzyKmeansModel: public TSModel<VectorBase> {
	protected:
	 	SimpleData<VectorBase> *gammadata; /**< this problem is solved during assembly  */
	 	SimpleData<VectorBase> *thetadata; /**< this problem is solved during assembly  */

		Fem<VectorBase> *fem;	/** instance of FEM used for reduction/prolongation */

		GeneralVector<VectorBase> *residuum; /**< temp vector for residuum computation */

		double fuzzifier; /**< fuzzy coeficient */
	public:

		/** @brief constructor from data
		 *
		 * @param tsdata time-series data on which model operates
		 */
		FuzzyKmeansModel(TSData<VectorBase> &tsdata, double fuzzifier, Fem<VectorBase> *new_fem = NULL);

		/** @brief destructor
		 */
		~FuzzyKmeansModel();

		/** @brief print info about model
		 *
		 * @param output where to print
		 */
		void print(ConsoleOutput &output) const;

		/** @brief print info about model
		 *
		 * @param output_global where to print global info
		 * @param output_local where to print local info
		 */
		void print(ConsoleOutput &output_global, ConsoleOutput &output_local) const;

		/** @brief print solution of the model
		 *
		 * @param output_global where to print global part
		 * @param output_local where to print local part
		 */
		void printsolution(ConsoleOutput &output_global, ConsoleOutput &output_local) const;

		std::string get_name() const;

		void gammasolver_initialize(GeneralSolver **gamma_solver);
		void gammasolver_updatebeforesolve(GeneralSolver *gamma_solver);
		void gammasolver_updateaftersolve(GeneralSolver *gamma_solver);
		void gammasolver_finalize(GeneralSolver **gamma_solver);

		void thetasolver_initialize(GeneralSolver **theta_solver);
		void thetasolver_updatebeforesolve(GeneralSolver *theta_solver);
		void thetasolver_updateaftersolve(GeneralSolver *theta_solver);
		void thetasolver_finalize(GeneralSolver **theta_solver);

		double get_L(GeneralSolver *gammasolver, GeneralSolver *thetasolver);

		SimpleData<VectorBase> *get_gammadata() const;
		SimpleData<VectorBase> *get_thetadata() const;
		BGMGraph<VectorBase> *get_graph() const;

		double get_aic(double L) const;

		int get_T_reduced() const;
		int get_T() const;

		Decomposition<VectorBase> *get_decomposition_reduced() const;
		double get_fem_reduce() const;

};


}
} /* end of namespace */


/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {
namespace model {

/* constructor */
template<class VectorBase>
FuzzyKmeansModel<VectorBase>::FuzzyKmeansModel(TSData<VectorBase> &tsdata, double fuzzifier, Fem<VectorBase> *new_fem) {
	LOG_FUNC_BEGIN

	//TODO

	LOG_FUNC_END
}

/* destructor */
template<class VectorBase>
FuzzyKmeansModel<VectorBase>::~FuzzyKmeansModel(){
	LOG_FUNC_BEGIN

	/* destroy auxiliary vectors */
//	delete this->decomposition_reduced;

	LOG_FUNC_END
}


/* print info about model */
template<class VectorBase>
void FuzzyKmeansModel<VectorBase>::print(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output << this->get_name() << std::endl;

	/* give information about presence of the data */
	output <<  " - T                 : " << this->tsdata->get_T() << std::endl;
	output <<  " - xdim              : " << this->tsdata->get_xdim() << std::endl;
	output <<  " - fuzzifier         : " << this->fuzzifier << std::endl;

	/* information of reduced problem */
	output <<  " - fem_reduce        : " << this->fem->get_fem_reduce() << std::endl;
	output.push();
	this->fem->print(output,output);
	output.pop();

	output <<  " - K                 : " << this->tsdata->get_K() << std::endl;
	output <<  " - R                 : " << this->tsdata->get_R() << std::endl;

	output <<  " - Graph             : " << std::endl;
	output.push();
	this->get_graph()->print(output);
	output.pop();

	output <<  " - thetalength: " << this->thetavectorlength_global << std::endl;

	output.synchronize();

	LOG_FUNC_END
}

/* print info about model */
template<class VectorBase>
void FuzzyKmeansModel<VectorBase>::print(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN

	output_global <<  this->get_name() << std::endl;

	/* give information about presence of the data */
	output_global <<  " - global info" << std::endl;
	output_global <<  "  - T                 : " << this->tsdata->get_T() << std::endl;
	output_global <<  "  - xdim              : " << this->tsdata->get_xdim() << std::endl;
	output_global <<  "  - fuzzifier         : " << this->fuzzifier << std::endl;

	/* information of reduced problem */
	output_global <<  "  - fem_reduce        : " << this->fem->get_fem_reduce() << std::endl;
	output_global.push();
	this->fem->print(output_global, output_local);
	output_global.pop();

	output_global <<  "  - K                 : " << this->tsdata->get_K() << std::endl;
	output_global <<  "  - R                 : " << this->tsdata->get_R() << std::endl;

	output_global.push();
	this->get_graph()->print(output_global);
	output_global.pop();

	output_global <<  "  - thetalength       : " << this->thetavectorlength_global << std::endl;

	/* give local info */
	output_global <<  " - local variables" << std::endl;
	output_global.push();
	output_local << "Tlocal =" << std::setw(6) << this->tsdata->get_Tlocal() << " (" << this->tsdata->get_Tbegin() << "," << this->tsdata->get_Tend() << "), ";
	output_local << "Rlocal =" << std::setw(6) << this->tsdata->get_Rlocal() << " (" << this->tsdata->get_Rbegin() << "," << this->tsdata->get_Rend() << "), ";
	output_local << "thetalength=" << std::setw(6) << this->thetavectorlength_local << std::endl;

	output_global.pop();
	output_local.synchronize();

	output_global.synchronize();

	LOG_FUNC_END
}

/* print model solution */
template<class VectorBase>
void FuzzyKmeansModel<VectorBase>::printsolution(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN

	output_global <<  this->get_name() << std::endl;

	//TODO

	LOG_FUNC_END
}

/* get name of the model */
template<class VectorBase>
std::string FuzzyKmeansModel<VectorBase>::get_name() const {
	return "FuzzyKmeans Model";
}

/* prepare gamma solver */
template<class VectorBase>
void FuzzyKmeansModel<VectorBase>::gammasolver_initialize(GeneralSolver **gammasolver){
	LOG_FUNC_BEGIN

	//TODO

	LOG_FUNC_END
}

/* prepare theta solver */
template<class VectorBase>
void FuzzyKmeansModel<VectorBase>::thetasolver_initialize(GeneralSolver **thetasolver){
	LOG_FUNC_BEGIN

	//TODO

	LOG_FUNC_END
}

/* destroy gamma solver */
template<class VectorBase>
void FuzzyKmeansModel<VectorBase>::gammasolver_finalize(GeneralSolver **gammasolver){
	LOG_FUNC_BEGIN

	/* I created this objects, I should destroy them */
	if(fem->is_reduced()){
		free(gammadata->get_x());
	}

	/* destroy data */
	free(gammadata);

	/* destroy solver */
	free(*gammasolver);

	LOG_FUNC_END
}

/* destroy theta solver */
template<class VectorBase>
void FuzzyKmeansModel<VectorBase>::thetasolver_finalize(GeneralSolver **thetasolver){
	LOG_FUNC_BEGIN

	/* I created this objects, I should destroy them */

	/* destroy data */
	free(thetadata);

	/* destroy solver */
	free(*thetasolver);

	LOG_FUNC_END
}

template<class VectorBase>
double FuzzyKmeansModel<VectorBase>::get_L(GeneralSolver *gammasolver, GeneralSolver *thetasolver){

	// TODO

	return -1;
}

template<class VectorBase>
SimpleData<VectorBase>* FuzzyKmeansModel<VectorBase>::get_gammadata() const {
	return gammadata;
}

template<class VectorBase>
SimpleData<VectorBase>* FuzzyKmeansModel<VectorBase>::get_thetadata() const {
	return thetadata;
}

template<class VectorBase>
BGMGraph<VectorBase> *FuzzyKmeansModel<VectorBase>::get_graph() const {
	return this->tsdata->get_decomposition()->get_graph();
}

template<class VectorBase>
void FuzzyKmeansModel<VectorBase>::gammasolver_updatebeforesolve(GeneralSolver *gammasolver){
	LOG_FUNC_BEGIN

	//TODO

	LOG_FUNC_END
}

template<class VectorBase>
void FuzzyKmeansModel<VectorBase>::gammasolver_updateaftersolve(GeneralSolver *gammasolver){
	LOG_FUNC_BEGIN

	//TODO

	LOG_FUNC_END
}


/* update theta solver */
template<class VectorBase>
void FuzzyKmeansModel<VectorBase>::thetasolver_updatebeforesolve(GeneralSolver *thetasolver){
	LOG_FUNC_BEGIN

	//TODO

	LOG_FUNC_END
}

template<class VectorBase>
void FuzzyKmeansModel<VectorBase>::thetasolver_updateaftersolve(GeneralSolver *thetasolver){
	LOG_FUNC_BEGIN

	//TODO

	LOG_FUNC_END
}


template<class VectorBase>
double FuzzyKmeansModel<VectorBase>::get_aic(double L) const{

	//TODO

	return 2*log(L) + this->tsdata->get_K();
}

template<class VectorBase>
int FuzzyKmeansModel<VectorBase>::get_T_reduced() const {
	return fem->get_decomposition_reduced()->get_T();
}

template<class VectorBase>
int FuzzyKmeansModel<VectorBase>::get_T() const {
	return this->tsdata->get_T();
}

template<class VectorBase>
double FuzzyKmeansModel<VectorBase>::get_fem_reduce() const {
	return this->fem_reduce;
}

template<class VectorBase>
Decomposition<VectorBase> *FuzzyKmeansModel<VectorBase>::get_decomposition_reduced() const {
	return this->fem->get_decomposition_reduced();
}



}
} /* end namespace */

#endif
