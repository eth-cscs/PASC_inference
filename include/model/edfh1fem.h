#ifndef PASC_EDFH1FEMMODEL_H
#define	PASC_EDFH1FEMMODEL_H

#ifndef USE_PETSCVECTOR
 #error 'EDFH1FEMMODEL is for PETSCVECTOR'
#endif

typedef petscvector::PetscVector PetscVector;

/* for debugging, if >= 100, then print info about ach called function */
extern int DEBUG_MODE;

#include "pascinference.h"

/* gamma problem */
#include "matrix/blockgraph.h"

#include "feasibleset/simplex_local.h"
#include "solver/qpsolver.h"
#include "data/qpdata.h"

/* theta problem */
#include "solver/simplesolver.h"
#include "data/simpledata.h"

#include "data/edfdata.h"


namespace pascinference {

template<class VectorBase>
class EdfH1FEMModel: public TSModel<VectorBase> {
	protected:
		QPData<VectorBase> *gammadata; /**< QP with simplex, will be solved by SPG-QP */
	 	SimpleData<VectorBase> *thetadata; /**< this problem is solved during assembly  */

		/* model specific variables */
		BGM_Graph *graph; /**< graph with stucture of the matrix */
		int R; /**< number of nodes of the graph */
		double epssqr; /**< penalty coeficient */
		
	public:
	
		EdfH1FEMModel(EdfData<VectorBase> &tsdata, BGM_Graph &new_graph, int K, double epssqr);
		~EdfH1FEMModel();

		void print(ConsoleOutput &output) const;
		void print(ConsoleOutput &output_global, ConsoleOutput &output_local) const;
		void printsolution(ConsoleOutput &output_global, ConsoleOutput &output_local) const;
		
		std::string get_name() const;
		
		void initialize_gammasolver(GeneralSolver **gamma_solver, const TSData<VectorBase> *tsdata);
		void initialize_thetasolver(GeneralSolver **theta_solver, const TSData<VectorBase> *tsdata);
		
		void finalize_gammasolver(GeneralSolver **gamma_solver, const TSData<VectorBase> *tsdata);
		void finalize_thetasolver(GeneralSolver **theta_solver, const TSData<VectorBase> *tsdata);

		void update_gammasolver(GeneralSolver *gamma_solver, const TSData<VectorBase> *tsdata);
		void update_thetasolver(GeneralSolver *theta_solver, const TSData<VectorBase> *tsdata);
	
		double get_L(GeneralSolver *gammasolver, GeneralSolver *thetasolver, const TSData<VectorBase> *tsdata);

		QPData<VectorBase> *get_gammadata() const;
		SimpleData<VectorBase> *get_thetadata() const;

};

} // end of namespace


/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {

/* constructor */
template<>
EdfH1FEMModel<PetscVector>::EdfH1FEMModel(EdfData<PetscVector> &tsdata, BGM_Graph &new_graph, int K, double epssqr) {
	LOG_FUNC_BEGIN

	this->xdim = 1; // TODO: now I can compute only 1D problem, sorry
	this->R = tsdata.get_R();

	/* set graph */
	this->graph = &new_graph;	

	/* get original PETSc Vec from data vector */
	Vec data_Vec = tsdata.get_datavector()->get_vector();

	/* get T from data */
	TRY( VecGetSize(data_Vec, &(this->datavectorlength_global)) );
	this->T = this->datavectorlength_global/(double)(R*xdim);

	/* get Tlocal from data */
	TRY( VecGetLocalSize(data_Vec, &(this->datavectorlength_local)) );
	this->Tlocal = this->datavectorlength_local/(double)(R*xdim);

	/* get T ranges from data */
	int low, high;
	TRY( VecGetOwnershipRange(data_Vec, &low, &high) );
	this->Tbegin = low/(double)(R*xdim);
	this->Tend = high/(double)(R*xdim);

	/* set given parameters */
	this->K = K;
	this->epssqr = epssqr;

	int my_rank = GlobalManager.get_rank();

	/* compute vector lengths */

	/* prepage global vectors with gamma */
	this->gammavectorlength_local = K*R*this->Tlocal;
	this->gammavectorlength_global = K*R*this->T;

	/* prepare sequential vector with Theta - yes, all procesors will have the same information */
	this->thetavectorlength_local = K;
	this->thetavectorlength_global = GlobalManager.get_size()*(this->thetavectorlength_local);

	/* set this model to data */
	tsdata.set_model(*this);
	
	LOG_FUNC_END
}

/* destructor */
template<class VectorBase>
EdfH1FEMModel<VectorBase>::~EdfH1FEMModel(){
	LOG_FUNC_BEGIN
	
	/* destroy auxiliary vectors */
	
	LOG_FUNC_END
}


/* print info about model */
template<class VectorBase>
void EdfH1FEMModel<VectorBase>::print(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output <<  this->get_name() << std::endl;
	
	/* give information about presence of the data */
	output <<  " - T:       " << this->T << std::endl;
	output <<  " - xdim:    " << this->xdim << std::endl;

	output <<  " - K:       " << this->K << std::endl;
	output <<  " - R:       " << this->R << std::endl;
	output <<  " - epssqr:  " << this->epssqr << std::endl;

	output <<  " - Graph:       " << std::endl;
	output.push();
	graph->print(output);
	output.pop();

	output <<  " - datalength:  " << this->datavectorlength_global << std::endl;
	output <<  " - gammalength: " << this->gammavectorlength_global << std::endl;
	output <<  " - thetalength: " << this->thetavectorlength_global << std::endl;
	
	output.synchronize();	

	LOG_FUNC_END
}

/* print info about model */
template<class VectorBase>
void EdfH1FEMModel<VectorBase>::print(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN

	output_global <<  this->get_name() << std::endl;
	
	/* give information about presence of the data */
	output_global <<  " - global info:  " << std::endl;
	output_global <<  "  - T:       " << this->T << std::endl;
	output_global <<  "  - xdim:    " << this->xdim << std::endl;

	output_global <<  " - K:       " << this->K << std::endl;
	output_global <<  " - R:       " << this->R << std::endl;
	output_global <<  " - epssqr:  " << this->epssqr << std::endl;

	output_global.push();
	graph->print(output_global);
	output_global.pop();

	output_global <<  " - datalength:  " << this->datavectorlength_global << std::endl;
	output_global <<  " - gammalength: " << this->gammavectorlength_global << std::endl;
	output_global <<  " - thetalength: " << this->thetavectorlength_global << std::endl;

	/* give local info */
	output_global <<  " - local variables:  " << std::endl;
	output_global.push();
	output_local << "Tlocal =" << std::setw(6) << this->Tlocal << " (" << this->Tbegin << "," << this->Tend << "), ";
	output_local << "datalength=" << std::setw(6) << this->datavectorlength_local << ", ";
	output_local << "gammalength=" << std::setw(6) << this->gammavectorlength_local << ", ";
	output_local << "thetalength=" << std::setw(6) << this->thetavectorlength_local << std::endl;

	output_global.pop();
	output_local.synchronize();

	output_global.synchronize();

	LOG_FUNC_END
}

/* print model solution */
template<>
void EdfH1FEMModel<PetscVector>::printsolution(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN

	output_global <<  this->get_name() << std::endl;
	
	/* give information about presence of the data */
	output_global <<  "theta:" << std::endl;

	std::ostringstream temp;
	
	double *theta;
	TRY( VecGetArray(thetadata->get_x()->get_vector(),&theta) );
	
	int k,i,n;

	output_local.push();
	output_local << "- proc: " << GlobalManager.get_rank() << std::endl;
	output_local.push();
	for(k=0;k<this->K;k++){
		output_local <<  "- k = " << k << std::endl;

		/* mu */
		output_local.push();
		output_local <<  "- theta = [";
		for(n=0;n<xdim;n++){
			temp << theta[k*xdim + n];
			output_local << temp.str();
			if(n < xdim-1){
				output_local << ", ";
			}
			temp.str("");
		}
		output_local <<  "]" << std::endl;
		output_local.pop();

	}
	output_local.pop();
	output_local.pop();

	output_local.synchronize();
	output_global.synchronize();

	LOG_FUNC_END
}

/* get name of the model */
template<class VectorBase>
std::string EdfH1FEMModel<VectorBase>::get_name() const {
	return "EDF-H1-FEM Time-Series Model";	
}

/* prepare gamma solver */
template<>
void EdfH1FEMModel<PetscVector>::initialize_gammasolver(GeneralSolver **gammasolver, const TSData<PetscVector> *tsdata){
	LOG_FUNC_BEGIN

	/* in this case, gamma problem is QP with simplex feasible set */
	
	/* create data */
	gammadata = new QPData<PetscVector>();
	
	gammadata->set_x0(tsdata->get_gammavector()); /* the initial approximation of QP problem is gammavector */
	gammadata->set_x(tsdata->get_gammavector()); /* the solution of QP problem is gamma */
	gammadata->set_b(new GeneralVector<PetscVector>(*gammadata->get_x0())); /* create new linear term of QP problem */

	gammadata->set_A(new BlockGraphMatrix<PetscVector>(*(gammadata->get_x0()), *(this->graph), this->K, this->epssqr*this->epssqr)); 
	gammadata->set_feasibleset(new SimplexFeasibleSet_Local(this->Tlocal,this->K)); /* the feasible set of QP is simplex */ 	

	/* create solver */
	*gammasolver = new QPSolver<PetscVector>(*gammadata);

	/* generate random data to gamma */
	gammadata->get_x0()->set_random();

	/* project random values to feasible set to be sure that initial approximation is feasible */
	gammadata->get_feasibleset()->project(*gammadata->get_x0());

	LOG_FUNC_END
}

/* prepare theta solver */
template<>
void EdfH1FEMModel<PetscVector>::initialize_thetasolver(GeneralSolver **thetasolver, const TSData<PetscVector> *tsdata){
	LOG_FUNC_BEGIN
	
	/* create data */
	thetadata = new SimpleData<PetscVector>();
	thetadata->set_x(tsdata->get_thetavector()); /* the solution of QP problem is gamma */

	/* create solver */
	*thetasolver = new SimpleSolver<PetscVector>(*thetadata);
	
	LOG_FUNC_END
}

/* destroy gamma solver */
template<class VectorBase>
void EdfH1FEMModel<VectorBase>::finalize_gammasolver(GeneralSolver **gammasolver, const TSData<VectorBase> *tsdata){
	LOG_FUNC_BEGIN

	/* I created this objects, I should destroy them */

	/* destroy data */
	free(gammadata->get_b());
	free(gammadata->get_A());
	free(gammadata->get_feasibleset());
	free(gammadata);

	/* destroy solver */
	free(*gammasolver);
	
	LOG_FUNC_END
}

/* destroy theta solver */
template<class VectorBase>
void EdfH1FEMModel<VectorBase>::finalize_thetasolver(GeneralSolver **thetasolver, const TSData<VectorBase> *tsdata){
	LOG_FUNC_BEGIN

	/* I created this objects, I should destroy them */

	/* destroy data */
	free(thetadata);

	/* destroy solver */
	free(*thetasolver);

	LOG_FUNC_END
}

template<class VectorBase>
double EdfH1FEMModel<VectorBase>::get_L(GeneralSolver *gammasolver, GeneralSolver *thetasolver, const TSData<VectorBase> *tsdata){
	// TODO: not suitable in every situation - I suppose that g was computed from actual x,b */
	return ((QPSolver<VectorBase> *)gammasolver)->get_fx();
}

template<class VectorBase>
QPData<VectorBase>* EdfH1FEMModel<VectorBase>::get_gammadata() const {
	return gammadata;
}

template<class VectorBase>
SimpleData<VectorBase>* EdfH1FEMModel<VectorBase>::get_thetadata() const {
	return thetadata;
}

template<>
void EdfH1FEMModel<PetscVector>::update_gammasolver(GeneralSolver *gammasolver, const TSData<PetscVector> *tsdata){
	LOG_FUNC_BEGIN

	/* update gamma_solver data - prepare new linear term */
	const double *theta_arr;
	TRY( VecGetArrayRead(tsdata->get_thetavector()->get_vector(), &theta_arr) );
	
	const double *data_arr;
	TRY( VecGetArrayRead(tsdata->get_datavector()->get_vector(), &data_arr) );
	
	double *b_arr;
	TRY( VecGetArray(gammadata->get_b()->get_vector(), &b_arr) );

	int k,t,n;
	double mydot;
	for(t=0;t<Tlocal;t++){
		for(k=0;k<K;k++){
			mydot = 0;
			for(n=0;n<xdim;n++){
				mydot += (data_arr[t*xdim+n] - theta_arr[k*xdim+n])*(data_arr[t*xdim+n] - theta_arr[k*xdim+n]); 
			}
			b_arr[k*Tlocal+t] = -mydot;
		}
	}

	/* restore arrays */
	TRY( VecRestoreArray(gammadata->get_b()->get_vector(), &b_arr) );
	TRY( VecRestoreArrayRead(tsdata->get_datavector()->get_vector(), &data_arr) );
	TRY( VecRestoreArrayRead(tsdata->get_thetavector()->get_vector(), &theta_arr) );
	
	LOG_FUNC_END
}

/* update theta solver */
template<>
void EdfH1FEMModel<PetscVector>::update_thetasolver(GeneralSolver *thetasolver, const TSData<PetscVector> *tsdata){
	LOG_FUNC_BEGIN


	LOG_FUNC_END
}



} /* end namespace */

#endif
