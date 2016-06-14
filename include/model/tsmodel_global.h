/** @file tsmodel_global.h
 *  @brief class for manipulation with global time-series models
 *
 *  @author Lukas Pospisil
 */

#ifndef PASC_TSMODEL_GLOBAL_H
#define	PASC_TSMODEL_GLOBAL_H

#define DEFAULT_T_SCATTER 100

#ifndef USE_PETSCVECTOR
 #error 'TSMODEL_GLOBAL is for PETSCVECTOR'
#endif

/* for debugging, if >= 100, then print info about ach called function */
extern int DEBUG_MODE;

#include "pascinference.h"
#include "data/tsdata_global.h"

namespace pascinference {

/* Maybe these classes are not defined yet */ 
class TSData_Global;

/** \class TSModel_Global
 *  \brief General class for manipulation with global time-series models.
 *
*/
class TSModel_Global: public GeneralModel {
	protected:
		int T; /**< length of time-series, including xmem */
		int Tlocal; /**< local part of time-series */
		
		int xdim; /**< number of components in each time-step */
		int xmem;

		int K; /**< number of clusters */

		int datavectorlength_global; /**< global length of datavector (time-series values) */
		int gammavectorlength_global; /**< global length of gammavector (switching functions) */
		int thetavectorlength_global; /**< global length of thetavector (model parameters) */

		int datavectorlength_local; /**< local length of datavector (time-series values) */
		int gammavectorlength_local; /**< local length of gammavector (switching functions) */
		int thetavectorlength_local; /**< local length of thetavector (model parameters) */

		int t_scatter; /**< size of scatter > xmem, how long time-series to scatter to all processors */
		Vec x_scatter;
				
	public:
		TSModel_Global();
		TSModel_Global(int T, int xdim, int K);
		~TSModel_Global();

		virtual void print(ConsoleOutput &output) const;
		virtual std::string get_name() const;

		/* global length */
		virtual int get_datavectorlength_global();
		virtual int get_gammavectorlength_global();
		virtual int get_thetavectorlength_global();

		/* local length */
		virtual int get_datavectorlength_local();
		virtual int get_gammavectorlength_local();
		virtual int get_thetavectorlength_local();

		/* get common variables */
		int get_T() const;
		int get_Tlocal() const;
		int get_xdim() const;
		int get_xmem() const;

		int get_K() const;
		
		/** @brief alloc memory for gamma solver
		 *  
		 *  Allocate memory for all data of gamma problem.
		 * 
		 */ 
		virtual void initialize_gammasolver(GeneralSolver **gammasolver, const TSData_Global *tsdata){
				*gammasolver = NULL;
		};

		/** @brief alloc memory for Theta solver
		 *  
		 *  Allocate memory for all data of Theta problem.
		 * 
		 */ 
		virtual void initialize_thetasolver(GeneralSolver **thetasolver, const TSData_Global *tsdata){
				*thetasolver = NULL;
		};

		/** @brief free memory of gamma solver
		 *  
		 *  Deallocate memory for all data of gamma problem.
		 * 
		 */ 
		virtual void finalize_gammasolver(GeneralSolver **gammasolver, const TSData_Global *tsdata){
		};

		/** @brief free memory of Theta solver
		 *  
		 *  Deallocate memory for all data of Theta problem.
		 * 
		 */ 
		virtual void finalize_thetasolver(GeneralSolver **thetasolver, const TSData_Global *tsdata){
		};

		/** @brief update data values of gamma solver
		 *  
		 *  Update data values of gamma solver, prepare it to the solving.
		 * 
		 */ 
		virtual void update_gammasolver(GeneralSolver *gammasolver, const TSData_Global *tsdata){
		};

		/** @brief update data values of Theta solver
		 *  
		 *  Update data values of Theta solver, prepare it to the solving.
		 * 
		 */ 
		virtual void update_thetasolver(GeneralSolver *thetasolver, const TSData_Global *tsdata){
		};

		/** @brief get the value of object function L(gamma,Theta,data)
		 *  
		 */ 
		virtual double get_L(GeneralSolver *gammasolver, GeneralSolver *thetasolver, const TSData_Global *tsdata){
			return std::numeric_limits<double>::max();
		};

		int get_t_scatter(){
			return t_scatter;
		};

		Vec get_x_scatter(){
			return x_scatter;
		};
		
		void scatter_xn(Vec &x_global, Vec &xn_local, int xdim, int n);
		void scatter_part(Vec &x_global, int is_begin, int is_end, int xdim);					
		virtual void compute_x_model(double *x_model, const double *x_arr, int t_x_arr, const double *theta_arr, const double *gamma_arr, int t_gamma_arr) {
			// TODO: here write something funny
		};

};

} // end of namespace

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {

/* constructor */
TSModel_Global::TSModel_Global(){
	if(DEBUG_MODE >= 100) coutMaster << "(TSModel_Global)CONSTRUCTOR" << std::endl;

	/* set initial content */
	this->T = 0;
	this->Tlocal = 0;
	this->xdim = 0;
	this->xmem = 0;
	this->K = 0;

	/* set other variables */
	consoleArg.set_option_value("tsmodel_global_t_scatter", &this->t_scatter, DEFAULT_T_SCATTER); /* > xmem, how long time-series to scatter to all processors */
	TRY( VecCreateSeq(PETSC_COMM_SELF, t_scatter*xdim, &x_scatter) );
	
}


TSModel_Global::TSModel_Global(int newT, int newxdim, int newK){
	if(DEBUG_MODE >= 100) coutMaster << "(TSModel_Global)CONSTRUCTOR" << std::endl;

	/* set initial content */
	this->T = newT;
	this->xdim = newxdim;
	this->xmem = 0;
	this->K = newK;

	/* set other variables */
	consoleArg.set_option_value("tsmodel_global_t_scatter", &this->t_scatter, DEFAULT_T_SCATTER); /* > xmem, how long time-series to scatter to all processors */
	TRY( VecCreateSeq(PETSC_COMM_SELF, t_scatter*xdim, &x_scatter) );

}

/* destructor */
TSModel_Global::~TSModel_Global(){
	if(DEBUG_MODE >= 100) coutMaster << "(TSModel_Global)DESTRUCTOR" << std::endl;
	
}


/* print info about model */
void TSModel_Global::print(ConsoleOutput &output) const {
	output <<  this->get_name() << std::endl;
	
	/* give information about presence of the data */
	output <<  " - nproc:  " << GlobalManager.get_size() << std::endl;
	
	output <<  " - T:      " << this->T << std::endl;
	output <<  " - xdim:    " << this->xdim << std::endl;

	output <<  " - K: " << this->K << std::endl;
		
}

/* get name of the model */
std::string TSModel_Global::get_name() const {
	return "General Time-Series-Global Model";	
}

/* ---- GET FUNCTIONS ---- */
int TSModel_Global::get_T() const{
	return T; 
}

int TSModel_Global::get_Tlocal() const{
	return Tlocal; 
}

int TSModel_Global::get_xdim() const{
	return xdim; 
}

int TSModel_Global::get_xmem() const{
	return xmem; 
}

int TSModel_Global::get_K() const{
	return K; 
}


int TSModel_Global::get_datavectorlength_global(){
	return this->datavectorlength_global;
}

int TSModel_Global::get_gammavectorlength_global(){
	return this->gammavectorlength_global;
}

int TSModel_Global::get_thetavectorlength_global(){
	return this->thetavectorlength_global;
}

int TSModel_Global::get_datavectorlength_local(){
	return this->datavectorlength_local;
}

int TSModel_Global::get_gammavectorlength_local(){
	return this->gammavectorlength_local;
}

int TSModel_Global::get_thetavectorlength_local(){
	return this->thetavectorlength_local;
}

void TSModel_Global::scatter_part(Vec &x_global, int is_begin, int is_end, int xdim){
	LOG_FUNC_BEGIN

	VecScatter ctx; 
	IS scatter_is;
	IS scatter_is_to;

	/* scatter part of time serie */
	TRY( ISCreateStride(PETSC_COMM_WORLD, (is_end-is_begin)*xdim, is_begin*xdim, 1, &scatter_is) );
	TRY( ISCreateStride(PETSC_COMM_SELF, (is_end-is_begin)*xdim, 0, 1, &scatter_is_to) );

	TRY( VecScatterCreate(x_global,scatter_is, x_scatter,scatter_is_to,&ctx) );
	TRY( VecScatterBegin(ctx,x_global,x_scatter,INSERT_VALUES,SCATTER_FORWARD) );
	TRY( VecScatterEnd(ctx,x_global,x_scatter,INSERT_VALUES,SCATTER_FORWARD) );
	TRY( VecScatterDestroy(&ctx) );

	TRY( ISDestroy(&scatter_is_to) );
	TRY( ISDestroy(&scatter_is) );

	LOG_FUNC_END
}

void TSModel_Global::scatter_xn(Vec &x_global, Vec &xn_local, int xdim, int n){
	LOG_FUNC_BEGIN
	
	/* compute T */
	int x_global_size;
	TRY( VecGetSize(x_global,&x_global_size) );
	int T = x_global_size/(double)xdim;

	VecScatter ctx; 
	IS xn_local_is;
	IS scatter_is_to;

	TRY( ISCreateStride(PETSC_COMM_SELF, T, n, xdim, &xn_local_is) );
	TRY( ISCreateStride(PETSC_COMM_SELF, T, 0, 1, &scatter_is_to) );

	TRY( VecScatterCreate(x_global, xn_local_is, xn_local,scatter_is_to,&ctx) );
	TRY( VecScatterBegin(ctx,x_global, xn_local, INSERT_VALUES,SCATTER_FORWARD) );
	TRY( VecScatterEnd(ctx,x_global, xn_local, INSERT_VALUES,SCATTER_FORWARD) );
	TRY( VecScatterDestroy(&ctx) );

	TRY( ISDestroy(&xn_local_is) );
	TRY( ISDestroy(&scatter_is_to) );

	LOG_FUNC_END
}

} /* end namespace */

#endif
