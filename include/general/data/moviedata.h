/** @file moviedata.h
 *  @brief class for manipulation with movie data
 *
 *  @author Lukas Pospisil
 */

#ifndef PASC_MOVIEDATA_H
#define	PASC_MOVIEDATA_H

#include <iostream>
#include "general/common/common.h"
#include "general/algebra/graph/bgmgraph.h"
#include "general/model/tsmodel.h"
#include "general/data/tsdata.h"

namespace pascinference {
namespace data {

template<class VectorBase>
class MovieData: public TSData<VectorBase> {
	protected:
		int width;		/**< width of image */
		int height;		/**< height of image */
	public:
		MovieData(Decomposition<VectorBase> &decomposition, int width, int height, std::string filename_data, int type = 1);
		MovieData(Decomposition<VectorBase> &decomposition, int width, int height);
		~MovieData();

		virtual void print(ConsoleOutput &output) const;
		virtual void print(ConsoleOutput &output_global, ConsoleOutput &output_local) const;

		virtual void printcontent(ConsoleOutput &output) const;
		virtual void printcontent(ConsoleOutput &output_global, ConsoleOutput &output_local) const;
		virtual std::string get_name() const;

		void saveMovie_datavector(std::string filename, int type = 1) const;
		void saveMovie_gammavector(std::string filename) const;
		void saveMovie_reconstructed(std::string filename, int type = 1) const;

		double compute_abserr_reconstructed(GeneralVector<VectorBase> &solution) const;

		int get_width() const;
		int get_height() const;
		int get_nvalues() const;

		static std::string get_type_name(int type);

};


}
} /* end of namespace */

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {
namespace data {

/* from filename */
template<class VectorBase>
MovieData<VectorBase>::MovieData(Decomposition<VectorBase> &new_decomposition, int width, int height, std::string filename_data, int type){
	LOG_FUNC_BEGIN

	//TODO

	LOG_FUNC_END
}

/* with blank datavector */
template<class VectorBase>
MovieData<VectorBase>::MovieData(Decomposition<VectorBase> &new_decomposition, int width, int height){
	LOG_FUNC_BEGIN

    //TODO

	LOG_FUNC_END
}

/* destructor */
template<class VectorBase>
MovieData<VectorBase>::~MovieData(){
	LOG_FUNC_BEGIN


	LOG_FUNC_END
}


/* print info about data */
template<class VectorBase>
void MovieData<VectorBase>::print(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output <<  this->get_name() << std::endl;

	/* give information about presence of the data */
	if(this->tsmodel){
		output <<  " - T:           " << this->get_T() << std::endl;
		output <<  " - xdim:        " << this->get_xdim() << std::endl;
		output <<  " - K:           " << this->get_K() << std::endl;
		output <<  " - model:       " << this->tsmodel->get_name() << std::endl;
	} else {
		output <<  " - model:       NO" << std::endl;
	}
	output <<  " - R:           " << this->get_R() << std::endl;

	output <<  " - datavector:  ";
	if(this->datavector){
		output << "YES (size: " << this->datavector->size() << ")" << std::endl;
	} else {
		output << "NO" << std::endl;
	}
	output <<   " - gammavector: ";
	if(this->gammavector){
		output << "YES (size: " << this->gammavector->size() << ")" << std::endl;
	} else {
		output << "NO" << std::endl;
	}
	output <<   " - thetavector: ";
	if(this->thetavector){
		output << "YES (size: " << this->thetavector->size() << ")" << std::endl;
	} else {
		output << "NO" << std::endl;
	}

	output.synchronize();

	LOG_FUNC_END
}

/* print info about data */
template<class VectorBase>
void MovieData<VectorBase>::print(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN

	output_global <<  this->get_name() << std::endl;

	/* give information about presence of the data */
	if(this->tsmodel){
		output_global <<  " - T:           " << this->get_T() << std::endl;
		output_local  <<  "  - Tlocal:     " << this->get_Tlocal() << std::endl;
		output_local.synchronize();

		output_global <<  " - xdim:        " << this->get_xdim() << std::endl;
		output_global <<  " - K:           " << this->get_K() << std::endl;

		output_global <<  " - model:       " << this->tsmodel->get_name() << std::endl;
	} else {
		output_global <<  " - model:       NO" << std::endl;
	}
	output_global <<  " - R:           " << this->get_R() << std::endl;

	output_global <<  " - datavector:  ";
	if(this->datavector){
		output_global << "YES (size: " << this->datavector->size() << ")" << std::endl;
		output_local  <<  "  - local size: " << this->datavector->local_size() << std::endl;
		output_local.synchronize();
	} else {
		output_global << "NO" << std::endl;
	}

	output_global <<   " - gammavector: ";
	if(this->gammavector){
		output_global << "YES (size: " << this->gammavector->size() << ")" << std::endl;
		output_local  <<  "  - local size: " << this->gammavector->local_size() << std::endl;
		output_local.synchronize();
	} else {
		output_global << "NO" << std::endl;
	}

	output_global <<   " - thetavector: ";
	if(this->thetavector){
		output_global << "YES (size: " << this->thetavector->size() << ")" << std::endl;
		output_local  <<  "  - local size: " << this->thetavector->local_size() << std::endl;
		output_local.synchronize();
	} else {
		output_global << "NO" << std::endl;
	}

	output_global.synchronize();

	LOG_FUNC_END
}

/* print content of all data */
template<class VectorBase>
void MovieData<VectorBase>::printcontent(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output <<  this->get_name() << std::endl;

	/* print the content of the data */
	output <<  " - datavector: ";
	if(this->datavector){
		output << *this->datavector << std::endl;
	} else {
		output << "not set" << std::endl;
	}

	output <<  " - gammavector: ";
	if(this->gammavector){
		output << *this->gammavector << std::endl;
	} else {
		output << "not set" << std::endl;
	}

	output <<  " - thetavector: ";
	if(this->thetavector){
		output << *this->thetavector << std::endl;
	} else {
		output << "not set" << std::endl;
	}

	LOG_FUNC_END
}

/* print content of all data */
template<class VectorBase>
void MovieData<VectorBase>::printcontent(ConsoleOutput &output_global,ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN

	output_global <<  this->get_name() << std::endl;

	/* print the content of the data */
	output_local <<  " - datavector: ";
	if(this->datavector){
		output_local << *this->datavector << std::endl;
	} else {
		output_local << "not set" << std::endl;
	}
	output_local.synchronize();

	output_local <<  " - gammavector: ";
	if(this->gammavector){
		output_local << *this->gammavector << std::endl;
	} else {
		output_local << "not set" << std::endl;
	}
	output_local.synchronize();

	output_local <<  " - thetavector: ";
	if(this->thetavector){
		output_local << *this->thetavector << std::endl;
	} else {
		output_local << "not set" << std::endl;
	}
	output_local.synchronize();

	output_global.synchronize();

	LOG_FUNC_END
}

template<class VectorBase>
std::string MovieData<VectorBase>::get_name() const {
	return "Movie Time-series Data";
}

template<class VectorBase>
void MovieData<VectorBase>::saveMovie_datavector(std::string filename, int type) const {
	LOG_FUNC_BEGIN

	//TODO

	LOG_FUNC_END
}

template<class VectorBase>
void MovieData<VectorBase>::saveMovie_gammavector(std::string filename) const {
	LOG_FUNC_BEGIN

	//TODO

	LOG_FUNC_END
}

template<class VectorBase>
void MovieData<VectorBase>::saveMovie_reconstructed(std::string filename, int type) const {
	LOG_FUNC_BEGIN

	//TODO

	LOG_FUNC_END
}

template<class VectorBase>
double MovieData<VectorBase>::compute_abserr_reconstructed(GeneralVector<VectorBase> &solution) const {
	LOG_FUNC_BEGIN

	//TODO

	LOG_FUNC_END

	return -1.0;
}

template<class VectorBase>
int MovieData<VectorBase>::get_width() const {
	return this->width;
}

template<class VectorBase>
int MovieData<VectorBase>::get_height() const {
	return this->height;
}

template<class VectorBase>
int MovieData<VectorBase>::get_nvalues() const {
	return this->get_xdim() * this->width * this->height;
}

template<class VectorBase>
std::string MovieData<VectorBase>::get_type_name(int type){
	std::ostringstream sout;

	if(type == 0) sout << "TRn";
	if(type == 1) sout << "TnR";
	if(type == 2) sout << "nTR";
	if(type < 0 | type > 2) sout << "error";

	return sout.str();
}


}
} /* end namespace */

#endif
