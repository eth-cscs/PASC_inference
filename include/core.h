#ifndef CORE_H
#define	CORE_H

#include "common.h"

class Core {
//  class Construction;
  class Datatypes;
//  class Gamma;	
  
//  class ModelFit {
	  
//  };
};

class Core::Datatypes {
	class CGamma {
		double discrete;
		double T;
		double K;
		double Ct;
		double nbins;
		Vec gamma;

//		femModelParameter

		public:
			CGamma(Vec, double, double, double, double); /*! constructor */
	};
};

#endif
