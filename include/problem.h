#ifndef PROBLEM_H
#define	PROBLEM_H

#include "common.h"
#include "data.h"
#include "model.h"

class Problem {
	protected:
		Timer timer_total; /* from init to finalize */
		
		Data data;
		Model model;
		
	public:
		void init();
		void finalize();
	
		void set_data(Data new_data);
		void set_model(Model new_model);

		Data get_data(); // TODO: only for testing
		Model get_model(); // TODO: only for testing

};


#endif
