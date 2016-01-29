#ifndef PROBLEM_H
#define	PROBLEM_H

#include "common.h"
#include "data.h"

void generate_problem(Data *data_out, int dataT);

class Problem {
	protected:
		Timer timer_total; /* from init to finalize */
		
		Data data;
		
	public:
		void init();
		void finalize();
	
		void set_data(Data new_data);
		Data get_data(); // TODO: only for testing
};


#endif
