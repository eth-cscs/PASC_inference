#ifndef PROBLEM_H
#define	PROBLEM_H

#include "common.h"
#include "data.h"
#include "model.h"

#include "inputoutput.h"


class Problem {
	protected:
		int it; /* outer iterations */
	
		Timer timer_total; /* from init to finalize */
		
		Data data;
		Model model;
		
	public:
		void init();
		void finalize();
	
		void set_data(Data new_data);
		void set_model(Model new_model);

		void solve(int max_s_steps, Scalar deltaL_eps);
		void print();
		void print_timers();
		void saveVTK(std::string name_of_file);

		Data get_data(); // TODO: only for testing
		Model get_model(); // TODO: only for testing

};


#endif
