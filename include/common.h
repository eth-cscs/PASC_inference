/*! This class should include all common functions and includes */

#ifndef PASC_COMMON_H
#define	PASC_COMMON_H

/* include common c++ header files */
#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <string>
#include <sstream>
#include <stack>
#include <limits>

/* include MINLIN */
#include <minlin/minlin.h>
#include <minlin/modules/threx/threx.h>
//#include <qpopt/smalbe.h>

using namespace minlin::threx;
MINLIN_INIT

/* include settings */
#include "settings.h"

/* for global time management */
class Timer;
extern Timer timer;

/* general utils */
void Initialize(int, char**);
void Finalize();

void Message(std::string text);
void Message_info(std::string text);
void Message_info_values(std::string text, std::string values);
void Message_info_value(std::string text, int value);
void Message_info_value(std::string text, double value);
void Message_info_time(std::string text, double value);
//void Message_error(string text);


/* structure for time management (measure computing time) */
class Timer {
		std::stack<double> time_stack;
		double time;

		double getUnixTime(void);
	public:
		void start();
		double stop();
	
};


#endif
