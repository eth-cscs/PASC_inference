/*! This class should include all common functions and includes */

#ifndef COMMON_H
#define	COMMON_H

#define datan 2
#define gammaK 3
#define max_s_steps 1000
#define deltaL_eps 0.0001
#define deltaL_eps_exact 0.01

#define DUALIZE 0

/* include common c++ header files */
#include <iostream>

/* include MINLIN */
#include <minlin/minlin.h>
#include <minlin/modules/threx/threx.h>
#include <qpopt/smalbe.h>

using namespace minlin::threx;

MINLIN_INIT

/* general utils */
void Initialize(int, char**);
void Finalize();

void Message(string text);
//void Message_error(string text);

#endif
