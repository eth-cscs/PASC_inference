#ifndef PROBLEM_H
#define	PROBLEM_H

#include "common.h"
#include "data.h"

void my_mvnrnd_D2(Scalar *mu, Scalar *covariance, Scalar *value1, Scalar *value2);

void get_problem_value1(Scalar *value1_out, Scalar *value2_out);
void get_problem_value2(Scalar *value1_out, Scalar *value2_out);
void get_problem_value3(Scalar *value1_out, Scalar *value2_out);

void generate_problem(Data *data_out, int dataT);

#endif
