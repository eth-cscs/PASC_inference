#ifndef PROBLEM_H
#define	PROBLEM_H

#include "common.h"
#include "data.h"


void my_mvnrnd_D2(PetscScalar *mu, PetscScalar *covariance, PetscScalar *value1, PetscScalar *value2);

void get_problem_value1(PetscScalar *value1_out, PetscScalar *value2_out);
void get_problem_value2(PetscScalar *value1_out, PetscScalar *value2_out);
void get_problem_value3(PetscScalar *value1_out, PetscScalar *value2_out);

void generate_problem(Data *data_out, PetscInt dataN);

#endif
