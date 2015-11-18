#ifndef PROBLEM_H
#define	PROBLEM_H

#include "common.h"
#include "data.h"


PetscErrorCode my_mvnrnd_D2(PetscScalar *mu, PetscScalar *covariance, PetscScalar *value1, PetscScalar *value2);

PetscErrorCode get_problem_value1(PetscScalar *value1_out, PetscScalar *value2_out);
PetscErrorCode get_problem_value2(PetscScalar *value1_out, PetscScalar *value2_out);
PetscErrorCode get_problem_value3(PetscScalar *value1_out, PetscScalar *value2_out);

PetscErrorCode generate_problem(Data *data_out, PetscInt dataN);

#endif
