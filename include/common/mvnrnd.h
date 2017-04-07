#ifndef PASC_COMMON_MVNRND_H
#define	PASC_COMMON_MVNRND_H

#include <math.h>
#include <stdlib.h>

namespace pascinference {
namespace common {

void my_mvnrnd_D2(double mu1, double mu2, double diag_covariance1, double diag_covariance2, double *value1, double *value2);
void my_mvnrnd_Dn(int n, double *mu, double *diag_covariance, double *value);

}
} /* end of namespace */

#endif
