#ifndef PASC_COMMON_MVNRND_H
#define	PASC_COMMON_MVNRND_H

namespace pascinference {
	
void my_mvnrnd_D2(double mu1, double mu2, double diag_covariance1, double diag_covariance2, double *value1, double *value2){
	double r1, r2, r1n, r2n; 
	
	double R, c, s;

	/* Compute normally distributed random numbers via Box-Muller */
	r1 = rand()/(double)(RAND_MAX);
	r1 = 1.-r1; /* to change from [0,1) to (0,1], which we need for the log */

	r2 = rand()/(double)(RAND_MAX);
	R = sqrt(-2.*log(r1));
	c = cos(2.*M_PI*r2);
	s = sin(2.*M_PI*r2);

	/* compute normal distributed random values */
	r1n = R*c;
	r2n = R*s;
	
	/* compute output values */
	/* y = L*randn(2,1) + mean */
	*value1 = sqrt(diag_covariance1)*r1n + mu1;
	*value2 = sqrt(diag_covariance2)*r2n + mu2;
}	

void my_mvnrnd_Dn(int n, double *mu, double *diag_covariance, double *value){
	int i;
	for(i=0;i<n/2.0;i++){
		my_mvnrnd_D2(mu[2*i], mu[2*i+1], diag_covariance[2*i], diag_covariance[2*i+1], &(value[2*i]), &(value[2*i+1]));
	}

	if(n%2==1){
		double trash;
		my_mvnrnd_D2(mu[n-1], 0, diag_covariance[n-1], 1.0, &(value[2*i]), &trash);
	}
}	




} /* end of namespace */

#endif
