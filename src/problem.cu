#include "problem.h"

void get_problem_value1(Scalar *value1_out, Scalar *value2_out)
{
	Scalar covariance[4] = {0.001, 0.0, 0.0, 0.1};
	Scalar mu[2] = {0.25, 0};

	my_mvnrnd_D2(mu, covariance, value1_out, value2_out);
}

void get_problem_value2(Scalar *value1_out, Scalar *value2_out)
{
	Scalar covariance[4] = {0.005, 0.0, 0.0, 0.05};
	Scalar mu[2] = {0.0, -0.5};

	my_mvnrnd_D2(mu, covariance, value1_out, value2_out);
}

void get_problem_value3(Scalar *value1_out, Scalar *value2_out)
{
	Scalar covariance[4] = {0.005, 0.0, 0.0, 0.05};
	Scalar mu[2] = {0.0, 0.5};

	my_mvnrnd_D2(mu, covariance, value1_out, value2_out);
}

void generate_problem(Data *data_out, int dataT)
{
	int t;
	Data data;
	Scalar random_value1, random_value2; 
	
	/* prepare the space for data */
	data.init(datan,dataT);

	/* generate random data */
	// TODO: this operation could be performed fully independently, no necessary to use for
	for(t=0;t<data.get_T();t++){ /* go through time serie */
		if(t >= 0 && t < data.get_T()/3.0){ 
			get_problem_value1(&random_value1, &random_value2);
		}
		if(t >= data.get_T()/3.0 && t < 2.0*data.get_T()/3.0){ 
			get_problem_value2(&random_value1, &random_value2);
		}
		if(t >= 2.0*data.get_T()/3.0 && t <= data.get_T()){ 
			get_problem_value3(&random_value1, &random_value2);
		}
		
		data.data_vec(t) = random_value1;
		data.data_vec(data.get_T()+t) = random_value2;
	}

	*data_out = data;
}

void my_mvnrnd_D2(Scalar *mu, Scalar *covariance, Scalar *value1, Scalar *value2)
{
	Scalar L[4];
	Scalar r1, r2, r1n, r2n; 
	
	Scalar R, c, s;

    /* Compute normally distributed random numbers via Box-Muller */
	r1 = rand()/(double)(RAND_MAX);
    r1   = 1.-r1; /* to change from [0,1) to (0,1], which we need for the log */

    r2 = rand()/(double)(RAND_MAX);
    R = sqrt(-2.*log(r1));
    c = cos(2.*M_PI*r2);
    s = sin(2.*M_PI*r2);

    /* compute normal distributed random values */
    r1n = R*c;
    r2n = R*s;
	
	/* choleski decomposition of SPD covariance matrix */
	L[0] = sqrt(covariance[0]);
	L[1] = 0;
	L[2] = covariance[1]/L[0];
	L[3] = sqrt(covariance[3] - L[2]*L[2]);

	/* compute output values */
	/* y = L*randn(2,1) + mean */
	*value1 = L[0]*r1n + L[1]*r2n + mu[0];
	*value2 = L[2]*r1n + L[3]*r2n + mu[1];
}

