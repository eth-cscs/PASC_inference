/** @file test_integral.cu
 *  @brief test numerical integration
 *
 *  @author Lukas Pospisil
 */


#include <iostream>
#include "dlib/matrix.h"
#include "dlib/numeric_constants.h"
#include "dlib/numerical_integration.h"

#include <time.h>

using namespace dlib;

/* random number from 0,1 */
double rand01() { return (double)std::rand() / RAND_MAX; }

typedef matrix<double,0,1> column_vector;
#define STOP_TOLERANCE 1e-06

double gg(double y, int order, column_vector& LM)
{
    long  x_size = LM.size();
    long  num_moments = x_size;
    column_vector z(num_moments);
    
    z = 0;
    for (int i = 0; i < num_moments; ++i)
        z(i) = pow(y,i+1);
    
    return pow(y,order)*(exp(-trans(LM)*z));
}

double mymc(column_vector& LM, int order, int N)
{
	double mysum = 0.0;
	double x;

    int num_moments = LM.size();
    column_vector z(num_moments);
	
	double f;
	
	for(int i=0; i < N; i++){
		x = -1.0 + 2.0*rand01();

		z = 0;
		for (int k = 0; k < num_moments; ++k)
			z(k) = pow(x,k+1);

		f = pow(x,order)*(exp(-trans(LM)*z));
		
//		std::cout << x << ": " << f << std::endl;
		
		mysum+=f;

	}
	
	return 2*mysum/((double)(N));
}

int main(int argc, const char * argv[]) {
	/* set the precision of cout */
	std::cout.precision(22); // to be honest, I don't know the largest possible number :)
    
    /* initialize random seed: */
	srand(time(NULL));
    
    int k = 5;
    column_vector lambda(k);
    lambda = 0.0;
    
	lambda(0) = -1;
	lambda(1) = 0.3;
	lambda(2) = -0.1;
	lambda(3) = 2;
	lambda(4) = 0.3;
    
    auto mom_function = [&](double x)->double { return gg(x, 3, lambda);};
    double F_Dlib = integrate_function_adapt_simp(mom_function, -1.0, 1.0, 1e-10);
    std::cout << " - Dlib  : " << F_Dlib << std::endl;


	double F_mc = mymc(lambda, 3, 1000000);
    std::cout << " - myMC  : " << F_mc << std::endl;
    
    return 0;
    
}

