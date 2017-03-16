//
//  Theta_problem
//
//  Created by Anya Marchenko on 28.02.17.
//  Copyright © 2017 Anna Marchenko. All rights reserved.
//

#include <iostream>
#include "dlib/matrix.h"
#include "dlib/numeric_constants.h"
#include "dlib/numerical_integration.h"
#include "dlib/optimization.h"

using namespace std;
using namespace dlib;

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

double get_functions_obj (const column_vector& LM, const column_vector& Mom, double eps)
//compute objective function value
{
    //compute normalization
    column_vector Vec = LM;
    auto mom_function = [&](double x)->double { return gg(x, 0, Vec);};//std::bind(gg, _1,  1, 2);
    double F_ = integrate_function_adapt_simp(mom_function, -1.0, 1.0, 1e-10);
    
    return trans(Mom)*LM + log(F_);// + eps*sum(LM);
    
}

column_vector get_functions_grad (const column_vector& LM, const column_vector& Mom, int k)
//compute gradient
{
    column_vector grad(k);
    column_vector I(k);
    
    //compute normalization
    column_vector Vec = LM;
    auto mom_function = [&](double x)->double { return gg(x, 0, Vec);};//std::bind(gg, _1,  1, 2);
    double F_ = integrate_function_adapt_simp(mom_function, -1.0, 1.0, 1e-10);
    
    //theoretical moments
    int i = 0;
    while (i < k)
    {
        auto mom_function = [&](double x)->double { return gg(x, i+1, Vec);};
        I(i) = integrate_function_adapt_simp(mom_function, -1.0, 1.0, 1e-10);
        i++;
    }
    
    for (int i = 0; i < k; ++i)
        grad(i) = Mom(i) - I(i)/F_;
    
//    double L1 = grad(0);
//    double L2 = grad(1);
    return grad;
}

matrix<double> get_functions_hess (const column_vector& LM, const column_vector& Mom, int k)
{
    matrix<double> hess(k, k);
    
    column_vector I(2*k);
    
    //compute normalization
    column_vector Vec = LM;
    auto mom_function = [&](double x)->double { return gg(x, 0, Vec);};//std::bind(gg, _1,  1, 2);
    double F_ = integrate_function_adapt_simp(mom_function, -1.0, 1.0, 1e-10);
    
    //theoretical moments
    int i = 0;
    while (i < 2*k)
    {
        auto mom_function = [&](double x)->double { return gg(x, i+1, Vec);};
        I(i) = integrate_function_adapt_simp(mom_function, -1.0, 1.0, 1e-10);
        i++;
    }
    
    for (int i = 0; i < k; ++i)
        for (int j = 0; j < k; ++j)
            hess(i,j) = I(i+j+1)/F_ - I(i)*I(j)/(F_*F_);
    
//    double L1 = hess(0,0);
//    double L2 = hess(0,1);
//    double L3 = hess(1,0);
//    double L4 = hess(1,1);
    return hess;
}




int main(int argc, const char * argv[]) {
    int k = 14;
    double eps = 0;
    column_vector Mom(14);
    Mom = 0.0;
    column_vector starting_point(k);
    starting_point = 0.0;
    
    Mom(0) = -0.059475865319633;
    Mom(1) = 0.129300495193170;
    Mom(2) = -0.010921752705306;
    Mom(3) = 0.052561399966191;
    Mom(4) = -0.000939888193712;
    Mom(5) = 0.033982495322740;
    Mom(6) = 0.000998053697573;
    Mom(7) = 0.027200750830201;
    Mom(8) = 0.001132732366870;
    Mom(9) = 0.024077478196280;
    Mom(10) = 0.000913083305865;
    Mom(11) = 0.022436802669071;
    Mom(12) = 0.000674858503998;
    Mom(13) = 0.021506588139508;
//    Mom(14) = 0.000484341797861;
//    Mom(15) = 0.020953408533696;
//    Mom(16) = 0.000343884350430;
//    Mom(17) = 0.020613645620372;
//    Mom(18) = 0.000243192736410;
//    Mom(19) = 0.020400057227090;


	/* set the precision of cout */
	std::cout.precision(22); // to be honest, I don't know the largest possible number :)
    
    auto get_functions_obj_lambda = [&](const column_vector& x)->double { return get_functions_obj(x, Mom, eps);};
    auto get_functions_grad_lambda = [&](const column_vector& x)->column_vector { return get_functions_grad(x, Mom, k);};
    auto get_functions_hess_lambda = [&](const column_vector& x)->matrix<double> { return get_functions_hess(x, Mom, k);};
    
    find_min_box_constrained(newton_search_strategy(get_functions_hess_lambda),
                             objective_delta_stop_strategy(STOP_TOLERANCE).be_verbose(),
                             get_functions_obj_lambda, get_functions_grad_lambda, starting_point, -1e12, 1e12 );
    
    double L1 = starting_point(0);
    double L2 = starting_point(1);
    double L3 = starting_point(2);
    double L4 = starting_point(3);
    double L5 = starting_point(4);
    double L6 = starting_point(5);
    
    double L7 = starting_point(6);
    double L8 = starting_point(7);
    double L9 = starting_point(8);
    double L10 = starting_point(9);
    double L11 = starting_point(10);
    double L12 = starting_point(11);

	for(int i=0;i<k;i++){
		std::cout << "solution[" << i+1 << "] = " << starting_point(i) << std::endl;
    }
    std::cout << " - function value  : " << get_functions_obj (starting_point, Mom, eps) << std::endl;
    
    return 0;
    
}

