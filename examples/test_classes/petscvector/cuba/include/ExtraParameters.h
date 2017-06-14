//
//  ExtraParameters.h
//  Theta_problem_MV
//
//  Created by Ganna Marchenko on 08/06/17.
//
//

#ifndef __Theta_problem_MV__ExtraParameters__
#define __Theta_problem_MV__ExtraParameters__

#include <stdio.h>
#include "dlib/matrix.h"

using namespace dlib;

typedef matrix<double,0,1> column_vector;
typedef matrix<double,1,0> row_vector;

class ExtraParameters
{
public:
    int k;
    column_vector Mom;
    column_vector LM;
    double L0;
    double eps;
    matrix<double> D;
    int type; /* type of integrant =0,1,2,3 */
    int order; /* row of D matrix */
    int order2;
    
    ExtraParameters();
    ExtraParameters(int _k, column_vector _Mom, column_vector _LM, double _L0, double _eps, dlib::matrix<double> _D, int _type, int _order);
    ~ExtraParameters();
    
    void Copy(ExtraParameters& _ExtraParameters);
};

#endif /* defined(__Theta_problem_MV__ExtraParameters__) */
