//
//  ExtraParameters.cpp
//  Theta_problem_MV
//
//  Created by Ganna Marchenko on 08/06/17.
//
//

#include "ExtraParameters.h"

ExtraParameters::ExtraParameters()
{
    k = 0;
    D = 0.0;
    eps= 0.0;
    Mom = 0.0;
    type = 0;
    LM = 0.0;
    L0 = 0.0;
    order = 0;
    order2 = 0;
}

ExtraParameters::ExtraParameters(int _k, column_vector _Mom, column_vector _LM, double _L0, double _eps, dlib::matrix<double> _D, int _type, int _order)
{
    k = _k;
    Mom = _Mom;
    LM = _LM;
    L0 = _L0;
    eps = _eps;
    D = _D;
    type = _type;
    order = _order;
    order2 = _order;
}

void ExtraParameters::Copy(ExtraParameters &_ExtraParameters)
{
    k = _ExtraParameters.k;
    Mom = _ExtraParameters.Mom;
    eps = _ExtraParameters.eps;
    D = _ExtraParameters.D;
    type = _ExtraParameters.type;
    LM = _ExtraParameters.LM;
    L0 = _ExtraParameters.L0;
    order = _ExtraParameters.order;
    order2 = _ExtraParameters.order2;
}

ExtraParameters::~ExtraParameters()
{
}