//
//  main.cpp
//  Theta_problem_MV
//
//  Created by Ganna Marchenko on 11/05/17.
//
//

#ifndef USE_DLIB
 #error 'This example is for DLIB'
#endif

#ifndef USE_CUBA
 #error 'This example is for CUBA'
#endif


#include <iostream>

/* DLIB */
#include "dlib/matrix.h"
#include "dlib/numeric_constants.h"
#include "dlib/numerical_integration.h"
#include "dlib/optimization.h"

/* CUBA */
#include "cuba.h"

/* user-defined */
#include "ExtraParameters.h"
#include "Integrator.h"

using namespace std;
using namespace dlib;

typedef matrix<double,0,1> column_vector;
typedef matrix<double,1,0> row_vector;

#define STOP_TOLERANCE 1e-06

#define KEY1 47
#define KEY2 1
#define KEY3 1
#define MAXPASS 5
#define BORDER 0.
#define MAXCHISQ 10.
#define MINDEVIATION .25
#define NGIVEN 0
#define LDXGIVEN NDIM
#define NEXTRA 0

#define KEY 0

static column_vector cLM;
static column_vector cgrad;
static matrix<double> chess;

double get_functions_obj (const column_vector& _LM, const column_vector& _Mom, double eps, const matrix<double>& mom_powers, int k)
//compute objective function value
{
//    //compute normalization
//    column_vector LM = _LM;
//    column_vector Mom = _Mom;
//    matrix<double> D = mom_powers;
//    
//    Integrator integrator;
//    ExtraParameters xp (k, Mom, LM, 0.0, 0.0, mom_powers, 0, 0);
//    integrator.USERDATA = &xp;
//    double F_ = integrator.computeVegas();
//    return trans(Mom)*LM + log(F_);// + eps*sum(LM);
    
    column_vector LM = _LM;
    cLM = _LM;
    column_vector Mom = _Mom;
    matrix<double> D = mom_powers;
    //# of variables
    long n = D.nr();
    
    column_vector I(n); //theoretical moments
    column_vector grad(n); //gradient
    matrix<double> hess; //hessian
    hess.set_size(n,n);
    
    Integrator integrator;
    //setting to compute normalization constant
    ExtraParameters xp(k, Mom, LM, 0.0, 0.0, D, 0, 0);
    integrator.USERDATA = &xp;
    
    double F_ = integrator.computeVegas();
    
    //modify settings to compute theoretical moments
    xp.type = 2;
    for (int j = 0; j < n; j++)
    {
        xp.order = j;
        I(j) = integrator.computeVegas();
        grad(j) = Mom(j) - I(j)/F_;
    }

    cout <<grad<< endl;
    cout <<LM<< endl;
    cout <<Mom<< endl;
    //modify settings to compute extra theoretical moments for hess
    xp.type = 3;
    double temp = 0.0;
    
    for (int i = 0; i < n; i++)
        for (int j = i; j < n; j++)
        {
            xp.order = i;
            xp.order2 = j;
            temp = integrator.computeVegas()/F_ - I(i)*I(j)/(F_*F_);
            
            hess(i,j) = temp;
            hess(j,i) = temp;
        }
    
    cgrad = grad;
    chess = hess;
    return trans(Mom)*LM + log(F_);// + eps*sum(LM);
}

column_vector get_functions_grad (const column_vector& _LM, const column_vector& _Mom, int k, const matrix<double>& mom_powers)
//compute gradient
{
    
    column_vector LM = _LM;
    
    if (LM != cLM)
    {
		cout << "--- GRADIENT: lagrange multipliers changed ----" << endl;

        column_vector Mom = _Mom;
        matrix<double> D = mom_powers;
        long n = D.nr();
        column_vector grad(n);
    
        Integrator integrator;
        ExtraParameters xp(k, Mom, LM, 0.0, 0.0, D, 0, 0);
        integrator.USERDATA = &xp;

        double F_ = integrator.computeVegas();
        xp.type = 2;
    
        for (int j = 0; j < n; j++)
        {
            xp.order = j;
            grad(j) = Mom(j) - integrator.computeVegas()/F_;
        }
        cgrad = grad;
        
    }
    
    return cgrad;
}

matrix<double> get_functions_hess (const column_vector& _LM, const column_vector& _Mom, int k, const matrix<double>& mom_powers)
{
    column_vector LM = _LM;
    
    if (LM != cLM)
    {
		cout << "--- HESSIAN: lagrange multipliers changed ----" << endl;
    
        column_vector Mom = _Mom;
        matrix<double> D = mom_powers;
        //# of variables
        long n = D.nr();
    
        column_vector I(n); //theoretical moments
        column_vector grad(n); //gradient
        matrix<double> hess; //hessian
        hess.set_size(n,n);
    
        Integrator integrator;
        //setting to compute normalization constant
        ExtraParameters xp(k, Mom, LM, 0.0, 0.0, D, 0, 0);
        integrator.USERDATA = &xp;
    
        double F_ = integrator.computeVegas();
    
        //modify settings to compute theoretical moments
        xp.type = 2;
        for (int j = 0; j < n; j++)
        {
            xp.order = j;
            I(j) = integrator.computeVegas();
            grad(j) = Mom(j) - I(j)/F_;
        }
  
        //modify settings to compute extra theoretical moments for hess
        xp.type = 3;
        double temp = 0.0;
    
        for (int i = 0; i < n; i++)
            for (int j = i; j < n; j++)
            {
                xp.order = i;
                xp.order2 = j;
                temp = integrator.computeVegas()/F_ - I(i)*I(j)/(F_*F_);
            
                hess(i,j) = temp;
                hess(j,i) = temp;
            }
        chess = hess;
    }

    return chess;

    }

matrix <double> ReadFile(std::string s_Path, long T, int d)
{
    matrix <double> m_Data;
    m_Data.set_size(T,d);
    string strFullName = s_Path;
    ifstream dataFile(strFullName.c_str());

    if (!dataFile)
    {
        std::cerr << "Error opening file!\n";
        return m_Data;
    }
    
    string line;
    getline(dataFile, line);
    dataFile.close();
   
    stringstream iss(line);
    
    int i = 0; int j = 0;
    while(true)
    {
        string val;
        double data;
        
        getline(iss, val, ',');
        
        if(!iss.good())
            break;
        
        stringstream convertor(trim(val));
        convertor >> data;
        m_Data(i,j) = data;
        if (i == T-1)
        {
            i = 0; j++;
        }
        else
            i++;
    }
   return m_Data;
}

int main(int argc, const char * argv[]) {
    
    //number of moments
    int k = 2;
    //length of the time-series
    long T = 800;
    //dimension of the data
    int d = 2;
    //read data
    matrix <double> xt;
    xt.set_size(T,d);
    xt = ReadFile("xt.txt", T, d);
 
 //   cout << xt << endl;
    
    //range of the moments p = [0,1,2,...k]
    matrix<long> p = range(0,k);
    
    //number of rows in the matrix of powers
    int mp_row = pow(k+1,d); /* upper bound on row number */
    
    //matrix of powers for calculations on joint moments
    matrix<double> mom_powers; /* = D from THE paper */
    mom_powers.set_size(mp_row, d);
    
    //steps for each dimension
    int step;
    int s;
    int ind = 0;
    //compute powers matrix
    for (int i = d-1; i >= 0; i--)
    {
        step = pow(k+1,ind);
        ind++;
        s = 0;
        for (int j = 0; j< mp_row; j = j + step)
        {
            set_subm(mom_powers,range(j,j+step-1), range(i,i)) = p(s);
            if (s == k)
                s = 0;
            else
                s++;
        }
    }
    
    //remove all rows where the sum of the elements is 0 or > k
    for (long j = 0; j< mp_row; j++)
        if (sum(subm(mom_powers, range(j,j), range(0,d-1))) > k || sum(subm(mom_powers, range(j,j), range(0,d-1))) == 0){
            mom_powers = remove_row(mom_powers, j);
            mp_row--;
            j--;
        }
    
    cout << "mom_powers ------------------------" << endl;
    cout << mom_powers << endl;
    cout << "-----------------------------------" << endl;
   
    /* Mom - sample moments */
    column_vector Mom(mp_row);
    cLM.set_size(mp_row,1);
    
    column_vector temp(T);
    matrix<double> data_powers; /* powers of data x_1^j * x_2^i * etc */
    data_powers.set_size(T, mp_row);
    
    //create data power matrix
    for (int i = 0; i < mp_row; i++)
    {
        temp = 1;
        for (int j=0; j<d; j++)
        {
            column_vector temp2(T);
            temp2 = pow(colm(xt,j),mom_powers(i,j));
//            cout << temp2 << endl;
            temp = pointwise_multiply(temp, temp2);
 //           cout << temp << endl;
        }
        set_colm(data_powers,i) = temp;
        Mom(i) = sum(temp)/T;
    }
    
//    cout << data_powers << endl;
    cout << "Moments ---------------------------" << endl;
    cout << Mom << endl;
    cout << "-----------------------------------" << endl;

    
    column_vector starting_point(mp_row);
    starting_point = 0.0;
    double eps = 0.0;
    
    auto get_functions_obj_lambda = [&](const column_vector& c)->double { return get_functions_obj(c, Mom, eps, mom_powers, k);};
    auto get_functions_grad_lambda = [&](const column_vector& c)->column_vector { return get_functions_grad(c, Mom, k, mom_powers);};
    
    auto get_functions_hess_lambda = [&](const column_vector& c)->matrix<double> { return get_functions_hess(c, Mom, k,mom_powers);};

    
    find_min_box_constrained(newton_search_strategy(get_functions_hess_lambda),
                             objective_delta_stop_strategy(STOP_TOLERANCE).be_verbose(),
                             get_functions_obj_lambda, get_functions_grad_lambda, starting_point, -1e12, 1e12 );

//    find_min_box_constrained(bfgs_search_strategy(),
//                             objective_delta_stop_strategy(STOP_TOLERANCE).be_verbose(),
//                            get_functions_obj_lambda, derivative(get_functions_obj_lambda), starting_point, -1e12, 1e12 );
////
//
//    starting_point = 0.0;
//    find_min_box_constrained(bfgs_search_strategy(),
//                             objective_delta_stop_strategy(STOP_TOLERANCE).be_verbose(), get_functions_obj_lambda, get_functions_grad_lambda, starting_point, -1e12, 1e12 );
////
    
    cout << starting_point << endl;
    
    
    matrix<double> residuum;
    residuum.set_size(T, 1);
    
    residuum = data_powers*starting_point;
    
    cout << residuum << endl;
//
//   
//
//    cout << LM << endl;
//    cout << mom_powers << endl;
//    
//    Integrator integrator;
//    ExtraParameters xp (k, Mom, LM, 0.5, 0.0, mom_powers, 0, 0);
//    integrator.USERDATA = &xp;
//    integrator.computeVegas();
//    
//    double fff = get_functions_obj (LM, Mom, 0.0, mom_powers, k);
//    
//    column_vector grad = get_functions_grad (LM, Mom, k, mom_powers);
//    cout << grad << endl;
//
//    matrix<double> hess = get_functions_hess (LM, Mom, k, mom_powers);
//    cout << hess << endl;
    ////    integrator.computeSuave();
////    integrator.computeDivonne();
////    integrator.computeCuhre();

    return 0;
}


