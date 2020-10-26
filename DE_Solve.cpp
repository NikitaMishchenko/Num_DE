
/*
 * van_der_pol_stiff.cpp
 *
 * Created on: Dec 12, 2011
 *
 * Copyright 2012 Karsten Ahnert
 * Copyright 2012-2013 Rajeev Singh
 * Copyright 2012-2013 Mario Mulansky
 *
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or
 * copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#include <iostream>
#include <fstream>
#include <utility>

#include <boost/numeric/odeint.hpp>

#include <boost/phoenix/core.hpp>
#include <boost/phoenix/operator.hpp>

#include "DE_Load.h"
#include "DE_Solve.h"

using namespace boost::numeric::odeint;
namespace phoenix = boost::phoenix;


typedef boost::numeric::ublas::vector< double > vector_type;
typedef boost::numeric::ublas::matrix< double > matrix_type;


model_coefficients D;///{1.0, 0, 0, 0.0000005, -0.05, 1.0};//{1, 0,0, 1.0, -4.0, 1.0};
//{1.0, 0, 0, 0.0000005, -0.05, 1.0}; outer 0.316 inner 0.05

void vdp_stiff::operator()( const vector_type &x , vector_type &dxdt , double t)
{
    double A, B, C, e0, e2, e4;
    A = D.A; B = D.B; C = D.C;
    e0 = D.e0; e2 = D.e2; e4 = D.e4;

    dxdt[0] = x[1];
    dxdt[1] = -(A*x[0] + B*x[0]*x[0]*x[0] + C*x[0]*x[0]*x[0]*x[0]*x[0])
                - x[1]*(e0+e2*x[0]*x[0] + e4*x[0]*x[0]*x[0]*x[0]);
    //dxdt[1] = -x[0] - mu * x[1] * (x[0]*x[0]-1.0);
}

void vdp_stiff_jacobi::operator()( const vector_type &x , matrix_type &J , const double &t , vector_type &dfdt)
{
    double A, B, C, e0, e2, e4;
    A = D.A; B = D.B; C = D.C;
    e0 = D.e0; e2 = D.e2; e4 = D.e4;

    J(0, 0) = 0.0;
    J(0, 1) = 1.0;
    J(1, 0) = - (A + 3.0*B*x[0]*x[0] + 5.0*C*x[0]*x[0]*x[0]*x[0])
                - x[1]*(2.0*e2*x[0] + 4.0*e4*x[0]*x[0]*x[0]);
    J(1, 1) = - (e0 + e2*x[0]*x[0] + e4*x[0]*x[0]*x[0]*x[0]);
    //J(1, 0) = -1.0 - 2.0*mu * x[0] * x[1];
    //J(1, 1) = -mu * ( x[0] * x[0] - 1.0);

    dfdt[0] = 0.0;
    dfdt[1] = 0.0;
}

void DE_Solve(model_coefficients &nD, const double *initial, const double *time )
{
    D = nD;
    //[ integrate_stiff_system
    vector_type x( 2 );

    /// initial conditions
        x[0] = initial[0];//0.7;
        x[1] = initial[1];//0.0;
    double time_start, time_end, time_d;
        time_start = time[0];
        time_end = time[1];
        time_d = time[2];


    std::ofstream fout("Rosenbrock4.txt");
    size_t num_of_steps = integrate_const( make_dense_output< rosenbrock4< double > >( 1.0e-6 , 1.0e-6 ) ,
            std::make_pair( vdp_stiff() , vdp_stiff_jacobi() ) ,
            x , time_start, time_end, time_d
            , fout << phoenix::arg_names::arg2 << " " << phoenix::arg_names::arg1[0] << " " << phoenix::arg_names::arg1[1] << "\n"
            );
    //]
    fout.close();
    std::clog << num_of_steps << std::endl;

    //[ integrate_stiff_system_alternative
    vector_type x2( 2 );
    // initial conditions
    for (int i=0; i<2; i++)
        x2[i] = 1.0; //(1.0 * rand()) / RAND_MAX;

    //size_t num_of_steps2 = integrate_const( make_dense_output< runge_kutta_dopri5< vector_type > >( 1.0e-6 , 1.0e-6 ) ,
    //        vdp_stiff() , x2 , 0.0 , 1000.0 , 1.0
    //        , cout << phoenix::arg_names::arg2 << " " << phoenix::arg_names::arg1[0] << " " << phoenix::arg_names::arg1[1] << "\n"
    //        );
    //]
    //clog << num_of_steps2 << endl;
}
