#ifndef DE_SOLVE_H_INCLUDED
#define DE_SOLVE_H_INCLUDED

#include <utility>

#include <boost/numeric/odeint.hpp>
#include <boost/phoenix/core.hpp>
#include <boost/phoenix/operator.hpp>


//using namespace boost::numeric::odeint;
//namespace phoenix = boost::phoenix;

typedef boost::numeric::ublas::vector< double > vector_type;
typedef boost::numeric::ublas::matrix< double > matrix_type;

struct model_coefficients
{
    double A, B, C, e0, e2, e4;
//        A = 1.0; B = 0.0; C = 0.0;
        //e0 = 0.0005; e2 = -0.05; e4 = 1.0;
  //      e0 = 1.0; e2 = -4.0; e4 = 1.0;
};

struct vdp_stiff
{
    void operator()( const vector_type &x , vector_type &dxdt , double t );
};

struct vdp_stiff_jacobi
{
    void operator()( const vector_type &x , matrix_type &J , const double &t , vector_type &dfdt );
};

void DE_Solve(model_coefficients &, const double*, const double*);

#endif // DE_SOLVE_H_INCLUDED
