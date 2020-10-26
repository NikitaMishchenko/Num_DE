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

struct vdp_stiff
{
    void operator()( const vector_type &x , vector_type &dxdt , double t );
};

struct vdp_stiff_jacobi
{
    void operator()( const vector_type &x , matrix_type &J , const double &t , vector_type &dfdt );
};

void general();

#endif // DE_SOLVE_H_INCLUDED
