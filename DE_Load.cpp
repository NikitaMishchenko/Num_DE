#include <fstream>
#include <string>

#include <boost/numeric/odeint.hpp>
#include <boost/phoenix/core.hpp>
#include <boost/phoenix/operator.hpp>

typedef boost::numeric::ublas::vector<double> vector_type;
typedef boost::numeric::ublas::matrix<double> matrix_type;

void DE_Load_model_coefficients(std::string s){
    std::ifstream fin(s);

    fin.close();
}

void DE_Load_Jacobi(std::string s, matrix_type &J){
    std::ifstream fin(s);

    fin.close();
}

void DE_Load_EQ(std::string s, vector_type &dxdt){
    std::ifstream fin(s);

    fin.close();
}

void DE_Save_Results(std::string s){
    std::ofstream fout(s);

    fout.close();
}
