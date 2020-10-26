#include <iostream>
#include <fstream>

#include "DE_Load.h"
#include "DE_Solve.h" ///DE_Solve

int main()
{
    std::cout << "SATART\n";
    std::string EQ_name_input_Jacobi, EQ_name_input_EQ, EQ_name_output_Results;

    ///repeat according manager
        //DE_Load_EQ(EQ_name_input_EQname_EQ);
       // DE_Load_Jacobi(EQ_name_input_Jacobi,name_Jacobi);

        ///A,B,C,e0,e2,e4 - coefficients of equation y'' = y'(e0y^2+e2y^4+e4y^4)+Ay+By^2+Cy^3
        model_coefficients D{1.0, 0, 0, 0.0000005, -0.05, 1.0};
        ///initial conditions
        double initial[2];
            initial[0] = 0.7;
            initial[1] = 0.0;
        ///range and time step
        double time[3];
            time[0] = 0.0;
            time[1] = 1000.0;
            time[2] = 0.01;

        DE_Solve(D, initial, time);

       // DE_Save_Results(EQ_name_output_Results);
}
