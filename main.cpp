#include <iostream>
#include <fstream>

#include "DE_Load.h"
#include "DE_Solve.h" ///DE_Solve

int main()
{
    std::cout << "SATART\n";
    std::string EQ_name_input_Jacobi, EQ_name_input_coeff, EQ_name_input_Init, EQ_name_output_Results;

    ///repeat according manager
    std::string instructions = "instructions.txt";
    std::ifstream fin(instructions);
        double initial[2];
        double time[3];
        model_coefficients D;
    std::string s_buff = "";
    double d_buff = 0.0;


    if(fin.is_open())
    while(!fin.eof()){        ///A,B,C,e0,e2,e4 - coefficients of equation y'' = -y'(e0y^2+e2y^4+e4y^4)-(Ay+By^2+Cy^3)

        ///name/id
        fin >> s_buff;
        ///range and time step;
        fin >> time[1];// = 1000.0;
        fin >> time[2];///time_step// = 0.01;

        ///initial conditions ///x0y0z0
        fin >> time[0];
        fin >> initial[0];// = 0.7;
        fin >> initial[1];// = 0.0;

        ///DE coefficients
        fin >> D.e0;// = 0.0000005;
        fin >> D.e2;// = -0.05;
        fin >> D.e4;// = 1.0;
        fin >> d_buff;///t_lag
            D.e0 *= d_buff;
            D.e2 *= d_buff;
            D.e4 *= d_buff;

        fin >> D.A;// = 1.0;
        fin >> D.B;// = 0.0
        fin >> D.C;// = 0.0;

        ///SOLVING DE
        std::cout << s_buff << "\t"
                    << time[0] << "\t" << time[1] << "\t" << time[2] << "\t"
                    << initial[0] << "\t" << initial[1] << "\t"
                        << D.e0 << "\t"
                        << D.e2 << "\t"
                        << D.e4 << "\t"
                        << D.A << "\t"
                        << D.B << "\t"
                        << D.C << "\n";

        DE_Solve(D, initial, time, s_buff + ".txt");
    }

       // DE_Save_Results(EQ_name_output_Results);
}
