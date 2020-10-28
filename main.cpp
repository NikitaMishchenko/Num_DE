#include <iostream>
#include <fstream>

#include "DE_Load.h"
#include "DE_Solve.h" ///DE_Solve

void equation_2order(const double a, const double b, const double c,
                        double& root1Re, double& root1Im, double& root2Re, double& root2Im){
    double D = b*b - 4.0*a*c;
    ///root1 +sqrt(D)
    if(D >= 0.0)
    {///onlyRe
        root1Re = (-1.0*b+sqrt(D))/2.0;
            root1Im = 0.0;
        root2Re = (-1.0*b-sqrt(D))/2.0;
            root2Im = 0.0;
    }else{///Re+Im
        root1Re = -1.0*b/2.0;
            root1Im = sqrt(-1.0*(D))/2.0;
        root2Re = -1.0*b/2.0;
            root2Im = -1.0*sqrt(-1.0*(D))/2.0;
    }
};

void Find_Dynamic_depend_on_LimitCycle(const double l1, const double l2, double &e0, double& e2, double &e4){
    e0 = l1*l1*l2*l2;
    e2 = -(l1*l1+l2*l2);
    e4 = 1.0;
}

void Dyn_check(std::string output_file_name, double l1,double l2) ///returns the dynamic coeff
{
    double e0, e2, e4;
    Find_Dynamic_depend_on_LimitCycle(l1, l2, e0,e2,e4);
        std::cout << e0 << "\t" << e2 << "\t" << e4 << std::endl;
    std::ofstream fout(output_file_name);
        fout << "limit cycles: \n\tl1 = " << l1 << "\n\tl2 = " << l2 << std::endl;
        fout << "e0\t" << "e2\t" << "e4\n";
        fout << e0 << "\t" << e2 << "\t" << e4 << std::endl;
    fout.close();

    double rx1,ix1,rx2,ix2;
        equation_2order(e4, e2, e0, rx1,ix1,rx2,ix2);
        std::cout << "rho1 = " << rx1 << " " << ix1 << std::endl
                    << "rho2 = " << rx2 << " " << ix2 << std::endl << std::endl;
        std::cout << "x1/2 " << sqrt(fabs(rx1)) << " " <<  sqrt(fabs(rx2)) << std::endl;
}

void Solver_Loop(std::string instructions){
    ///repeat according manager
    //std::string instructions = "instructions.txt";
    std::ifstream fin(instructions);
        double initial[2];
        double time[3];
        model_coefficients D;
    std::string s_buff = "";
    double d_buff = 0.0;

    double rx1,ix1,rx2,ix2;

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
        equation_2order(D.e4, D.e2, D.e4, rx1,ix1,rx2,ix2);
        std::cout << "rho1 = " << rx1 << " " << ix1 << std::endl
                    << "rho2 = " << rx2 << " " << ix2 << std::endl << std::endl;
        std::cout << "x1/2 " << sqrt(fabs(rx1)) << " " <<  sqrt(fabs(rx2)) << std::endl;

        double abs_err, real_err;
            abs_err = real_err = 1.0e-10;
        DE_Solve(D, initial, time, abs_err, real_err, s_buff);
        std::cout << s_buff << " solved\n\n";
    }
}

int main()
{
    std::string mode;
    std::cout << "Enter mode:\n\teq - for solving DE\n\telse Dyn_check\n";
    std::cin >> mode;
    //if( mode == "eq")
     //   Solver_Loop("instructions4.txt");
   // else
      //  Dyn_check("limit_cycles.txt", 0.39, 0.215);

    double rx1,ix1,rx2,ix2;
    equation_2order(6.0,-0.1*2,0.0009, rx1,ix1,rx2,ix2);
    std::cout << rx1 << "\t"
                << ix1 << "\t"
                 << rx2 << "\t"
                  << ix2 << "\n"
                << sqrt(rx1) << "\t"
                << sqrt(rx2) << "\t";


}
