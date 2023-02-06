#include <iostream>
#include <cmath>
#include <array>
#include <algorithm>

//array sized 10 for storing simulation variables
typedef std::array<double,14> Arrayofdouble;

class eos
{
protected:

    std::vector<std::vector<double> > eos_array_1;

public:
    eos();
    
    double compute_pressure(double alpha_a, double alpha_rho_a, double alpha_rho_b, double energy, double v_x, double v_y, double gamma_a, double pinf_a, double gamma_b, double pinf_b, double p_max, double epsinf_a, double epsinf_b, double cv_a);
    double compute_energy(double alpha_a, double alpha_rho_a, double alpha_rho_b, double& pressure, double v_x, double v_y, double gamma_a, double pinf_a, double gamma_b, double pinf_b, double p_max, double epsinf_a, double epsinf_b);
    double compute_soundspeed(double alpha_a, double alpha_rho_a, double alpha_rho_b, double& pressure, double gamma_a, double pinf_a, double gamma_b, double pinf_b, double p_max, double epsinf_a, double epsinf_b, double cv_a, double cv_b);
    
    Arrayofdouble prim_to_u(Arrayofdouble prim, double gamma_a, double pinf_a, double gamma_b, double pinf_b, double p_max, double epsinf_a, double epsinf_b);  
    Arrayofdouble u_to_prim(Arrayofdouble u, double gamma_a, double pinf_a, double gamma_b, double pinf_b, double p_max, double epsinf_a, double epsinf_b, double cv_a);
    
    Arrayofdouble fluxf(Arrayofdouble u, Arrayofdouble prim);
    Arrayofdouble yfluxf(Arrayofdouble u, Arrayofdouble prim);

    Arrayofdouble eos_interpolation(double pressure_input);
    Arrayofdouble eos_interpolation2(double vl_input);
    Arrayofdouble eos_interpolation3(double vg_input);

    double compute_massfraction(double alpha_rho_a, double alpha_rho_b, double pressure);
    double pressure_bisectionsearch(double p_max, double p_min, double alpha_rho_a, double alpha_rho_b, double gamma_a, double pinf_a, double gamma_b, double pinf_b, double epsinf_a, double epsinf_b, double eps_real);

};


