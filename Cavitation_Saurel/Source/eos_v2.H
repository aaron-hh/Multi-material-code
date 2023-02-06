#include <iostream>
#include <cmath>
#include <array>
#include <algorithm>

//array sized 10 for storing simulation variables
typedef std::array<double,14> Arrayofdouble;

class eos
{
protected:

public:
    eos();
    
    //double compute_pressure(double alpha_a, double rho_a, double eps_a, double gamma_a, double pinf_a, double rho_b, double eps_b, double gamma_b, double pinf_b);
    double epstopressure(double rho, double eps, double pinf, double gamma);
    double pressuretoeps(double rho, double p, double pinf, double gamma);

    double compute_pressure(double alpha_a, double rho_a, double gamma_a, double pinf_a, double rho_b, double gamma_b, double pinf_b, double energy, double v_x, double v_y);
    double compute_energy(double alpha_a, double rho_a, double eps_a, double rho_b, double eps_b, double v_x, double v_y);
    
    Arrayofdouble prim_to_u(Arrayofdouble prim, double gamma_a, double pinf_a, double gamma_b, double pinf_b);  
    Arrayofdouble u_to_prim(Arrayofdouble u, double gamma_a, double pinf_a, double gamma_b, double pinf_b);
    Arrayofdouble fluxf(Arrayofdouble u, Arrayofdouble prim);
    Arrayofdouble yfluxf(Arrayofdouble u, Arrayofdouble prim);

    Arrayofdouble secondprim_to_u(Arrayofdouble prim, double gamma_a, double pinf_a, double gamma_b, double pinf_b);
    Arrayofdouble secondu_to_prim(Arrayofdouble u, double gamma_a, double pinf_a, double gamma_b, double pinf_b);

    double energyflux(Arrayofdouble recon_var, double gamma_a, double pinf_a, double gamma_b, double pinf_b);
    double yenergyflux(Arrayofdouble yrecon_var, double gamma_a, double pinf_a, double gamma_b, double pinf_b);
};
