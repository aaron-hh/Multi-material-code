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
    
    double compute_pressure(double alpha_a, double alpha_rho_a, double alpha_rho_b, double energy, double v_x, double v_y, double gamma_a, double pinf_a, double gamma_b, double pinf_b);
    double compute_energy(double alpha_a, double alpha_rho_a, double alpha_rho_b, double pressure, double v_x, double v_y, double gamma_a, double pinf_a, double gamma_b, double pinf_b);
    double compute_soundspeed(double alpha_a, double alpha_rho_a, double alpha_rho_b, double pressure, double gamma_a, double pinf_a, double gamma_b, double pinf_b);
      
    Arrayofdouble prim_to_u(Arrayofdouble prim, double gamma_a, double pinf_a, double gamma_b, double pinf_b);  
    Arrayofdouble u_to_prim(Arrayofdouble u, double gamma_a, double pinf_a, double gamma_b, double pinf_b);
    
    Arrayofdouble fluxf(Arrayofdouble u, Arrayofdouble prim);
    Arrayofdouble yfluxf(Arrayofdouble u, Arrayofdouble prim);

};
