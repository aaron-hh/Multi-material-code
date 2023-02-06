// Defining the functions for variable conversion between primitive and conservative, flux functions for 1-D Multi-material solver with diffuse interface approach
// one-D conservative variable: alpha_a, alpha_rho_a, alpha_rho_b, RHO_v_x, RHO_E, RHO_v_y
// one-D primitive variable: alpha_a, rho_a, rho_b, v_x, P, v_y
// a - material 1; b - material 2

#include "eos_v2.H"

//constructor 
eos::eos()
{}

//defining compute mixture pressure function
double eos::compute_pressure(double alpha_a, double alpha_rho_a, double alpha_rho_b, double energy, double v_x, double v_y, double gamma_a, double pinf_a, double gamma_b, double pinf_b)
{
    double alpha_b = 1 - alpha_a;

    double rho = alpha_rho_a + alpha_rho_b;

    double rho_eps = rho * (energy - 0.5 * (v_x * v_x + v_y * v_y));

    double gamma_minus1 = 1 / ((alpha_a / (gamma_a - 1)) + (alpha_b / (gamma_b - 1)));

    double gamma_pinf = gamma_minus1 * ((alpha_a * pinf_a * gamma_a) / (gamma_a - 1) + (alpha_b * pinf_b * gamma_b) / (gamma_b - 1));

    double pressure = gamma_minus1 * rho_eps - gamma_pinf;

    return pressure;
}

//defining compute total energy function
double eos::compute_energy(double alpha_a, double alpha_rho_a, double alpha_rho_b, double pressure, double v_x, double v_y, double gamma_a, double pinf_a, double gamma_b, double pinf_b)
{
    double alpha_b = 1 - alpha_a;

    double rho = alpha_rho_a + alpha_rho_b;

    double sum1 = (alpha_a / (gamma_a - 1)) + (alpha_b / (gamma_b - 1));

    double sum2 = ((alpha_a * pinf_a * gamma_a) / (gamma_a - 1) + (alpha_b * pinf_b * gamma_b) / (gamma_b - 1));

    double rho_eps = pressure * sum1 + sum2;

    double energy = rho_eps / rho + 0.5 * (v_x * v_x + v_y * v_y);

    return energy;
}

//defining primitive to conservative function
Arrayofdouble eos::prim_to_u(Arrayofdouble prim, double gamma_a, double pinf_a, double gamma_b, double pinf_b) 
{
    Arrayofdouble u;

    double alpha_b = 1 - prim[0];
    double rho = prim[0] * prim[1] + alpha_b * prim[2];
    double energy = compute_energy(prim[0], prim[0]*prim[1], alpha_b*prim[2], prim[4], prim[3], prim[5], gamma_a, pinf_a, gamma_b, pinf_b);

    u[0] = prim[0]; //alpha_a
    u[1] = prim[0] * prim[1]; //alpha_a * rho_a
    u[2] = alpha_b * prim[2]; //alpha_b * rho_b
    u[3] = rho * prim[3]; //rho * v_x
    u[4] = rho * energy; // rho * energy
    u[5] = rho * prim[5]; // rho * v_y

    return u;
}

//defining conservative to primitive function
Arrayofdouble eos::u_to_prim(Arrayofdouble u, double gamma_a, double pinf_a, double gamma_b, double pinf_b)
{  
    Arrayofdouble prim;

    double alpha_b = 1 - u[0];
    double rho = u[1] + u[2];

    prim[0] = u[0]; //alpha_a
    prim[1] = u[1] / u[0]; //rho_a
    prim[2] = u[2] / alpha_b; //rho_b
    prim[3] = u[3] / rho; // v_x
    prim[5] = u[5] / rho; // v_y

    double energy = u[4] / rho;
    prim[4] = compute_pressure(u[0], u[1], u[2], energy, prim[3], prim[5], gamma_a, pinf_a, gamma_b, pinf_b); // mixture pressure

    return prim;
}

//defining flux function
Arrayofdouble eos::fluxf(Arrayofdouble u, Arrayofdouble prim)
{
    Arrayofdouble x_fluxfunction;

    x_fluxfunction[0] = prim[0] * prim[3]; // alpha_a * v_x 
    x_fluxfunction[1] = u[1] * prim[3]; // alpha_a * rho_a * v_x
    x_fluxfunction[2] = u[2] * prim[3]; // alpha_b * rho_b * v_x
    x_fluxfunction[3] = u[3] * prim[3] + prim[4]; // density * v_x * v_x + pressure
    x_fluxfunction[4] = (u[4] + prim[4]) * prim[3]; // (density * energy + pressure) * v_x
    x_fluxfunction[5] = u[3] * prim[5]; // density * v_x * v_y
    
    return x_fluxfunction;
}

//defining flux function
Arrayofdouble eos::yfluxf(Arrayofdouble u, Arrayofdouble prim)
{
    Arrayofdouble y_fluxfunction;

    y_fluxfunction[0] = prim[0] * prim[5]; // alpha_a * v_y 
    y_fluxfunction[1] = u[1] * prim[5]; // alpha_a * rho_a * v_y
    y_fluxfunction[2] = u[2] * prim[5]; // alpha_b * rho_b * v_y
    y_fluxfunction[3] = u[5] * prim[3]; // density * v_y * v_x
    y_fluxfunction[4] = (u[4] + prim[4]) * prim[5]; // (density * energy + pressure) * v_y
    y_fluxfunction[5] = u[5] * prim[5] + prim[4]; // density * v_y * v_y + pressure
    
    return y_fluxfunction;
}

//defining compute sound speed function
double eos::compute_soundspeed(double alpha_a, double alpha_rho_a, double alpha_rho_b, double pressure, double gamma_a, double pinf_a, double gamma_b, double pinf_b)
{
    double alpha_b = 1 - alpha_a;

    double rho = alpha_rho_a + alpha_rho_b;

    double rho_a = alpha_rho_a / alpha_a;
    double rho_b = alpha_rho_b / alpha_b;

    double y_a = alpha_rho_a / rho;
    double y_b = alpha_rho_b / rho;

    double xi_a = 1 / (gamma_a - 1);
    double xi_b = 1 / (gamma_b - 1);

    double xi = alpha_a * xi_a + alpha_b * xi_b; // alpha_a * xi_a + alpha_b * xi_b

    double c_a = gamma_a * (pressure + pinf_a) / rho_a;
    double c_b = gamma_b * (pressure + pinf_b) / rho_b;

    double c = sqrt((y_a * xi_a * c_a + y_b * xi_b * c_b) / xi);

    return c;
}