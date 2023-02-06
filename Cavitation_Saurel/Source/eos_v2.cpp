// Defining the functions for variable conversion between primitive and conservative, flux functions for 1-D Multi-material solver with diffuse interface approach
// one-D primitive variable: alpha_xl, rho_xl, rho_xr, eps_xl, eps_xr, v_x, p_x
// a - material 1; b - material 2

#include "eos_v2.H"

//constructor 
eos::eos()
{}

//defining compute single component pressure function 
double eos::epstopressure(double rho, double eps, double pinf, double gamma)
{
    double p = (gamma - 1) * rho * eps - gamma * pinf;

    p = std::max(p, 1.0e-5);

    return p;
}

//defining compute single component specific internal energy function 
double eos::pressuretoeps(double rho, double p, double pinf, double gamma)
{
    double eps = (p + gamma * pinf) / ((gamma - 1) * rho);

    return eps;
}

//defining compute mixture pressure function
double eos::compute_pressure(double alpha_a, double rho_a, double gamma_a, double pinf_a, double rho_b, double gamma_b, double pinf_b, double energy, double v_x, double v_y)
{
    double alpha_b = 1 - alpha_a;

    double alpha_pinf_sum = (alpha_a * gamma_a * pinf_a) / (gamma_a - 1) + (alpha_b * gamma_b * pinf_b) / (gamma_b - 1);

    double alpha_sum = alpha_a / (gamma_a - 1) + alpha_b / (gamma_b - 1);

    double density = alpha_a * rho_a + alpha_b * rho_b;

    double p_mix_new = (density * energy - 0.5 * density * (v_x * v_x + v_y * v_y) - alpha_pinf_sum) / alpha_sum;

    return p_mix_new;
}

//defining compute total energy function
double eos::compute_energy(double alpha_a, double rho_a, double eps_a, double rho_b, double eps_b, double v_x, double v_y)
{
    double alpha_b = 1 - alpha_a;

    double rho = alpha_a * rho_a + alpha_b * rho_b; // total mixture density

    double y_a = alpha_a * rho_a / rho; // mass fraction for material a

    double y_b = alpha_b * rho_b / rho; // mass fraction for material b

    double energy = y_a * eps_a + y_b * eps_b + 0.5 * (v_x * v_x + v_y * v_y); // total mixture energy

    return energy;
}

//defining primitive to conservative function
Arrayofdouble eos::prim_to_u(Arrayofdouble prim, double gamma_a, double pinf_a, double gamma_b, double pinf_b) 
{
    Arrayofdouble u;

    double alpha_b = 1 - prim[0];
    double density = prim[0] * prim[1] + alpha_b * prim[2];
    double energy = compute_energy(prim[0], prim[1], prim[3], prim[2], prim[4], prim[5], prim[7]);

    u[0] = prim[0]; //alpha_a
    u[1] = prim[0] * prim[1]; //alpha_a * rho_a
    u[2] = alpha_b * prim[2]; //alpha_b * rho_b
    u[3] = prim[0] * prim[1] * prim[3]; //alpha_a * rho_a * eps_xa
    u[4] = alpha_b * prim[2] * prim[4]; //alpha_b * rho_b * eps_xb
    u[5] = density * prim[5]; //rho * v_x
    u[6] = density * energy; // rho * energy
    u[7] = density * prim[7]; // rho * v_y

    return u;
}

//defining conservative to primitive function
Arrayofdouble eos::u_to_prim(Arrayofdouble u, double gamma_a, double pinf_a, double gamma_b, double pinf_b)
{  
    Arrayofdouble prim;

    double alpha_b = 1 - u[0];

    prim[0] = u[0]; //alpha_a
    prim[1] = u[1] / u[0]; //rho_a
    prim[2] = u[2] / alpha_b; //rho_b
    prim[3] = u[3] / u[1]; //eps_a
    prim[4] = u[4] / u[2]; //eps_b

    double density = prim[0] * prim[1] + alpha_b * prim[2]; 
    prim[5] = u[5] / density; // v_x
    prim[7] = u[7] / density; // v_y

    double energy = u[6] / density;
    prim[6] = compute_pressure(prim[0], prim[1], gamma_a, pinf_a, prim[2], gamma_b, pinf_b, energy, prim[5], prim[7]); // mixture pressure

    return prim;
}

//defining flux function
Arrayofdouble eos::fluxf(Arrayofdouble u, Arrayofdouble prim)
{
    Arrayofdouble x_fluxfunction;

    x_fluxfunction[0] = prim[0] * prim[5]; // alpha_a * v_x 
    x_fluxfunction[1] = u[1] * prim[5]; // alpha_a * rho_a * v_x
    x_fluxfunction[2] = u[2] * prim[5]; // alpha_b * rho_b * v_x
    x_fluxfunction[3] = u[3] * prim[5]; // alpha_a * rho_a * eps_a * v_x
    x_fluxfunction[4] = u[4] * prim[5]; // alpha_b * rho_b * eps_b * v_x
    x_fluxfunction[5] = u[5] * prim[5] + prim[6]; // density * v_x * v_x + pressure
    x_fluxfunction[6] = (u[6] + prim[6]) * prim[5]; // (density * energy + pressure) * v_x
    x_fluxfunction[7] = u[5] * prim[7]; // density * v_x * v_y
    
    return x_fluxfunction;
}

//defining flux function
Arrayofdouble eos::yfluxf(Arrayofdouble u, Arrayofdouble prim)
{
    Arrayofdouble y_fluxfunction;

    y_fluxfunction[0] = prim[0] * prim[7]; // alpha_a * v_y 
    y_fluxfunction[1] = u[1] * prim[7]; // alpha_a * rho_a * v_y
    y_fluxfunction[2] = u[2] * prim[7]; // alpha_b * rho_b * v_y
    y_fluxfunction[3] = u[3] * prim[7]; // alpha_a * rho_a * eps_a * v_y
    y_fluxfunction[4] = u[4] * prim[7]; // alpha_b * rho_b * eps_b * v_y
    y_fluxfunction[5] = u[5] * prim[7]; // density * v_y * v_x
    y_fluxfunction[6] = (u[6] + prim[6]) * prim[7]; // (density * energy + pressure) * v_y
    y_fluxfunction[7] = u[7] * prim[7] + prim[6]; // density * v_y * v_y + pressure
    
    return y_fluxfunction;
}

//defining primitive to conservative function for second order extension
Arrayofdouble eos::secondprim_to_u(Arrayofdouble prim, double gamma_a, double pinf_a, double gamma_b, double pinf_b) 
{
    Arrayofdouble u;

    double alpha_b = 1 - prim[0];
    double density = prim[0] * prim[1] + alpha_b * prim[2];
    double energy = compute_energy(prim[0], prim[1], prim[3], prim[2], prim[4], prim[5], prim[6]);

    u[0] = prim[0]; //alpha_a
    u[1] = prim[0] * prim[1]; //alpha_a * rho_a
    u[2] = alpha_b * prim[2]; //alpha_b * rho_b
    // pressure is used as primitive variables for second order extension
    double eps_a = pressuretoeps(prim[1], prim[7], pinf_a, gamma_a);
    double eps_b = pressuretoeps(prim[2], prim[7], pinf_b, gamma_b);
    u[3] = prim[0] * prim[1] * eps_a; //alpha_a * rho_a * eps_xa
    u[4] = alpha_b * prim[2] * eps_b; //alpha_b * rho_b * eps_xb
    u[5] = density * prim[5]; //rho * v_x
    u[6] = density * prim[6]; // rho * v_y
    u[7] = density * energy; // rho * energy

    return u;
}

//defining conservative to primitive function for second order extension
Arrayofdouble eos::secondu_to_prim(Arrayofdouble u, double gamma_a, double pinf_a, double gamma_b, double pinf_b)
{  
    Arrayofdouble prim;

    double alpha_b = 1 - u[0];
    double density = u[1] + u[2];

    prim[0] = u[0]; //alpha_a
    prim[1] = u[1] / u[0]; //rho_a
    prim[2] = u[2] / alpha_b; //rho_b
    // use pressure as primitive variables for second order extension
    prim[3] = epstopressure(prim[1], u[3]/u[1], pinf_a, gamma_a); //p_a
    prim[4] = epstopressure(prim[2], u[4]/u[2], pinf_b, gamma_b); //p_b
    prim[5] = u[5] / density; // v_x
    prim[6] = u[6] / density; // v_y

    double energy = u[7] / density;
    prim[7] = compute_pressure(prim[0], prim[1], gamma_a, pinf_a, prim[2], gamma_b, pinf_b, energy, prim[5], prim[6]); // pressure

    return prim;
}

// defining compute energy xflux function for second extension
double eos::energyflux(Arrayofdouble recon_var, double gamma_a, double pinf_a, double gamma_b, double pinf_b)
{
    double density = recon_var[0] * recon_var[1] + (1.0 - recon_var[0]) * recon_var[2];

    double recon_p = compute_pressure(recon_var[0], recon_var[1], gamma_a, pinf_a, recon_var[2], gamma_b, pinf_b, recon_var[7]/density, recon_var[5], recon_var[6]);

    double energyflux = (recon_p + recon_var[7]) * recon_var[5];

    return energyflux;
}

// defining compute energy yflux function for second extension
double eos::yenergyflux(Arrayofdouble yrecon_var, double gamma_a, double pinf_a, double gamma_b, double pinf_b)
{
    double density = yrecon_var[0] * yrecon_var[1] + (1.0 - yrecon_var[0]) * yrecon_var[2];
    
    double recon_p = compute_pressure(yrecon_var[0], yrecon_var[1], gamma_a, pinf_a, yrecon_var[2], gamma_b, pinf_b, yrecon_var[7]/density, yrecon_var[5], yrecon_var[6]);

    double yenergyflux = (recon_p + yrecon_var[7]) * yrecon_var[6];

    return yenergyflux;
}