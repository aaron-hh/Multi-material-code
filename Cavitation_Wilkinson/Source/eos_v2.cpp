// Defining the functions for variable conversion between primitive and conservative, flux functions for 1-D Multi-material solver with diffuse interface approach
// one-D conservative variable: alpha_a, alpha_rho_a, alpha_rho_b, RHO_v_x, RHO_E, RHO_v_y
// one-D primitive variable: alpha_a, rho_a, rho_b, v_x, P, v_y
// a - material 1; b - material 2
#include <iostream>
#include <cmath>
#include <vector>
#include <array>
#include <algorithm>

#include "eos_v2.H"
#include <fstream>

typedef std::vector<double> Vectorofdouble;

//constructor 
eos::eos()
{
    std::ifstream eos_data("water_eos.txt");

    if (!eos_data.is_open()) 
    std::cout<<"Error opening file" ;

    int row = 193;
    int col = 4;

    for(size_t i=0; i < row; i++)
    {
        std::vector<double> vec;

        for(size_t j=0; j < col; j++)
        {
            double temp;
            eos_data >> temp;
            vec.push_back(temp);   
        }
        eos_array_1.push_back(vec);
    }

    eos_data.close();
};


//defining primitive to conservative function
Arrayofdouble eos::prim_to_u(Arrayofdouble prim, double gamma_a, double pinf_a, double gamma_b, double pinf_b, double p_max, double epsinf_a, double epsinf_b) 
{
    Arrayofdouble u;

    double alpha_b = 1.0 - prim[0];
    double rho = prim[0] * prim[1] + alpha_b * prim[2];
    double energy = compute_energy(prim[0], prim[0]*prim[1], alpha_b*prim[2], prim[4], prim[3], prim[5], gamma_a, pinf_a, gamma_b, pinf_b, p_max, epsinf_a, epsinf_b);

    u[0] = prim[0]; //alpha_a
    u[1] = prim[0] * prim[1]; //alpha_a * rho_a
    u[2] = alpha_b * prim[2]; //alpha_b * rho_b
    u[3] = rho * prim[3]; //rho * v_x
    u[4] = rho * energy; // rho * energy
    u[5] = rho * prim[5]; // rho * v_y

    return u;
}

//defining conservative to primitive function
Arrayofdouble eos::u_to_prim(Arrayofdouble u, double gamma_a, double pinf_a, double gamma_b, double pinf_b, double p_max, double epsinf_a, double epsinf_b, double cv_a)
{  
    Arrayofdouble prim;

    double alpha_b = 1.0 - u[0];
    double rho = u[1] + u[2];

    prim[0] = u[0]; //alpha_a
    prim[1] = u[1] / u[0]; //rho_a
    prim[2] = u[2] / alpha_b; //rho_b
    prim[3] = u[3] / rho; // v_x
    prim[5] = u[5] / rho; // v_y

    double energy = u[4] / rho;
    prim[4] = compute_pressure(u[0], u[1], u[2], energy, prim[3], prim[5], gamma_a, pinf_a, gamma_b, pinf_b, p_max, epsinf_a, epsinf_b, cv_a); // mixture pressure

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

//defining eos interpolation function
Arrayofdouble eos::eos_interpolation(double pressure_input)
{
    int row = 193;
    int col = 4;
    int index_i, index_iplus1;

    for(int i=0; i<row; i++)
    {
        if(eos_array_1[i][1] >= pressure_input)
        {
            index_iplus1 = i;
            index_i = i-1;
            break;
        }
    }

    double temp_i, temp_iplus1, pressure_i, pressure_iplus1, vl_i, vl_iplus1, vg_i, vg_iplus1;

    temp_i = eos_array_1[index_i][0];
    temp_iplus1 = eos_array_1[index_iplus1][0];

    pressure_i = eos_array_1[index_i][1];
    pressure_iplus1 = eos_array_1[index_iplus1][1];

    vl_i = 1.0/eos_array_1[index_i][2];
    vl_iplus1 = 1.0/eos_array_1[index_iplus1][2];

    vg_i = 1.0/eos_array_1[index_i][3];
    vg_iplus1 = 1.0/eos_array_1[index_iplus1][3];

    double temp_output, vl_output, vg_output;

    temp_output = (temp_iplus1 - temp_i) * (pressure_input - pressure_i) / (pressure_iplus1 - pressure_i) + temp_i;
    vl_output = (vl_iplus1 - vl_i) * (pressure_input - pressure_i) / (pressure_iplus1 - pressure_i) + vl_i;
    vg_output = (vg_iplus1 - vg_i) * (pressure_input - pressure_i) / (pressure_iplus1 - pressure_i) + vg_i;

    Arrayofdouble eos_output;

    eos_output[0] = temp_output;
    eos_output[1] = vl_output;
    eos_output[2] = vg_output;

    return eos_output;
}

//defining eos interpolation2 function which takes in specific volume for liquid
Arrayofdouble eos::eos_interpolation2(double vl_input)
{
    int row = 193;
    int col = 4;
    int index_i, index_iplus1;

    double rhol_input = 1.0/vl_input;

    for(int i=0; i<row; i++)
    {
        if(eos_array_1[i][2] <= rhol_input)
        {
            index_iplus1 = i;
            index_i = i-1;
            break;
        }
    }

    double temp_i, temp_iplus1, pressure_i, pressure_iplus1, vl_i, vl_iplus1, vg_i, vg_iplus1;

    temp_i = eos_array_1[index_i][0];
    temp_iplus1 = eos_array_1[index_iplus1][0];

    pressure_i = eos_array_1[index_i][1];
    pressure_iplus1 = eos_array_1[index_iplus1][1];

    vl_i = 1.0/eos_array_1[index_i][2];
    vl_iplus1 = 1.0/eos_array_1[index_iplus1][2];

    vg_i = 1.0/eos_array_1[index_i][3];
    vg_iplus1 = 1.0/eos_array_1[index_iplus1][3];

    double temp_output, pressure_output, vg_output;

    temp_output = (temp_iplus1 - temp_i) * (vl_input - vl_i) / (vl_iplus1 - vl_i) + temp_i;
    pressure_output = (pressure_iplus1 - pressure_i) * (vl_input - vl_i) / (vl_iplus1 - vl_i) + pressure_i;
    vg_output = (vg_iplus1 - vg_i) * (vl_input - vl_i) / (vl_iplus1 - vl_i) + vg_i;

    Arrayofdouble eos_output;

    eos_output[0] = temp_output;
    eos_output[1] = pressure_output;
    eos_output[2] = vg_output;

    return eos_output;
}

//defining eos interpolation3 function which takes in specific volume for vapour
Arrayofdouble eos::eos_interpolation3(double vg_input)
{
    int row = 193;
    int col = 4;
    int index_i, index_iplus1;

    for(int i=0; i<row; i++)
    {
        if(eos_array_1[i][3] >= 1.0/vg_input)
        {
            index_iplus1 = i;
            index_i = i-1;
            break;
        }
    }

    double temp_i, temp_iplus1, pressure_i, pressure_iplus1, vl_i, vl_iplus1, vg_i, vg_iplus1;

    temp_i = eos_array_1[index_i][0];
    temp_iplus1 = eos_array_1[index_iplus1][0];

    pressure_i = eos_array_1[index_i][1];
    pressure_iplus1 = eos_array_1[index_iplus1][1];

    vl_i = 1.0/eos_array_1[index_i][2];
    vl_iplus1 = 1.0/eos_array_1[index_iplus1][2];

    vg_i = 1.0/eos_array_1[index_i][3];
    vg_iplus1 = 1.0/eos_array_1[index_iplus1][3];

    double temp_output, pressure_output, vl_output;

    temp_output = (temp_iplus1 - temp_i) * (vg_input - vg_i) / (vg_iplus1 - vg_i) + temp_i;
    pressure_output = (pressure_iplus1 - pressure_i) * (vg_input - vg_i) / (vg_iplus1 - vg_i) + pressure_i;
    vl_output = (vl_iplus1 - vl_i) * (vg_input - vg_i) / (vg_iplus1 - vg_i) + vl_i;

    Arrayofdouble eos_output;

    eos_output[0] = temp_output;
    eos_output[1] = pressure_output;
    eos_output[2] = vl_output;

    return eos_output;
}

double eos::compute_massfraction(double alpha_rho_a, double alpha_rho_b, double pressure) // computing the mass fraction for liquid
{
    double v = 1 / (alpha_rho_a + alpha_rho_b);

    Arrayofdouble eos_output;

    double p_input = std::max(612.0, std::min(pressure, 22064000.0)); // making sure interpolation is possible

    eos_output = eos_interpolation(p_input);

    double vl_star = eos_output[1];
    double vg_star = eos_output[2];

    double lambda;

    lambda = (v - vg_star) / (vl_star - vg_star);

    return lambda;
}

double eos::pressure_bisectionsearch(double p_max, double p_min, double alpha_rho_a, double alpha_rho_b, double gamma_a, double pinf_a, double gamma_b, double pinf_b, double epsinf_a, double epsinf_b, double eps_real)
{
    double v = 1.0 / (alpha_rho_a + alpha_rho_b);

    double delta = 1e6;

    double p_guess = p_max;

    while (abs(delta)>1e-6)
    {

    // compute p_guess by taking the average between p_min and p_max 
    p_guess = (p_max + p_min) / 2.0;

    // compute eps_new
    double lambda = compute_massfraction(alpha_rho_a, alpha_rho_b, p_guess);

    Arrayofdouble eos_output;

    p_guess = std::max(612.0, std::min(p_guess, 22064000.0)); // making sure interpolation is possible

    eos_output = eos_interpolation(p_guess);

    double vl_star = eos_output[1];
    double vg_star = eos_output[2];

    double epsl_star = epsinf_a + vl_star * (p_guess + gamma_a * pinf_a) / (gamma_a - 1.0);

    double epsg_star = epsinf_b + vg_star * (p_guess + gamma_b * pinf_b) / (gamma_b - 1.0);

    double eps_new = lambda * epsl_star + (1.0 - lambda) * epsg_star;

    // iteration condition
    delta = eps_new - eps_real;

    if(eps_new * eps_real > 0.0) // On the same size of y domain
    {
        if(delta > 0.0)
        {
        p_max = p_guess;
        }
        else if(delta < 0.0)
        {
        p_min = p_guess;
        }
    }
    else if(eps_new * eps_real < 0.0) // On the opposite site of y domain
    {
        if(eps_new < 0.0)
        {
            p_min = p_guess;
        }
        else if(eps_new > 0.0)
        {
            p_max = p_guess;
        }
    }

    }
    return p_guess;
}

//defining compute mixture pressure function
double eos::compute_pressure(double alpha_a, double alpha_rho_a, double alpha_rho_b, double energy, double v_x, double v_y, double gamma_a, double pinf_a, double gamma_b, double pinf_b, double p_max, double epsinf_a, double epsinf_b, double cv_a)
{
    double pressure = 0.0;

    double alpha_b = 1.0 - alpha_a;

    // logic for checking simulation's location on phase diagram
    double rho = alpha_rho_a + alpha_rho_b;

    double v = 1.0 / rho;

    double v_critical = 1.0 / 322.0; 

    double p_sat = 0.0;

    double p_min = 612;

    // compute vlstar_min and vgstar_max (at minimum P and T)
    double vlstar_min = 1.0 / 999.793; // triple point
    double vgstar_max = 1.0 / 0.00485; // triple point

    // compute real specific internal energy 
    double eps_real = energy - 0.5 * (v_x * v_x + v_y * v_y);

    double eps_sat = 0.0;

    if(v >= vlstar_min && v <= vgstar_max)
    {

    if(v < v_critical) // on the left of saturation dome
    {
        Arrayofdouble eos_output2;

        eos_output2 = eos_interpolation2(v); // determine the saturated pressure point using v_a

        p_sat = eos_output2[1];

        eps_sat = epsinf_a + v * (p_sat + gamma_a * pinf_a) / (gamma_a - 1);
    }

    else if(v > v_critical) // on the right of saturation dome
    {
        Arrayofdouble eos_output3;

        eos_output3 = eos_interpolation3(v); // determine the saturated pressure point using v_b

        p_sat = eos_output3[1];

        eps_sat = epsinf_b + v * (p_sat + gamma_b * pinf_b) / (gamma_b - 1);
    }

    if(eps_real <= eps_sat)
    {
        pressure = pressure_bisectionsearch(p_max, p_min, alpha_rho_a, alpha_rho_b, gamma_a, pinf_a, gamma_b, pinf_b, epsinf_a, epsinf_b, eps_real);
    }

    else
    {
        // compute pressure outside saturation dome

        double rho_eps = rho * (energy - 0.5 * (v_x * v_x + v_y * v_y));

        double gamma_minus1 = 1 / ((alpha_a / (gamma_a - 1)) + (alpha_b / (gamma_b - 1)));

        double gamma_pinf = ((alpha_a * pinf_a * gamma_a) / (gamma_a - 1) + (alpha_b * pinf_b * gamma_b) / (gamma_b - 1));

        double alpha_epsinf = alpha_rho_a * epsinf_a + alpha_rho_b * epsinf_b; // reference internal energy

        pressure = gamma_minus1 * (rho_eps - gamma_pinf - alpha_epsinf);
    }
    }

    else
    {
        // compute pressure outside saturation dome

        double rho_eps = rho * (energy - 0.5 * (v_x * v_x + v_y * v_y));

        double gamma_minus1 = 1 / ((alpha_a / (gamma_a - 1)) + (alpha_b / (gamma_b - 1)));

        double gamma_pinf = ((alpha_a * pinf_a * gamma_a) / (gamma_a - 1) + (alpha_b * pinf_b * gamma_b) / (gamma_b - 1));

        double alpha_epsinf = alpha_rho_a * epsinf_a + alpha_rho_b * epsinf_b; // reference internal energy

        pressure = gamma_minus1 * (rho_eps - gamma_pinf - alpha_epsinf);
    }

    return pressure;
}

//defining compute total energy function
double eos::compute_energy(double alpha_a, double alpha_rho_a, double alpha_rho_b, double& pressure, double v_x, double v_y, double gamma_a, double pinf_a, double gamma_b, double pinf_b, double p_max, double epsinf_a, double epsinf_b)
{
    double energy;

    // logic for checking simulation's location on phase diagram
    double rho = alpha_rho_a + alpha_rho_b;

    double v = 1 / (alpha_rho_a + alpha_rho_b);

    Arrayofdouble eos_output;

    if(pressure<=p_max)
    {
        pressure = std::max(612.0, std::min(pressure, 22064000.0)); // making sure interpolation is possible

        eos_output = eos_interpolation(pressure);

        double vl_star = eos_output[1];
        double vg_star = eos_output[2];

        if(v>=vl_star && v<=vg_star)
        {
            // compute energy inside saturation dome 
            double lambda;

            lambda = compute_massfraction(alpha_rho_a, alpha_rho_b, pressure);

            double epsl_star = epsinf_a + vl_star * (pressure + gamma_a * pinf_a) / (gamma_a - 1);

            double epsg_star = epsinf_b + vg_star * (pressure + gamma_b * pinf_b) / (gamma_b - 1);

            double eps_pt = lambda * epsl_star + (1 - lambda) * epsg_star;

            energy = eps_pt + 0.5 * (v_x * v_x + v_y * v_y);
        }
        else
        {
            // compute energy outside saturation dome
            double alpha_b = 1.0 - alpha_a;

            double sum0 = alpha_rho_a * epsinf_a + alpha_rho_b * epsinf_b; // reference internal energy

            double sum1 = (alpha_a / (gamma_a - 1)) + (alpha_b / (gamma_b - 1));

            double sum2 = ((alpha_a * pinf_a * gamma_a) / (gamma_a - 1) + (alpha_b * pinf_b * gamma_b) / (gamma_b - 1));

            double rho_eps = sum0 + pressure * sum1 + sum2;

            energy = rho_eps / rho + 0.5 * (v_x * v_x + v_y * v_y);
        }
    }
    else
    {
        // compute energy outside saturation dome
        double alpha_b = 1.0 - alpha_a;

        double sum0 = alpha_rho_a * epsinf_a + alpha_rho_b * epsinf_b; // reference internal energy

        double sum1 = (alpha_a / (gamma_a - 1)) + (alpha_b / (gamma_b - 1));

        double sum2 = ((alpha_a * pinf_a * gamma_a) / (gamma_a - 1) + (alpha_b * pinf_b * gamma_b) / (gamma_b - 1));

        double rho_eps = sum0 + pressure * sum1 + sum2;

        energy = rho_eps / rho + 0.5 * (v_x * v_x + v_y * v_y);
    }

    return energy;
}

//defining compute sound speed function
double eos::compute_soundspeed(double alpha_a, double alpha_rho_a, double alpha_rho_b, double& pressure, double gamma_a, double pinf_a, double gamma_b, double pinf_b, double p_max, double epsinf_a, double epsinf_b, double cv_a, double cv_b)
{
    double c;

    // logic for checking simulation's location on phase diagram
    double rho = alpha_rho_a + alpha_rho_b;

    double v = 1 / (alpha_rho_a + alpha_rho_b);

    Arrayofdouble eos_output;

    if(pressure<=p_max)
    {
        pressure = std::max(612.0, std::min(pressure, 22064000.0)); // making sure interpolation is possible

        eos_output = eos_interpolation(pressure);

        double T_sat = eos_output[0];
        double vl_star = eos_output[1];
        double vg_star = eos_output[2];

        if (v >= vl_star && v <= vg_star)
        {
            // vl/p; vl/T; epsl/p; epsl/v
            double vl_p = - (gamma_a - 1) * cv_a * T_sat / pow((pressure + pinf_a), 2);
            double vl_T = (gamma_a - 1) * cv_a / (pressure + pinf_a);
            double epsl_p = vl_star / (gamma_a - 1);
            double epsl_v = (pressure + gamma_a * pinf_a) / (gamma_a - 1);

            // vg/p; vg/T; epsg/p; epsg/v
            double vg_p = - (gamma_b - 1) * cv_b * T_sat / pow((pressure + pinf_b), 2);
            double vg_T = (gamma_b - 1) * cv_b / (pressure + pinf_b);
            double epsg_p = vg_star / (gamma_b - 1);
            double epsg_v = (pressure + gamma_b * pinf_b) / (gamma_b - 1);

            // Tsatl/p
            double Tsatl_p = vl_star / ((gamma_a - 1) * cv_a);

            // Tsatg/p
            double Tsatg_p = vg_star / ((gamma_b - 1) * cv_b);

            // vl_star/p
            double vl_star_p = vl_T * Tsatl_p + vl_p;

            // vg_star/p
            double vg_star_p = vg_T * Tsatg_p + vg_p;

            // lambda
            double lambda;

            lambda = compute_massfraction(alpha_rho_a, alpha_rho_b, pressure);

            // epsl_star
            double epsl_star = epsinf_a + vl_star * (pressure + gamma_a * pinf_a) / (gamma_a - 1);

            // epsg_star
            double epsg_star = epsinf_b + vg_star * (pressure + gamma_b * pinf_b) / (gamma_b - 1);

            // lambda/p
            double lambda_p = (1/(vg_star - vl_star)) * (lambda * vl_star_p + (1 - lambda) * vg_star_p);

            // epsl_star/p
            double epsl_star_p = epsl_v * vl_star_p + epsl_p;

            // epsg_star/p
            double epsg_star_p = epsg_v * vg_star_p + epsg_p;

            // eps/p
            double eps_p = lambda * epsl_star_p + (1 - lambda) * epsg_star_p + lambda_p * (epsl_star - epsg_star);

            // eps/v
            double eps_v = (epsg_star - epsl_star) / (vg_star - vl_star);

            // sound speed
            c = sqrt(v*v*(eps_v + pressure) / eps_p);

        }
        else
        {
            double alpha_b = 1.0 - alpha_a;

            double rho_a = alpha_rho_a / alpha_a;
            double rho_b = alpha_rho_b / alpha_b;

            double y_a = alpha_rho_a / rho;
            double y_b = alpha_rho_b / rho;

            double xi_a = 1 / (gamma_a - 1); // reference internal energy no change as rho*epsinf is constant
            double xi_b = 1 / (gamma_b - 1); // reference internal energy no change as rho*epsinf is constant

            double xi = alpha_a * xi_a + alpha_b * xi_b; // alpha_a * xi_a + alpha_b * xi_b

            double c_a = gamma_a * (pressure + pinf_a) / rho_a; // reference internal energy no change as rho*epsinf is constant
            double c_b = gamma_b * (pressure + pinf_b) / rho_b; // reference internal energy no change as rho*epsinf is constant

            c = sqrt((y_a * xi_a * c_a + y_b * xi_b * c_b) / xi);

        }
    }
    else
    {
        double alpha_b = 1.0 - alpha_a;

        double rho_a = alpha_rho_a / alpha_a;
        double rho_b = alpha_rho_b / alpha_b;

        double y_a = alpha_rho_a / rho;
        double y_b = alpha_rho_b / rho;

        double xi_a = 1 / (gamma_a - 1); // reference internal energy no change as rho*epsinf is constant
        double xi_b = 1 / (gamma_b - 1); // reference internal energy no change as rho*epsinf is constant

        double xi = alpha_a * xi_a + alpha_b * xi_b; // alpha_a * xi_a + alpha_b * xi_b

        double c_a = gamma_a * (pressure + pinf_a) / rho_a; // reference internal energy no change as rho*epsinf is constant
        double c_b = gamma_b * (pressure + pinf_b) / rho_b; // reference internal energy no change as rho*epsinf is constant

        c = sqrt((y_a * xi_a * c_a + y_b * xi_b * c_b) / xi);
    }

    return c;
}



