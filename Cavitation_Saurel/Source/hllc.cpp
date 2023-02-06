// 2020H - File added which contains functions for HLLC solver
#include <iostream>
#include <cmath>
#include <vector>
#include <array>
#include <algorithm>

#include "eos_v2.H"
#include "hllc.H"
#include "AMReX_Array4.H"

eos* SGEOS = new eos();

using namespace amrex;

//constructor
numerical_method::numerical_method()
{}

void numerical_method::wavespeedestimate(Arrayofdouble const& conserv_l, Arrayofdouble const&  conserv_r, double gamma_a, double gamma_b, double pinf_a, double pinf_b, Arrayofdouble& wavespeed)
{
  // converting input data from conservative variable to primitive variable
  Arrayofdouble prim_l, prim_r;
  prim_l = SGEOS->u_to_prim(conserv_l, gamma_a, pinf_a, gamma_b, pinf_b);
  prim_r = SGEOS->u_to_prim(conserv_r, gamma_a, pinf_a, gamma_b, pinf_b);

  // defining primitive variable left
  double alpha_al = prim_l[0];
  double alpha_bl = 1 - prim_l[0];
  double rho_al = prim_l[1];
  double rho_bl = prim_l[2];
  double eps_al = prim_l[3];
  double eps_bl = prim_l[4];
  double p_al = SGEOS->epstopressure(rho_al, eps_al, pinf_a, gamma_a);
  double p_bl = SGEOS->epstopressure(rho_bl, eps_bl, pinf_b, gamma_b);
  double vx_l = prim_l[5];
  //double p_l = prim_l[6]; // mixture pressure left
  double p_l = alpha_al * p_al + alpha_bl * p_bl; // sm

  // defining primitive variable right 
  double alpha_ar = prim_r[0];
  double alpha_br = 1 - prim_r[0];
  double rho_ar = prim_r[1];
  double rho_br = prim_r[2];
  double eps_ar = prim_r[3];
  double eps_br = prim_r[4];
  double p_ar = SGEOS->epstopressure(rho_ar, eps_ar, pinf_a, gamma_a);
  double p_br = SGEOS->epstopressure(rho_br, eps_br, pinf_b, gamma_b);
  double vx_r = prim_r[5];
  //double p_r = prim_r[6]; // mixture pressure right
  double p_r = alpha_ar * p_ar + alpha_br * p_br; // sm

  // computing sound speed left from sound speed for material a and b
  double c_al = (gamma_a * (p_al + pinf_a) / rho_al); // negative density leads to negative c_l
  double c_bl = (gamma_b * (p_bl + pinf_b) / rho_bl);
  double rho_l = alpha_al * rho_al + alpha_bl * rho_bl; // total density for mass fraction calculation
  double y_al = alpha_al * rho_al / rho_l; // mass fraction for material a
  double y_bl = alpha_bl * rho_bl / rho_l; // mass fraction for material b
  double c_l = sqrt(y_al*c_al + y_bl*c_bl); // frozen mixture sound speed

  // computing sound speed right from sound speed for material a and b
  double c_ar = (gamma_a * (p_ar + pinf_a) / rho_ar);
  double c_br = (gamma_b * (p_br + pinf_b) / rho_br);
  double rho_r = alpha_ar * rho_ar + alpha_br * rho_br; // total density for mass fraction calculation
  double y_ar = alpha_ar * rho_ar / rho_r; // mass fraction for material a
  double y_br = alpha_br * rho_br / rho_r; // mass fraction for material b
  double c_r = sqrt(y_ar*c_ar + y_br*c_br); // frozen mixture sound speed

  // computing Sl and Sr
  double Sr = std::max(vx_l + c_l, vx_r + c_r);
  double Sl = std::min(vx_l - c_l, vx_r - c_r);

  // computing S_star
  double rho_vx_l = rho_l * vx_l;
  double rho_vx_r = rho_r * vx_r;
  double mom_flux_l = rho_l * vx_l * vx_l + p_l;
  double mom_flux_r = rho_r * vx_r * vx_r + p_r;
  double S_star = (mom_flux_l - mom_flux_r - Sl*rho_vx_l + Sr*rho_vx_r) / (rho_vx_l - rho_vx_r - Sl*rho_l + Sr*rho_r);

  wavespeed[0] = Sl;
  wavespeed[1] = Sr;
  wavespeed[2] = S_star;

  if(Sl > 0)
  {
    wavespeed[3] = vx_l;
  }
  else if(Sr < 0)
  {
    wavespeed[3] = vx_r;
  }
  else
  {
    wavespeed[3] = S_star;
  }
}

//function to find Minbee constant
Arrayofdouble numerical_method::Minbee(Arrayofdouble ui, Arrayofdouble uiMinus1, Arrayofdouble uiPlus1, int nVar)
{
  Arrayofdouble r;
  Arrayofdouble xi;
  Arrayofdouble Mb;

  double beta_plushalf = 1.0;//2.0 / (1.0 - 0.3);

  for(int i=0; i<nVar; i++)
  {
    if (fabs(uiPlus1[i] - ui[i])>0)
    {
      r[i] = (ui[i] - uiMinus1[i]) / (uiPlus1[i] - ui[i]);
    }
    else 
    {
      r[i] = 0;
    }
    // r[i] = 0;
  }

  for(int i=0; i<nVar; i++)
  {
    xi[i] = 2 * beta_plushalf / (1 + r[i]);
  }

  for(int i=0; i<nVar; i++)
  {
    if(r[i]<=0)
    {
      Mb[i] = 0;
    }
    else if(r[i]>0 && r[i]<=1)
    {
      Mb[i] = r[i];
    }
    else if(r[i] > 1)
    {
      if((2/(1+r[i])>1))
      {
        Mb[i] = 1;
      }
      else
      {
        Mb[i] = 2/(1+r[i]);
      }
    }
  }
  return Mb;
}

// function to find Minbee constant using mixture density
Arrayofdouble numerical_method::Minbee2(Arrayofdouble ui, Arrayofdouble uiMinus1, Arrayofdouble uiPlus1, int nVar)
{
  Arrayofdouble r; // output is still an array

  // defining mixture densities
  double u_i = ui[1] + ui[2];
  double u_iMinus1 = uiMinus1[1] + uiMinus1[2];
  double u_iPlus1 = uiPlus1[1] + uiPlus1[2];

  for(int i=0; i<nVar; i++)
  {
    if (fabs(uiPlus1[i] - ui[i])>0)
    {
      r[i] = (u_i - u_iMinus1) / (u_iPlus1 - u_i);
    }
    else 
    {
      r[i] = 0;
    }
    // r[i] = 0;
  }

  Arrayofdouble Mb;
  for(int i=0; i<nVar; i++)
  {
    if(r[i]<=0)
    {
      Mb[i] = 0;
    }
    else if(r[i]>0 && r[i]<=1)
    {
      Mb[i] = r[i];
    }
    else
    {
      if((2/(1+r[i])>1))
      {
        Mb[i] = 1;
      }
      else
      {
        Mb[i] = 2/(1+r[i]);
      }
    }
  }
  return Mb;
}

//function to find Vanleer constant
Arrayofdouble numerical_method::Vanleer(Arrayofdouble ui, Arrayofdouble uiMinus1, Arrayofdouble uiPlus1, int nVar)
{
  Arrayofdouble r;
  Arrayofdouble xi;
  Arrayofdouble xi_l;
  Arrayofdouble xi_r;
  Arrayofdouble Vl;

  double beta_plushalf = 1.0;

  for(int i=0; i<nVar; i++)
  {
    if (fabs(uiPlus1[i] - ui[i])>0)
    {
      r[i] = (ui[i] - uiMinus1[i]) / (uiPlus1[i] - ui[i]);
    }
    else 
    {
      r[i] = 0;
    }
    // r[i] = 0;
  }

  for(int i=0; i<nVar; i++)
  {
    {
      xi[i] = 2 * beta_plushalf / (1 + r[i]);
    }
  }

  for(int i=0; i<nVar; i++)
  {
    if(r[i]<=0)
    {
      Vl[i] = 0;
    }
    else if(r[i]>0)
    {
      if(xi[i]>(2*r[i]/(1+r[i])))
      {
        Vl[i] = 2*r[i]/(1+r[i]);
      }
      else 
      {
        Vl[i] = xi[i];
      }
    }
  }
  return Vl;
}

//function to find Vanleer constant
Arrayofdouble numerical_method::Superbee(Arrayofdouble ui, Arrayofdouble uiMinus1, Arrayofdouble uiPlus1, int nVar)
{
  Arrayofdouble r;
  Arrayofdouble xi;
  Arrayofdouble Sb_min1;
  Arrayofdouble Sb;

  double beta_plushalf = 1.0;//2.0 / (1.0 - 0.3);

  for(int i=0; i<nVar; i++)
  {
    if (fabs(uiPlus1[i] - ui[i])>0)
    {
      r[i] = (ui[i] - uiMinus1[i]) / (uiPlus1[i] - ui[i]);
    }
    else 
    {
      r[i] = 0;
    }
    // r[i] = 0;
  }

  for(int i=0; i<nVar; i++)
  {
    {
      xi[i] = 2 * beta_plushalf / (1 + r[i]);
    }
  }

  for(int i=0; i<nVar; i++)
  {
    if(r[i]<=0)
    {
      Sb[i] = 0;
    }
    else if(r[i]>0 && r[i]<0.5)
    {
      Sb[i] = 2 * r[i];
    }
    else if(r[i]>=0.5 && r[i]<=1.0)
    {
      Sb[i] = 1.0;
    }
    else if (r[i]>=1)
    {
      Sb_min1[i] = std::min(r[i], Sb[i]);
      Sb[i] = std::min(Sb_min1[i], 2.0);
    }
  }
  return Sb;
}

//computing nplushalf for xHLLC flux
void numerical_method::secondorder_extension(Arrayofdouble u_n, Arrayofdouble u, Arrayofdouble u_p, double gamma_a, double gamma_b, double pinf_a, double pinf_b, Arrayofdouble& primiL, Arrayofdouble& primiR, int nVar, double cell_size, double sl_const)
{
  Arrayofdouble primi, primi_minus1, primi_plus1, Vl_prim, Mb_prim, Sb_prim, Mb_prim2, delta_prim_plushalf, delta_prim_minushalf, delta_prim, Vl_cons, Mb_cons, delta_cons_plushalf, delta_cons_minushalf, delta_cons;

  // nVar passed in = 8 (NUM_STATE -7)
  // compute primitive variables from input 
  primi_minus1 = SGEOS->secondu_to_prim(u_n, gamma_a, pinf_a, gamma_b, pinf_b);
  primi = SGEOS->secondu_to_prim(u, gamma_a, pinf_a, gamma_b, pinf_b);
  primi_plus1 = SGEOS->secondu_to_prim(u_p, gamma_a, pinf_a, gamma_b, pinf_b);

  // gradient limitation
  // finding Mb using primitive variables
  Mb_prim = Minbee(primi, primi_minus1, primi_plus1, nVar);
  double Mb_prim_min = 1.0e5;
  for(int l=0; l<nVar; l++)
  {
    if(Mb_prim[l] < Mb_prim_min)
    {
      Mb_prim_min = Mb_prim[l];
    }
  }

  // finding Mb using conservative variables
  Mb_cons = Minbee(u, u_n, u_p, nVar);
  double Mb_cons_min = 1.0e5;
  for(int l=0; l<nVar; l++)
  {
    if(Mb_cons[l] < Mb_cons_min)
    {
      Mb_cons_min = Mb_cons[l];
    }
  }

  // finding Vl using primitive variables
  Vl_prim = Vanleer(primi, primi_minus1, primi_plus1, nVar);
  double Vl_prim_min = 1.0e5;
  for(int l=0; l<nVar; l++)
  {
    if(Vl_prim[l] < Vl_prim_min)
    {
      Vl_prim_min = Vl_prim[l];
    }
  }

  // finding Sb using primitive variables
  Sb_prim = Superbee(primi, primi_minus1, primi_plus1, nVar);
  double Sb_prim_min = 1e5;
  for(int l=0; l<nVar; l++)
  {
    if(Sb_prim[l] < Sb_prim_min)
    {
      Sb_prim_min = Sb_prim[l];
    }
  }

  // finding Vl using conservative variables
  Vl_cons = Vanleer(u, u_n, u_p, nVar);
  double Vl_cons_min = 1.0e5;
  for(int l=0; l<nVar; l++)
  {
    if(Vl_cons[l] < Vl_cons_min)
    {
      Vl_cons_min = Vl_cons[l];
    }
  }

  // define delta using primitive variables
  for(int l=0; l<nVar; l++)
  {
    delta_prim_minushalf[l] = (primi[l] - primi_minus1[l]) / cell_size;
    delta_prim_plushalf[l] = (primi_plus1[l] - primi[l]) / cell_size;
    delta_prim[l] = 0.5 * (1 + sl_const) * delta_prim_minushalf[l] + 0.5 * (1 - sl_const) * delta_prim_plushalf[l];
  }

  // define delta using conservative variables
  for(int l=0; l<nVar; l++)
  {
    delta_cons_minushalf[l] = (u[l] - u_n[l]) / cell_size;
    delta_cons_plushalf[l] = (u_p[l] - u[l]) / cell_size;
    delta_cons[l] = 0.5 * (1 + sl_const) * delta_cons_minushalf[l] + 0.5 * (1 - sl_const) * delta_cons_plushalf[l];
  }

  // variables extrapolation
  // define primiL, primiR
  for(int l=0; l<nVar; l++)
  {
    primiL[l] = primi[l] - 0.5 * Vl_prim_min * cell_size * delta_prim[l];
    primiR[l] = primi[l] + 0.5 * Vl_prim_min * cell_size * delta_prim[l];
  }

  // reconstruct convservative form for energy and insert into primiL[nVar-1], primiR[nVar-1]
    primiL[nVar-1] = u[nVar-1] - 0.5 * Vl_prim_min * cell_size * delta_cons[nVar-1];
    primiR[nVar-1] = u[nVar-1] + 0.5 * Vl_prim_min * cell_size * delta_cons[nVar-1];

}

void numerical_method::halftimestep_update(Arrayofdouble primi, Arrayofdouble primiL, Arrayofdouble primiR, Arrayofdouble yprimi, Arrayofdouble yprimiL, Arrayofdouble yprimiR, double dt, double dx, double dy, Arrayofdouble& uiL_nplushalf, Arrayofdouble& uiR_nplushalf, double gamma_a, double gamma_b, double pinf_a, double pinf_b)
{
  Arrayofdouble primiL_nplushalf, primiR_nplushalf;

  // calculate primiL_nplushalf
  primiL_nplushalf[0] = primiL[0] + 0.5 * dt / dx * (primi[5] * (primiL[0] - primiR[0])) + 0.5 * dt / dy * (yprimi[6] * (yprimiL[0] - yprimiR[0]));
  primiL_nplushalf[1] = primiL[1] + 0.5 * dt / dx * (primi[5] * (primiL[1] - primiR[1]) + primi[1] * (primiL[5] - primiR[5])) + 0.5 * dt / dy * (yprimi[6] * (yprimiL[1] - yprimiR[1]) + yprimi[1] * (yprimiL[6] - yprimiR[6]));
  primiL_nplushalf[2] = primiL[2] + 0.5 * dt / dx * (primi[5] * (primiL[2] - primiR[2]) + primi[2] * (primiL[5] - primiR[5])) + 0.5 * dt / dy * (yprimi[6] * (yprimiL[2] - yprimiR[2]) + yprimi[2] * (yprimiL[6] - yprimiR[6]));

  double c_al = (gamma_a * (primi[3] + pinf_a) / primi[1]); // sound speed for material a in x-direction
  double c_bl = (gamma_b * (primi[4] + pinf_b) / primi[2]); // sound speed for material b in x-direction
  double yc_al = (gamma_a * (yprimi[3] + pinf_a) / yprimi[1]); // sound speed for material a in y-direction
  double yc_bl = (gamma_b * (yprimi[4] + pinf_b) / yprimi[2]); // sound speed for material b in y-direction

  primiL_nplushalf[3] = primiL[3] + 0.5 * dt / dx * (primi[5] * (primiL[3] - primiR[3]) + primi[1] * c_al * (primiL[5] - primiR[5])) + 0.5 * dt / dy * (yprimi[6] * (yprimiL[3] - yprimiR[3]) + yprimi[1] * yc_al * (yprimiL[6] - yprimiR[6]));
  primiL_nplushalf[4] = primiL[4] + 0.5 * dt / dx * (primi[5] * (primiL[4] - primiR[4]) + primi[2] * c_bl * (primiL[5] - primiR[5])) + 0.5 * dt / dy * (yprimi[6] * (yprimiL[4] - yprimiR[4]) + yprimi[2] * yc_bl * (yprimiL[6] - yprimiR[6]));

  double densityL = primi[0] * primi[1] + (1 - primi[0]) * primi[2]; // computing mixture density in x-direction
  double alpha_oneL = (primi[3] - primi[4]) / densityL * (primiL[0] - primiR[0]); // terms for vx reconstruction in x-direction
  double p_oneL = primi[0] / densityL * (primiL[3] - primiR[3]); // terms for vx reconstruction in x-direction
  double p_twoL = (1 - primi[0]) / densityL * (primiL[4] - primiR[4]); // terms for vx reconstruction in x-direction
  primiL_nplushalf[5] = primiL[5] + 0.5 * dt / dx * (primi[5] * (primiL[5] - primiR[5]) + alpha_oneL + p_oneL + p_twoL) + 0.5 * dt / dy * (yprimi[6] * (yprimiL[5] - yprimiR[5]));
 
  double ydensityL = yprimi[0] * yprimi[1] + (1 - yprimi[0]) * yprimi[2]; // computing mixture density in y-direction
  double yalpha_oneL = (yprimi[3] - yprimi[4]) / ydensityL * (yprimiL[0] - yprimiR[0]); // terms for vx reconstruction in y-direction
  double yp_oneL = yprimi[0] / ydensityL * (yprimiL[3] - yprimiR[3]); // terms for vx reconstruction in y-direction
  double yp_twoL = (1 - yprimi[0]) / ydensityL * (yprimiL[4] - yprimiR[4]); // terms for vx reconstruction in y-direction
  primiL_nplushalf[6] = primiL[6] + 0.5 * dt / dx * primi[5] * (primiL[6] - primiR[6]) + 0.5 * dt / dy * (yprimi[6] * (yprimiL[6] - yprimiR[6]) + yalpha_oneL + yp_oneL + yp_twoL);
  primiL_nplushalf[7] = primi[7];

  // calculate primiR_nplushalf
  primiR_nplushalf[0] = primiR[0] + 0.5 * dt / dx * (primi[5] * (primiL[0] - primiR[0])) + 0.5 * dt / dy * (yprimi[6] * (yprimiL[0] - yprimiR[0]));
  primiR_nplushalf[1] = primiR[1] + 0.5 * dt / dx * (primi[5] * (primiL[1] - primiR[1]) + primi[1] * (primiL[5] - primiR[5])) + 0.5 * dt / dy * (yprimi[6] * (yprimiL[1] - yprimiR[1]) + yprimi[1] * (yprimiL[6] - yprimiR[6]));
  primiR_nplushalf[2] = primiR[2] + 0.5 * dt / dx * (primi[5] * (primiL[2] - primiR[2]) + primi[2] * (primiL[5] - primiR[5])) + 0.5 * dt / dy * (yprimi[6] * (yprimiL[2] - yprimiR[2]) + yprimi[2] * (yprimiL[6] - yprimiR[6]));

  double c_ar = (gamma_a * (primi[3] + pinf_a) / primi[1]);
  double c_br = (gamma_b * (primi[4] + pinf_b) / primi[2]);
  double yc_ar = (gamma_a * (yprimi[3] + pinf_a) / yprimi[1]);
  double yc_br = (gamma_b * (yprimi[4] + pinf_b) / yprimi[2]);

  primiR_nplushalf[3] = primiR[3] + 0.5 * dt / dx * (primi[5] * (primiL[3] - primiR[3]) + primi[1] * c_ar * (primiL[5] - primiR[5])) + 0.5 * dt / dy * (yprimi[6] * (yprimiL[3] - yprimiR[3]) + yprimi[1] * yc_ar * (yprimiL[6] - yprimiR[6]));
  primiR_nplushalf[4] = primiR[4] + 0.5 * dt / dx * (primi[5] * (primiL[4] - primiR[4]) + primi[2] * c_br * (primiL[5] - primiR[5])) + 0.5 * dt / dy * (yprimi[6] * (yprimiL[4] - yprimiR[4]) + yprimi[2] * yc_br * (yprimiL[6] - yprimiR[6]));

  double densityR = primi[0] * primi[1] + (1 - primi[0]) * primi[2]; // computing mixture density
  double alpha_oneR = (primi[3] - primi[4]) / densityR * (primiL[0] - primiR[0]); // terms for vx reconstruction
  double p_oneR = primi[0] / densityR * (primiL[3] - primiR[3]); // terms for vx reconstruction
  double p_twoR = (1 - primi[0]) / densityR * (primiL[4] - primiR[4]); // terms for vx reconstruction
  primiR_nplushalf[5] = primiR[5] + 0.5 * dt / dx * (primi[5] * (primiL[5] - primiR[5]) + alpha_oneR + p_oneR + p_twoR) + 0.5 * dt / dy * (yprimi[6] * (yprimiL[5] - yprimiR[5]));

  double ydensityR = yprimi[0] * yprimi[1] + (1 - yprimi[0]) * yprimi[2]; // computing mixture density
  double yalpha_oneR = (yprimi[3] - yprimi[4]) / ydensityR * (yprimiL[0] - yprimiR[0]); // terms for vx reconstruction
  double yp_oneR = yprimi[0] / ydensityR * (yprimiL[3] - yprimiR[3]); // terms for vx reconstruction
  double yp_twoR = (1 - yprimi[0]) / ydensityR * (yprimiL[4] - yprimiR[4]); // terms for vx reconstruction
  primiR_nplushalf[6] = primiR[6] + 0.5 * dt / dx * primi[5] * (primiL[6] - primiR[6]) + 0.5 * dt / dy * (yprimi[6] * (yprimiL[6] - yprimiR[6]) + yalpha_oneR + yp_oneR + yp_twoR);
  primiR_nplushalf[7] = primi[7];

  // compute conservative variables for output
  uiL_nplushalf = SGEOS->secondprim_to_u(primiL_nplushalf, gamma_a, pinf_a, gamma_b, pinf_b);
  uiR_nplushalf = SGEOS->secondprim_to_u(primiR_nplushalf, gamma_a, pinf_a, gamma_b, pinf_b);
}

void numerical_method::yhalftimestep_update(Arrayofdouble primi, Arrayofdouble primiL, Arrayofdouble primiR, Arrayofdouble yprimi, Arrayofdouble yprimiL, Arrayofdouble yprimiR, double dt, double dx, double dy, Arrayofdouble& uiL_nplushalf, Arrayofdouble& uiR_nplushalf, double gamma_a, double gamma_b, double pinf_a, double pinf_b)
{
  Arrayofdouble primiL_nplushalf, primiR_nplushalf;

  // calculate primiL_nplushalf
  primiL_nplushalf[0] = yprimiL[0] + 0.5 * dt / dx * (primi[5] * (primiL[0] - primiR[0])) + 0.5 * dt / dy * (yprimi[6] * (yprimiL[0] - yprimiR[0]));
  primiL_nplushalf[1] = yprimiL[1] + 0.5 * dt / dx * (primi[5] * (primiL[1] - primiR[1]) + primi[1] * (primiL[5] - primiR[5])) + 0.5 * dt / dy * (yprimi[6] * (yprimiL[1] - yprimiR[1]) + yprimi[1] * (yprimiL[6] - yprimiR[6]));
  primiL_nplushalf[2] = yprimiL[2] + 0.5 * dt / dx * (primi[5] * (primiL[2] - primiR[2]) + primi[2] * (primiL[5] - primiR[5])) + 0.5 * dt / dy * (yprimi[6] * (yprimiL[2] - yprimiR[2]) + yprimi[2] * (yprimiL[6] - yprimiR[6]));

  double c_al = (gamma_a * (primi[3] + pinf_a) / primi[1]); // sound speed for material a in x-direction
  double c_bl = (gamma_b * (primi[4] + pinf_b) / primi[2]); // sound speed for material b in x-direction
  double yc_al = (gamma_a * (yprimi[3] + pinf_a) / yprimi[1]); // sound speed for material a in y-direction
  double yc_bl = (gamma_b * (yprimi[4] + pinf_b) / yprimi[2]); // sound speed for material b in y-direction

  primiL_nplushalf[3] = yprimiL[3] + 0.5 * dt / dx * (primi[5] * (primiL[3] - primiR[3]) + primi[1] * c_al * (primiL[5] - primiR[5])) + 0.5 * dt / dy * (yprimi[6] * (yprimiL[3] - yprimiR[3]) + yprimi[1] * yc_al * (yprimiL[6] - yprimiR[6]));
  primiL_nplushalf[4] = yprimiL[4] + 0.5 * dt / dx * (primi[5] * (primiL[4] - primiR[4]) + primi[2] * c_bl * (primiL[5] - primiR[5])) + 0.5 * dt / dy * (yprimi[6] * (yprimiL[4] - yprimiR[4]) + yprimi[2] * yc_bl * (yprimiL[6] - yprimiR[6]));

  double densityL = primi[0] * primi[1] + (1 - primi[0]) * primi[2]; // computing mixture density in x-direction
  double alpha_oneL = (primi[3] - primi[4]) / densityL * (primiL[0] - primiR[0]); // terms for vx reconstruction in x-direction
  double p_oneL = primi[0] / densityL * (primiL[3] - primiR[3]); // terms for vx reconstruction in x-direction
  double p_twoL = (1 - primi[0]) / densityL * (primiL[4] - primiR[4]); // terms for vx reconstruction in x-direction
  primiL_nplushalf[5] = yprimiL[5] + 0.5 * dt / dx * (primi[5] * (primiL[5] - primiR[5]) + alpha_oneL + p_oneL + p_twoL) + 0.5 * dt / dy * (yprimi[6] * (yprimiL[5] - yprimiR[5]));
 
  double ydensityL = yprimi[0] * yprimi[1] + (1 - yprimi[0]) * yprimi[2]; // computing mixture density in y-direction
  double yalpha_oneL = (yprimi[3] - yprimi[4]) / ydensityL * (yprimiL[0] - yprimiR[0]); // terms for vx reconstruction in y-direction
  double yp_oneL = yprimi[0] / ydensityL * (yprimiL[3] - yprimiR[3]); // terms for vx reconstruction in y-direction
  double yp_twoL = (1 - yprimi[0]) / ydensityL * (yprimiL[4] - yprimiR[4]); // terms for vx reconstruction in y-direction
  primiL_nplushalf[6] = yprimiL[6] + 0.5 * dt / dx * primi[5] * (primiL[6] - primiR[6]) + 0.5 * dt / dy * (yprimi[6] * (yprimiL[6] - yprimiR[6]) + yalpha_oneL + yp_oneL + yp_twoL);
  primiL_nplushalf[7] = primi[7];

  // calculate primiR_nplushalf
  primiR_nplushalf[0] = yprimiR[0] + 0.5 * dt / dx * (primi[5] * (primiL[0] - primiR[0])) + 0.5 * dt / dy * (yprimi[6] * (yprimiL[0] - yprimiR[0]));
  primiR_nplushalf[1] = yprimiR[1] + 0.5 * dt / dx * (primi[5] * (primiL[1] - primiR[1]) + primi[1] * (primiL[5] - primiR[5])) + 0.5 * dt / dy * (yprimi[6] * (yprimiL[1] - yprimiR[1]) + yprimi[1] * (yprimiL[6] - yprimiR[6]));
  primiR_nplushalf[2] = yprimiR[2] + 0.5 * dt / dx * (primi[5] * (primiL[2] - primiR[2]) + primi[2] * (primiL[5] - primiR[5])) + 0.5 * dt / dy * (yprimi[6] * (yprimiL[2] - yprimiR[2]) + yprimi[2] * (yprimiL[6] - yprimiR[6]));

  double c_ar = (gamma_a * (primi[3] + pinf_a) / primi[1]);
  double c_br = (gamma_b * (primi[4] + pinf_b) / primi[2]);
  double yc_ar = (gamma_a * (yprimi[3] + pinf_a) / yprimi[1]);
  double yc_br = (gamma_b * (yprimi[4] + pinf_b) / yprimi[2]);

  primiR_nplushalf[3] = yprimiR[3] + 0.5 * dt / dx * (primi[5] * (primiL[3] - primiR[3]) + primi[1] * c_ar * (primiL[5] - primiR[5])) + 0.5 * dt / dy * (yprimi[6] * (yprimiL[3] - yprimiR[3]) + yprimi[1] * yc_ar * (yprimiL[6] - yprimiR[6]));
  primiR_nplushalf[4] = yprimiR[4] + 0.5 * dt / dx * (primi[5] * (primiL[4] - primiR[4]) + primi[2] * c_br * (primiL[5] - primiR[5])) + 0.5 * dt / dy * (yprimi[6] * (yprimiL[4] - yprimiR[4]) + yprimi[2] * yc_br * (yprimiL[6] - yprimiR[6]));

  double densityR = primi[0] * primi[1] + (1 - primi[0]) * primi[2]; // computing mixture density
  double alpha_oneR = (primi[3] - primi[4]) / densityR * (primiL[0] - primiR[0]); // terms for vx reconstruction
  double p_oneR = primi[0] / densityR * (primiL[3] - primiR[3]); // terms for vx reconstruction
  double p_twoR = (1 - primi[0]) / densityR * (primiL[4] - primiR[4]); // terms for vx reconstruction
  primiR_nplushalf[5] = yprimiR[5] + 0.5 * dt / dx * (primi[5] * (primiL[5] - primiR[5]) + alpha_oneR + p_oneR + p_twoR) + 0.5 * dt / dy * (yprimi[6] * (yprimiL[5] - yprimiR[5]));

  double ydensityR = yprimi[0] * yprimi[1] + (1 - yprimi[0]) * yprimi[2]; // computing mixture density
  double yalpha_oneR = (yprimi[3] - yprimi[4]) / ydensityR * (yprimiL[0] - yprimiR[0]); // terms for vx reconstruction
  double yp_oneR = yprimi[0] / ydensityR * (yprimiL[3] - yprimiR[3]); // terms for vx reconstruction
  double yp_twoR = (1 - yprimi[0]) / ydensityR * (yprimiL[4] - yprimiR[4]); // terms for vx reconstruction
  primiR_nplushalf[6] = yprimiR[6] + 0.5 * dt / dx * primi[5] * (primiL[6] - primiR[6]) + 0.5 * dt / dy * (yprimi[6] * (yprimiL[6] - yprimiR[6]) + yalpha_oneR + yp_oneR + yp_twoR);
  primiR_nplushalf[7] = primi[7];

  // compute conservative variables for output
  uiL_nplushalf = SGEOS->secondprim_to_u(primiL_nplushalf, gamma_a, pinf_a, gamma_b, pinf_b);
  uiR_nplushalf = SGEOS->secondprim_to_u(primiR_nplushalf, gamma_a, pinf_a, gamma_b, pinf_b);
}

void numerical_method::compute_HLLCflux(amrex::Array4<amrex::Real> const& arr, amrex::Array4<amrex::Real> const& arr_xrecon_lht, amrex::Array4<amrex::Real> const& arr_xrecon_rht, Arrayofdouble& Fhllc, Arrayofdouble& wavespeed, int nVar, double gamma_a, double pinf_a, double gamma_b, double pinf_b, int i, int j, int k, double dx, double dt)
{
  //defining left and right reconstructed conservative variables for HLLC solver
  Arrayofdouble conserv_l, conserv_r;
  for(int l=0; l<nVar; l++)
  {
    // conserv_l[l] = arr_xrecon_rht(i-1,j,k,l);
    // conserv_r[l] = arr_xrecon_lht(i,j,k,l);

    conserv_l[l] = arr(i-1,j,k,l);
    conserv_r[l] = arr(i,j,k,l);
  }

  // calling the compute wavespeed function
  wavespeedestimate(conserv_l, conserv_r, gamma_a, gamma_b, pinf_a, pinf_b, wavespeed);

  // define input variables 
  double Sl = wavespeed[0];
  double Sr = wavespeed[1];
  double S_star = wavespeed[2];

  // converting input data from conservative variable to primitive variable
  Arrayofdouble prim_l, prim_r;
  prim_l = SGEOS->u_to_prim(conserv_l, gamma_a, pinf_a, gamma_b, pinf_b);
  prim_r = SGEOS->u_to_prim(conserv_r, gamma_a, pinf_a, gamma_b, pinf_b);

  // defining primitive variable left
  double alpha_al = prim_l[0];
  double alpha_bl = 1 - prim_l[0];
  double rho_al = prim_l[1];
  double rho_bl = prim_l[2];
  double eps_al = prim_l[3];
  double eps_bl = prim_l[4];
  double p_al = SGEOS->epstopressure(rho_al, eps_al, pinf_a, gamma_a);
  double p_bl = SGEOS->epstopressure(rho_bl, eps_bl, pinf_b, gamma_b);
  double vx_l = prim_l[5];
  //double p_l = prim_l[6];
  double p_l = alpha_al * p_al + alpha_bl * p_bl; // sm
  double vy_l = prim_l[7];

  // defining primitive variable right 
  double alpha_ar = prim_r[0];
  double alpha_br = 1 - prim_r[0];
  double rho_ar = prim_r[1];
  double rho_br = prim_r[2];
  double eps_ar = prim_r[3];
  double eps_br = prim_r[4];
  double p_ar = SGEOS->epstopressure(rho_ar, eps_ar, pinf_a, gamma_a);
  double p_br = SGEOS->epstopressure(rho_br, eps_br, pinf_b, gamma_b);
  double vx_r = prim_r[5];
  //double p_r = prim_r[6];
  double p_r = alpha_ar * p_ar + alpha_br * p_br; // sm
  double vy_r = prim_r[7];

  // computing flux left and right 
  Arrayofdouble Fl,Fr;
  Fl = SGEOS->fluxf(conserv_l, prim_l);
  Fr = SGEOS->fluxf(conserv_r, prim_r);

  // computing Ql*
  Arrayofdouble Ql_star, Fl_star;

  //alpha_a
  Ql_star[0] = alpha_al; 

  //alpha_a * rho_a
  Ql_star[1] = alpha_al * rho_al * (Sl - vx_l) / (Sl - S_star); 

  //alpha_b * rho_b
  Ql_star[2] = alpha_bl * rho_bl * (Sl - vx_l) / (Sl - S_star); 

  // computing intermediate mixture density
  double rho_alstar = rho_al * (vx_l - Sl) / (S_star - Sl); // fluid density for material a
  double rho_blstar = rho_bl * (vx_l - Sl) / (S_star - Sl); // fluid density for material b 
  double rho_lstar = alpha_al * rho_alstar + alpha_bl * rho_blstar; // mixture fluid density

  //alpha_a * rho_a * eps_xa 
  double gamma_p_rho_al = (gamma_a + 1) * rho_al;
  double gamma_m_rho_al = (gamma_a - 1) * rho_al;
  double gamma_p_rho_alstar = (gamma_a + 1) * rho_alstar;
  double gamma_m_rho_alstar = (gamma_a - 1) * rho_alstar;
  //double p_alstar = (p_al + pinf_a) * (gamma_m_rho_al - gamma_p_rho_alstar) / (gamma_m_rho_alstar - gamma_p_rho_al) - pinf_a ; // phasic pressure for material a
  double p_alstar = p_al + rho_al * (Sl - vx_l) * (S_star - vx_l); // sm
  double eps_alstar = (p_alstar + gamma_a * pinf_a) / ((gamma_a - 1) * rho_alstar); // computing intermediate internal energy for material a left from eos
  Ql_star[3] = alpha_al * rho_alstar * eps_alstar;

  //alpha_b * rho_b * eps_xb 
  double gamma_p_rho_bl = (gamma_b + 1) * rho_bl;
  double gamma_m_rho_bl = (gamma_b - 1) * rho_bl;
  double gamma_p_rho_blstar = (gamma_b + 1) * rho_blstar;
  double gamma_m_rho_blstar = (gamma_b - 1) * rho_blstar;
  //double p_blstar = (p_bl + pinf_b) * (gamma_m_rho_bl - gamma_p_rho_blstar) / (gamma_m_rho_blstar - gamma_p_rho_bl) - pinf_b; // phasic pressure for material b
  double p_blstar = p_bl + rho_bl * (Sl - vx_l) * (S_star - vx_l); // sm
  double eps_blstar = (p_blstar + gamma_b * pinf_b) / ((gamma_b - 1) * rho_blstar); // computing intermediate internal energy for material b left from eos
  Ql_star[4] = alpha_bl * rho_blstar * eps_blstar;

  //rho * v_x
  double rho_l = alpha_al * rho_al + alpha_bl * rho_bl; // computing mixture density left
  Ql_star[5] = rho_lstar * S_star;

  // rho * energy
  double p_lstar = p_l + rho_l * vx_l * (vx_l - Sl) - rho_lstar * S_star * (S_star - Sl); // computing intermediate pressure left
  double y_al = alpha_al * rho_al / rho_l; // mass fraction for material a left
  double y_bl = alpha_bl * rho_bl / rho_l; // mass fraction for material b left
  //double energy_l = y_al * eps_al + y_bl * eps_bl + 0.5 * vx_l * vx_l; // computing total energy left 
  double energy_l = SGEOS->compute_energy(alpha_al, rho_al, eps_al, rho_bl, eps_bl, vx_l, vy_l);
  double energy_lstar = (rho_l * energy_l * (vx_l - Sl) + p_l * vx_l - p_lstar * S_star) / (rho_lstar * (S_star - Sl)); // computing total intermediate energy left
  //Ql_star[6] = rho_lstar * energy_lstar;
  Ql_star[6] = (Sl - vx_l) / (Sl - S_star) * (rho_l * energy_l + (S_star - vx_l) * (S_star * rho_l + p_l / (Sl - vx_l)));

  // rho * v_y
  Ql_star[7] = rho_lstar * vy_l;

  // computing Qr*
  Arrayofdouble Qr_star, Fr_star;

  //alpha_a
  Qr_star[0] = alpha_ar; 

  //alpha_a * rho_a
  Qr_star[1] = alpha_ar * rho_ar * (Sr - vx_r) / (Sr - S_star); 

  //alpha_b * rho_b
  Qr_star[2] = alpha_br * rho_br * (Sr - vx_r) / (Sr - S_star); 

  // computing intermediate mixture density
  double rho_arstar = rho_ar * (vx_r - Sr) / (S_star - Sr); // fluid density for material a
  double rho_brstar = rho_br * (vx_r - Sr) / (S_star - Sr); // fluid density for material b 
  double rho_rstar = alpha_ar * rho_arstar + alpha_br * rho_brstar; // mixture fluid density

  //alpha_a * rho_a * eps_xa 
  double gamma_p_rho_ar = (gamma_a + 1) * rho_ar;
  double gamma_m_rho_ar = (gamma_a - 1) * rho_ar;
  double gamma_p_rho_arstar = (gamma_a + 1) * rho_arstar;
  double gamma_m_rho_arstar = (gamma_a - 1) * rho_arstar;
  //double p_arstar = (p_ar + pinf_a) * (gamma_m_rho_ar - gamma_p_rho_arstar) / (gamma_m_rho_arstar - gamma_p_rho_ar) - pinf_a; // phasic pressure for material a
  double p_arstar = p_ar + rho_ar * (Sr - vx_r) * (S_star - vx_r); // sm
  double eps_arstar = (p_arstar + gamma_a * pinf_a) / ((gamma_a - 1) * rho_arstar); // computing intermediate internal energy for material a left from eos
  Qr_star[3] = alpha_ar * rho_arstar * eps_arstar; 

  //alpha_b * rho_b * eps_xb 
  double gamma_p_rho_br = (gamma_b + 1) * rho_br;
  double gamma_m_rho_br = (gamma_b - 1) * rho_br;
  double gamma_p_rho_brstar = (gamma_b + 1) * rho_brstar;
  double gamma_m_rho_brstar = (gamma_b - 1) * rho_brstar;
  //double p_brstar = (p_br + pinf_b) * (gamma_m_rho_br - gamma_p_rho_brstar) / (gamma_m_rho_brstar - gamma_p_rho_br) - pinf_b; // phasic pressure for material b
  double p_brstar = p_br + rho_br * (Sr - vx_r) * (S_star - vx_r); // sm
  double eps_brstar = (p_brstar + gamma_b * pinf_b) / ((gamma_b - 1) * rho_brstar); // computing intermediate internal energy for material b left from eos
  Qr_star[4] = alpha_br * rho_brstar * eps_brstar; 

  //rho * v_x
  double rho_r = alpha_ar * rho_ar + alpha_br * rho_br; // computing mixture density right
  Qr_star[5] = rho_rstar * S_star;

  // rho * energy
  double p_rstar = p_r + rho_r * vx_r * (vx_r - Sr) - rho_rstar * S_star * (S_star - Sr); // computing intermediate pressure right
  double y_ar = alpha_ar * rho_ar / rho_r; // mass fraction for material a left
  double y_br = alpha_br * rho_br / rho_r; // mass fraction for material b left
  //double energy_r = y_ar * eps_ar + y_br * eps_br + 0.5 * vx_r * vx_r; // computing total energy left 
  double energy_r = SGEOS->compute_energy(alpha_ar, rho_ar, eps_ar, rho_br, eps_br, vx_r, vy_r);
  double energy_rstar = (rho_r * energy_r * (vx_r - Sr) + p_r * vx_r - p_rstar * S_star) / (rho_rstar * (S_star - Sr)); // computing total intermediate energy left
  //Qr_star[6] = rho_rstar * energy_rstar;
  //Qr_star[6] = rho_rstar * (energy_r + (S_star - vx_r) * (S_star * rho_r + p_r / (Sr - vx_r)));
  Qr_star[6] = (Sr - vx_r) / (Sr - S_star) * (rho_r * energy_r + (S_star - vx_r) * (S_star * rho_r + p_r / (Sr - vx_r)));

  // rho * v_y
  Qr_star[7] = rho_rstar * vy_r;

  // computing Fl* and Fr*
  for(int l=0; l<nVar; l++)
  {
  Fl_star[l] = Fl[l]  + Sl*(Ql_star[l] - conserv_l[l]); 
  Fr_star[l] = Fr[l]  + Sr*(Qr_star[l] - conserv_r[l]);
  }

  // implementing logic for choosing flux
	if(Sl >= 0)
	{	
		Fhllc = Fl;
    Fhllc[0] = conserv_l[0];
    Fhllc[3] = conserv_l[3];
    Fhllc[4] = conserv_l[4];
		return; 
	}

	else if(S_star >= 0)
	{
		Fhllc = Fl_star;
    Fhllc[0] = Ql_star[0];
    Fhllc[3] = Ql_star[3];
    Fhllc[4] = Ql_star[4];
		return;
	}
	
	else if(Sr > 0)
	{
		Fhllc = Fr_star;
    Fhllc[0] = Qr_star[0];
    Fhllc[3] = Qr_star[3];
    Fhllc[4] = Qr_star[4];
		return;
	}

	else
	{	Fhllc = Fr;
    Fhllc[0] = conserv_r[0];
    Fhllc[3] = conserv_r[3];
    Fhllc[4] = conserv_r[4];
		return; 
	}		
}

void numerical_method::pressure_relaxation(amrex::Array4<amrex::Real> const arr, int i, int j, int k, double gamma_a, double pinf_a, double gamma_b, double pinf_b, Arrayofdouble& relaxed_terms)
{
  // defining some initial variables
  double alpha_a_init = arr(i,j,k,0); 
  double alpha_b_init = 1.0 - alpha_a_init;  
  double rho_a_init = arr(i,j,k,1) / alpha_a_init;
  double rho_b_init = arr(i,j,k,2) / alpha_b_init;
  double eps_a_init = arr(i,j,k,3) / (alpha_a_init * rho_a_init);
  double eps_b_init = arr(i,j,k,4) / (alpha_b_init * rho_b_init);

  double vol_a_init = 1 / rho_a_init;
  double vol_b_init = 1 / rho_b_init;

  double p_a_init = SGEOS->epstopressure(rho_a_init, eps_a_init, pinf_a, gamma_a);
  double p_b_init = SGEOS->epstopressure(rho_b_init, eps_b_init, pinf_b, gamma_b);

  double density = arr(i,j,k,1) + arr(i,j,k,2);

  // volume fraction for material a
  double A_a = arr(i,j,k,1) * vol_a_init * (p_a_init + gamma_a * pinf_a);
  double B_a = arr(i,j,k,1) * vol_a_init * (gamma_a - 1);
  double C_a = gamma_a * pinf_a;

  // volume fraction for material b 
  double A_b = arr(i,j,k,2) * vol_b_init * (p_b_init + gamma_b * pinf_b);
  double B_b = arr(i,j,k,2) * vol_b_init * (gamma_b - 1);
  double C_b = gamma_b * pinf_b;

  // solving the quadratic equation to find relaxed pressure 
  double D = B_a * gamma_b + B_b * gamma_a - gamma_a * gamma_b; // 'a' for quadratic formula
  double E = A_a * gamma_b + B_a * C_b + A_b * gamma_a + B_b * C_a - C_a * gamma_b - C_b * gamma_a; // 'b' for quadratic formula
  double F = A_a * C_b + A_b * C_a - C_a * C_b; // 'c' for quadratic formula

  double p_relaxed = (- E - sqrt(E*E - 4*D*F)) / (2*D); 
  // p_relaxed = std::max(p_relaxed, 1.0e-5);

  double v_x = arr(i,j,k,5) / density;
  double energy = arr(i,j,k,6) / density; 
  double v_y = arr(i,j,k,7) / density;

  // correcting the volume fractions for material a using relaxed pressure
  double vol_a_new = (A_a + B_a * p_relaxed) / ((C_a + gamma_a * p_relaxed) * arr(i,j,k,1));
  double alpha_a_new = arr(i,j,k,1) * vol_a_new; // new volume fraction for material a using specific volume computed from relaxed pressure

  // correcting the volume fractions for material b using relaxed pressure
  double vol_b_new = (A_b + B_b * p_relaxed) / ((C_b + gamma_b * p_relaxed) * arr(i,j,k,2)); 
  //SM
  double alpha_b_new = 1 - alpha_a_new;// new volume fraction for material b using specific volume computed from relaxed pressure

  // computing the mixture pressure 
  double alpha_pinf_sum = (alpha_a_new * gamma_a * pinf_a) / (gamma_a - 1) + (alpha_b_new * gamma_b * pinf_b) / (gamma_b - 1);
  double alpha_sum = alpha_a_new / (gamma_a - 1) + alpha_b_new / (gamma_b - 1);
  double density_new = alpha_a_new / vol_a_new + alpha_b_new / vol_b_new;
  double p_mix_new = (density * energy - 0.5 * density * (v_x * v_x + v_y * v_y) - alpha_pinf_sum) / alpha_sum;
  // p_mix_new = std::max(p_mix_new, 1.0e-5);

  // correcting the volume fraction for material a
  // relaxed_terms[0] = alpha_a_new;
  relaxed_terms[0] = std::min(alpha_a_new, 0.999999);

  // correcting the volume fraction * density * internal energy for material a using mixture pressure and SG eos
  double eps_a_new = (p_mix_new + gamma_a * pinf_a) * vol_a_new / (gamma_a - 1);
  relaxed_terms[1] = arr(i,j,k,1) * eps_a_new;

  // correcting the volume fraction * density * internal energy for material b using mixture pressure and SG eos
  double eps_b_new = (p_mix_new + gamma_b * pinf_b) * vol_b_new / (gamma_b - 1);
  relaxed_terms[2] = arr(i,j,k,2) * eps_b_new;
}

void numerical_method::pressure_relaxation2(amrex::Array4<amrex::Real> const& arr, int i, int j, int k, double gamma_a, double pinf_a, double gamma_b, double pinf_b, Arrayofdouble& relaxed_terms, int nVar)
{
  // defining some initial variables
  double alpha_a = arr(i,j,k,0); 
  double alpha_b = 1.0 - alpha_a;  
  double rho_a = arr(i,j,k,1) / alpha_a;
  double rho_b = arr(i,j,k,2) / alpha_b;
  double eps_a = arr(i,j,k,3) / (alpha_a * rho_a);
  double eps_b = arr(i,j,k,4) / (alpha_b * rho_b);

  double vol_a = 1 / rho_a;
  double vol_b = 1 / rho_b;

  double p_a = SGEOS->epstopressure(rho_a, eps_a, pinf_a, gamma_a);
  double p_b = SGEOS->epstopressure(rho_b, eps_b, pinf_b, gamma_b);

  double density = arr(i,j,k,1) + arr(i,j,k,2);

  // computing interface pressure
  double c_a = gamma_a * (p_a + pinf_a) / rho_a;
  double c_b = gamma_b * (p_b + pinf_b) / rho_b;

  double p_i = (rho_b * c_a * p_a + rho_a * c_b * p_b) / (rho_b * c_a + rho_a * c_b);

  // finding mixture pressure by solving the quadratic equation 
  double C_a = 2 * gamma_a * pinf_a + (gamma_a - 1) * p_i;
  double C_b = 2 * gamma_b * pinf_b + (gamma_b - 1) * p_i;
  double A = 1 + gamma_b * alpha_a + gamma_a * alpha_b;
  double B = C_a * alpha_b + C_b * alpha_a - (1 + gamma_b) * alpha_a * p_a - (1 + gamma_a) * alpha_b * p_b;
  double D = - C_b * alpha_a * p_a - C_a * alpha_b * p_b;

  double p_mix = (-B + sqrt(B*B - 4*A*D)) / (2*A);

  // finding alpha mixture 
  double alpha_mix = alpha_a * ((gamma_a - 1) * p_mix + 2 * p_a + C_a) / ((gamma_a + 1) * p_mix + C_a);

  // defining an array of primitive and conservative variables
  Arrayofdouble prim_relaxed, cons_relaxed;

  prim_relaxed[0] = alpha_mix;
  prim_relaxed[1] = alpha_a * rho_a / alpha_mix;
  prim_relaxed[2] = alpha_b * rho_b / (1 - alpha_mix);
  prim_relaxed[3] = (p_mix + gamma_a * pinf_a) / ((gamma_a - 1) * prim_relaxed[1]);
  prim_relaxed[4] = (p_mix + gamma_b * pinf_b) / ((gamma_b - 1) * prim_relaxed[2]);
  double density_relaxed = alpha_mix * prim_relaxed[1] + (1 - alpha_mix) * prim_relaxed[2];
  prim_relaxed[5] = arr(i,j,k,5) / density_relaxed;
  prim_relaxed[6] = p_mix;
  prim_relaxed[7] = arr(i,j,k,7) / density_relaxed;

  cons_relaxed = SGEOS->prim_to_u(prim_relaxed, gamma_a, pinf_a, gamma_b, pinf_b);

  // primitive variable correction
  for(int i=0; i<nVar; i++)
  {
    relaxed_terms[i] = cons_relaxed[i];
  }
}

void numerical_method::ywavespeedestimate(Arrayofdouble const& conserv_l, Arrayofdouble const&  conserv_r, double gamma_a, double gamma_b, double pinf_a, double pinf_b, Arrayofdouble& wavespeed)
{
  // converting input data from conservative variable to primitive variable
  Arrayofdouble prim_l, prim_r;
  prim_l = SGEOS->u_to_prim(conserv_l, gamma_a, pinf_a, gamma_b, pinf_b);
  prim_r = SGEOS->u_to_prim(conserv_r, gamma_a, pinf_a, gamma_b, pinf_b);

  // defining primitive variable left
  double alpha_al = prim_l[0];
  double alpha_bl = 1 - prim_l[0];
  double rho_al = prim_l[1];
  double rho_bl = prim_l[2];
  double eps_al = prim_l[3];
  double eps_bl = prim_l[4];
  double p_al = SGEOS->epstopressure(rho_al, eps_al, pinf_a, gamma_a);
  double p_bl = SGEOS->epstopressure(rho_bl, eps_bl, pinf_b, gamma_b);
  double vx_l = prim_l[5];
  //double p_l = prim_l[6]; // mixture pressure left
  double p_l = alpha_al * p_al + alpha_bl * p_bl; // sm
  double vy_l = prim_l[7]; // vy left

  // defining primitive variable right 
  double alpha_ar = prim_r[0];
  double alpha_br = 1 - prim_r[0];
  double rho_ar = prim_r[1];
  double rho_br = prim_r[2];
  double eps_ar = prim_r[3];
  double eps_br = prim_r[4];

  double p_ar = SGEOS->epstopressure(rho_ar, eps_ar, pinf_a, gamma_a);
  double p_br = SGEOS->epstopressure(rho_br, eps_br, pinf_b, gamma_b);
  double vx_r = prim_r[5];
  //double p_r = prim_r[6]; // mixture pressure right
  double p_r = alpha_ar * p_ar + alpha_br * p_br; // sm
  double vy_r = prim_r[7]; // vy right

  // computing sound speed left from sound speed for material a and b
  double c_al = sqrt(gamma_a * (p_al + pinf_a) / rho_al);
  double c_bl = sqrt(gamma_b * (p_bl + pinf_b) / rho_bl);
  double rho_l = alpha_al * rho_al + alpha_bl * rho_bl; // total density for mass fraction calculation
  double y_al = alpha_al * rho_al / rho_l; // mass fraction for material a
  double y_bl = alpha_bl * rho_bl / rho_l; // mass fraction for material b
  double c_l = sqrt(y_al*c_al*c_al + y_bl*c_bl*c_bl); // frozen mixture sound speed

  // computing sound speed right from sound speed for material a and b
  double c_ar = sqrt(gamma_a * (p_ar + pinf_a) / rho_ar);
  double c_br = sqrt(gamma_b * (p_br + pinf_b) / rho_br);
  double rho_r = alpha_ar * rho_ar + alpha_br * rho_br; // total density for mass fraction calculation
  double y_ar = alpha_ar * rho_ar / rho_r; // mass fraction for material a
  double y_br = alpha_br * rho_br / rho_r; // mass fraction for material b
  double c_r = sqrt(y_ar*c_ar*c_ar + y_br*c_br*c_br); // frozen mixture sound speed

  // computing Sl and Sr
  double Sr = std::max(vy_l + c_l, vy_r + c_r);
  double Sl = std::min(vy_l - c_l, vy_r - c_r);

  // computing S_star
  double rho_vy_l = rho_l * vy_l; 
  double rho_vy_r = rho_r * vy_r; 
  double mom_flux_l = rho_l * vy_l * vy_l + p_l;
  double mom_flux_r = rho_r * vy_r * vy_r + p_r;
  double S_star = (mom_flux_l - mom_flux_r - Sl*rho_vy_l + Sr*rho_vy_r) / (rho_vy_l - rho_vy_r - Sl*rho_l + Sr*rho_r);

  wavespeed[0] = Sl;
  wavespeed[1] = Sr;
  wavespeed[2] = S_star;

  if(Sl > 0)
  {
    wavespeed[3] = vy_l;
  }
  else if(Sr < 0)
  {
    wavespeed[3] = vy_r;
  }
  else
  {
    wavespeed[3] = S_star;
  }
}

void numerical_method::ycompute_HLLCflux(amrex::Array4<amrex::Real> const& arr, amrex::Array4<amrex::Real> const& arr_yrecon_lht, amrex::Array4<amrex::Real> const& arr_yrecon_rht, Arrayofdouble& yFhllc, Arrayofdouble& ywavespeed, int nVar, double gamma_a, double pinf_a, double gamma_b, double pinf_b, int i, int j, int k, double dy, double dt)
{
  //defining left and right reconstructed conservative variables for HLLC solver
  Arrayofdouble conserv_l, conserv_r;
  for(int l=0; l<nVar; l++)
  {
    // conserv_l[l] = arr_yrecon_rht(i,j-1,k,l);
    // conserv_r[l] = arr_yrecon_lht(i,j,k,l);

    conserv_l[l] = arr(i,j-1,k,l);
    conserv_r[l] = arr(i,j,k,l);
  }

  ywavespeedestimate(conserv_l, conserv_r, gamma_a, gamma_b, pinf_a, pinf_b, ywavespeed);

  // define input variables 
  double Sl = ywavespeed[0];
  double Sr = ywavespeed[1];
  double S_star = ywavespeed[2];

  // converting input data from conservative variable to primitive variable
  Arrayofdouble prim_l, prim_r;
  prim_l = SGEOS->u_to_prim(conserv_l, gamma_a, pinf_a, gamma_b, pinf_b);
  prim_r = SGEOS->u_to_prim(conserv_r, gamma_a, pinf_a, gamma_b, pinf_b);

  // defining primitive variable left
  double alpha_al = prim_l[0];
  double alpha_bl = 1 - prim_l[0];
  double rho_al = prim_l[1];
  double rho_bl = prim_l[2];
  double eps_al = prim_l[3];
  double eps_bl = prim_l[4];
  double p_al = SGEOS->epstopressure(rho_al, eps_al, pinf_a, gamma_a);
  double p_bl = SGEOS->epstopressure(rho_bl, eps_bl, pinf_b, gamma_b);
  double vx_l = prim_l[5];
  //double p_l = prim_l[6];
  double p_l = alpha_al * p_al + alpha_bl * p_bl; // sm
  double vy_l = prim_l[7];

  // defining primitive variable right 
  double alpha_ar = prim_r[0];
  double alpha_br = 1 - prim_r[0];
  double rho_ar = prim_r[1];
  double rho_br = prim_r[2];
  double eps_ar = prim_r[3];
  double eps_br = prim_r[4];
  double p_ar = SGEOS->epstopressure(rho_ar, eps_ar, pinf_a, gamma_a);
  double p_br = SGEOS->epstopressure(rho_br, eps_br, pinf_b, gamma_b);
  double vx_r = prim_r[5];
  //double p_r = prim_r[6];
  double p_r = alpha_ar * p_ar + alpha_br * p_br; // sm
  double vy_r = prim_r[7];

  // computing flux left and right 
  Arrayofdouble Fl,Fr;
  Fl = SGEOS->yfluxf(conserv_l, prim_l);
  Fr = SGEOS->yfluxf(conserv_r, prim_r);

  // computing Ql*
  Arrayofdouble Ql_star, Fl_star;

  //alpha_a
  Ql_star[0] = alpha_al; 

  //alpha_a * rho_a
  Ql_star[1] = alpha_al * rho_al * (Sl - vy_l) / (Sl - S_star); 

  //alpha_b * rho_b
  Ql_star[2] = alpha_bl * rho_bl * (Sl - vy_l) / (Sl - S_star); 

  // computing intermediate mixture density
  double rho_alstar = rho_al * (vy_l - Sl) / (S_star - Sl); // fluid density for material a
  double rho_blstar = rho_bl * (vy_l - Sl) / (S_star - Sl); // fluid density for material b 
  double rho_lstar = alpha_al * rho_alstar + alpha_bl * rho_blstar; // mixture fluid density

  //alpha_a * rho_a * eps_xa 
  double gamma_p_rho_al = (gamma_a + 1) * rho_al;
  double gamma_m_rho_al = (gamma_a - 1) * rho_al;
  double gamma_p_rho_alstar = (gamma_a + 1) * rho_alstar;
  double gamma_m_rho_alstar = (gamma_a - 1) * rho_alstar;
  //double p_alstar = (p_al + pinf_a) * (gamma_m_rho_al - gamma_p_rho_alstar) / (gamma_m_rho_alstar - gamma_p_rho_al) - pinf_a; // phasic pressure for material a
  double p_alstar = p_al + rho_al * (Sl - vy_l) * (S_star - vy_l); // sm
  double eps_alstar = (p_alstar + gamma_a * pinf_a) / ((gamma_a - 1) * rho_alstar); // computing intermediate internal energy for material a left from eos
  Ql_star[3] = alpha_al * rho_alstar * eps_alstar;

  //alpha_b * rho_b * eps_xb 
  double gamma_p_rho_bl = (gamma_b + 1) * rho_bl;
  double gamma_m_rho_bl = (gamma_b - 1) * rho_bl;
  double gamma_p_rho_blstar = (gamma_b + 1) * rho_blstar;
  double gamma_m_rho_blstar = (gamma_b - 1) * rho_blstar;
  //double p_blstar = (p_bl + pinf_b) * (gamma_m_rho_bl - gamma_p_rho_blstar) / (gamma_m_rho_blstar - gamma_p_rho_bl) - pinf_b; // phasic pressure for material b
  double p_blstar = p_bl + rho_bl * (Sl - vy_l) * (S_star - vy_l); // sm
  double eps_blstar = (p_blstar + gamma_b * pinf_b) / ((gamma_b - 1) * rho_blstar); // computing intermediate internal energy for material b left from eos
  Ql_star[4] = alpha_bl * rho_blstar * eps_blstar;

  //rho * v_x
  double rho_l = alpha_al * rho_al + alpha_bl * rho_bl; // computing mixture density left
  Ql_star[5] = rho_lstar * vx_l;

  // rho * energy
  double p_lstar = p_l + rho_l * vy_l * (vy_l - Sl) - rho_lstar * S_star * (S_star - Sl); // computing intermediate pressure left
  double y_al = alpha_al * rho_al / rho_l; // mass fraction for material a left
  double y_bl = alpha_bl * rho_bl / rho_l; // mass fraction for material b left
  double energy_l = SGEOS->compute_energy(alpha_al, rho_al, eps_al, rho_bl, eps_bl, vx_l, vy_l); // computing total energy left
  double energy_lstar = (rho_l * energy_l * (vy_l - Sl) + p_l * vy_l - p_lstar * S_star) / (rho_lstar * (S_star - Sl)); // computing total intermediate energy left
  //Ql_star[6] = rho_lstar * energy_lstar;
  //Ql_star[6] = rho_lstar * (energy_l + (S_star - vy_l) * (S_star * rho_l + p_l / (Sl - vy_l)));
  Ql_star[6] = (Sl - vy_l) / (Sl - S_star) * (rho_l * energy_l + (S_star - vy_l) * (S_star * rho_l + p_l / (Sl - vy_l)));

  //rho * v_y
  Ql_star[7] = rho_lstar * S_star;

  // computing Qr*
  Arrayofdouble Qr_star, Fr_star;

  //alpha_a
  Qr_star[0] = alpha_ar; 

  //alpha_a * rho_a
  Qr_star[1] = alpha_ar * rho_ar * (Sr - vy_r) / (Sr - S_star); 

  //alpha_b * rho_b
  Qr_star[2] = alpha_br * rho_br * (Sr - vy_r) / (Sr - S_star); 

  // computing intermediate mixture density
  double rho_arstar = rho_ar * (vy_r - Sr) / (S_star - Sr); // fluid density for material a
  double rho_brstar = rho_br * (vy_r - Sr) / (S_star - Sr); // fluid density for material b 
  double rho_rstar = alpha_ar * rho_arstar + alpha_br * rho_brstar; // mixture fluid density

  //alpha_a * rho_a * eps_xa 
  double gamma_p_rho_ar = (gamma_a + 1) * rho_ar;
  double gamma_m_rho_ar = (gamma_a - 1) * rho_ar;
  double gamma_p_rho_arstar = (gamma_a + 1) * rho_arstar;
  double gamma_m_rho_arstar = (gamma_a - 1) * rho_arstar;
  //double p_arstar = (p_ar + pinf_a) * (gamma_m_rho_ar - gamma_p_rho_arstar) / (gamma_m_rho_arstar - gamma_p_rho_ar) - pinf_a; // phasic pressure for material a
  double p_arstar = p_ar + rho_ar * (Sr - vy_r) * (S_star - vy_r); // sm
  double eps_arstar = (p_arstar + gamma_a * pinf_a) / ((gamma_a - 1) * rho_arstar); // computing intermediate internal energy for material a left from eos
  Qr_star[3] = alpha_ar * rho_arstar * eps_arstar; 

  //alpha_b * rho_b * eps_xb 
  double gamma_p_rho_br = (gamma_b + 1) * rho_br;
  double gamma_m_rho_br = (gamma_b - 1) * rho_br;
  double gamma_p_rho_brstar = (gamma_b + 1) * rho_brstar;
  double gamma_m_rho_brstar = (gamma_b - 1) * rho_brstar;
  //double p_brstar = (p_br + pinf_b) * (gamma_m_rho_br - gamma_p_rho_brstar) / (gamma_m_rho_brstar - gamma_p_rho_br) - pinf_b; // phasic pressure for material b
  double p_brstar = p_br + rho_br * (Sr - vy_r) * (S_star - vy_r); // sm
  double eps_brstar = (p_brstar + gamma_b * pinf_b) / ((gamma_b - 1) * rho_brstar); // computing intermediate internal energy for material b left from eos
  Qr_star[4] = alpha_br * rho_brstar * eps_brstar; 

  //rho * v_x
  double rho_r = alpha_ar * rho_ar + alpha_br * rho_br; // computing mixture density right
  Qr_star[5] = rho_rstar * vx_r;

  // rho * energy
  double p_rstar = p_r + rho_r * vy_r * (vy_r - Sr) - rho_rstar * S_star * (S_star - Sr); // computing intermediate pressure right
  double y_ar = alpha_ar * rho_ar / rho_r; // mass fraction for material a left
  double y_br = alpha_br * rho_br / rho_r; // mass fraction for material b left
  double energy_r = SGEOS->compute_energy(alpha_ar, rho_ar, eps_ar, rho_br, eps_br, vx_r, vy_r); // computing total energy left 
  double energy_rstar = (rho_r * energy_r * (vy_r - Sr) + p_r * vy_r - p_rstar * S_star) / (rho_rstar * (S_star - Sr)); // computing total intermediate energy left
  //Qr_star[6] = rho_rstar * energy_rstar;
  //Qr_star[6] = rho_rstar * (energy_r + (S_star - vy_r) * (S_star * rho_r + p_r / (Sr - vy_r)));
  Qr_star[6] = (Sr - vy_r) / (Sr - S_star) * (rho_r * energy_r + (S_star - vy_r) * (S_star * rho_r + p_r / (Sr - vy_r)));

  // rho * v_v
  Qr_star[7] = rho_rstar * S_star;

  // computing Fl* and Fr*
  for(int l=0; l<nVar; l++)
  {
  Fl_star[l] = Fl[l]  + Sl*(Ql_star[l] - conserv_l[l]);
  Fr_star[l] = Fr[l]  + Sr*(Qr_star[l] - conserv_r[l]);
  }

  // implementing logic for choosing flux
	if(Sl >= 0)
	{	
		yFhllc = Fl;
    yFhllc[0] = conserv_l[0];
    yFhllc[3] = conserv_l[3];
    yFhllc[4] = conserv_l[4];
		return; 
	}

	if(Sl < 0 && S_star >= 0)
	{
		yFhllc = Fl_star;
    yFhllc[0] = Ql_star[0];
    yFhllc[3] = Ql_star[3];
    yFhllc[4] = Ql_star[4];
		return;
	}
	
	if(S_star < 0 && Sr > 0)
	{
		yFhllc = Fr_star;
    yFhllc[0] = Qr_star[0];
    yFhllc[3] = Qr_star[3];
    yFhllc[4] = Qr_star[4];
		return;
	}

	if(Sr <= 0)
	{	
    yFhllc = Fr;
    yFhllc[0] = conserv_r[0];
    yFhllc[3] = conserv_r[3];
    yFhllc[4] = conserv_r[4];
		return; 
	}		
}