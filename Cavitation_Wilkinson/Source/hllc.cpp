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

void numerical_method::wavespeedestimate(Arrayofdouble const& conserv_l, Arrayofdouble const&  conserv_r, double gamma_a, double gamma_b, double pinf_a, double pinf_b, Arrayofdouble& wavespeed, double p_max, double epsinf_a, double epsinf_b, double cv_a, double cv_b)
{
  // converting input data from conservative variable to primitive variable
  Arrayofdouble prim_l, prim_r;
  prim_l = SGEOS->u_to_prim(conserv_l, gamma_a, pinf_a, gamma_b, pinf_b, p_max, epsinf_a, epsinf_b, cv_a);
  prim_r = SGEOS->u_to_prim(conserv_r, gamma_a, pinf_a, gamma_b, pinf_b, p_max, epsinf_a, epsinf_b, cv_a);

  // defining primitive variable left
  double alpha_al = prim_l[0];
  double alpha_bl = 1.0 - prim_l[0];
  double rho_al = prim_l[1];
  double rho_bl = prim_l[2];
  double vx_l = prim_l[3];
  double p_l = prim_l[4]; // mixture pressure left

  // defining primitive variable right 
  double alpha_ar = prim_r[0];
  double alpha_br = 1.0 - prim_r[0];
  double rho_ar = prim_r[1];
  double rho_br = prim_r[2];
  double vx_r = prim_r[3];
  double p_r = prim_r[4]; // mixture pressure right

  // computing sound speed left from sound speed for material a and b
  double c_l = SGEOS->compute_soundspeed(conserv_l[0], conserv_l[1], conserv_l[2], prim_l[4], gamma_a, pinf_a, gamma_b, pinf_b, p_max, epsinf_a, epsinf_b, cv_a, cv_b);

  // computing sound speed right from sound speed for material a and b
  double c_r = SGEOS->compute_soundspeed(conserv_r[0], conserv_r[1], conserv_r[2], prim_r[4], gamma_a, pinf_a, gamma_b, pinf_b, p_max, epsinf_a, epsinf_b, cv_a, cv_b);

  // computing Sl and Sr
  double Sr = std::max(vx_l + c_l, vx_r + c_r);
  double Sl = std::min(vx_l - c_l, vx_r - c_r);

  // computing S_star
  double rho_l = alpha_al * rho_al + alpha_bl * rho_bl; // total density for mass fraction calculation
  double rho_r = alpha_ar * rho_ar + alpha_br * rho_br; // total density for mass fraction calculation
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

void numerical_method::compute_HLLCflux(amrex::Array4<amrex::Real> const& arr, amrex::Array4<amrex::Real> const& arr_xrecon_lht, amrex::Array4<amrex::Real> const& arr_xrecon_rht, Arrayofdouble& Fhllc, Arrayofdouble& wavespeed, int nVar, double gamma_a, double pinf_a, double gamma_b, double pinf_b, int i, int j, int k, double dx, double dt, double p_max, double epsinf_a, double epsinf_b, double cv_a, double cv_b)
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
  wavespeedestimate(conserv_l, conserv_r, gamma_a, gamma_b, pinf_a, pinf_b, wavespeed, p_max, epsinf_a, epsinf_b, cv_a, cv_b);

  // define input variables 
  double Sl = wavespeed[0];
  double Sr = wavespeed[1];
  double S_star = wavespeed[2];

  // converting input data from conservative variable to primitive variable
  Arrayofdouble prim_l, prim_r;
  prim_l = SGEOS->u_to_prim(conserv_l, gamma_a, pinf_a, gamma_b, pinf_b, p_max, epsinf_a, epsinf_b, cv_a);
  prim_r = SGEOS->u_to_prim(conserv_r, gamma_a, pinf_a, gamma_b, pinf_b, p_max, epsinf_a, epsinf_b, cv_a);

  // defining primitive variable left
  double alpha_al = prim_l[0];
  double alpha_bl = 1.0 - prim_l[0];
  double rho_al = prim_l[1];
  double rho_bl = prim_l[2];
  double vx_l = prim_l[3];
  double p_l = prim_l[4];
  double vy_l = prim_l[5];

  // defining primitive variable right 
  double alpha_ar = prim_r[0];
  double alpha_br = 1.0 - prim_r[0];
  double rho_ar = prim_r[1];
  double rho_br = prim_r[2];
  double vx_r = prim_r[3];
  double p_r = prim_r[4];
  double vy_r = prim_r[5];

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

  //rho * v_x
  double rho_l = alpha_al * rho_al + alpha_bl * rho_bl; // computing mixture density left
  Ql_star[3] = rho_lstar * S_star;

  // rho * energy
  double p_lstar = p_l + rho_l * vx_l * (vx_l - Sl) - rho_lstar * S_star * (S_star - Sl); // computing intermediate pressure left
  double energy_l = SGEOS->compute_energy(conserv_l[0], conserv_l[1], conserv_l[2], prim_l[4], prim_l[3], prim_l[5], gamma_a, pinf_a, gamma_b, pinf_b, p_max, epsinf_a, epsinf_b);
  double energy_lstar = (rho_l * energy_l * (vx_l - Sl) + p_l * vx_l - p_lstar * S_star) / (rho_lstar * (S_star - Sl)); // computing total intermediate energy left
  Ql_star[4] = (Sl - vx_l) / (Sl - S_star) * (rho_l * energy_l + (S_star - vx_l) * (S_star * rho_l + p_l / (Sl - vx_l)));

  // rho * v_y
  Ql_star[5] = rho_lstar * vy_l;

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

  //rho * v_x
  double rho_r = alpha_ar * rho_ar + alpha_br * rho_br; // computing mixture density right
  Qr_star[3] = rho_rstar * S_star;

  // rho * energy
  double p_rstar = p_r + rho_r * vx_r * (vx_r - Sr) - rho_rstar * S_star * (S_star - Sr); // computing intermediate pressure right
  double energy_r = SGEOS->compute_energy(conserv_r[0], conserv_r[1], conserv_r[2], prim_r[4], prim_r[3], prim_r[5], gamma_a, pinf_a, gamma_b, pinf_b, p_max, epsinf_a, epsinf_b);
  double energy_rstar = (rho_r * energy_r * (vx_r - Sr) + p_r * vx_r - p_rstar * S_star) / (rho_rstar * (S_star - Sr)); // computing total intermediate energy left
  //Qr_star[6] = rho_rstar * energy_rstar;
  //Qr_star[6] = rho_rstar * (energy_r + (S_star - vx_r) * (S_star * rho_r + p_r / (Sr - vx_r)));
  Qr_star[4] = (Sr - vx_r) / (Sr - S_star) * (rho_r * energy_r + (S_star - vx_r) * (S_star * rho_r + p_r / (Sr - vx_r)));

  // rho * v_y
  Qr_star[5] = rho_rstar * vy_r;

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
		return; 
	}

	else if(S_star >= 0)
	{
		Fhllc = Fl_star;
    Fhllc[0] = Ql_star[0];
		return;
	}
	
	else if(Sr > 0)
	{
		Fhllc = Fr_star;
    Fhllc[0] = Qr_star[0];
		return;
	}

	else
	{	Fhllc = Fr;
    Fhllc[0] = conserv_r[0];
		return; 
	}		
}

void numerical_method::ywavespeedestimate(Arrayofdouble const& conserv_l, Arrayofdouble const&  conserv_r, double gamma_a, double gamma_b, double pinf_a, double pinf_b, Arrayofdouble& wavespeed, double p_max, double epsinf_a, double epsinf_b, double cv_a, double cv_b)
{
  // converting input data from conservative variable to primitive variable
  Arrayofdouble prim_l, prim_r;
  prim_l = SGEOS->u_to_prim(conserv_l, gamma_a, pinf_a, gamma_b, pinf_b, p_max, epsinf_a, epsinf_b, cv_a);
  prim_r = SGEOS->u_to_prim(conserv_r, gamma_a, pinf_a, gamma_b, pinf_b, p_max, epsinf_a, epsinf_b, cv_a);

  // defining primitive variable left
  double alpha_al = prim_l[0];
  double alpha_bl = 1.0 - prim_l[0];
  double rho_al = prim_l[1];
  double rho_bl = prim_l[2];
  double vx_l = prim_l[3];
  double p_l = prim_l[4]; // mixture pressure left
  double vy_l = prim_l[5]; // vy left

  // defining primitive variable right 
  double alpha_ar = prim_r[0];
  double alpha_br = 1.0 - prim_r[0];
  double rho_ar = prim_r[1];
  double rho_br = prim_r[2];
  double vx_r = prim_r[3];
  double p_r = prim_r[4]; // mixture pressure right
  double vy_r = prim_r[5]; // vy right

  // computing sound speed left from sound speed for material a and b
  double c_l = SGEOS->compute_soundspeed(conserv_l[0], conserv_l[1], conserv_l[2], prim_l[4], gamma_a, pinf_a, gamma_b, pinf_b, p_max, epsinf_a, epsinf_b, cv_a, cv_b);

  // computing sound speed right from sound speed for material a and b
  double c_r = SGEOS->compute_soundspeed(conserv_r[0], conserv_r[1], conserv_r[2], prim_r[4], gamma_a, pinf_a, gamma_b, pinf_b, p_max, epsinf_a, epsinf_b, cv_a, cv_b);

  // computing Sl and Sr
  double Sr = std::max(vy_l + c_l, vy_r + c_r);
  double Sl = std::min(vy_l - c_l, vy_r - c_r);

  // computing S_star
  double rho_l = alpha_al * rho_al + alpha_bl * rho_bl; // total density for mass fraction calculation
  double rho_r = alpha_ar * rho_ar + alpha_br * rho_br; // total density for mass fraction calculation
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

void numerical_method::ycompute_HLLCflux(amrex::Array4<amrex::Real> const& arr, amrex::Array4<amrex::Real> const& arr_yrecon_lht, amrex::Array4<amrex::Real> const& arr_yrecon_rht, Arrayofdouble& yFhllc, Arrayofdouble& ywavespeed, int nVar, double gamma_a, double pinf_a, double gamma_b, double pinf_b, int i, int j, int k, double dy, double dt, double p_max, double epsinf_a, double epsinf_b, double cv_a, double cv_b)
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

  ywavespeedestimate(conserv_l, conserv_r, gamma_a, gamma_b, pinf_a, pinf_b, ywavespeed, p_max, epsinf_a, epsinf_b, cv_a, cv_b);

  // define input variables 
  double Sl = ywavespeed[0];
  double Sr = ywavespeed[1];
  double S_star = ywavespeed[2];

  // converting input data from conservative variable to primitive variable
  Arrayofdouble prim_l, prim_r;
  prim_l = SGEOS->u_to_prim(conserv_l, gamma_a, pinf_a, gamma_b, pinf_b, p_max, epsinf_a, epsinf_b, cv_a);
  prim_r = SGEOS->u_to_prim(conserv_r, gamma_a, pinf_a, gamma_b, pinf_b, p_max, epsinf_a, epsinf_b, cv_a);

  // defining primitive variable left
  double alpha_al = prim_l[0];
  double alpha_bl = 1.0 - prim_l[0];
  double rho_al = prim_l[1];
  double rho_bl = prim_l[2];
  double vx_l = prim_l[3];
  double p_l = prim_l[4];
  double vy_l = prim_l[5];

  // defining primitive variable right 
  double alpha_ar = prim_r[0];
  double alpha_br = 1.0 - prim_r[0];
  double rho_ar = prim_r[1];
  double rho_br = prim_r[2];
  double vx_r = prim_r[3];
  double p_r = prim_r[4];
  double vy_r = prim_r[5];

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

  //rho * v_x
  double rho_l = alpha_al * rho_al + alpha_bl * rho_bl; // computing mixture density left
  Ql_star[3] = rho_lstar * vx_l;

  // rho * energy
  double p_lstar = p_l + rho_l * vy_l * (vy_l - Sl) - rho_lstar * S_star * (S_star - Sl); // computing intermediate pressure left
  double energy_l = SGEOS->compute_energy(conserv_l[0], conserv_l[1], conserv_l[2], prim_l[4], prim_l[3], prim_l[5], gamma_a, pinf_a, gamma_b, pinf_b, p_max, epsinf_a, epsinf_b);
  double energy_lstar = (rho_l * energy_l * (vy_l - Sl) + p_l * vy_l - p_lstar * S_star) / (rho_lstar * (S_star - Sl)); // computing total intermediate energy left
  //Ql_star[6] = rho_lstar * energy_lstar;
  //Ql_star[6] = rho_lstar * (energy_l + (S_star - vy_l) * (S_star * rho_l + p_l / (Sl - vy_l)));
  Ql_star[4] = (Sl - vy_l) / (Sl - S_star) * (rho_l * energy_l + (S_star - vy_l) * (S_star * rho_l + p_l / (Sl - vy_l)));

  //rho * v_y
  Ql_star[5] = rho_lstar * S_star;

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

  //rho * v_x
  double rho_r = alpha_ar * rho_ar + alpha_br * rho_br; // computing mixture density right
  Qr_star[3] = rho_rstar * vx_r;

  // rho * energy
  double p_rstar = p_r + rho_r * vy_r * (vy_r - Sr) - rho_rstar * S_star * (S_star - Sr); // computing intermediate pressure right
  double energy_r = SGEOS->compute_energy(conserv_r[0], conserv_r[1], conserv_r[2], prim_r[4], prim_r[3], prim_r[5], gamma_a, pinf_a, gamma_b, pinf_b, p_max, epsinf_a, epsinf_b);
  double energy_rstar = (rho_r * energy_r * (vy_r - Sr) + p_r * vy_r - p_rstar * S_star) / (rho_rstar * (S_star - Sr)); // computing total intermediate energy left
  //Qr_star[6] = rho_rstar * energy_rstar;
  //Qr_star[6] = rho_rstar * (energy_r + (S_star - vy_r) * (S_star * rho_r + p_r / (Sr - vy_r)));
  Qr_star[4] = (Sr - vy_r) / (Sr - S_star) * (rho_r * energy_r + (S_star - vy_r) * (S_star * rho_r + p_r / (Sr - vy_r)));

  // rho * v_v
  Qr_star[5] = rho_rstar * S_star;

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
		return; 
	}

	if(Sl < 0 && S_star >= 0)
	{
		yFhllc = Fl_star;
    yFhllc[0] = Ql_star[0];
		return;
	}
	
	if(S_star < 0 && Sr > 0)
	{
		yFhllc = Fr_star;
    yFhllc[0] = Qr_star[0];
		return;
	}

	if(Sr <= 0)
	{	
    yFhllc = Fr;
    yFhllc[0] = conserv_r[0];
		return; 
	}		
}

//function to find Vanleer constant
Arrayofdouble numerical_method::Vanleer(Arrayofdouble ui, Arrayofdouble uiMinus1, Arrayofdouble uiPlus1, int nVar)
{
  Arrayofdouble r;
  Arrayofdouble xi;
  Arrayofdouble Vl;

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
    xi[i] = 2/(1+r[i]);
    // r[i] = 0;
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

//computing nplushalf for xHLLC flux
void numerical_method::secondorder_extension(Arrayofdouble u_n, Arrayofdouble u, Arrayofdouble u_p, double gamma_a, double gamma_b, double pinf_a, double pinf_b, Arrayofdouble& primiL, Arrayofdouble& primiR, int nVar, double cell_size, double sl_const, double p_max, double epsinf_a, double epsinf_b, double cv_a)
{
  Arrayofdouble primi, primi_minus1, primi_plus1, Vl_prim, delta_prim_plushalf, delta_prim_minushalf, delta_prim;

  // compute primitive variables from input 
  primi_minus1 = SGEOS->u_to_prim(u_n, gamma_a, pinf_a, gamma_b, pinf_b, p_max, epsinf_a, epsinf_b, cv_a);
  primi = SGEOS->u_to_prim(u, gamma_a, pinf_a, gamma_b, pinf_b, p_max, epsinf_a, epsinf_b, cv_a);
  primi_plus1 = SGEOS->u_to_prim(u_p, gamma_a, pinf_a, gamma_b, pinf_b, p_max, epsinf_a, epsinf_b, cv_a);

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

  // define delta using primitive variables
  for(int l=0; l<nVar; l++)
  {
    delta_prim_minushalf[l] = (primi[l] - primi_minus1[l]) / cell_size;
    delta_prim_plushalf[l] = (primi_plus1[l] - primi[l]) / cell_size;
    delta_prim[l] = 0.5 * (1 + sl_const) * delta_prim_minushalf[l] + 0.5 * (1 - sl_const) * delta_prim_plushalf[l];
  }

  // variables extrapolation
  // define primiL, primiR
  for(int l=0; l<nVar; l++)
  {
    primiL[l] = primi[l] - 0.5 * Vl_prim_min * cell_size * delta_prim[l];
    primiR[l] = primi[l] + 0.5 * Vl_prim_min * cell_size * delta_prim[l];
  }
}

void numerical_method::halftimestep_update(Arrayofdouble primi, Arrayofdouble primiL, Arrayofdouble primiR, Arrayofdouble yprimi, Arrayofdouble yprimiL, Arrayofdouble yprimiR, double dt, double dx, double dy, Arrayofdouble& uiL_nplushalf, Arrayofdouble& uiR_nplushalf, double gamma_a, double gamma_b, double pinf_a, double pinf_b, double p_max, double epsinf_a, double epsinf_b, double cv_a, double cv_b)
{
  Arrayofdouble primiL_nplushalf, primiR_nplushalf;

  // defining centered variables
  double density = primi[0] * primi[1] + (1 - primi[0]) * primi[2]; // computing mixture density in x-direction
  double ydensity = yprimi[0] * yprimi[1] + (1 - yprimi[0]) * yprimi[2]; // computing mixture density in y-direction
  double c = SGEOS->compute_soundspeed(primi[0], primi[0]*primi[1], (1-primi[0])*primi[2], primi[4], gamma_a, pinf_a, gamma_b, pinf_b, p_max, epsinf_a, epsinf_b, cv_a, cv_b);
  double yc = SGEOS->compute_soundspeed(yprimi[0], yprimi[0]*yprimi[1], (1-yprimi[0])*yprimi[2], yprimi[4], gamma_a, pinf_a, gamma_b, pinf_b, p_max, epsinf_a, epsinf_b, cv_a, cv_b);

  // calculate primiL_nplushalf
  primiL_nplushalf[0] = primiL[0] + 0.5 * dt / dx * (primi[3] * (primiL[0] - primiR[0])) + 0.5 * dt / dy * (yprimi[5] * (yprimiL[0] - yprimiR[0]));
  primiL_nplushalf[1] = primiL[1] + 0.5 * dt / dx * (primi[3] * (primiL[1] - primiR[1]) + primi[1] * (primiL[3] - primiR[3])) + 0.5 * dt / dy * (yprimi[5] * (yprimiL[1] - yprimiR[1]) + yprimi[1] * (yprimiL[5] - yprimiR[5]));
  primiL_nplushalf[2] = primiL[2] + 0.5 * dt / dx * (primi[3] * (primiL[2] - primiR[2]) + primi[2] * (primiL[3] - primiR[3])) + 0.5 * dt / dy * (yprimi[5] * (yprimiL[2] - yprimiR[2]) + yprimi[2] * (yprimiL[5] - yprimiR[5]));
  primiL_nplushalf[3] = primiL[3] + 0.5 * dt / dx * (primi[3] * (primiL[3] - primiR[3]) + 1/density * (primiL[4] - primiR[4])) + 0.5 * dt / dy * (yprimi[5] * (yprimiL[3] - yprimiR[3]));
  primiL_nplushalf[4] = primiL[4] + 0.5 * dt / dx * (primi[3] * (primiL[4] - primiR[4]) + density * c * c * (primiL[3] - primiR[3])) + 0.5 * dt / dy * (yprimi[5] * (yprimiL[4] - yprimiR[4]) + ydensity * yc * yc * (yprimiL[5] - yprimiR[5]));
  primiL_nplushalf[5] = primiL[5] + 0.5 * dt / dx * (primi[3] * (primiL[5] - primiR[5])) + 0.5 * dt / dy * (yprimi[5] * (yprimiL[5] - yprimiR[5]) + 1/ydensity * (yprimiL[4] - yprimiR[4]));

  // calculate primiR_nplushalf
  primiR_nplushalf[0] = primiR[0] + 0.5 * dt / dx * (primi[3] * (primiL[0] - primiR[0])) + 0.5 * dt / dy * (yprimi[5] * (yprimiL[0] - yprimiR[0]));
  primiR_nplushalf[1] = primiR[1] + 0.5 * dt / dx * (primi[3] * (primiL[1] - primiR[1]) + primi[1] * (primiL[3] - primiR[3])) + 0.5 * dt / dy * (yprimi[5] * (yprimiL[1] - yprimiR[1]) + yprimi[1] * (yprimiL[5] - yprimiR[5]));
  primiR_nplushalf[2] = primiR[2] + 0.5 * dt / dx * (primi[3] * (primiL[2] - primiR[2]) + primi[2] * (primiL[3] - primiR[3])) + 0.5 * dt / dy * (yprimi[5] * (yprimiL[2] - yprimiR[2]) + yprimi[2] * (yprimiL[5] - yprimiR[5]));
  primiR_nplushalf[3] = primiR[3] + 0.5 * dt / dx * (primi[3] * (primiL[3] - primiR[3]) + 1/density * (primiL[4] - primiR[4])) + 0.5 * dt / dy * (yprimi[5] * (yprimiL[3] - yprimiR[3]));
  primiR_nplushalf[4] = primiR[4] + 0.5 * dt / dx * (primi[3] * (primiL[4] - primiR[4]) + density * c * c * (primiL[3] - primiR[3])) + 0.5 * dt / dy * (yprimi[5] * (yprimiL[4] - yprimiR[4]) + ydensity * yc * yc * (yprimiL[5] - yprimiR[5]));
  primiR_nplushalf[5] = primiR[5] + 0.5 * dt / dx * (primi[3] * (primiL[5] - primiR[5])) + 0.5 * dt / dy * (yprimi[5] * (yprimiL[5] - yprimiR[5]) + 1/ydensity * (yprimiL[4] - yprimiR[4]));

  // compute conservative variables for output
  uiL_nplushalf = SGEOS->prim_to_u(primiL_nplushalf, gamma_a, pinf_a, gamma_b, pinf_b, p_max, epsinf_a, epsinf_b);
  uiR_nplushalf = SGEOS->prim_to_u(primiR_nplushalf, gamma_a, pinf_a, gamma_b, pinf_b, p_max, epsinf_a, epsinf_b);
}

void numerical_method::yhalftimestep_update(Arrayofdouble primi, Arrayofdouble primiL, Arrayofdouble primiR, Arrayofdouble yprimi, Arrayofdouble yprimiL, Arrayofdouble yprimiR, double dt, double dx, double dy, Arrayofdouble& uiL_nplushalf, Arrayofdouble& uiR_nplushalf, double gamma_a, double gamma_b, double pinf_a, double pinf_b, double p_max, double epsinf_a, double epsinf_b, double cv_a, double cv_b)
{
  Arrayofdouble primiL_nplushalf, primiR_nplushalf;

  // defining centered variables
  double density = primi[0] * primi[1] + (1 - primi[0]) * primi[2]; // computing mixture density in x-direction
  double ydensity = yprimi[0] * yprimi[1] + (1 - yprimi[0]) * yprimi[2]; // computing mixture density in y-direction
  double c = SGEOS->compute_soundspeed(primi[0], primi[0]*primi[1], (1-primi[0])*primi[2], primi[4], gamma_a, pinf_a, gamma_b, pinf_b, p_max, epsinf_a, epsinf_b, cv_a, cv_b);
  double yc = SGEOS->compute_soundspeed(yprimi[0], yprimi[0]*yprimi[1], (1-yprimi[0])*yprimi[2], yprimi[4], gamma_a, pinf_a, gamma_b, pinf_b, p_max, epsinf_a, epsinf_b, cv_a, cv_b);

  // calculate primiL_nplushalf
  primiL_nplushalf[0] = yprimiL[0] + 0.5 * dt / dx * (primi[3] * (primiL[0] - primiR[0])) + 0.5 * dt / dy * (yprimi[5] * (yprimiL[0] - yprimiR[0]));
  primiL_nplushalf[1] = yprimiL[1] + 0.5 * dt / dx * (primi[3] * (primiL[1] - primiR[1]) + primi[1] * (primiL[3] - primiR[3])) + 0.5 * dt / dy * (yprimi[5] * (yprimiL[1] - yprimiR[1]) + yprimi[1] * (yprimiL[5] - yprimiR[5]));
  primiL_nplushalf[2] = yprimiL[2] + 0.5 * dt / dx * (primi[3] * (primiL[2] - primiR[2]) + primi[2] * (primiL[3] - primiR[3])) + 0.5 * dt / dy * (yprimi[5] * (yprimiL[2] - yprimiR[2]) + yprimi[2] * (yprimiL[5] - yprimiR[5]));
  primiL_nplushalf[3] = yprimiL[3] + 0.5 * dt / dx * (primi[3] * (primiL[3] - primiR[3]) + 1/density * (primiL[4] - primiR[4])) + 0.5 * dt / dy * (yprimi[5] * (yprimiL[3] - yprimiR[3]));
  primiL_nplushalf[4] = yprimiL[4] + 0.5 * dt / dx * (primi[3] * (primiL[4] - primiR[4]) + density * c * c * (primiL[3] - primiR[3])) + 0.5 * dt / dy * (yprimi[5] * (yprimiL[4] - yprimiR[4]) + ydensity * yc * yc * (yprimiL[5] - yprimiR[5]));
  primiL_nplushalf[5] = yprimiL[5] + 0.5 * dt / dx * (primi[3] * (primiL[5] - primiR[5])) + 0.5 * dt / dy * (yprimi[5] * (yprimiL[5] - yprimiR[5]) + 1/ydensity * (yprimiL[4] - yprimiR[4]));

  // calculate primiR_nplushalf
  primiR_nplushalf[0] = yprimiR[0] + 0.5 * dt / dx * (primi[3] * (primiL[0] - primiR[0])) + 0.5 * dt / dy * (yprimi[5] * (yprimiL[0] - yprimiR[0]));
  primiR_nplushalf[1] = yprimiR[1] + 0.5 * dt / dx * (primi[3] * (primiL[1] - primiR[1]) + primi[1] * (primiL[3] - primiR[3])) + 0.5 * dt / dy * (yprimi[5] * (yprimiL[1] - yprimiR[1]) + yprimi[1] * (yprimiL[5] - yprimiR[5]));
  primiR_nplushalf[2] = yprimiR[2] + 0.5 * dt / dx * (primi[3] * (primiL[2] - primiR[2]) + primi[2] * (primiL[3] - primiR[3])) + 0.5 * dt / dy * (yprimi[5] * (yprimiL[2] - yprimiR[2]) + yprimi[2] * (yprimiL[5] - yprimiR[5]));
  primiR_nplushalf[3] = yprimiR[3] + 0.5 * dt / dx * (primi[3] * (primiL[3] - primiR[3]) + 1/density * (primiL[4] - primiR[4])) + 0.5 * dt / dy * (yprimi[5] * (yprimiL[3] - yprimiR[3]));
  primiR_nplushalf[4] = yprimiR[4] + 0.5 * dt / dx * (primi[3] * (primiL[4] - primiR[4]) + density * c * c * (primiL[3] - primiR[3])) + 0.5 * dt / dy * (yprimi[5] * (yprimiL[4] - yprimiR[4]) + ydensity * yc * yc * (yprimiL[5] - yprimiR[5]));
  primiR_nplushalf[5] = yprimiR[5] + 0.5 * dt / dx * (primi[3] * (primiL[5] - primiR[5])) + 0.5 * dt / dy * (yprimi[5] * (yprimiL[5] - yprimiR[5]) + 1/ydensity * (yprimiL[4] - yprimiR[4]));

  // compute conservative variables for output
  uiL_nplushalf = SGEOS->prim_to_u(primiL_nplushalf, gamma_a, pinf_a, gamma_b, pinf_b, p_max, epsinf_a, epsinf_b);
  uiR_nplushalf = SGEOS->prim_to_u(primiR_nplushalf, gamma_a, pinf_a, gamma_b, pinf_b, p_max, epsinf_a, epsinf_b);
}