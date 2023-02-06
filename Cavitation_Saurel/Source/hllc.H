#include <array>
#include <iostream>
#include <fstream>
#include <string>

#include "AMReX_Array4.H"

//array sized 10 for storing simulation variables
typedef std::array<double,14> Arrayofdouble;

class numerical_method
{
private:

public:
    numerical_method();
    
    void wavespeedestimate(Arrayofdouble const& conserv_l, Arrayofdouble const&  conserv_r, double gamma_a, double gamma_b, double pinf_a, double pinf_b, Arrayofdouble& wavespeed);
    void compute_HLLCflux(amrex::Array4<amrex::Real> const& arr, amrex::Array4<amrex::Real> const& arr_xrecon_lht, amrex::Array4<amrex::Real> const& arr_xrecon_rht, Arrayofdouble& Fhllc, Arrayofdouble& wavespeed, int nVar, double gamma_a, double pinf_a, double gamma_b, double pinf_b, int i, int j, int k, double dx, double dt);
    void pressure_relaxation(amrex::Array4<amrex::Real> const arr, int i, int j, int k, double gamma_a, double pinf_a, double gamma_b, double pinf_b, Arrayofdouble& relaxed_terms);
    void pressure_relaxation2(amrex::Array4<amrex::Real> const& arr, int i, int j, int k, double gamma_a, double pinf_a, double gamma_b, double pinf_b, Arrayofdouble& relaxed_terms, int nVar); // primitive variable pressure relaxation
    
    void ywavespeedestimate(Arrayofdouble const& conserv_l, Arrayofdouble const&  conserv_r, double gamma_a, double gamma_b, double pinf_a, double pinf_b, Arrayofdouble& wavespeed);
    void ycompute_HLLCflux(amrex::Array4<amrex::Real> const& arr, amrex::Array4<amrex::Real> const& arr_yrecon_lht, amrex::Array4<amrex::Real> const& arr_yrecon_rht, Arrayofdouble& yFhllc, Arrayofdouble& ywavespeed, int nVar, double gamma_a, double pinf_a, double gamma_b, double pinf_b, int i, int j, int k, double dy, double dt);

    Arrayofdouble Minbee(Arrayofdouble ui, Arrayofdouble uiMinus1, Arrayofdouble uiPlus1, int nVar);
    Arrayofdouble Minbee2(Arrayofdouble ui, Arrayofdouble uiMinus1, Arrayofdouble uiPlus1, int nVar);
    Arrayofdouble Vanleer(Arrayofdouble ui, Arrayofdouble uiMinus1, Arrayofdouble uiPlus1, int nVar);
    Arrayofdouble Superbee(Arrayofdouble ui, Arrayofdouble uiMinus1, Arrayofdouble uiPlus1, int nVar);
    void secondorder_extension(Arrayofdouble u_n, Arrayofdouble u, Arrayofdouble u_p, double gamma_a, double gamma_b, double pinf_a, double pinf_b, Arrayofdouble& primiL, Arrayofdouble& primiR, int nVar, double cell_size, double sl_const);
    void halftimestep_update(Arrayofdouble primi, Arrayofdouble primiL, Arrayofdouble primiR, Arrayofdouble yprimi, Arrayofdouble yprimiL, Arrayofdouble yprimiR, double dt, double dx, double dy, Arrayofdouble& uiL_nplushalf, Arrayofdouble& uiR_nplushalf, double gamma_a, double gamma_b, double pinf_a, double pinf_b);
    void yhalftimestep_update(Arrayofdouble primi, Arrayofdouble primiL, Arrayofdouble primiR, Arrayofdouble yprimi, Arrayofdouble yprimiL, Arrayofdouble yprimiR, double dt, double dx, double dy, Arrayofdouble& uiL_nplushalf, Arrayofdouble& uiR_nplushalf, double gamma_a, double gamma_b, double pinf_a, double pinf_b);
};



