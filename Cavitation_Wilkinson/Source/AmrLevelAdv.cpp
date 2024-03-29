
#include <AmrLevelAdv.H>
#include <Adv_F.H>
#include <AMReX_VisMF.H>
#include <AMReX_TagBox.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BCUtil.H>

#include "eos_v2.H"
#include "hllc.H"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <iomanip>

typedef std::array<double,14> Arrayofdouble;

using namespace amrex;

//initiating classes
eos* EOS = new eos();
numerical_method* NM = new numerical_method();

int      AmrLevelAdv::verbose         = 0;//# VERBOSITY - controls the number of messages output to screen
Real     AmrLevelAdv::cfl             = 0.9; // Default value - can be overwritten in settings file
int      AmrLevelAdv::do_reflux       = 1;  

int      AmrLevelAdv::NUM_STATE       = 11;  // One variable in the state
int      AmrLevelAdv::NUM_GROW        = 2;  // number of ghost cells = 2 for transmissive boundary condition

// material specific data 
double gamma_a = 0.0;
double pinf_a = 0.0;
double cv_a = 0.0;
double epsinf_a = 0.0;

double gamma_b = 0.0;
double pinf_b = 0.0;
double cv_b = 0.0;
double epsinf_b = 0.0;

double p_max = 0.0;
double sl_const = 0.0;

//
//Default constructor.  Builds invalid object.
//
AmrLevelAdv::AmrLevelAdv ()
{
  // Flux registers store fluxes at patch boundaries to ensure fluxes are conservative between AMR levels
  flux_reg = 0;
}

//
//The basic constructor.
//papa - parent; lev - level
AmrLevelAdv::AmrLevelAdv (Amr&            papa,
     	                  int             lev,
                          const Geometry& level_geom,
                          const BoxArray& bl,
                          const DistributionMapping& dm,
                          Real            time)
  :
  AmrLevel(papa,lev,level_geom,bl,dm,time) 
{
  // Flux registers are only required if AMR is actually used, and if flux fix up is being done (recommended)
  flux_reg = 0;
  if (level > 0 && do_reflux)
  {
    flux_reg = new FluxRegister(grids,dmap,crse_ratio,level,NUM_STATE);
  }
}

//
//The destructor.
//
AmrLevelAdv::~AmrLevelAdv () 
{
    delete flux_reg;
}

//
//Restart from a checkpoint file.
//
// AMReX can save simultion state such
// that if the code crashes, it can be restarted, with different
// settings files parameters if necessary (e.g. to output about the
// point of the crash).
//
void
AmrLevelAdv::restart (Amr&          papa,
	              std::istream& is,
                      bool          bReadSpecial)//istream means input; ostream means output
{
  AmrLevel::restart(papa,is,bReadSpecial);
  
  BL_ASSERT(flux_reg == 0);//for debugging, refer to Assert header file for details
  if (level > 0 && do_reflux)
  {
    flux_reg = new FluxRegister(grids,dmap,crse_ratio,level,NUM_STATE);
  }
  
}

//
// Write a checkpoint file - format is handled automatically by AMReX
void 
AmrLevelAdv::checkPoint (const std::string& dir,
		         std::ostream&      os,
                         VisMF::How         how,
                         bool               dump_old) 
{
  AmrLevel::checkPoint(dir, os, how, dump_old);
}

//
//Write a plotfile to specified directory - format is handled automatically by AMReX.
//
void
AmrLevelAdv::writePlotFile (const std::string& dir,
	 	            std::ostream&      os,
                            VisMF::How         how)//dir - directory
{
  AmrLevel::writePlotFile (dir,os,how);//ensure base class function is still called regardless of being overwridden

  //outputing data
  std::ofstream output;
  output.precision(15);
  std::ofstream dc;
  output.open(dir + "amrex_adv.dat", std::ios_base::app);
  dc.open(dir + "conv.txt", std::ios_base::app);   
  const Real* prob_lo = geom.ProbLo();
  const Real probLoX = prob_lo[0]; 
  const Real probLoY = (amrex::SpaceDim > 1 ? prob_lo[1] : 0.0);
  const Real* dx  = geom.CellSize();
  const Real dX = dx[0];
  const Real dY = (amrex::SpaceDim > 1 ? dx[1] : 0.0);//if 1D, dY = 0

  MultiFab& S_new = get_new_data(Phi_Type);

  // 2020H: amended to output data for 1-D plot
  // Loop over all the patches at this level
  for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi)
  {
    const Box& bx = mfi.tilebox();
    const Dim3 lo = lbound(bx);
    const Dim3 hi = ubound(bx);

    const auto& arr = S_new.array(mfi);

    for(int k = lo.z; k <= hi.z; k++)
    {
      for(int j = lo.y; j <= hi.y; j++)
      {
        for(int i = lo.x; i <= hi.x; i++)
        {
          double alpha_a = arr(i,j,k,0); // volume fraction of material a
          double density = arr(i,j,k,1) + arr(i,j,k,2); // mixture density
          double v = 1.0 / density; // specific volume of the mixture
          
          double vx = arr(i,j,k,3) / density; // normal velocity vx
          double vy = arr(i,j,k,5) / density; // normal velocity vy

          double energy = arr(i,j,k,4) / density;//*EOS->compute_energy(alpha_a, rho_a, eps_a, rho_b, eps_b, vx, vy);
          double pressure = EOS->compute_pressure(arr(i,j,k,0), arr(i,j,k,1), arr(i,j,k,2), energy, vx, vy, gamma_a, pinf_a, gamma_b, pinf_b, p_max, epsinf_a, epsinf_b, cv_a);

          double lambda_liquid;
          if(pressure < 22.064e6)
          {
            Arrayofdouble eos_output;

            eos_output = EOS->eos_interpolation(pressure);

            double vl_star = eos_output[1];
            double vg_star = eos_output[2];

            if(v >= vl_star && v <= vg_star)
            {
              lambda_liquid = EOS->compute_massfraction(arr(i,j,k,1), arr(i,j,k,2), pressure);
            }
            else
            {
              lambda_liquid = 1.0; // specific to rarefaction advection test
            }
          }
          else
          {
            lambda_liquid = 1.0; // specific to rarefaction advection test
          }

          const Real x = probLoX + (double(i)+0.5) * dX;
          const Real y = probLoY + (double(j)+0.5) * dY;
          output<<x<<" "<<y<<" "<<arr(i,j,k,1)<<" "<<arr(i,j,k,2)<<" "<<alpha_a<<" "<<density<<" "<<vx<<" "<<pressure<<" "<<energy<<" "<<lambda_liquid<<" "<<vy<<std::endl;
        }
        output<<std::endl;
      }
    }  
  }
  output<<std::endl;
}

//
//Define data descriptors.
//
// This is how the variables in a simulation are defined.  In the case
// of the advection equation, a single variable, phi, is defined.
//
void
AmrLevelAdv::variableSetUp ()
{
  BL_ASSERT(desc_lst.size() == 0);
  //desc_lst - descriptor list

  // A function which contains all processing of the settings file,
  // setting up initial data, choice of numerical methods and
  // boundary conditions

  read_params();
  
  const int storedGhostZones = 0;
    
  // Setting up a container for a variable, or vector of variables:
  // Phi_Type: Enumerator for this variable type
  // IndexType::TheCellType(): AMReX can support cell-centred and vertex-centred variables (cell centred here)
  // StateDescriptor::Point: Data can be a point in time, or an interval over time (point here)
  // storedGhostZones: Ghost zones can be stored (e.g. for output).  Generally set to zero.
  // NUM_STATE: Number of variables in the variable vector (1 in the case of advection equation)
  // cell_cons_interp: Controls interpolation between levels - cons_interp is good for finite volume
  desc_lst.addDescriptor(Phi_Type,IndexType::TheCellType(),
			 StateDescriptor::Point,storedGhostZones,NUM_STATE,
			 &cell_cons_interp);

  //Set up boundary conditions, all boundaries can be set
  //independently, including for individual variables, but lo (left) and hi (right) are useful ways to
  //store them, for consistent access notation for the boundary
  //locations

// Boundary condition for 1D test 1,2,3
  int lo_bc[amrex::SpaceDim];//SpaceDim - number of spatial dimensions
  int hi_bc[amrex::SpaceDim];
  // AMReX has pre-set BCs, including periodic (int_dir) and transmissive (foextrap)
  for (int i = 0; i < amrex::SpaceDim; ++i) {
    lo_bc[i] = hi_bc[i] = BCType::foextrap;   // AN transmissive boundaries
  }

  // Object for storing all the boundary conditions
  BCRec bc(lo_bc, hi_bc);

  // Set up variable-specific information; needs to be done for each variable in NUM_STATE
  // Phi_Type: Enumerator for the variable type being set
  // 0: Position of the variable in the variable vector.  Single variable for advection.
  // phi: Name of the variable - appears in output to identify what is being plotted
  // bc: Boundary condition object for this variable (defined above)
  // BndryFunc: Function for setting boundary conditions.  For basic BCs, AMReX can handle these automatically

  // 2020H: Added new conservative variables for multi material solver using diffuse interface approach
  StateDescriptor::BndryFunc bndryfunc2 (cns_bcfill);

  desc_lst.setComponent(Phi_Type, 0, "alpha_a", bc, bndryfunc2);
  desc_lst.setComponent(Phi_Type, 1, "alpha_rho_a", bc, bndryfunc2);
  desc_lst.setComponent(Phi_Type, 2, "alpha_rho_b", bc, bndryfunc2);
  desc_lst.setComponent(Phi_Type, 3, "rho_vx", bc, bndryfunc2);
  desc_lst.setComponent(Phi_Type, 4, "rho_energy", bc, bndryfunc2);
  desc_lst.setComponent(Phi_Type, 5, "rho_vy", bc, bndryfunc2);
  desc_lst.setComponent(Phi_Type, 6, "S_star", bc, bndryfunc2);
  desc_lst.setComponent(Phi_Type, 7, "yS_star", bc, bndryfunc2);
  desc_lst.setComponent(Phi_Type, 8, "rho", bc, bndryfunc2);
  desc_lst.setComponent(Phi_Type, 9, "vx", bc, bndryfunc2);
  desc_lst.setComponent(Phi_Type, 10, "Pressure", bc, bndryfunc2);
}

//
//Cleanup data descriptors at end of run.
//
void
AmrLevelAdv::variableCleanUp () 
{
    desc_lst.clear();
}

//
//Initialize grid data at problem start-up.
//
void
AmrLevelAdv::initData ()
{
  //
  // Loop over grids, call FORTRAN function to init with data.
  //
  const Real* dx  = geom.CellSize();//geom is an object within Amrlevel; Amrleveladv is a child class of amr level
  //geom its own is a class which cell size is a member function of geom class
  // Position of the bottom left corner of the domain
  const Real* prob_lo = geom.ProbLo();

  // Create a multifab which can store the initial data
  MultiFab& S_new = get_new_data(Phi_Type);
  Real cur_time   = state[Phi_Type].curTime();

  // amrex::Print works like std::cout, but in parallel only prints from the root processor
  if (verbose) {
    amrex::Print() << "Initializing the data at level " << level << std::endl;
  }

  // Slightly messy way to ensure uninitialised data is not used.
  // AMReX has an XDim3 object, but a function needs to be written to
  // convert Real* to XDim3
  const Real dX = dx[0];
  const Real dY = (amrex::SpaceDim > 1 ? dx[1] : 0.0);//if 1D, dY = 0
  const Real dZ = (amrex::SpaceDim > 2 ? dx[2] : 0.0);//if 1D, dZ = 0

  const Real probLoX = prob_lo[0];
  const Real probLoY = (amrex::SpaceDim > 1 ? prob_lo[1] : 0.0);
  const Real probLoZ = (amrex::SpaceDim > 2 ? prob_lo[2] : 0.0);

  // 2020H: definition of initial data
  // primitive variable left
  double alpha_al = 0.0; // primitive variable 0
  double rho_al = 0.0; // primitive variable 1
  double rho_bl = 0.0; // primitive variable 2
  double vx_l = 0.0; // primitive variable 5
  double p_l = 0.0; // primitive variable 6
  double vy_l = 0.0; // primitive variable 7

  // primitive variable right
  double alpha_ar = 0.0; // primitive variable 0
  double rho_ar = 0.0;// primitive variable 1
  double rho_br = 0.0; // primitive variable 2
  double vx_r = 0.0; // primitive variable 5
  double p_r = 0.0; // primitive variable 6
  double vy_r = 0.0; // primitive variable 7

  ParmParse pp;

  pp.query("alpha_al", alpha_al);
  pp.query("rho_al", rho_al);
  pp.query("rho_bl", rho_bl);
  pp.query("vx_l", vx_l);
  pp.query("p_l", p_l);
  pp.query("vy_l", vy_l);

  pp.query("alpha_ar", alpha_ar);
  pp.query("rho_ar", rho_ar);
  pp.query("rho_br", rho_br);
  pp.query("vx_r", vx_r);
  pp.query("p_r", p_r);
  pp.query("vy_r", vy_r);

  double alpha_bl = 1.0 - alpha_al;
  double alpha_br = 1.0 - alpha_ar;

  // defining initial condition in conservative variable
  // Loop over all the patches at this level
  for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
  {
    Box bx = mfi.tilebox();
    const Dim3 lo = lbound(bx);
    const Dim3 hi = ubound(bx);

    const auto& arr = S_new.array(mfi);

    for(int k = lo.z; k <= hi.z; k++)
    {
      const Real z = probLoZ + (double(k)+0.5) * dZ;
      for(int j = lo.y; j <= hi.y; j++)
      {
        const Real y = probLoY + (double(j)+0.5) * dY;
        for(int i = lo.x; i <= hi.x; i++)
        {
          const Real x = probLoX + (double(i)+0.5) * dX;

          // Initial condition for test 1, 2, 3, 3.2
          if(x<1.0)//(sqrt(pow((x-1),2)+pow((y-1),2))<0.4)
          {
            arr(i,j,k,0) = alpha_al; // alpha_a
            arr(i,j,k,1) = alpha_al * rho_al; // alpha_a * rho_a
            arr(i,j,k,2) = alpha_bl * rho_bl; // alpha_b * rho_b

            // rho * v_x
            double rho_l = alpha_al * rho_al + alpha_bl * rho_bl;
            arr(i,j,k,3) = rho_l * vx_l;

            // rho * energy
            double energy_l = EOS->compute_energy(arr(i,j,k,0), arr(i,j,k,1), arr(i,j,k,2), p_l, vx_l, vy_l, gamma_a, pinf_a, gamma_b, pinf_b, p_max, epsinf_a, epsinf_b);
            arr(i,j,k,4) = rho_l * energy_l;

            // rho * v_y
            arr(i,j,k,5) = rho_l * vy_l;
          }
          else
          {
            arr(i,j,k,0) = alpha_ar; // alpha_a
            arr(i,j,k,1) = alpha_ar * rho_ar; // alpha_a * rho_a
            arr(i,j,k,2) = alpha_br * rho_br; // alpha_b * rho_b

            // rho * v_x
            double rho_r = alpha_ar * rho_ar + alpha_br * rho_br;
            arr(i,j,k,3) = rho_r * vx_r;

            // rho * energy
            double energy_r = EOS->compute_energy(arr(i,j,k,0), arr(i,j,k,1), arr(i,j,k,2), p_r, vx_r, vy_r, gamma_a, pinf_a, gamma_b, pinf_b, p_max, epsinf_a, epsinf_b);
            arr(i,j,k,4) = rho_r * energy_r; 

            // rho * v_y
            arr(i,j,k,5) = rho_r * vy_r;
          }

            // initialising rest of the NUM_STATE
            arr(i,j,k,6) = -1e5;
            arr(i,j,k,7) = -1e5;
            arr(i,j,k,8) = -1e5;
            arr(i,j,k,9) = -1e5;
            arr(i,j,k,10) = -1e5;

        }
      }
    }
  }

  if (verbose) {
    amrex::Print() << "Done initializing the level " << level 
		   << " data " << std::endl;
  }
}

//
//Initialize data on this level from another AmrLevelAdv (during regrid).
// These are standard AMReX commands which are unlikely to need altering
//
void
AmrLevelAdv::init (AmrLevel &old)
{
  
  AmrLevelAdv* oldlev = (AmrLevelAdv*) &old;
  //
  // Create new grid data by fillpatching from old.
  //
  Real dt_new    = parent->dtLevel(level);
  Real cur_time  = oldlev->state[Phi_Type].curTime();
  Real prev_time = oldlev->state[Phi_Type].prevTime();
  Real dt_old    = cur_time - prev_time;
  setTimeLevel(cur_time,dt_old,dt_new);
  
  MultiFab& S_new = get_new_data(Phi_Type);

  const int zeroGhosts = 0;
  // FillPatch takes the data from the first argument (which contains
  // all patches at a refinement level) and fills (copies) the
  // appropriate data onto the patch specified by the second argument:
  // old: Source data
  // S_new: destination data
  // zeroGhosts: If this is non-zero, ghost zones could be filled too - not needed for init routines
  // cur_time: AMReX can attempt interpolation if a different time is specified - not recommended for advection eq.
  // Phi_Type: Specify the type of data being set
  // 0: This is the first data index that is to be copied
  // NUM_STATE: This is the number of states to be copied
  FillPatch(old, S_new, zeroGhosts, cur_time, Phi_Type, 0, NUM_STATE);

  // Note: In this example above, the all states in Phi_Type (which is
  // only 1 to start with) are being copied.  However, the FillPatch
  // command could be used to create a velocity vector from a
  // primitive variable vector.  In this case, the `0' argument is
  // replaced with the position of the first velocity component in the
  // primitive variable vector, and the NUM_STATE arguement with the
  // dimensionality - this argument is the number of variables that
  // are being filled/copied, and NOT the position of the final
  // component in e.g. the primitive variable vector.
}

//
//Initialize data on this level after regridding if old level did not previously exist
// These are standard AMReX commands which are unlikely to need altering
//
void
AmrLevelAdv::init ()
{
    Real dt        = parent->dtLevel(level);
    Real cur_time  = getLevel(level-1).state[Phi_Type].curTime();
    Real prev_time = getLevel(level-1).state[Phi_Type].prevTime();

    Real dt_old = (cur_time - prev_time)/(Real)parent->MaxRefRatio(level-1);

    setTimeLevel(cur_time,dt_old,dt);
    MultiFab& S_new = get_new_data(Phi_Type);

    // See first init function for documentation
    FillCoarsePatch(S_new, 0, cur_time, Phi_Type, 0, NUM_STATE);
}

//
//Advance grids at this level in time.
//  This function is the one that actually calls the flux functions.
//
Real
AmrLevelAdv::advance (Real time,
                      Real dt,
                      int  iteration,
                      int  ncycle)
{

  MultiFab& S_mm = get_new_data(Phi_Type);//S max and min

  // Note that some useful commands exist - the maximum and minumum
  // values on the current level can be computed directly - here the
  // max and min of variable 0 are being calculated, and output.
  Real maxval = S_mm.max(0);
  Real minval = S_mm.min(0);
  amrex::Print() << "phi max = " << maxval << ", min = " << minval  << std::endl;

  // This ensures that all data computed last time step is moved from
  // `new' data to `old data' - this should not need changing If more
  // than one type of data were declared in variableSetUp(), then the
  // loop ensures that all of it is updated appropriately
  for (int k = 0; k < NUM_STATE_TYPE; k++) {
    state[k].allocOldData();
    state[k].swapTimeLevels(dt);//different dt needed for level
  }

  // S_new is the MultiFab that will be operated upon to update the data
  MultiFab& S_new = get_new_data(Phi_Type);

  const Real prev_time = state[Phi_Type].prevTime();
  const Real cur_time = state[Phi_Type].curTime();
  const Real ctr_time = 0.5*(prev_time + cur_time);

  const Real* dx = geom.CellSize();
  const Real* prob_lo = geom.ProbLo();

  //
  // Get pointers to Flux registers, or set pointer to zero if not there.
  //
  FluxRegister *fine    = 0;
  FluxRegister *current = 0;
    
  int finest_level = parent->finestLevel();

  // If we are not on the finest level, fluxes may need correcting
  // from those from finer levels.  To start this process, we set the
  // flux register values to zero
  if (do_reflux && level < finest_level) {
    fine = &getFluxReg(level+1);
    fine->setVal(0.0);
  }

  // If we are not on the coarsest level, the fluxes are going to be
  // used to correct those on coarser levels.  We get the appropriate
  // flux level to include our fluxes within
  if (do_reflux && level > 0)
  {
    current = &getFluxReg(level);
  }

  // Set up a dimensional multifab that will contain the fluxes
  MultiFab fluxes[amrex::SpaceDim];

  // Define the appropriate size for the flux MultiFab.
  // Fluxes are defined at cell faces - this is taken care of by the
  // surroundingNodes(j) command, ensuring the size of the flux
  // storage is increased by 1 cell in the direction of the flux.
  // This is only needed if refluxing is happening, otherwise fluxes
  // don't need to be stored, just used
  for (int j = 0; j < amrex::SpaceDim; j++)
  {
    BoxArray ba = S_new.boxArray();
    ba.surroundingNodes(j);
    fluxes[j].define(ba, dmap, NUM_STATE, 0);
  }

  // Advection velocity - AMReX allows the defintion of a vector
  // object (similar functionality to C++ std::array<N>, since its size must
  // be known, but was implemented before array was added to C++)
  //const Vector<Real> vel{1.0,1.0,0.0};
  //const Vector<Real> vel{0.0,0.0,0.0};

  // State with ghost cells - this is used to compute fluxes and perform the update.
  MultiFab Sborder(grids, dmap, NUM_STATE, NUM_GROW);
  // See init function for details about the FillPatch function
  FillPatch(*this, Sborder, NUM_GROW, time, Phi_Type, 0, NUM_STATE);//*this is source data, Sborder is the destination
  // Fill periodic boundaries where they exist.  More accurately, the
  // FillBoundary call will fill overlapping boundaries (with periodic
  // domains effectively being overlapping).  It also takes care of
  // AMR patch and CPU boundaries.

  // 2020H: Defining additional multifabs for storing the x reconstructed values and y reconstructed values
  MultiFab Sborder_xrecon_l(grids, dmap, NUM_STATE, NUM_GROW);
  FillPatch(*this, Sborder_xrecon_l, NUM_GROW, time, Phi_Type, 0, NUM_STATE);
  Sborder_xrecon_l = 0.0;

  MultiFab Sborder_xrecon_r(grids, dmap, NUM_STATE, NUM_GROW);
  FillPatch(*this, Sborder_xrecon_r, NUM_GROW, time, Phi_Type, 0, NUM_STATE);
  Sborder_xrecon_r = 0.0;


  MultiFab Sborder_yrecon_l(grids, dmap, NUM_STATE, NUM_GROW);
  FillPatch(*this, Sborder_yrecon_l, NUM_GROW, time, Phi_Type, 0, NUM_STATE);
  Sborder_yrecon_l = 0.0;

  MultiFab Sborder_yrecon_r(grids, dmap, NUM_STATE, NUM_GROW);
  FillPatch(*this, Sborder_yrecon_r, NUM_GROW, time, Phi_Type, 0, NUM_STATE);
  Sborder_yrecon_r = 0.0;


  MultiFab Sborder_xrecon_lht(grids, dmap, NUM_STATE, NUM_GROW);
  FillPatch(*this, Sborder_xrecon_lht, NUM_GROW, time, Phi_Type, 0, NUM_STATE);
  Sborder_xrecon_lht = 0.0;

  MultiFab Sborder_yrecon_lht(grids, dmap, NUM_STATE, NUM_GROW);
  FillPatch(*this, Sborder_yrecon_lht, NUM_GROW, time, Phi_Type, 0, NUM_STATE);
  Sborder_yrecon_lht = 0.0;


  MultiFab Sborder_xrecon_rht(grids, dmap, NUM_STATE, NUM_GROW);
  FillPatch(*this, Sborder_xrecon_rht, NUM_GROW, time, Phi_Type, 0, NUM_STATE);
  Sborder_xrecon_rht = 0.0;

  MultiFab Sborder_yrecon_rht(grids, dmap, NUM_STATE, NUM_GROW);
  FillPatch(*this, Sborder_yrecon_rht, NUM_GROW, time, Phi_Type, 0, NUM_STATE);
  Sborder_yrecon_rht = 0.0;

  Vector<BCRec> ArrayofBC(NUM_STATE);

  // 2020H: defining boundary condition for test 1,2,3
  for(int l=0; l<NUM_STATE; ++l)//++i and i++ are the same in a for loop
  {
    for(int d=0; d<amrex::SpaceDim; ++d)
    {
      ArrayofBC[l].setLo(d, BCType::foextrap);
      ArrayofBC[l].setHi(d, BCType::foextrap);
    }
  }

  Sborder.FillBoundary(geom.periodicity());
  FillDomainBoundary(Sborder, geom, ArrayofBC);

  // 2020H: Introduced new function calls for extending system to second order accurate
  // Loop over all the patches at this level
  for (MFIter mfi(Sborder, true); mfi.isValid(); ++mfi)
  {
    const Box& bx = mfi.tilebox();
    const Dim3 lo = lbound(bx);
    const Dim3 hi = ubound(bx);

    // Indexable arrays for the data, and the directional flux
    // Based on the vertex-centred definition of the flux array, the
    // data array runs from e.g. [0,N] and the flux array from [0,N+1]
    const auto& arr = Sborder.array(mfi);
    const auto& fluxArr = fluxes[0].array(mfi);
    const auto& yfluxArr = fluxes[1].array(mfi);

    const auto& arr_xrecon_l = Sborder_xrecon_l.array(mfi);
    const auto& arr_xrecon_r = Sborder_xrecon_r.array(mfi);

    const auto& arr_yrecon_l = Sborder_yrecon_l.array(mfi);
    const auto& arr_yrecon_r = Sborder_yrecon_r.array(mfi);

    // computing xflux
    Arrayofdouble flux, wavespeed;

    double nVar = 6;

    // second order reconstruction in x-direction
    for(int k = lo.z; k <= hi.z; k++)
    {
      for(int j = lo.y; j <= hi.y+2; j++)
      {
        for(int i = lo.x; i <= hi.x+2; i++)
        {
          //defining u_n, u, u_p for input into secondorder_extension function
          Arrayofdouble ul_n, ul, ul_p;
          
          for(int l=0; l<nVar; l++) // skip energy, vy and only pass in 6 variables into second order extension function
          {
            // reconstructing in x-direction from (-1,-1) to (8,8)
            //ui-1
            ul_n[l] = arr(i-2,j-1,k,l);

            //ui
            ul[l] = arr(i-1,j-1,k,l);

            //ui+1
            ul_p[l] = arr(i,j-1,k,l);
          }

          // reconstructing i-1 cell
          Arrayofdouble primiL, primiR;
          NM->secondorder_extension(ul_n, ul, ul_p, gamma_a, gamma_b, pinf_a, pinf_b, primiL, primiR, nVar, dx[0], sl_const, p_max, epsinf_a, epsinf_b, cv_a);
        
          for(int m=0; m<nVar; m++) // storing reconstructed values as the same index as cell as number of cells equal to number of left and right wall
          {
          arr_xrecon_l(i-1,j-1,k,m) = primiL[m]; // left wall of cell (i-1, j-1)
          arr_xrecon_r(i-1,j-1,k,m) = primiR[m]; // right wall of cell (i-1, j-1)
          }
        }
      }
    }

    // second order reconstruction in y-direction
    for(int k = lo.z; k <= hi.z; k++)
    {
      for(int j = lo.y; j <= hi.y+2; j++)
      {
        for(int i = lo.x; i <= hi.x+2; i++)
        {
          //defining u_n, u, u_p for input into secondorder_extension function
          Arrayofdouble yul_n, yul, yul_p;

          for(int l=0; l<nVar; l++) // skip energy, vy and only pass in 6 variables into second order extension function
          {
            // reconstructing in y-direction from (-1,-1) to (8,8)
            //ui-1
            yul_n[l] = arr(i-1,j-2,k,l);

            //ui
            yul[l] = arr(i-1,j-1,k,l);

            //ui+1
            yul_p[l] = arr(i-1,j,k,l);
          }

          // reconstructing j-1 cell
          Arrayofdouble yprimiL, yprimiR;
          NM->secondorder_extension(yul_n, yul, yul_p, gamma_a, gamma_b, pinf_a, pinf_b, yprimiL, yprimiR, nVar, dx[1], sl_const, p_max, epsinf_a, epsinf_b, cv_a);
        
          for(int m=0; m<nVar; m++)
          {
          arr_yrecon_l(i-1,j-1,k,m) = yprimiL[m]; // left (bottom) wall of cell (i-1, j-1)
          arr_yrecon_r(i-1,j-1,k,m) = yprimiR[m]; // right (top) wall of cell (i-1, j-1)
          }
        
        }
      }
    }
  } // close MFIter loop here

  // 2020H: second order update
  // Loop over all the patches at this level
  for (MFIter mfi(Sborder, true); mfi.isValid(); ++mfi)
  {
    const Box& bx = mfi.tilebox();

    const Dim3 lo = lbound(bx);
    const Dim3 hi = ubound(bx);

    // Indexable arrays for the data, and the directional flux
    // Based on the vertex-centred definition of the flux array, the
    // data array runs from e.g. [0,N] and the flux array from [0,N+1]
    const auto& arr = Sborder.array(mfi);

    // HHN21: PVMH
    const auto& arr_xrecon_l = Sborder_xrecon_l.array(mfi);
    const auto& arr_xrecon_r = Sborder_xrecon_r.array(mfi);

    const auto& arr_yrecon_l = Sborder_yrecon_l.array(mfi);
    const auto& arr_yrecon_r = Sborder_yrecon_r.array(mfi);

    const auto& arr_xrecon_lht = Sborder_xrecon_lht.array(mfi);
    const auto& arr_xrecon_rht = Sborder_xrecon_rht.array(mfi);

    const auto& arr_yrecon_lht = Sborder_yrecon_lht.array(mfi);
    const auto& arr_yrecon_rht = Sborder_yrecon_rht.array(mfi);

    double nVar = 6;

      // half time step evolution
      for(int k = lo.z; k <= hi.z; k++)
      {
        for(int j = lo.y; j <= hi.y+2; j++) // previous dimensional loop already closed, still 1 for half time step evolution
        {
          for(int i = lo.x; i <= hi.x+2; i++) // previous dimensional loop already closed, still 1 for half time step evolution
          {

            Arrayofdouble ui, primi, primiL, primiR, yprimi, yprimiL, yprimiR, uiL_nplushalf, uiR_nplushalf, yuiL_nplushalf, yuiR_nplushalf;

            for(int l=0; l<nVar; l++)
            {
              ui[l] = arr(i-1,j-1,k,l);
            }

            primi = EOS->u_to_prim(ui, gamma_a, pinf_a, gamma_b, pinf_b, p_max, epsinf_a, epsinf_b, cv_a);

            for(int l = 0; l<nVar; l++)
            {
              primiL[l] = arr_xrecon_l(i-1,j-1,k,l);
              primiR[l] = arr_xrecon_r(i-1,j-1,k,l);

              yprimi[l] = primi[l];
              yprimiL[l] = arr_yrecon_l(i-1,j-1,k,l);
              yprimiR[l] = arr_yrecon_r(i-1,j-1,k,l);
            }

            // half time step evolution for left, right reconstructed variables (all variables except energy)
            NM->halftimestep_update(primi, primiL, primiR, yprimi, yprimiL, yprimiR, dt, dx[0], dx[1], uiL_nplushalf, uiR_nplushalf, gamma_a, gamma_b, pinf_a, pinf_b, p_max, epsinf_a, epsinf_b, cv_a, cv_b);

            for(int l=0; l<nVar; l++)
            {
              arr_xrecon_lht(i-1,j-1,k,l) = uiL_nplushalf[l];
              arr_xrecon_rht(i-1,j-1,k,l) = uiR_nplushalf[l];
            }

            // half time step evolution for top, bottom reconstructed variables (all variables except energy)
            NM->yhalftimestep_update(primi, primiL, primiR, yprimi, yprimiL, yprimiR, dt, dx[0], dx[1], yuiL_nplushalf, yuiR_nplushalf, gamma_a, gamma_b, pinf_a, pinf_b, p_max, epsinf_a, epsinf_b, cv_a, cv_b);

            for(int l=0; l<nVar; l++)
            {
              arr_yrecon_lht(i-1,j-1,k,l) = yuiL_nplushalf[l];
              arr_yrecon_rht(i-1,j-1,k,l) = yuiR_nplushalf[l];
            }
          }
        }
      }
  }

  // 2020H: Introduced new function calls for computing flux and conservative update
  //conservative update from one dimension to another
  for (int d = 0; d < amrex::SpaceDim ; d++)
  {
    const int iOffset = ( d == 0 ? 1 : 0);//iOffset = 1 if d is 0 (one dimensional), else iOffset = 0
    const int jOffset = ( d == 1 ? 1 : 0);
    const int kOffset = ( d == 2 ? 1 : 0);
  
    // Loop over all the patches at this level
    for (MFIter mfi(Sborder, true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();

      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);

      // Indexable arrays for the data, and the directional flux
      // Based on the vertex-centred definition of the flux array, the
      // data array runs from e.g. [0,N] and the flux array from [0,N+1]
      const auto& arr = Sborder.array(mfi);
      const auto& arr_xrecon_lht = Sborder_xrecon_lht.array(mfi);
      const auto& arr_xrecon_rht = Sborder_xrecon_rht.array(mfi);
      const auto& arr_yrecon_lht = Sborder_yrecon_lht.array(mfi);
      const auto& arr_yrecon_rht = Sborder_yrecon_rht.array(mfi);
      const auto& fluxArr = fluxes[0].array(mfi);
      const auto& yfluxArr = fluxes[1].array(mfi);

      // computing xflux
      Arrayofdouble flux, wavespeed;

      if(d == 0)
      {
      double nVar = 6;
      // x flux computation
      for(int k = lo.z; k <= hi.z+kOffset; k++)
      {
        for(int j = lo.y; j <= hi.y+jOffset; j++)
        {
          for(int i = lo.x; i <= hi.x+iOffset; i++)
          {

            // x flux computation
            NM->compute_HLLCflux(arr, arr_xrecon_lht, arr_xrecon_rht, flux, wavespeed, nVar, gamma_a, pinf_a, gamma_b, pinf_b, i, j ,k, dx[d], dt, p_max, epsinf_a, epsinf_b, cv_a, cv_b);
            
            for(int l = 0; l<nVar; l++)
            {
              fluxArr(i,j,k,l) = flux[l];
            }

            fluxArr(i,j,k,6) = wavespeed[2]; // intermediate wave speed     
          }
        }
      }
      } 

      else if(d == 1)
      {
      double nVar = 6;
      // y flux computation
      for(int k = lo.z; k <= hi.z+kOffset; k++)
      {
        for(int j = lo.y; j <= hi.y+jOffset; j++)
        {
          for(int i = lo.x; i <= hi.x+iOffset; i++)
          {
            // computing yflux
            Arrayofdouble yflux, ywavespeed;
            
            NM->ycompute_HLLCflux(arr, arr_yrecon_lht, arr_yrecon_rht, yflux, ywavespeed, nVar, gamma_a, pinf_a, gamma_b, pinf_b, i, j ,k, dx[d], dt, p_max, epsinf_a, epsinf_b, cv_a, cv_b);

            for(int l = 0; l<nVar; l++)
            {
              yfluxArr(i,j,k,l) = yflux[l];
            }

            yfluxArr(i,j,k,7) = ywavespeed[2]; // intermediate wave speed   
          }
        }
      }
      }
    } // close MFIter loop here
  } // close dimensional loop here for computing xflux and yflux

  // 2020H: Godunov update
  // Loop over all the patches at this level
  for (MFIter mfi(Sborder, true); mfi.isValid(); ++mfi)
  {
    const Box& bx = mfi.tilebox();

    const Dim3 lo = lbound(bx);
    const Dim3 hi = ubound(bx);

    // Indexable arrays for the data, and the directional flux
    // Based on the vertex-centred definition of the flux array, the
    // data array runs from e.g. [0,N] and the flux array from [0,N+1]
    const auto& arr = Sborder.array(mfi);
    const auto& fluxArr = fluxes[0].array(mfi);
    const auto& yfluxArr = fluxes[1].array(mfi);
    const auto& arr_xrecon_lht = Sborder_xrecon_lht.array(mfi);
    const auto& arr_xrecon_rht = Sborder_xrecon_rht.array(mfi);

    double nVar = 6;

      //conservative variable update 
      for(int k = lo.z; k <= hi.z; k++)
      {
        for(int j = lo.y; j <= hi.y; j++)
        {
          for(int i = lo.x; i <= hi.x; i++)
          {

            // defining terms to store variables needed to do conservative update
            double alpha_a = arr(i,j,k,0);
            double alpha_b = 1.0 - alpha_a;
            double rho_a = arr(i,j,k,1) / alpha_a;
            double rho_b = arr(i,j,k,2) / alpha_b;

            double S_star = fluxArr(i,j,k,6);
            double S_star_p = fluxArr(i+1,j,k,6);
            double yS_star = yfluxArr(i,j,k,7);
            double yS_star_p = yfluxArr(i,j+1,k,7);

            double rho = arr(i,j,k,1) + arr(i,j,k,2);
            double vx = arr(i,j,k,3) / rho;

            for(int l = 0; l<nVar; l++)
            {

              if(l == 0)
              {
                double vx_rp = arr_xrecon_rht(i,j,k,3) / (arr_xrecon_rht(i,j,k,1) + arr_xrecon_rht(i,j,k,2));
                double vx_lp = arr_xrecon_lht(i+1,j,k,3) / (arr_xrecon_lht(i+1,j,k,1) + arr_xrecon_lht(i+1,j,k,2));

                double vx_r = arr_xrecon_rht(i-1,j,k,3) / (arr_xrecon_rht(i-1,j,k,1) + arr_xrecon_rht(i-1,j,k,2));
                double vx_l = arr_xrecon_lht(i,j,k,3) / (arr_xrecon_lht(i,j,k,1) + arr_xrecon_lht(i,j,k,2));

                double avgvx_p = (vx_rp + vx_lp) / 2.0;
                double avgvx = (vx_r + vx_l) / 2.0;

                double q_p = 0.5 * (arr_xrecon_rht(i,j,k,0) * vx_rp + arr_xrecon_lht(i+1,j,k,0) * vx_lp) - 0.5 * std::abs(avgvx_p) * (arr_xrecon_rht(i,j,k,0) - arr_xrecon_lht(i+1,j,k,0));
                double q = 0.5 * (arr_xrecon_rht(i-1,j,k,0) * vx_r + arr_xrecon_lht(i,j,k,0) * vx_l) - 0.5 * std::abs(avgvx) * (arr_xrecon_rht(i-1,j,k,0) - arr_xrecon_lht(i,j,k,0));

                arr(i,j,k,l) = arr(i,j,k,l) - (dt / dx[0]) * (q_p - q) + (dt / dx[0]) * alpha_a * (avgvx_p - avgvx);
              }
              else
              {
                arr(i,j,k,l) = arr(i,j,k,l) - (dt / dx[0]) * (fluxArr(i+1,j,k,l) - fluxArr(i,j,k,l)) - (dt / dx[1]) * (yfluxArr(i,j+1,k,l) - yfluxArr(i,j,k,l)); 
              }
            }

            arr(i,j,k,8) = arr(i,j,k,1) + arr(i,j,k,2); // mixture density for visit
            arr(i,j,k,9) = arr(i,j,k,3)/arr(i,j,k,8); // vx for visit

            double alpha_a2 = arr(i,j,k,0);
            double alpha_b2 = 1.0 - arr(i,j,k,0);
            double rho2 = arr(i,j,k,1) + arr(i,j,k,2);

            arr(i,j,k,10) = EOS->compute_pressure(alpha_a2, arr(i,j,k,1), arr(i,j,k,2), arr(i,j,k,4)/rho2, arr(i,j,k,3)/rho2, arr(i,j,k,5)/rho2, gamma_a, pinf_a, gamma_b, pinf_b, p_max, epsinf_a, epsinf_b, cv_a);
          }
        }
      }
    }// close MFIter loop here

    // // We need to compute boundary conditions again after each update
    Sborder.FillBoundary(geom.periodicity());
    FillDomainBoundary(Sborder, geom, ArrayofBC);

    //scaling needed before reflux as variables such as density are point values, need to times it by area
    // The fluxes now need scaling for the reflux command.
    // This scaling is by the size of the boundary through which the flux passes, e.g. the x-flux needs scaling by the dy, dz and dt
    for(int d=0; d<amrex::SpaceDim; ++d)
    {
    if(do_reflux)
    {
      Real scaleFactor = dt;
      for(int scaledir = 0; scaledir < amrex::SpaceDim; ++scaledir)
      {
        // Fluxes don't need scaling by dx[d]
        if(scaledir == d)
        {
          continue;
        }
        scaleFactor *= dx[scaledir];
      }
      // The mult function automatically multiplies entries in a multifab by a scalar
      // scaleFactor: The scalar to multiply by
      // 0: The first data index in the multifab to multiply
      // NUM_STATE:  The total number of data indices that will be multiplied
      fluxes[d].mult(scaleFactor, 0, NUM_STATE);
    }
    }//close 2nd dimension for loop here

  // The updated data is now copied to the S_new multifab.  This means
  // it is now accessible through the get_new_data command, and AMReX
  // can automatically interpolate or extrapolate between layers etc.
  // S_new: Destination
  // Sborder: Source
  // Third entry: Starting variable in the source array to be copied (the zeroth variable in this case)
  // Fourth entry: Starting variable in the destination array to receive the copy (again zeroth here)
  // NUM_STATE: Total number of variables being copied
  // Sixth entry: Number of ghost cells to be included in the copy (zero in this case, since only real
  //              data is needed for S_new)
  MultiFab::Copy(S_new, Sborder, 0, 0, NUM_STATE, 0);

  // Refluxing at patch boundaries.  Amrex automatically does this
  // where needed, but you need to state a few things to make sure it
  // happens correctly:
  // FineAdd: If we are not on the coarsest level, the fluxes at this level will form part of the correction
  //          to a coarse level
  // CrseInit:  If we are not the finest level, the fluxes at patch boundaries need correcting.  Since we
  //            know that the coarse level happens first, we initialise the boundary fluxes through this
  //            function, and subsequently FineAdd will modify things ready for the correction
  // Both functions have the same arguments:
  // First: Name of the flux MultiFab (this is done dimension-by-dimension
  // Second: Direction, to ensure the correct vertices are being corrected
  // Third: Source component - the first entry of the flux MultiFab that is to be copied (it is possible that
  //        some variables will not need refluxing, or will be computed elsewhere (not in this example though)
  // Fourth: Destinatinon component - the first entry of the flux register that this call to FineAdd sends to
  // Fifth: NUM_STATE - number of states being added to the flux register
  // Sixth: Multiplier - in general, the least accurate (coarsest) flux is subtracted (-1) and the most
  //        accurate (finest) flux is added (+1)
  if (do_reflux) {
    if (current) {
      for (int i = 0; i < amrex::SpaceDim ; i++)
	current->FineAdd(fluxes[i],i,0,0,NUM_STATE,1.);
    }
    if (fine) {
      for (int i = 0; i < amrex::SpaceDim ; i++)
	fine->CrseInit(fluxes[i],i,0,0,NUM_STATE,-1.);
    }
  }  
  return dt;
}

//
//Estimate time step.
// This function is called by all of the other time step functions in AMReX, and is the only one that should
// need modifying
//
int counter = 0;

Real
AmrLevelAdv::estTimeStep (Real)
{
  // This is just a dummy value to start with 
  Real dt_est  = 1.0e+20;

  const Real* dx = geom.CellSize();
  const Real* prob_lo = geom.ProbLo();
  const Real cur_time = state[Phi_Type].curTime();
  const MultiFab& S_new = get_new_data(Phi_Type);

  // This should not really be hard coded
  //const Real velMag = sqrt(2.);

  //initiate new multifab to store cs
  //FillPatch(*this, Cs, NUM_GROW, time, Phi_Type, 0, 1);
  double velMag = 0.0;

  // 2020H: Amended to compute time step for Euler system
  //finding amax
  // Loop over all the patches at this level
  for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi)
  {
    const Box& bx = mfi.tilebox();
    const Dim3 lo = lbound(bx);
    const Dim3 hi = ubound(bx);

    const auto& arr = S_new.array(mfi);

    for(int k = lo.z; k <= hi.z; k++)
    {
      for(int j = lo.y; j <= hi.y; j++)
      {
        for(int i = lo.x; i <= hi.x; i++)
        {
          // computing vx
          double alpha_a = arr(i,j,k,0);
          double alpha_b = 1.0 - alpha_a;
          double rho = arr(i,j,k,1) + arr(i,j,k,2);
          double vx = arr(i,j,k,3) / rho;
          double vy = arr(i,j,k,5) / rho;
          double vm = sqrt(vx*vx + vy*vy);

          // computing sound speed
          double p = EOS->compute_pressure(arr(i,j,k,0), arr(i,j,k,1), arr(i,j,k,2), arr(i,j,k,4)/rho, arr(i,j,k,3)/rho, arr(i,j,k,5)/rho, gamma_a, pinf_a, gamma_b, pinf_b, p_max, epsinf_a, epsinf_b, cv_a);
          double cs = EOS->compute_soundspeed(arr(i,j,k,0), arr(i,j,k,1), arr(i,j,k,2), p, gamma_a, pinf_a, gamma_b, pinf_b, p_max, epsinf_a, epsinf_b, cv_a, cv_b);

          if(cs+fabs(vm)>velMag)
          {
            velMag = (cs + fabs(vm));
          }
        }
      }
    }  
  }

  for(unsigned int d = 0; d < amrex::SpaceDim; ++d)
  {
    if(dx[d]/velMag < dt_est)
    {
      dt_est = dx[d]/velMag; // getting the minimum of dx[d] 
    }
  }
  
  // Ensure that we really do have the minimum across all processors
  ParallelDescriptor::ReduceRealMin(dt_est);
  dt_est *= cfl;

  counter ++;

  if(counter<5)
  {
    dt_est = dt_est * 0.2;
  }

  if (verbose) {
    amrex::Print() << "AmrLevelAdv::estTimeStep at level " << level 
		   << ":  dt_est = " << dt_est << std::endl;
  }
  
  return dt_est;
}

//
//Compute initial time step.
//
Real
AmrLevelAdv::initialTimeStep ()
{
  return estTimeStep(0.0);
}

//
//Compute initial `dt'.
//
void
AmrLevelAdv::computeInitialDt (int                   finest_level,
	  	               int                   sub_cycle,
                               Vector<int>&           n_cycle,
                               const Vector<IntVect>& ref_ratio,
                               Vector<Real>&          dt_level,
                               Real                  stop_time)
{
  //
  // Grids have been constructed, compute dt for all levels.
  //
  // AMReX's AMR Level mode assumes that the time step only needs
  // calculating on the coarsest level - all subsequent time steps are
  // reduced by the refinement factor
  if (level > 0)
    return;

  // Initial guess
  Real dt_0 = 1.0e+100;
  int n_factor = 1;
  //computing time step for each fine levels using the coarse level
  for (int i = 0; i <= finest_level; i++)
  {
    dt_level[i] = getLevel(i).initialTimeStep();//initialTimeStep - estTimeStep
    n_factor   *= n_cycle[i];//n_factor = n_factor * n_cycle[i]
    dt_0 = std::min(dt_0,n_factor*dt_level[i]);
  }

  //
  // Limit dt's by the value of stop_time.
  //
  const Real eps = 0.001*dt_0;
  Real cur_time  = state[Phi_Type].curTime();
  if (stop_time >= 0.0) {
    if ((cur_time + dt_0) > (stop_time - eps))
      dt_0 = stop_time - cur_time;
  }

  n_factor = 1;
  for (int i = 0; i <= finest_level; i++)
  {
    n_factor *= n_cycle[i];
    dt_level[i] = dt_0/n_factor;
  }
}

//
//Compute new `dt'.
//
void
AmrLevelAdv::computeNewDt (int                   finest_level,
		           int                   sub_cycle,
                           Vector<int>&           n_cycle,
                           const Vector<IntVect>& ref_ratio,
                           Vector<Real>&          dt_min,
                           Vector<Real>&          dt_level,
                           Real                  stop_time,
                           int                   post_regrid_flag)
{
  //
  // We are at the end of a coarse grid timecycle.
  // Compute the timesteps for the next iteration.
  //
  if (level > 0)
    return;

  // Although we only compute the time step on the finest level, we
  // need to take information from all levels into account.  The
  // sharpest features may be smeared out on coarse levels, so not
  // using finer levels could cause instability
  for (int i = 0; i <= finest_level; i++)
  {
    AmrLevelAdv& adv_level = getLevel(i);
    dt_min[i] = adv_level.estTimeStep(dt_level[i]);
  }

  // A couple of things are implemented to ensure that time step's
  // don't suddenly grow by a lot, as this could lead to errors - for
  // sensible mesh refinement choices, these shouldn't really change
  // anything
  if (post_regrid_flag == 1) 
  {
    //
    // Limit dt's by pre-regrid dt
    //
    for (int i = 0; i <= finest_level; i++)
    {
      dt_min[i] = std::min(dt_min[i],dt_level[i]);
    }
  }
  else 
  {
    //
    // Limit dt's by change_max * old dt
    //
    static Real change_max = 1.1;
    for (int i = 0; i <= finest_level; i++)
    {
      dt_min[i] = std::min(dt_min[i],change_max*dt_level[i]);
    }
  }
    
  //
  // Find the minimum over all levels
  //
  Real dt_0 = 1.0e+100;
  int n_factor = 1;
  for (int i = 0; i <= finest_level; i++)
  {
    n_factor *= n_cycle[i];
    dt_0 = std::min(dt_0,n_factor*dt_min[i]);
  }

  //
  // Limit dt's by the value of stop_time.
  //
  const Real eps = 0.001*dt_0;
  Real cur_time  = state[Phi_Type].curTime();
  if (stop_time >= 0.0) {
    if ((cur_time + dt_0) > (stop_time - eps))
      dt_0 = stop_time - cur_time;
  }
  
  n_factor = 1;
  for (int i = 0; i <= finest_level; i++)
  {
    n_factor *= n_cycle[i];
    dt_level[i] = dt_0/n_factor;
  }
}

//
//Do work after timestep().
// If something has to wait until all processors have done their advance function, the post_timestep function
// is the place to put it.  Refluxing and averaging down are two standard examples for AMR
//
void
AmrLevelAdv::post_timestep (int iteration)
{
  //
  // Integration cycle on fine level grids is complete
  // do post_timestep stuff here.
  //
  int finest_level = parent->finestLevel();
  
  if (do_reflux && level < finest_level)
    reflux();
  
  if (level < finest_level)
    avgDown();
}

//
//Do work after regrid().
// Nothing normally needs doing here, but if something was calculated on a per-patch basis, new patches might
// this to be calcuated immediately
//
void
AmrLevelAdv::post_regrid (int lbase, int new_finest)
{

}

//
//Do work after a restart().
// Similar to post_regrid, nothing normally needs doing here
//
void
AmrLevelAdv::post_restart() 
{

}

//
//Do work after init().
// Once new patches have been initialised, work may need to be done to ensure consistency, for example,
// averaging down - though for linear interpolation, this probably won't change anything
//
void
AmrLevelAdv::post_init (Real stop_time)
{
  if (level > 0)
    return;
  //
  // Average data down from finer levels
  // so that conserved data is consistent between levels.
  //
  int finest_level = parent->finestLevel();
  for (int k = finest_level-1; k>= 0; k--)
    getLevel(k).avgDown();
}

//
//Error estimation for regridding.
//  Determine which parts of the domain need refinement
//
void
AmrLevelAdv::errorEst (TagBoxArray& tags,
	               int          clearval,
                       int          tagval,
                       Real         time,
                       int          n_error_buf,
                       int          ngrow)
{
  const Real* dx        = geom.CellSize();
  const Real* prob_lo   = geom.ProbLo();

  MultiFab& S_new = get_new_data(Phi_Type);

  Vector<int> itags;
	
  for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
  {
    const Box&  tilebx  = mfi.tilebox();

    // An AMReX construction, effectively a boolean array which is true in positions that are valid for refinement
    TagBox&     tagfab  = tags[mfi];

    // Traditionally, a lot of the array-based operations in AMReX happened in Fortran.  The standard template
    // for these is short and easy to read, flagging on values or gradients (first order calculation)
    // We cannot pass tagfab to Fortran becuase it is BaseFab<char>.
    // So we are going to get a temporary integer array.
    tagfab.get_itags(itags, tilebx);
	    
    // data pointer and index space
    int*        tptr    = itags.dataPtr();
    const int*  tlo     = tilebx.loVect();
    const int*  thi     = tilebx.hiVect();

    // Various macros exist to convert the C++ data structures to Fortran
    state_error(tptr,  AMREX_ARLIM_3D(tlo), AMREX_ARLIM_3D(thi),
		BL_TO_FORTRAN_3D(S_new[mfi]),
		&tagval, &clearval, 
		AMREX_ARLIM_3D(tilebx.loVect()), AMREX_ARLIM_3D(tilebx.hiVect()), 
		AMREX_ZFILL(dx), AMREX_ZFILL(prob_lo), &time, &level);
    //
    // Now update the tags in the TagBox.
    //
    tagfab.tags_and_untags(itags, tilebx);
  }
}

//
// This function reads the settings file
//
void
AmrLevelAdv::read_params ()
{
  // Make sure that this is only done once
  static bool done = false;

  if (done) return;

  done = true;

  // A ParmParse object allows settings, with the correct prefix, to be read in from the settings file
  // The prefix can help identify what a settings parameter is used for
  // AMReX has some default ParmParse names, amr and geometry are two commonly needed ones
  ParmParse pp("adv");   

  // ParmParse has two options; query and get.  Query will only alter
  // a parameter if it can be found (if these aren't in the settings
  // file, then the values at the top of this file will be used).  Get
  // will throw an error if the parameter is not found in the settings
  // file.

  pp.query("v",verbose);
  pp.query("cfl",cfl);
  pp.query("do_reflux",do_reflux);

  pp.query("gamma_a",gamma_a);
  pp.query("pinf_a",pinf_a);
  pp.query("cv_a",cv_a);
  pp.query("epsinf_a",epsinf_a);

  pp.query("gamma_b",gamma_b);
  pp.query("pinf_b",pinf_b);
  pp.query("cv_b",cv_b);
  pp.query("epsinf_b",epsinf_b);

  pp.query("p_max",p_max);
  pp.query("sl_const",sl_const);

  // Vector variables can be read in; these require e.g.\ pp.queryarr
  // and pp.getarr, so that the ParmParse object knows to look for
  // more than one variable

  // Geometries can be Cartesian, cylindrical or spherical - some
  // functions (e.g. divergence in linear solvers) are coded with this
  // geometric dependency
  Geometry const* gg = AMReX::top()->getDefaultGeometry();

  // This tutorial code only supports Cartesian coordinates.
  if (! gg->IsCartesian()) {
    amrex::Abort("Please set geom.coord_sys = 0");
  }

  // This tutorial code only supports periodic boundaries.
  // The periodicity is read from the settings file in AMReX source code, but can be accessed here
  // if (! gg->isAllPeriodic()) {
  //   amrex::Abort("Please set geometry.is_periodic = 1 1 1");
  // }

  //
  // read tagging parameters from probin file
  //
  // Tradtionally, the inputs file with ParmParse functionality is handled by C++.  However, a Fortran settings
  // file, by default named probin, can also supply variables.  Mostly used for mesh refinement (tagging) critera
  std::string probin_file("probin");

  ParmParse ppa("amr");
  ppa.query("probin_file",probin_file);

  int probin_file_length = probin_file.length();
  Vector<int> probin_file_name(probin_file_length);

  for (int i = 0; i < probin_file_length; i++)
    probin_file_name[i] = probin_file[i];

  // use a fortran routine to
  // read in tagging parameters from probin file
  get_tagging_params(probin_file_name.dataPtr(), &probin_file_length);

}

//
// AMReX has an inbuilt reflux command, but we still have the freedom
// to decide what goes into it (for example, which variables are
// actually refluxed).  This also gives a little flexibility as to
// where flux registers are stored.  In this example, they are stored
// on levels [1,fine] but not level 0.  
//
void
AmrLevelAdv::reflux ()
{
  BL_ASSERT(level<parent->finestLevel());

  const Real strt = amrex::second();

  // Call the reflux command with the appropriate data.  Because there
  // are no flux registers on the coarse level, they start from the
  // first level.  But the coarse level to the (n-1)^th are the ones
  // that need refluxing, hence the `level+1'.  
  getFluxReg(level+1).Reflux(get_new_data(Phi_Type),1.0,0,0,NUM_STATE,geom);
    
  if (verbose)
  {
    const int IOProc = ParallelDescriptor::IOProcessorNumber();
    Real      end    = amrex::second() - strt;
    
    ParallelDescriptor::ReduceRealMax(end,IOProc);
    
    amrex::Print() << "AmrLevelAdv::reflux() at level " << level 
		   << " : time = " << end << std::endl;
  }
}

//
// Generic function for averaging down - in this case it just makes sure it doesn't happen on the finest level
//
void
AmrLevelAdv::avgDown ()
{
  if (level == parent->finestLevel())
  {
    return;
  }
  // Can select which variables averaging down will happen on - only one to choose from in this case!
  avgDown(Phi_Type);
}

//
// Setting up the call to the AMReX-implemented average down function
//
void
AmrLevelAdv::avgDown (int state_indx)
{
  // For safety, again make sure this only happens if a finer level exists
  if (level == parent->finestLevel()) return;

  // You can access data at other refinement levels, use this to
  // specify your current data, and the finer data that is to be
  // averaged down
  AmrLevelAdv& fine_lev = getLevel(level+1);
  MultiFab&  S_fine   = fine_lev.get_new_data(state_indx);
  MultiFab&  S_crse   = get_new_data(state_indx);

  // Call the AMReX average down function:
  // S_fine: Multifab with the fine data to be averaged down
  // S_crse: Multifab with the coarse data to receive the fine data where necessary
  // fine_lev.geom:  Geometric information (cell size etc.) for the fine level
  // geom: Geometric information for the coarse level (i.e. this level)
  // 0: First variable to be averaged (as not all variables need averaging down
  // S_fine.nComp(): Number of variables to average - this can be computed automatically from a multifab
  // refRatio: The refinement ratio between this level and the finer level
  amrex::average_down(S_fine,S_crse,
		      fine_lev.geom,geom,
		      0,S_fine.nComp(),parent->refRatio(level));
}