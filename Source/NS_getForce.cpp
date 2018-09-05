
#include <NavierStokesBase.H>
#include <AMReX_BLFort.H>
#include <PROB_NS_F.H>

using namespace amrex;

//
// Virtual access function for getting the forcing terms for the
// velocities and scalars.  The base version computes a buoyancy.
//
// As NavierStokesBase is currently implemented.  Velocities are integrated
// according to the equation
//
//     ui_t + uj ui_j = S_ui        ===> tforces = rho S_ui
//
// and scalars psi where (psi = rho q) as
//
//     psi_t + (uj psi)_j = S_psi   ===> tforces = S_psi = rho S_q
//
// q is a concentration.  This function returns a rho weighted
// source term, which requires a division by rho in the predict_velocity
// and velocity_advection routines.
//

#ifdef BOUSSINESQ
void
NavierStokesBase::getForce (FArrayBox&       force,
			    int              gridno,
			    int              ngrow,
			    int              scomp,
			    int              ncomp,
			    const Real       time,
			    const FArrayBox& Scal)
{
    force.resize(amrex::grow(grids[gridno],ngrow),ncomp);

    if (scomp == Xvel && ncomp == BL_SPACEDIM)
    {
       if (Scal.nComp() > 1) 
       {
	 amrex::Print() << "OOPS -- ONLY SUPPOSED TO BE ONE COMPONENT IN SCALAR " << std::endl;
	 exit(0);
       }

       const Real* dx       = geom.CellSize();
       const int*  f_lo     = force.loVect();
       const int*  f_hi     = force.hiVect();
       const int*  s_lo     = Scal.loVect();
       const int*  s_hi     = Scal.hiVect();
   
       RealBox gridloc = RealBox(grids[gridno],geom.CellSize(),geom.ProbLo());

       const Real* ScalDataPtr = Scal.dataPtr(0);
       FORT_MAKEFORCE (force.dataPtr(), ScalDataPtr,
   		       ARLIM(f_lo), ARLIM(f_hi),
   		       ARLIM(s_lo), ARLIM(s_hi),
   		       dx, gridloc.lo(), gridloc.hi(),
   		       &scomp,&ncomp);
    } else {
       force.setVal(0.0);
    }
}
#else
void
NavierStokesBase::getForce (FArrayBox&       force,
			    int              gridno,
			    int              ngrow,
			    int              scomp,
			    int              ncomp,
			    const FArrayBox& Rho,
			    int              RComp)
{
    BL_ASSERT(Rho.nComp() > RComp);

    force.resize(amrex::grow(grids[gridno],ngrow),ncomp);

    BL_ASSERT(Rho.box().contains(force.box()));

    const Real grav = gravity;

    for (int dc = 0; dc < ncomp; dc++)
    {
        const int sc = scomp + dc;
#if (BL_SPACEDIM == 2)
        if (sc == Yvel && std::fabs(grav) > 0.001) 
#endif
#if (BL_SPACEDIM == 3)
        if (sc == Zvel && std::fabs(grav) > 0.001) 
#endif
        {
            //
            // Set force to -rho*g.
            //
     	    force.copy(Rho,RComp,dc,1);
            force.mult(grav,dc,1);
        }
        else
        {
            force.setVal(0,dc);
        }
    }    
}
#endif
