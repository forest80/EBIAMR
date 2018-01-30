
#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Amr.H>
#include <AMReX_EBTower.H>

#include <NavierStokesBase.H>

using namespace amrex;

void initialize_EBIS(const int max_level);

int
main (int   argc,
      char* argv[])
{
    amrex::Initialize(argc,argv);

    BL_PROFILE_REGION_START("main()");
    BL_PROFILE_VAR("main()", pmain);

    Real timer_tot = ParallelDescriptor::second();
    Real timer_init = 0.;
    Real timer_advance = 0.;

    int  max_step;
    int  num_steps;
    Real strt_time;
    Real stop_time;

    ParmParse pp;

    max_step  = -1; 
    num_steps = -1; 
    strt_time =  0.0;
    stop_time = -1.0;

    pp.query("max_step",  max_step);
    pp.query("num_steps", num_steps);
    pp.query("strt_time", strt_time);
    pp.query("stop_time", stop_time);

    if (strt_time < 0.0)
    {
        amrex::Abort("MUST SPECIFY a non-negative strt_time");
    }

    if (max_step < 0 && stop_time < 0)
    {
        amrex::Abort("Exiting because neither max_step nor stop_time is non-negative.");
    }

    {
        timer_init = ParallelDescriptor::second();

        Amr amr;

        initialize_EBIS(amr.maxLevel());
        EBTower::Build();
        AMReX_EBIS::reset();  // We no longer needs the EBIndexSpace singleton.
        AmrLevel::SetEBSupportLevel(EBSupport::full);
        AmrLevel::SetEBMaxGrowCells(5,4,2);

        amr.init(strt_time,stop_time);

        if (num_steps > 0)
        {
            if (max_step < 0)
            {
                max_step = num_steps + amr.levelSteps(0);
            }
            else
            {
                max_step = std::min(max_step, num_steps + amr.levelSteps(0));
            }

            amrex::Print() << "Using effective max_step = " << max_step << '\n';
        }
        //
        // If we set the regrid_on_restart flag and if we are *not* going to take
        // a time step then we want to go ahead and regrid here.
        //
        if (amr.RegridOnRestart())
        {
            if (    (amr.levelSteps(0) >= max_step ) ||
                    ( (stop_time >= 0.0) &&
                      (amr.cumTime() >= stop_time)  )    )
            {
                //
                // Regrid only!
                //
                amr.RegridOnly(amr.cumTime());
            }
        }

        timer_init = ParallelDescriptor::second() - timer_init;
        timer_advance = ParallelDescriptor::second();

        while ( amr.okToContinue()                            &&
                (amr.levelSteps(0) < max_step || max_step < 0) &&
                (amr.cumTime() < stop_time || stop_time < 0.0) )
        {
            amr.coarseTimeStep(stop_time);
        }

        timer_advance = ParallelDescriptor::second() - timer_advance;

        if (amr.stepOfLastCheckPoint() < amr.levelSteps(0))
        {
            amr.checkPoint();
        }

        if (amr.stepOfLastPlotFile() < amr.levelSteps(0))
        {
            amr.writePlotFile();
        }

        EBTower::Destroy();
    }

    timer_tot = ParallelDescriptor::second() - timer_tot;

    ParallelDescriptor::ReduceRealMax({timer_tot, timer_init, timer_advance},
                                      ParallelDescriptor::IOProcessorNumber());

    amrex::Print() << "Run Time total        = " << timer_tot     << "\n"
                   << "Run Time init         = " << timer_init    << "\n"
                   << "Run Time advance      = " << timer_advance << "\n";

    BL_PROFILE_VAR_STOP(pmain);
    BL_PROFILE_REGION_STOP("main()");
    BL_PROFILE_SET_RUN_TIME(run_stop);
    BL_PROFILE_FINALIZE();


    amrex::Finalize();

    return 0;
}
