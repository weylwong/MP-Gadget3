/*! \file allvars.h
 *  \brief declares the All structure.
 *
 */

#ifndef ALLVARS_H
#define ALLVARS_H

#include <mpi.h>
#include <omp.h>

#include "cosmology.h"
#include "gravity.h"
#include "physconst.h"
#include "types.h"

/*! This structure contains data which is the SAME for all tasks (mostly code parameters read from the
 * parameter file).  Holding this data in a structure is convenient for writing/reading the restart file, and
 * it allows the introduction of new global variables in a simple way. The only thing to do is to introduce
 * them into this structure.
 */
extern struct global_data_all_processes
{
    /* The following variables are set by petaio_read_header */
    int64_t TotNumPartInit; /* The initial total number of particles; we probably want to get rid of all references to this. */
    int64_t NTotalInit[6]; /* The initial number of total particles in the IC. */
    double TimeInit;		/* time of simulation start: if restarting from a snapshot this holds snapshot time.*/
    double TimeIC;       /* Time when the simulation ICs were generated*/
    double BoxSize;   /* Boxsize in case periodic boundary conditions are used */
    double MassTable[6]; /* Initial mass of particles */
    double UnitMass_in_g;		/*!< factor to convert internal mass unit to grams/h */
    double UnitVelocity_in_cm_per_s;	/*!< factor to convert intqernal velocity unit to cm/sec */
    double UnitLength_in_cm;		/*!< factor to convert internal length unit to cm/h */


/* end of read_header parameters */

    struct {
        size_t BytesPerFile;   /* Number of bytes per physical file; this decides how many files bigfile creates each block */
        int WritersPerFile;    /* Number of concurrent writers per file; this decides number of writers */
        int NumWriters;        /* Number of concurrent writers, this caps number of writers */
        int MinNumWriters;        /* Min Number of concurrent writers, this caps number of writers */
        int EnableAggregatedIO;  /* Enable aggregated IO policy for small files.*/
        size_t AggregatedIOThreshold; /* bytes per writer above which to use non-aggregated IO (avoid OOM)*/
        /* Changes the comoving factors of the snapshot outputs. Set in the ICs.
         * If UsePeculiarVelocity = 1 then snapshots save to the velocity field the physical peculiar velocity, v = a dx/dt (where x is comoving distance).
         * If UsePeculiarVelocity = 0 then the velocity field is a * v = a^2 dx/dt in snapshots
         * and v / sqrt(a) = sqrt(a) dx/dt in the ICs. Note that snapshots never match Gadget-2, which
         * saves physical peculiar velocity / sqrt(a) in both ICs and snapshots. */
        int UsePeculiarVelocity;
    } IO;

    double PartAllocFactor;	/*!< in order to maintain work-load balance, the particle load will usually
                              NOT be balanced.  Each processor allocates memory for PartAllocFactor times
                              the average number of particles to allow for that */

    double SlotsIncreaseFactor; /* !< What percentage to increase the slot allocation by when requested*/
    int OutputPotential;        /*!< Flag whether to include the potential in snapshots*/
    int OutputDebugFields;      /* Flag whether to include a lot of debug output in snapshots*/
    int ShowBacktrace;          /* Flag to enable or disable the backtrace printing code*/

    double RandomParticleOffset; /* If > 0, a random shift of max RandomParticleOffset * BoxSize is applied to every particle
                                  * every time a full domain decomposition is done. The box is periodic and the offset
                                  * is subtracted on output, so this only affects the internal gravity solver.
                                  * The purpose of this is to avoid correlated errors in the tree code, which occur when
                                  * the tree opening conditions are similar in every timestep and accumulate over a
                                  * long period of time. Upstream Arepo says this substantially improves momentum conservation,
                                  * and it has the side-effect of guarding against periodicity bugs.
                                  */
    /* Random shift applied to the box. This is changed
     * every domain decomposition to prevent correlated
     * errors building up in the tree force. */
    double CurrentParticleOffset[3];

    /* some SPH parameters */

    double InitGasTemp;		/*!< may be used to set the temperature in the IC's */
    double MinEgySpec; /* Minimum internal energy for timestepping, converted from MinGasTemp*/

    /* system of units  */

    double UnitTime_in_s,		/*!< factor to convert internal time unit to seconds/h */
           UnitDensity_in_cgs,		/*!< factor to convert internal length unit to g/cm^3*h^2 */
           UnitEnergy_in_cgs,		/*!< factor to convert internal energy to cgs units */
           UnitTime_in_Megayears,	/*!< factor to convert internal time to megayears/h */
           G;				/*!< Gravity-constant in internal units */
    /* Cosmology */
    Cosmology CP;

    /* Code options */
    int CoolingOn;  /* if cooling is enabled */
    int HydroOn;  /*  if hydro force is enabled */
    int DensityOn;  /*  if SPH density computation is enabled */
    int TreeGravOn;     /* tree gravity force is enabled*/

    int BlackHoleOn;  /* if black holes are enabled */
    int StarformationOn;  /* if star formation is enabled */
    int WindOn; /* if Wind is enabled */

    int MassiveNuLinRespOn; /*!< flags that massive neutrinos using the linear
                                 response code of Ali-Haimoud & Bird 2013.*/
    int HybridNeutrinosOn; /*!< Flags that hybrid neutrinos are enabled */
    double HybridVcrit; /*!< Critical velocity switching between particle
                          and analytic solvers when hybrid neutrinos are on*/
    double HybridNuPartTime; /*!< Redshift at which hybrid neutrinos switch on*/

    int FastParticleType; /*!< flags a particle species to exclude timestep calculations.*/
    /* parameters determining output frequency */

    int SnapshotFileCount;	/*!< number of snapshot that is written next */
    int InitSnapshotCount;  /*!< Number of first snapshot written this run*/
    double AutoSnapshotTime;    /*!< cpu-time between regularly generated snapshots. */
    double TimeBetweenSeedingSearch; /*Factor to multiply TimeInit by to find the next seeding check.*/

    /* Current time of the simulation, global step, and end of simulation */

    double Time,			/*!< current time of the simulation */
           TimeStep,			/*!< difference between current times of previous and current timestep */
           TimeMax;			/*!< marks the point of time until the simulation is to be evolved */

    struct {
        double a;
        double a3inv;
        double a2inv;
        double hubble;
    } cf;

    /* variables for organizing discrete timeline */

    inttime_t Ti_Current;		/*!< current time on integer timeline */

    int Nmesh;

    /* variables that keep track of cumulative CPU consumption */

    double TimeLimitCPU;

    /*! The scale of the short-range/long-range force split in units of FFT-mesh cells */
    double Asmth;
    enum ShortRangeForceWindowType ShortRangeForceWindowType;	/*!< method of the feedback*/

    double MaxMemSizePerNode;

    double HydroCostFactor; /* cost factor for hydro in load balancing. */

    double GravitySoftening; /* Softening as a fraction of DM mean separation. */
    double GravitySofteningGas;  /* if 0, enable adaptive gravitational softening for gas particles, which uses the Hsml as ForceSoftening */

    double MeanSeparation[6]; /* mean separation between particles. 0 if the species doesn't exist. */

    /* some filenames */
    char InitCondFile[100],
         OutputDir[100],
         SnapshotFileBase[100],
         FOFFileBase[100],
         EnergyFile[100],
         CpuFile[100];

    /*Should we store the energy to EnergyFile on PM timesteps.*/
    int OutputEnergyDebug;

    double OutputListTimes[1024];
    int OutputListLength;

    int SnapshotWithFOF; /*Flag that doing FOF for snapshot outputs is on*/

    int RandomSeed; /*Initial seed for the random number table*/
}
All;

#endif
