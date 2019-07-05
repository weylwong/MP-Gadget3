#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>

#include "allvars.h"
#include "cooling.h"
#include "densitykernel.h"
#include "treewalk.h"
#include "timefac.h"
#include "slotsmanager.h"
#include "timestep.h"
#include "winds.h"
#include "utils.h"

#define MAXITER 400

/* The evolved entropy at drift time: evolved dlog a.
 * Used to predict pressure and entropy for SPH */
static MyFloat
SPH_EntVarPred(int i)
{
        double dloga = dloga_from_dti(P[i].Ti_drift - P[i].Ti_kick);
        double EntVarPred = SPHP(i).Entropy + SPHP(i).DtEntropy * dloga;
        /*Entropy limiter for the predicted entropy: makes sure entropy stays positive. */
        if(dloga > 0 && EntVarPred < 0.5*SPHP(i).Entropy)
            EntVarPred = 0.5 * SPHP(i).Entropy;
        EntVarPred = pow(EntVarPred, 1/GAMMA);
        return EntVarPred;
}

/* Get the predicted velocity for a particle
 * at the Force computation time, which always coincides with the Drift inttime.
 * for gravity and hydro forces.*/
static void
SPH_VelPred(int i, MyFloat * VelPred)
{
    const int ti = P[i].Ti_drift;
    const double Fgravkick2 = get_gravkick_factor(P[i].Ti_kick, ti);
    const double Fhydrokick2 = get_hydrokick_factor(P[i].Ti_kick, ti);
    const double FgravkickB = get_gravkick_factor(PM.Ti_kick, ti);
    int j;
    for(j = 0; j < 3; j++) {
        VelPred[j] = P[i].Vel[j] + Fgravkick2 * P[i].GravAccel[j]
            + P[i].GravPM[j] * FgravkickB + Fhydrokick2 * SPHP(i).HydroAccel[j];
    }
}

/*! Structure for communication during the density computation. Holds data that is sent to other processors.
*/
typedef struct {
    TreeWalkNgbIterBase base;
    DensityKernel kernel;
    double kernel_volume;
} TreeWalkNgbIterDensity;

typedef struct
{
    TreeWalkQueryBase base;
    MyFloat Hsml;
} TreeWalkQueryBHDensity;

typedef struct
{
    TreeWalkQueryBHDensity base;
    double Vel[3];
    MyFloat DelayTime;
} TreeWalkQueryDensity;

typedef struct {
    TreeWalkResultBase base;
    MyFloat Ngb;
    int Ninteractions;
} TreeWalkResultBHDensity;

typedef struct {
    TreeWalkResultBHDensity base;

    /*These are only used for density independent SPH*/
    MyFloat EgyRho;
    MyFloat DhsmlEgyDensity;

    MyFloat Rho;
    MyFloat DhsmlDensity;
    MyFloat Div;
    MyFloat Rot[3];

    /*Only used if sfr_need_to_compute_sph_grad_rho is true*/
    MyFloat GradRho[3];
} TreeWalkResultDensity;

struct BHDensityPriv {
    MyFloat *Left, *Right, *NumNgb;
    int NIteration;
    int *NPLeft;
    double desnumngb;
};

struct DensityPriv {
    struct BHDensityPriv dp;
    MyFloat (*Rot)[3];
    /* This is the DhsmlDensityFactor for the pure density,
     * not the entropy weighted density.
     * If DensityIndependentSphOn = 0 then DhsmlEgyDensityFactor and DhsmlDensityFactor
     * are the same and this is not used.
     * If DensityIndependentSphOn = 1 then this is used to set DhsmlEgyDensityFactor.*/
    MyFloat * DhsmlDensityFactor;
    int update_hsml;
    int DoEgyDensity;
};

#define DENSITY_GET_PRIV(tw) ((struct DensityPriv*) ((tw)->priv))
#define BHDENSITY_GET_PRIV(tw) ((struct BHDensityPriv*) ((tw)->priv))

static void
density_ngbiter(
        TreeWalkQueryDensity * I,
        TreeWalkResultDensity * O,
        TreeWalkNgbIterDensity * iter,
        LocalTreeWalk * lv);

static void
density_bh_ngbiter(
        TreeWalkQueryBHDensity * I,
        TreeWalkResultBHDensity * O,
        TreeWalkNgbIterDensity * iter,
        LocalTreeWalk * lv);

static int density_haswork(int n, TreeWalk * tw);

static void density_postprocess(int i, TreeWalk * tw);
static void density_check_neighbours(int i, TreeWalk * tw);

static void density_reduce_bh(int place, TreeWalkResultBHDensity * remote, enum TreeWalkReduceMode mode, TreeWalk * tw);
static void density_copy_bh(int place, TreeWalkQueryBHDensity * I, TreeWalk * tw);

static void density_reduce_sph(int place, TreeWalkResultDensity * remote, enum TreeWalkReduceMode mode, TreeWalk * tw);
static void density_copy_sph(int place, TreeWalkQueryDensity * I, TreeWalk * tw);

/*! \file density.c
 *  \brief SPH density computation and smoothing length determination
 *
 *  This file contains the "first SPH loop", where the SPH densities and some
 *  auxiliary quantities are computed.  There is also functionality that
 *  corrects the smoothing length if needed.
 */

/* Do initialisation common to BH and gas*/
static void
density_init_bhdensity(TreeWalk * tw, int size)
{
    int i;
    BHDENSITY_GET_PRIV(tw)->Left = (MyFloat *) mymalloc("DENSITY->Left", size * sizeof(MyFloat));
    memset(BHDENSITY_GET_PRIV(tw)->Left, 0, size * sizeof(MyFloat));
    BHDENSITY_GET_PRIV(tw)->Right = (MyFloat *) mymalloc("DENSITY->Right", size * sizeof(MyFloat));
    BHDENSITY_GET_PRIV(tw)->NumNgb = (MyFloat *) mymalloc("DENSITY->NumNgb", size * sizeof(MyFloat));
    BHDENSITY_GET_PRIV(tw)->NIteration = 0;

    #pragma omp parallel for
    for(i = 0; i < size; i++)
    {
        BHDENSITY_GET_PRIV(tw)->Right[i] = All.BoxSize;
    }
}

struct density_queues
{
    int * SPH_Queue;
    int nsphqueue;
    int * BH_Queue;
    int nbhqueue;
};

/* This builds a queue for the active SPH particles and another for the active BH particles,
 * saving doing this twice.*/
static struct density_queues
density_init_queues(int * active_set, const int size)
{
    int i;
    int NumThreads = omp_get_max_threads();
    /* Since we use a static schedule below we only need size / tw->NThread elements per thread.
     * Add 2 for non-integer parts.*/
    int sph_tsize = SlotsManager->info[0].size / NumThreads + 2;
    int bh_tsize = SlotsManager->info[5].size / NumThreads + 2;

    struct density_queues denque;
    denque.BH_Queue = (int *) mymalloc2("BH_Queue", bh_tsize * sizeof(int) * NumThreads);
    denque.SPH_Queue = (int *) mymalloc("SPH_Queue", sph_tsize * sizeof(int) * NumThreads);

    /*We want a lockless algorithm which preserves the ordering of the particle list.*/
    size_t *nqthr_sph = ta_malloc("nqthr_s", size_t, NumThreads);
    int **thrqueue_sph = ta_malloc("thrqueue_s", int *, NumThreads);
    /*We want a lockless algorithm which preserves the ordering of the particle list.*/
    size_t *nqthr_bh = ta_malloc("nqthr_b", size_t, NumThreads);
    int **thrqueue_bh = ta_malloc("thrqueue_b", int *, NumThreads);

    gadget_setup_thread_arrays(denque.SPH_Queue, thrqueue_sph, nqthr_sph, sph_tsize, NumThreads);
    gadget_setup_thread_arrays(denque.BH_Queue, thrqueue_bh, nqthr_bh, bh_tsize, NumThreads);

    /* We enforce schedule static to ensure that each thread executes on contiguous particles.*/
    #pragma omp parallel for schedule(static)
    for(i=0; i < size; i++)
    {
        const int tid = omp_get_thread_num();
        /*Use raw particle number if active_set is null, otherwise use active_set*/
        const int p_i = active_set ? active_set[i] : i;

        /* Skip the garbage particles */
        if(P[p_i].IsGarbage) continue;

        /* this has to be done before treewalk so that
         * all particles are ran for the first loop.
         * The iteration will gradually turn DensityIterationDone on more particles.
         * */
        P[p_i].DensityIterationDone = 0;

        if(P[p_i].Type == 0) {
            thrqueue_sph[tid][nqthr_sph[tid]] = p_i;
            nqthr_sph[tid]++;
        }
        if(P[p_i].Type == 5) {
            thrqueue_bh[tid][nqthr_bh[tid]] = p_i;
            nqthr_bh[tid]++;
        }
    }
    /*Merge step for the queue.*/
    denque.nbhqueue = gadget_compact_thread_arrays(denque.BH_Queue, thrqueue_bh, nqthr_bh, NumThreads);
    ta_free(thrqueue_bh);
    ta_free(nqthr_bh);
    /*Shrink memory*/
    denque.BH_Queue = myrealloc(denque.BH_Queue, sizeof(int) * denque.nbhqueue);
    denque.nsphqueue = gadget_compact_thread_arrays(denque.SPH_Queue, thrqueue_sph, nqthr_sph, NumThreads);
    ta_free(thrqueue_sph);
    ta_free(nqthr_sph);
    /*Shrink memory*/
    denque.SPH_Queue = myrealloc(denque.SPH_Queue, sizeof(int) * denque.nsphqueue);
    return denque;
}

static void
density_do_iterations(TreeWalk * tw, int * active_set, int size)
{
    int64_t ntot = 0;
    int i;
    BHDENSITY_GET_PRIV(tw)->NPLeft = ta_malloc("NPLeft", int, All.NumThreads);

    /* we will repeat the whole thing for those particles where we didn't find enough neighbours */
    do {
        memset(BHDENSITY_GET_PRIV(tw)->NPLeft, 0, sizeof(int)*All.NumThreads);

        treewalk_run(tw, active_set, size);

        /* Set the haswork function to checking DensityInterationDone*/
        tw->haswork = density_haswork;
        int Nleft = 0;

        for(i = 0; i< All.NumThreads; i++)
            Nleft += BHDENSITY_GET_PRIV(tw)->NPLeft[i];

        sumup_large_ints(1, &Nleft, &ntot);

        if(ntot == 0) break;

        BHDENSITY_GET_PRIV(tw)->NIteration ++;
        /*
        if(ntot < 1 ) {
            foreach(ActiveParticle)
            {
                if(density_haswork(i) && !P[i].DensityIterationDone) {
                    MyFloat Left = DENSITY_GET_PRIV(tw)->Left[PI];
                    MyFloat Right = DENSITY_GET_PRIV(tw)->Right[PI];
                    message (1, "i=%d task=%d ID=%llu type=%d, Hsml=%g Left=%g Right=%g Ngbs=%g Right-Left=%g\n   pos=(%g|%g|%g)\n",
                         i, ThisTask, P[i].ID, P[i].Type, P[i].Hsml, Left, Right,
                         (float) P[i].NumNgb, Right - Left, P[i].Pos[0], P[i].Pos[1], P[i].Pos[2]);
                }
            }

        }
        */

        if(BHDENSITY_GET_PRIV(tw)->NIteration > 0) {
            message(0, "ngb iteration %d: need to repeat for %ld particles.\n", BHDENSITY_GET_PRIV(tw)->NIteration, ntot);
        }

        if(BHDENSITY_GET_PRIV(tw)->NIteration > MAXITER) {
            endrun(1155, "failed to converge in neighbour iteration in density()\n");
        }
    } while(1);

    ta_free(BHDENSITY_GET_PRIV(tw)->NPLeft);
}

/*! This function computes the local density for each active SPH particle, the
 * number of neighbours in the current smoothing radius, and the divergence
 * and rotation of the velocity field.  The pressure is updated as well.  If a
 * particle with its smoothing region is fully inside the local domain, it is
 * not exported to the other processors. The function also detects particles
 * that have a number of neighbours outside the allowed tolerance range. For
 * these particles, the smoothing length is adjusted accordingly, and the
 * density() computation is called again.  Note that the smoothing length is
 * not allowed to fall below the lower bound set by MinGasHsml (this may mean
 * that one has to deal with substantially more than normal number of
 * neighbours.)
 */
void
density(int update_hsml, int DoEgyDensity, ForceTree * tree)
{
    if(!All.DensityOn)
        return;

    TreeWalk tw[1] = {{0}};
    struct DensityPriv priv[1];

    tw->ev_label = "DENSITY";
    tw->visit = (TreeWalkVisitFunction) treewalk_visit_ngbiter;
    tw->ngbiter_type_elsize = sizeof(TreeWalkNgbIterDensity);
    tw->ngbiter = (TreeWalkNgbIterFunction) density_ngbiter;
    tw->haswork = NULL;
    tw->fill = (TreeWalkFillQueryFunction) density_copy_sph;
    tw->reduce = (TreeWalkReduceResultFunction) density_reduce_sph;
    tw->postprocess = (TreeWalkProcessFunction) density_postprocess;
    tw->UseNodeList = 1;
    tw->query_type_elsize = sizeof(TreeWalkQueryDensity);
    tw->result_type_elsize = sizeof(TreeWalkResultDensity);
    tw->priv = priv;
    tw->tree = tree;

    int i;
    double timeall = 0;
    double timecomp, timecomm, timewait;

    walltime_measure("/Misc");

    struct density_queues denque = density_init_queues(ActiveParticle, NumActiveParticle);

    BHDENSITY_GET_PRIV(tw)->desnumngb = All.DesNumNgb;

    density_init_bhdensity(tw, SlotsManager->info[0].size);

    DENSITY_GET_PRIV(tw)->Rot = (MyFloat (*) [3]) mymalloc("DENSITY_GET_PRIV(tw)->Rot", SlotsManager->info[0].size * sizeof(priv->Rot[0]));
    if(DoEgyDensity)
        DENSITY_GET_PRIV(tw)->DhsmlDensityFactor = (MyFloat *) mymalloc("DENSITY_GET_PRIV(tw)->DhsmlDensity", SlotsManager->info[0].size * sizeof(MyFloat));
    else
        DENSITY_GET_PRIV(tw)->DhsmlDensityFactor = NULL;

    DENSITY_GET_PRIV(tw)->update_hsml = update_hsml;
    DENSITY_GET_PRIV(tw)->DoEgyDensity = DoEgyDensity;

    #pragma omp parallel for
    for(i = 0; i < PartManager->NumPart; i++)
    {
        if(P[i].Type == 0) {
            const int PI = P[i].PI;
            SphP_scratch->EntVarPred[PI] = SPH_EntVarPred(i);
            SPH_VelPred(i, SphP_scratch->VelPred + 3 * PI);
        }
    }

    walltime_measure("/SPH/Density/Init");

    density_do_iterations(tw, denque.SPH_Queue, denque.nsphqueue);

    if(DoEgyDensity)
        myfree(DENSITY_GET_PRIV(tw)->DhsmlDensityFactor);
    myfree(DENSITY_GET_PRIV(tw)->Rot);
    myfree(BHDENSITY_GET_PRIV(tw)->NumNgb);
    myfree(BHDENSITY_GET_PRIV(tw)->Right);
    myfree(BHDENSITY_GET_PRIV(tw)->Left);

    myfree(denque.SPH_Queue);
    timeall = walltime_measure(WALLTIME_IGNORE);

    /* Now do a treewalk to work out Hsml for the black holes*/
    if(All.BlackHoleOn && update_hsml && SlotsManager->info[5].size > 0) {
        struct BHDensityPriv bhpriv[1];

        tw->ev_label = "BHDENSITY";
        tw->haswork = NULL;
        tw->ngbiter = (TreeWalkNgbIterFunction) density_bh_ngbiter;
        tw->fill = (TreeWalkFillQueryFunction) density_copy_bh;
        tw->reduce = (TreeWalkReduceResultFunction) density_reduce_bh;
        tw->postprocess = (TreeWalkProcessFunction) density_check_neighbours;
        tw->query_type_elsize = sizeof(TreeWalkQueryBHDensity);
        tw->result_type_elsize = sizeof(TreeWalkResultBHDensity);
        tw->priv = bhpriv;

        BHDENSITY_GET_PRIV(tw)->desnumngb = All.DesNumNgb * All.BlackHoleNgbFactor;

        density_init_bhdensity(tw, SlotsManager->info[5].size);

        walltime_measure("/SPH/Density/Init");

        density_do_iterations(tw, denque.BH_Queue, denque.nbhqueue);

        myfree(BHDENSITY_GET_PRIV(tw)->NumNgb);
        myfree(BHDENSITY_GET_PRIV(tw)->Right);
        myfree(BHDENSITY_GET_PRIV(tw)->Left);
    }

    /* collect some timing information */

    myfree(denque.BH_Queue);

    timeall += walltime_measure(WALLTIME_IGNORE);

    timecomp = tw->timecomp3 + tw->timecomp1 + tw->timecomp2;
    timewait = tw->timewait1 + tw->timewait2;
    timecomm = tw->timecommsumm1 + tw->timecommsumm2;

    walltime_add("/SPH/Density/Compute", timecomp);
    walltime_add("/SPH/Density/Wait", timewait);
    walltime_add("/SPH/Density/Comm", timecomm);
    walltime_add("/SPH/Density/Misc", timeall - (timecomp + timewait + timecomm));
}

static void
density_copy_bh(int place, TreeWalkQueryBHDensity * I, TreeWalk * tw)
{
    I->Hsml = P[place].Hsml;
}

static void
density_copy_sph(int place, TreeWalkQueryDensity * I, TreeWalk * tw)
{
    density_copy_bh(place, &I->base, tw);
    I->Vel[0] = SphP_scratch->VelPred[3 * P[place].PI];
    I->Vel[1] = SphP_scratch->VelPred[3 * P[place].PI + 1];
    I->Vel[2] = SphP_scratch->VelPred[3 * P[place].PI + 2];
    I->DelayTime = SPHP(place).DelayTime;
}

static void
density_reduce_bh(int place, TreeWalkResultBHDensity * remote, enum TreeWalkReduceMode mode, TreeWalk * tw)
{
    int pi = P[place].PI;
    TREEWALK_REDUCE(BHDENSITY_GET_PRIV(tw)->NumNgb[pi], remote->Ngb);
    /* these will be added */
    P[place].GravCost += All.HydroCostFactor * All.cf.a * remote->Ninteractions;
}

static void
density_reduce_sph(int place, TreeWalkResultDensity * remote, enum TreeWalkReduceMode mode, TreeWalk * tw)
{
    density_reduce_bh(place, &remote->base, mode, tw);

    TREEWALK_REDUCE(SPHP(place).Density, remote->Rho);
    TREEWALK_REDUCE(SPHP(place).DivVel, remote->Div);
    int pi = P[place].PI;
    TREEWALK_REDUCE(DENSITY_GET_PRIV(tw)->Rot[pi][0], remote->Rot[0]);
    TREEWALK_REDUCE(DENSITY_GET_PRIV(tw)->Rot[pi][1], remote->Rot[1]);
    TREEWALK_REDUCE(DENSITY_GET_PRIV(tw)->Rot[pi][2], remote->Rot[2]);

    if(SphP_scratch->GradRho) {
        TREEWALK_REDUCE(SphP_scratch->GradRho[3*pi], remote->GradRho[0]);
        TREEWALK_REDUCE(SphP_scratch->GradRho[3*pi+1], remote->GradRho[1]);
        TREEWALK_REDUCE(SphP_scratch->GradRho[3*pi+2], remote->GradRho[2]);
    }
    /*Only used for density independent SPH*/
    if(DENSITY_GET_PRIV(tw)->DoEgyDensity) {
        TREEWALK_REDUCE(SPHP(place).EgyWtDensity, remote->EgyRho);
        TREEWALK_REDUCE(SPHP(place).DhsmlEgyDensityFactor, remote->DhsmlEgyDensity);
        TREEWALK_REDUCE(DENSITY_GET_PRIV(tw)->DhsmlDensityFactor[pi], remote->DhsmlDensity);
    }
    else
        TREEWALK_REDUCE(SPHP(place).DhsmlEgyDensityFactor, remote->DhsmlDensity);
}

/******
 *
 *  This function represents the core of the SPH density computation.
 *
 *  The neighbours of the particle in the Query are enumerated, and results
 *  are stored into the Result object.
 *
 *  Upon start-up we initialize the iterator with the density kernels used in
 *  the computation. The assumption is the density kernels are slow to
 *  initialize.
 *
 */

static void
density_bh_ngbiter(
        TreeWalkQueryBHDensity * I,
        TreeWalkResultBHDensity * O,
        TreeWalkNgbIterDensity * iter,
        LocalTreeWalk * lv)
{
    if(iter->base.other == -1) {
        const double h = I->Hsml;
        density_kernel_init(&iter->kernel, h);
        iter->kernel_volume = density_kernel_volume(&iter->kernel);

        iter->base.Hsml = h;
        iter->base.mask = 1; /* gas only */
        iter->base.symmetric = NGB_TREEFIND_ASYMMETRIC;
        return;
    }
    const double r = iter->base.r;
    const double r2 = iter->base.r2;

    if(r2 < iter->kernel.HH)
    {
        const double u = r * iter->kernel.Hinv;
        const double wk = density_kernel_wk(&iter->kernel, u);
        O->Ngb += wk * iter->kernel_volume;
    }

    /* some performance measures not currently used */
    O->Ninteractions ++;
}


/******
 *
 *  This function represents the core of the SPH density computation.
 *
 *  The neighbours of the particle in the Query are enumerated, and results
 *  are stored into the Result object.
 *
 *  Upon start-up we initialize the iterator with the density kernels used in
 *  the computation. The assumption is the density kernels are slow to
 *  initialize.
 *
 */

static void
density_ngbiter(
        TreeWalkQueryDensity * I,
        TreeWalkResultDensity * O,
        TreeWalkNgbIterDensity * iter,
        LocalTreeWalk * lv)
{
    density_bh_ngbiter(&I->base, &O->base, iter, lv);
    if(iter->base.other == -1)
            return;

    const int other = iter->base.other;
    const double r = iter->base.r;
    const double r2 = iter->base.r2;
    const double * dist = iter->base.dist;

    if(All.WindOn) {
        if(winds_is_particle_decoupled(other))
            if(!(I->DelayTime > 0))	/* if I'm not wind, then ignore the wind particle */
                return;
    }

    if(P[other].Mass == 0) {
        endrun(12, "Encountered zero mass particle during density;"
                  " We haven't implemented tracer particles and this shall not happen\n");
    }

    if(r2 < iter->kernel.HH)
    {
        const double u = r * iter->kernel.Hinv;
        const double wk = density_kernel_wk(&iter->kernel, u);

        const double dwk = density_kernel_dwk(&iter->kernel, u);

        const double mass_j = P[other].Mass;

        O->Rho += (mass_j * wk);

        /* Hinv is here because O->DhsmlDensity is drho / dH.
         * nothing to worry here */
        double density_dW = density_kernel_dW(&iter->kernel, u, wk, dwk);
        O->DhsmlDensity += mass_j * density_dW;

        if(DENSITY_GET_PRIV(lv->tw)->DoEgyDensity) {
            const double EntPred = SphP_scratch->EntVarPred[P[other].PI];
            O->EgyRho += mass_j * EntPred * wk;
            O->DhsmlEgyDensity += mass_j * EntPred * density_dW;
        }


        if(SphP_scratch->GradRho) {
            if(r > 0)
            {
                int d;
                for (d = 0; d < 3; d ++) {
                    O->GradRho[d] += mass_j * dwk * dist[d] / r;
                }
            }
        }

        if(r > 0)
        {
            double fac = mass_j * dwk / r;
            double dv[3];
            double rot[3];
            int d;
            for(d = 0; d < 3; d ++) {
                dv[d] = I->Vel[d] - SphP_scratch->VelPred[3 * P[other].PI + d];
            }
            O->Div += -fac * dotproduct(dist, dv);

            crossproduct(dv, dist, rot);
            for(d = 0; d < 3; d ++) {
                O->Rot[d] += fac * rot[d];
            }
        }
    }
}

static int
density_haswork(int n, TreeWalk * tw)
{
    if(P[n].DensityIterationDone) return 0;

    return 1;
}

static void density_check_neighbours_int (int i, TreeWalk * tw, MyFloat DensFac);

static void
density_postprocess(int i, TreeWalk * tw)
{
    int PI = P[i].PI;
    MyFloat * DhsmlDens;
    if(DENSITY_GET_PRIV(tw)->DoEgyDensity)
        DhsmlDens = &(DENSITY_GET_PRIV(tw)->DhsmlDensityFactor[PI]);
    else
        DhsmlDens = &(SPHP(i).DhsmlEgyDensityFactor);
    if(SPHP(i).Density <= 0) {
        if(BHDENSITY_GET_PRIV(tw)->NumNgb[PI] == 0) {
            SPHP(i).Density = 1;
            *DhsmlDens = 1;
        } else
            endrun(12, "Particle %d has bad density: %g\n", i, SPHP(i).Density);
    }

    *DhsmlDens *= P[i].Hsml / (NUMDIMS * SPHP(i).Density);
    *DhsmlDens = 1 / (1 + *DhsmlDens);

    /*Compute the EgyWeight factors, which are only useful for density independent SPH */
    if(DENSITY_GET_PRIV(tw)->DoEgyDensity) {
        const double EntPred = SphP_scratch->EntVarPred[P[i].PI];
        if(EntPred <= 0 || SPHP(i).EgyWtDensity <=0)
            endrun(12, "Particle %d has bad predicted entropy: %g or EgyWtDensity: %g\n", i, EntPred, SPHP(i).EgyWtDensity);
        SPHP(i).DhsmlEgyDensityFactor *= P[i].Hsml/ (NUMDIMS * SPHP(i).EgyWtDensity);
        SPHP(i).DhsmlEgyDensityFactor *= - (*DhsmlDens);
        SPHP(i).EgyWtDensity /= EntPred;
    }

    MyFloat * Rot = DENSITY_GET_PRIV(tw)->Rot[PI];
    SPHP(i).CurlVel = sqrt(Rot[0] * Rot[0] + Rot[1] * Rot[1] + Rot[2] * Rot[2]) / SPHP(i).Density;

    SPHP(i).DivVel /= SPHP(i).Density;

    /* This is slightly more complicated so we put it in a different function */
    if(DENSITY_GET_PRIV(tw)->update_hsml)
        density_check_neighbours_int(i, tw, *DhsmlDens);
}

static void
density_check_neighbours(int i, TreeWalk * tw)
{
    density_check_neighbours_int(i, tw, -1);
}

static void
density_check_neighbours_int (int i, TreeWalk * tw, MyFloat DensFac)
{
    /* now check whether we had enough neighbours */

    const double desnumngb = BHDENSITY_GET_PRIV(tw)->desnumngb;

    MyFloat * Left = BHDENSITY_GET_PRIV(tw)->Left;
    MyFloat * Right = BHDENSITY_GET_PRIV(tw)->Right;
    MyFloat * NumNgb = BHDENSITY_GET_PRIV(tw)->NumNgb;
    int PI = P[i].PI;

    if(NumNgb[PI] < (desnumngb - All.MaxNumNgbDeviation) ||
            (NumNgb[PI] > (desnumngb + All.MaxNumNgbDeviation)))
    {
        /* need to redo this particle */
        if(P[i].DensityIterationDone) {
            /* should have been 0*/
            endrun(999993, "Already has DensityIterationDone set, bad memory intialization.");
        }

        /* This condition is here to prevent the density code looping forever if it encounters
         * multiple particles at the same position. If this happens you likely have worse
         * problems anyway, so warn also. */
        if((Right[PI] - Left[PI]) < 1.0e-6 * Right[PI])
        {
            /* If this happens probably the exchange is screwed up and all your particles have moved to (0,0,0)*/
            message(1, "Very tight Hsml bounds for i=%d ID=%lu Hsml=%g Left=%g Right=%g Ngbs=%g Right-Left=%g pos=(%g|%g|%g)\n",
            i, P[i].ID, P[i].Hsml, Left[PI], Right[PI], NumNgb[PI], Right[PI] - Left[PI], P[i].Pos[0], P[i].Pos[1], P[i].Pos[2]);
            P[i].Hsml = Right[PI];
            P[i].DensityIterationDone = 1;
            return;
        }

        /* If we need more neighbours, move the lower bound up. If we need fewer, move the upper bound down.*/
        if(NumNgb[PI] < desnumngb) {
                Left[PI] = P[i].Hsml;
        } else {
                Right[PI] = P[i].Hsml;
        }

        /* Next step is geometric mean of previous. */
        if(Right[PI] < 0.99 * All.BoxSize && Left[PI] > 0)
            P[i].Hsml = pow(0.5 * (pow(Left[PI], 3) + pow(Right[PI], 3)), 1.0 / 3);
        else
        {
            if(Right[PI] > 0.99 * All.BoxSize && Left[PI] <= 0)
                endrun(8188, "Cannot occur. Check for memory corruption: L = %g R = %g N=%g.", Left[PI], Right[PI], NumNgb[PI]);

            double fac = 1.26;
            /* If this is the first step we can be faster by increasing or decreasing current Hsml by a constant factor*/
            if(Right[PI] > 0.99 * All.BoxSize && Left[PI] > 0)
                fac = 1.26;
            if(Right[PI] < 0.99*All.BoxSize && Left[PI] == 0)
                fac = 1/1.26;

            /* Check whether this actually helps. If it does, why not for BH as well? DensFac ~ 1.*/
            if(DensFac > 0 && fabs(NumNgb[PI] - desnumngb) < 0.5 * desnumngb) {
                fac = 1 - (NumNgb[PI] - desnumngb) / (NUMDIMS * NumNgb[PI]) * DensFac;
                if(fac > 1.26)
                    fac = 1.26;
                if(fac < 1/1.26)
                    fac = 1/1.26;
            }
            P[i].Hsml *= fac;
        }

        if(Right[PI] < All.MinGasHsml) {
            P[i].Hsml = All.MinGasHsml;
            P[i].DensityIterationDone = 1;
        }
    }
    else {
        /* We might have got here by serendipity, without bounding Right.*/
        if(P[i].Hsml < All.MinGasHsml)
            P[i].Hsml = All.MinGasHsml;
        P[i].DensityIterationDone = 1;
    }

    if(BHDENSITY_GET_PRIV(tw)->NIteration >= MAXITER - 10)
    {
         message(1, "i=%d ID=%lu Hsml=%g Left=%g Right=%g Ngbs=%g Right-Left=%g\n   pos=(%g|%g|%g)\n",
             i, P[i].ID, P[i].Hsml, Left[PI], Right[PI],
             NumNgb[PI], Right[PI] - Left[PI], P[i].Pos[0], P[i].Pos[1], P[i].Pos[2]);
    }

    if(!P[i].DensityIterationDone) {
        int tid = omp_get_thread_num();
        BHDENSITY_GET_PRIV(tw)->NPLeft[tid] ++;
    }
}

