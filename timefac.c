#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>

#include "allvars.h"
#include "timefac.h"
#include "timebinmgr.h"
#include "cosmology.h"
#include "endrun.h"

#define WORKSIZE 10000

/* Simple single-value cache for the drift table,
 * which is the same for all particles*/
static inttime_t df_last_ti0 = -1, df_last_ti1 = -1;
static double df_last_value;

/*inttimes for the kick cache*/
static inttime_t kk_last_ti0[TIMEBINS] = {-1}, kk_last_ti1[TIMEBINS] = {-1};
static inttime_t kk_pred_ti;

/*Values for the kick caches*/
static double hk_last_value[TIMEBINS];
static double gk_last_value[TIMEBINS];
/*This is to store the last predicted gravkick value for this PM step*/
static double gk_last_PM;

/* Integrand for the drift table*/
static double drift_integ(double a, void *param)
{
  double h;
  h = hubble_function(a);
  return 1 / (h * a * a * a);
}

/* Integrand for the gravkick table*/
static double gravkick_integ(double a, void *param)
{
  double h;

  h = hubble_function(a);

  return 1 / (h * a * a);
}

/* Integrand for the hydrokick table.
 * Note this is the same function as drift.*/
static double hydrokick_integ(double a, void *param)
{
  double h;

  h = hubble_function(a);

  return 1 / (h * pow(a, 3 * GAMMA_MINUS1) * a);
}

/*Do the integral required to get a factor.*/
static double
get_exact_factor(inttime_t t0, inttime_t t1, double (*factor) (double, void *))
{
    double result, abserr;
    if(t0 == t1)
        return 0;
    double a0 = exp(loga_from_ti(t0));
    double a1 = exp(loga_from_ti(t1));
    gsl_function F;
    gsl_integration_workspace *workspace;
    workspace = gsl_integration_workspace_alloc(WORKSIZE);
    F.function = factor;
    gsl_integration_qag(&F, a0, a1, 0, 1.0e-8, WORKSIZE, GSL_INTEG_GAUSS61, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);
    return result;
}

/*Update the computed drift factor for a new timestep*/
void
update_drift_factor(inttime_t ti0, inttime_t ti1, int reset)
{
    /* Rebuild the drift cache*/
    if(!reset && df_last_ti1 > 0 && df_last_ti1 != ti0)
        endrun(2,"Drift should start where last step finished! Before had: %d -> %d, now %d->%d\n",
            df_last_ti0, df_last_ti1, ti0, ti1);
    df_last_ti0 = ti0;
    df_last_ti1 = ti1;
    df_last_value = get_exact_factor(ti0, ti1, &drift_integ);
}

/* Update the computed kick factors for a new timestep.
 * We will need two sets of factor: one for the kick itself
 * and one for the velocity prediction (which has ti_drift as the end time).
 * We could one day alter this to only do the integral for active timebins.
 */
void
update_kick_factors(inttime_t Drift_Ti, int direction, int reset)
{
    message(1,"update_kick_factor: %d %d\n", Drift_Ti);

    /* Rebuild the kick cache for each timebin.*/
    int bin;
    for(bin = 0; bin < TIMEBINS; bin++)
    {
        inttime_t ti0 =  Drift_Ti - direction * dti_from_timebin(bin)/2;
        inttime_t ti1 = Drift_Ti + dti_from_timebin(bin)/2;
        if(!reset && kk_last_ti1[bin] > 0 && kk_last_ti1[bin] != ti0)
            endrun(2,"Kick should start where last step finished! Before had: %d -> %d, now %d->%d\n",
                kk_last_ti0, kk_last_ti1, ti0, ti1);
        kk_last_ti0[bin] = ti0;
        kk_last_ti1[bin] = ti1;
        kk_pred_ti = Drift_Ti;

        gk_last_value[bin] = get_exact_factor(ti0, ti1, &gravkick_integ);
        hk_last_value[bin] = get_exact_factor(Drift_Ti, ti1, &gravkick_integ);
        gk_last_PM = get_exact_factor(PM_Ti, Drift_Ti, &gravkick_integ);

        hk_todrift[bin] = get_exact_factor(ti0, Drift_Ti, &hydrokick_integ);
        hk_fromdrift[bin] = get_exact_factor(Drift_Ti, ti1, &hydrokick_integ);
    }
}

/*! This function integrates the cosmological prefactor for a drift
 *   step between ti0 and ti1. The value returned is
 *  \f[ \int_{a_0}^{a_1} \frac{{\rm d}a}{H(a) a^3}
 *  \f]
 *
 *  A cached value is used for reasons of speed.
 */
double
get_drift_factor(inttime_t t0, inttime_t t1)
{
    /* New step: cache is invalid in this case*/
    if(df_last_ti0 != t0 || df_last_ti1 != t1)
        endrun(2,"get_drift_factor called outside of cache! Cache is: %d -> %d. Want: %d->%d\n",
                df_last_ti0, df_last_ti1, t0, t1);
    return df_last_value;
}

/*! This function integrates the cosmological prefactor for a grav kick
 *   step between ti0 and ti1. The value returned is
 *  \f[ \int_{a_0}^{a_1} \frac{{\rm d}a}{H(a) a^2}
 *  \f]
 *  A per-bin cache is used for reasons of speed.
 *  If bin == -1, just do the integral (this is for a PM step)
 */
double
get_gravkick_factor(inttime_t ti0, inttime_t ti1, int bin)
{
    if(bin < 0) {
        /*This is a cache for the predicted kick from the PM step*/
        if(ti1 == kk_pred_ti)
            return gk_last_PM;
        else
            return get_exact_factor(ti0, ti1, &gravkick_integ);
    }
    if(kk_last_ti0[bin] == ti0 && kk_pred_ti == ti1)
        return gk_todrift[bin];
    /*This is predicting velocities at the drift time*/
    else if(kk_pred_ti == ti0 && kk_last_ti1[bin] == ti1)
        return gk_fromdrift[bin];
    else {
        endrun(2,"get_gravkick_factor bad! bin = %d. Cache is: %d -> %d (drift: %d). Want: %d->%d\n",
                bin, kk_last_ti0[bin], kk_last_ti1[bin], kk_pred_ti, ti0, ti1);
        return 0;
    }
}

/*! This function integrates the cosmological prefactor for a hydro kick
 *   step between ti0 and ti1. The value returned is
 *  \f[ \int_{a_0}^{a_1} \frac{{\rm d}a}{H(a) a^3}
 *  \f]
 *  A per-bin cache is used for reasons of speed.
 *  If bin == -1, just do the integral (this is for a PM step)
 */
double
get_hydrokick_factor(inttime_t ti0, inttime_t ti1, int bin)
{
    if(bin < 0)
        return get_exact_factor(ti0, ti1, &hydrokick_integ);
    if(kk_last_ti0[bin] == ti0 && kk_pred_ti == ti1)
        return hk_todrift[bin];
    /*This is predicting velocities at the drift time*/
    else if(kk_pred_ti == ti0 && kk_last_ti1[bin] == ti1)
        return hk_fromdrift[bin];
    else {
        endrun(2,"get_hydrokick_factor bad! bin = %d. Cache is: %d -> %d (drift: %d). Want: %d->%d\n",
                bin, kk_last_ti0[bin], kk_last_ti1[bin], kk_pred_ti, ti0, ti1);
        return 0;
    }
}
