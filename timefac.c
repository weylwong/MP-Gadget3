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
#define DRIFT_TABLE_LENGTH 2000

static double logTimeInit;
static double logTimeMax;

/*! table for the cosmological kick factor for gravitational forces */
static double GravKickTable[DRIFT_TABLE_LENGTH];

/*! table for the cosmological kick factor for hydrodynmical forces */
static double HydroKickTable[DRIFT_TABLE_LENGTH];

/* Simple single-value cache for the drift table,
 * which is the same for all particles*/
static inttime_t df_last_ti0 = -1, df_last_ti1 = -1;
static double df_last_value;

static inttime_t hk_last_ti0 = -1, hk_last_ti1 = -1;
static double hk_last_value;
#pragma omp threadprivate(hk_last_ti0, hk_last_ti1, hk_last_value)

static inttime_t gk_last_ti0 = -1, gk_last_ti1 = -1;
static double gk_last_value;
#pragma omp threadprivate(gk_last_ti0, gk_last_ti1, gk_last_value)

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
static double get_exact_factor(inttime_t t0, inttime_t t1, double (*factor) (double, void *))
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
void update_drift_factor(inttime_t t0, inttime_t t1, int reset)
{
    /* Rebuild the drift cache*/
    if(!reset && df_last_ti1 > 0 && df_last_ti1 != t0)
        endrun(2,"Drift should start where last step finished! Before had: %d -> %d, now %d->%d",
            df_last_ti0, df_last_ti1, t0, t1);
    df_last_ti0 = t0;
    df_last_ti1 = t1;
    df_last_value = get_exact_factor(t0, t1, &drift_integ);
}

/*! This function integrates the cosmological prefactor for a drift
 *   step between ti0 and ti1. The value returned is
 *  \f[ \int_{a_0}^{a_1} \frac{{\rm d}a}{H(a) a^3}
 *  \f]
 *
 *  A cached value is used for reasons of speed.
 */
double get_drift_factor(int t0, int t1)
{
    /* New step: cache is invalid in this case*/
    if(df_last_ti0 != t0 || df_last_ti1 != t1)
        endrun(2,"get_drift_factor called outside of cache! Cache is: %d -> %d. Want: %d->%d",
                df_last_ti0, df_last_ti1, t0, t1);
    return df_last_value;
}

/*
double get_gravkick_factor(int t1, int t2)
{
    cache_lookup(t1, t2, 1);
    return get_exact_factor(t1, t2, 1);
}

double get_hydrokick_factor(int t1, int t2)
{
    cache_lookup(t1, t2, 2);
    return get_exact_factor(t1, t2, 2);
}
*/

/*Variables for the cache*/
/* In practice we will always be asking for
 * drift factors from one timestep to the current one.
 * So a cache of size corresponding to the number of timebins suffices.
 * Kick factors are more complicated: we will be asking for kicks from
 * current timestep to a midpoint timestep, but again all particles will
 * be covered by the number of timebins, plus one for the PM timestep.*/
/*
#define CACHESIZE 32;
int NCache=0;
int last_t_drift = -1;
int last_t0[TIMEBINS];
double drift_factors[TIMEBINS];
double gravkick_factors[TIMEBINS];
double hydrokick_factors[TIMEBINS];

void refresh_cache(int Ti_Current, int *Ti_drift)
{
    int i;
    for(i=0; i< TIMEBINS; i++) {
        if(Ti_drift[i] > 0) {
            last_t0[i] = Ti_drift[i];
            drift_factors[i] = get_exact_factor(a1, a2, 0);
            gravkick_factors[i] = get_exact_factor(a1, a2, 1);
            hydrokick_factors[i] = get_exact_factor(a1, a2, 2);
        }
        else
            last_t0[i] = -1;
    }
    last_t11 = Ti_Current;
}
*/

void init_drift_table(double timeBegin, double timeMax)
{
  int i;
  double result, abserr;

  gsl_function F;
  gsl_integration_workspace *workspace;

  logTimeInit = log(timeBegin);
  logTimeMax = log(timeMax);
  if(logTimeMax <=logTimeInit)
      endrun(1,"Error: Invalid drift table range: (%d->%d)\n", timeBegin, timeMax);

  workspace = gsl_integration_workspace_alloc(WORKSIZE);

  for(i = 0; i < DRIFT_TABLE_LENGTH; i++)
    {
      F.function = &gravkick_integ;
      gsl_integration_qag(&F, exp(logTimeInit),
			  exp(logTimeInit + ((logTimeMax - logTimeInit) / DRIFT_TABLE_LENGTH) * (i + 1)), 0,
			  1.0e-8, WORKSIZE, GSL_INTEG_GAUSS41, workspace, &result, &abserr);
      GravKickTable[i] = result;


      F.function = &hydrokick_integ;
      gsl_integration_qag(&F, exp(logTimeInit),
			  exp(logTimeInit + ((logTimeMax - logTimeInit) / DRIFT_TABLE_LENGTH) * (i + 1)), 0,
			  1.0e-8, WORKSIZE, GSL_INTEG_GAUSS41, workspace, &result, &abserr);
      HydroKickTable[i] = result;

    }
  gsl_integration_workspace_free(workspace);
  df_last_ti0 = df_last_ti1 = gk_last_ti0 = gk_last_ti1 = hk_last_ti0 = hk_last_ti1 = -1;
}

/*Find which bin in the table we are looking up.
 * Pointer argument gives the full floating point value for interpolation.*/
int find_bin_number(inttime_t ti0, double *rem)
{
  double a1 = loga_from_ti(ti0);
  double u1;
  int i1;
  u1 = (a1 - logTimeInit) / (logTimeMax - logTimeInit) * DRIFT_TABLE_LENGTH;
  i1 = (int) u1;
  /*Bound u1*/
  if(i1 >= DRIFT_TABLE_LENGTH)
    i1 = DRIFT_TABLE_LENGTH - 1;
  if(i1 <=1)
      i1=1;
  *rem = u1;
  return i1;
}


/*! This function integrates the cosmological prefactor for a drift
 *   step between ti0 and ti1. The value returned is
 *  \f[ \int_{a_0}^{a_1} \frac{{\rm d}a}{H(a) a^3}
 *  \f]
 *  
 *  A lookup-table is used for reasons of speed. 
 */


double get_gravkick_factor(inttime_t ti0, inttime_t ti1)
{
  double df1, df2, u1, u2;
  int i1, i2;

  if(ti0 == gk_last_ti0 && ti1 == gk_last_ti1)
    return gk_last_value;

  /* note: will only be called for cosmological integration */
  i1 = find_bin_number(ti0, &u1);
  if(i1 <= 1)
    df1 = u1 * GravKickTable[0];
  else
    df1 = GravKickTable[i1 - 1] + (GravKickTable[i1] - GravKickTable[i1 - 1]) * (u1 - i1);

  i2 = find_bin_number(ti1, &u2);
  if(i2 <= 1)
    df2 = u2 * GravKickTable[0];
  else
    df2 = GravKickTable[i2 - 1] + (GravKickTable[i2] - GravKickTable[i2 - 1]) * (u2 - i2);

  gk_last_ti0 = ti0;
  gk_last_ti1 = ti1;

  return gk_last_value = (df2 - df1);
}

double get_hydrokick_factor(inttime_t ti0, inttime_t ti1)
{
  double df1, df2,u1,u2;
  int i1, i2;

  if(ti0 == hk_last_ti0 && ti1 == hk_last_ti1)
    return hk_last_value;

  /* note: will only be called for cosmological integration */

  i1 = find_bin_number(ti0, &u1);
  if(i1 <= 1)
    df1 = u1 * HydroKickTable[0];
  else
    df1 = HydroKickTable[i1 - 1] + (HydroKickTable[i1] - HydroKickTable[i1 - 1]) * (u1 - i1);

  i2 = find_bin_number(ti1, &u2);
  if(i2 <= 1)
    df2 = u2 * HydroKickTable[0];
  else
    df2 = HydroKickTable[i2 - 1] + (HydroKickTable[i2] - HydroKickTable[i2 - 1]) * (u2 - i2);

  hk_last_ti0 = ti0;
  hk_last_ti1 = ti1;

  return hk_last_value = (df2 - df1);
}
