#ifndef __TIMEFAC_H
#define __TIMEFAC_H

/* Update the kick cache. This stores a value for every timebin,
 * both to current drift time and for the kick timestep.*/
void update_kick_factors(inttime_t PM_Ti, inttime_t Drift_Ti, int reset);

/*Get the kick factors for a given bin. The time variables are passed as a check.
 * If the bin is -1, the value is recomputed (unless it matches the PM kick step).
 */
double get_hydrokick_factor(inttime_t ti0, inttime_t ti1, int bin);
double get_gravkick_factor(inttime_t ti0, inttime_t ti1, int bin);

/* Update the last used drift cache.
 * If reset is true, the new value is not checked to be just after
 * the cached value. */
void update_drift_factor(inttime_t ti0, inttime_t ti1, int reset);
/* Get the drift factor at given time. Run is terminated if
 * ti0 and ti1 are not in the cache*/
double get_drift_factor(inttime_t ti0, inttime_t ti1);

#endif
