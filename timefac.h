#ifndef __TIMEFAC_H
#define __TIMEFAC_H

void init_drift_table(double timeBegin, double timeMax);
double get_hydrokick_factor(inttime_t ti0, inttime_t ti1);
double get_gravkick_factor(inttime_t ti0, inttime_t ti1);

/* Update the last used drift cache.
 * If reset is true, the new value is not checked to be just after
 * the cached value. */
void update_drift_factor(inttime_t ti0, inttime_t ti1, int reset);
/* Get the drift factor at given time. Run is terminated if
 * ti0 and ti1 are not in the cache*/
double get_drift_factor(inttime_t ti0, inttime_t ti1);

#endif
