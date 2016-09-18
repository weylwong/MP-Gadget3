#ifndef _COOLING_H_
#define _COOLING_H_
struct UVBG {
    double J_UV;
    double gJH0;
    double gJHep;
    double gJHe0;
    double epsH0;
    double epsHep;
    double epsHe0;
} ;

void GetParticleUVBG(const double * pos, struct UVBG * uvbg, const double Time);
double AbundanceRatios(double u, double rho, struct UVBG * uvbg, double *ne_guess, double *nH0_pointer, double *nHeII_pointer);
double GetCoolingTime(double u_old, double rho, struct UVBG * uvbg,  double *ne_guess, double Z, const double Time);
double DoCooling(double u_old, double rho, double dt, struct UVBG * uvbg, double *ne_guess, double Z, const double Time);
double ConvertInternalEnergy2Temperature(double u, double ne);

void InitCool(int CoolingOn, const double TimeBegin, const char * TreeCoolFile, const char * MetalCoolFile, const char * UVFluctuationFile, const double UnitDensity_in_cgs, const double HubbleParam, const double UnitTime_in_s, const double UnitPressure_in_cgs, const double MinGasTemp);
void   IonizeParams(const double Time);
void   MakeCoolingTable(const double MinGasTemp);

#ifdef _COOLING_PRIVATE
/*These functions are private to cooling.c. Only here for testing.*/
struct abundance {
    double ne;
    double nH0;
    double nHp;
    double nHe0;
    double nHep;
    double nHepp;
};

struct rates {
    double aHp;
    double aHep;
    double aHepp;
    double ad;
    double geH0;
    double geHe0;
    double geHep;
    double bH0;
    double bHep;
    double bff;
};

void find_abundances_and_rates(double logT, double nHcgs, struct UVBG * uvbg, struct abundance * y, struct rates * r);
double solve_equilibrium_temp(double u, double nHcgs, struct UVBG * uvbg, struct abundance * y);

double PrimordialCoolingRate(double logT, double nHcgs, struct UVBG * uvbg, double *nelec, const double redshift);
double CoolingRateFromU(double u, double nHcgs, struct UVBG * uvbg, double *ne_guess, double Z, const double Time);
void InitCoolMemory(void);

#endif

#endif
