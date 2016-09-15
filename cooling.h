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

void GetParticleUVBG(int i, struct UVBG * uvbg, const double Time);
void GetGlobalUVBG(struct UVBG * uvbg);
double AbundanceRatios(double u, double rho, struct UVBG * uvbg, double *ne_guess, double *nH0_pointer, double *nHeII_pointer);
double GetCoolingTime(double u_old, double rho, struct UVBG * uvbg,  double *ne_guess, double Z, const double Time);
double DoCooling(double u_old, double rho, double dt, struct UVBG * uvbg, double *ne_guess, double Z, const double Time);
double ConvertInternalEnergy2Temperature(double u, double ne);

void InitCool(int CoolingOn, const double TimeBegin, const char * TreeCoolFile, const char * MetalCoolFile, const char * UVFluctuationFile, const double UnitDensity_in_cgs, const double HubbleParam, const double UnitTime_in_s, const double UnitPressure_in_cgs, const double MinGasTemp);
void   IonizeParams(const double Time);
void   MakeCoolingTable(const double MinGasTemp);
void   SetZeroIonization(void);

#endif
