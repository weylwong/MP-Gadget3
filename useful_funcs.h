#ifndef _USEFUL_FUNCS_
#define _USEFUL_FUNCS_
/*Various functions that should really be in the C standard library!*/

/*Arguably this belongs in physconst.h, but it is so unbelievably common in the code I am putting it in here!*/
#ifndef  GAMMA
#define  GAMMA         (5.0/3.0)	/*!< adiabatic index of simulated gas */
#endif

#define  GAMMA_MINUS1  (GAMMA-1)

static inline double
dotproduct(double v1[3], double v2[3])
{
    double r =0;
    int d;
    for(d = 0; d < 3; d ++) {
        r += v1[d] * v2[d];
    }
    return r;
}

static inline void crossproduct(double v1[3], double v2[3], double out[3])
{
    static int D2[3] = {1, 2, 0};
    static int D3[3] = {2, 0, 1};

    int d1, d2, d3;

    for(d1 = 0; d1 < 3; d1++)
    {
        d2 = D2[d1];
        d3 = D3[d1];

        out[d1] = (v1[d2] * v2[d3] -  v2[d2] * v1[d3]);
    }
}

static inline double DMAX(double a, double b) {
    if(a > b) return a;
    return b;
}
static inline double DMIN(double a, double b) {
    if(a < b) return a;
    return b;
}
static inline int IMAX(int a, int b) {
    if(a > b) return a;
    return b;
}
static inline int IMIN(int a, int b) {
    if(a < b) return a;
    return b;
}

#endif
