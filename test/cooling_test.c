#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <cmocka.h>
#include <math.h>
#include <stdio.h>
#define _COOLING_PRIVATE
#include "../cooling.h"
#include "../physconst.h"

#define YHELIUM ((1 - HYDROGEN_MASSFRAC) / (4 * HYDROGEN_MASSFRAC))

static int setup_cooling(void ** state)
{
    InitCoolMemory();
    MakeCoolingTable(200);
    return 0;
}

static void test_eq_temp(void ** state)
{
    double uin, rhoin, muin, nein;

    //tempin = 34.0025;
    uin = 6.01329e+09;
    rhoin = 7.85767e-29;
    muin = 0.691955;
    static struct UVBG uvbg = {0};
    nein = (1 + 4 * YHELIUM) / muin - (1 + YHELIUM);
    struct abundance y = {0};
    y.ne = nein;
    double nHcgs = rhoin * HYDROGEN_MASSFRAC / PROTONMASS;
    assert_true(fabs(solve_equilibrium_temp(uin, nHcgs, &uvbg, &y) - 59.1275) < 1e5);
}

int main(void) {
    const struct CMUnitTest tests[] = {
        cmocka_unit_test(test_eq_temp),
    };
    return cmocka_run_group_tests(tests, setup_cooling, NULL);
}
