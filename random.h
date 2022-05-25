#ifndef _RANDOM_H_
#define _RANDOM_H_

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>

extern void   rnd_init();
extern double rnd_get();    //(0, 1]
extern void rnd_seed_set();
extern void rnd_free();

inline int rnd_poisson(const double &lamda)
{
    const double L = exp(-lamda);
    int k = 0;
    double p = double(1.0);
    do {
        k++;
        p *= rnd_get();
    } while (p > L);
    return k - 1;
}

inline double rnd_gauss() {
    static double rnd_gauss2;
    static bool rnd_gauss2_exist = false;
    
    //---
    if (rnd_gauss2_exist) {
        rnd_gauss2_exist = false;
        return rnd_gauss2;
    } else {
        double x, y, s;
        do {
            x = rnd_get() * double (2.0) - double(1.0);
            y = rnd_get() * double (2.0) - double(1.0);
            s = x * x + y * y;
        } while (s == double(0.0) || s > double(1.0));
        s = sqrt((-double(2.0) * log(s)) / s);
        rnd_gauss2 = y * s;
        rnd_gauss2_exist = true;
        return x * s;
    }
}
    //!!!BE CAREFFULY FOR MULTITHREADS APPLICATIONS!!!
    //-->NOT BE USE FOR IT!!!

/*algorithm poisson random number (Knuth):
    init:
         Let L ← exp(−λ), k ← 0 and p ← 1.
    do:
         k ← k + 1.
         Generate uniform random number u in [0,1] and let p ← p × u.
    while p > L.
    return k − 1.
 */

#endif
