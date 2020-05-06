#include <stdlib.h>
#include <math.h>
#include <iostream>
// #include "../TYPES/structs.h"

#ifndef INTERPOLATION_H
#define	INTERPOLATION_H

namespace wmm {

    template<class T> void solve_tridiag_equation(T* a, T* b, T* c, T* d, int n) {
        n--;
        c[0] /= b[0];
        d[0] /= b[0];

        for (int i = 1; i < n; i++) {
            c[i] /= b[i] - a[i]*c[i-1];
            d[i] = (d[i] - a[i]*d[i-1]) / (b[i] - a[i]*c[i-1]);
        }

        d[n] = (d[n] - a[n]*d[n-1]) / (b[n] - a[n]*c[n-1]);

        for (int i = n; i-- > 0;) {
            d[i] -= c[i]*d[i+1];
        }
    }



}
#endif	/* INTERPOLATION_H */

