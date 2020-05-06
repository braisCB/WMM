#include <stdlib.h>
#include <math.h>
#include <iostream>
// #include "../TYPES/structs.h"

#ifndef INTERPOLATION_H
#define	INTERPOLATION_H

namespace wmm {

    template<typename T> T checkPchipSlope(const T d0, const T d1, int isborder) {
        if (isborder) {
            T v = (T(3) * d0 - d1) / T(2);
            if ((v * d0) <= 0)
                return T(0);
            else if (((d0 * d1) <= T(0)) && (fabs(v) > fabs(T(3) * d0)))
                return T(3) * d0;
        }
        else {
            return (d1 * d0 <= T(0)) ? T(0) : T(2) * d1 * d0 / (d1 + d0);
        }
    }

    /*template<typename T> Array_<T> checkPchipSlope(const Array_<T> &d0, const Array_<T> &d1, int isborder) {
        Array_<T> result(d0.shape, d0.ndim, d0.mode);
        for (int i=0; i <= d0.size; i++) {
            result.data[i] = checkPchipSlope(d0.data[i], d1.data[i], isborder);
        }
        return result;
    }*/

    template<typename T> T checkHermiteSlope(const T d0, const T d1) {
        T m = (T(3)*d0 - d1)/T(2);
        return (fabs(m) > fabs(T(3)*d0)) ? T(3)*d0 : m;
    }

    /*template<typename T> Array_<T> checkSlope(Array_<T> &d0, Array_<T> &d1) {
        Array_<T> result(d0.shape, d0.ndim, d0.mode);
        for (int i=0; i <= d0.size; i++) {
            result.data[i] = checkHermiteSlope(d0.data[i], d1.data[i]);
        }
        return result;
    }*/

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

    template<class T> void get_spline_coeffs(T *y, T *m, int gamma, int init, int end) {
        T a[gamma], b[gamma], c[gamma], d[gamma];
        a[0] = T(0); b[0] = T(2); c[0] = T(1); d[0] = T(3) * (y[1] - y[0]);
        for (int i=1; i < gamma - 1; i++) {
            a[i] = T(1); c[i] = T(1); b[i] = T(4); d[i] = T(3) * (y[i + 1] - y[i - 1]);
        }
        a[gamma-1] = T(1); b[gamma-1] = T(2); c[gamma-1] = T(0); d[gamma-1] = T(3) * (y[gamma-1] - y[gamma - 2]);
        solve_tridiag_equation(a, b, c, d, gamma);
        for (int i=init; i<end; i++) {
            m[2*(i - init)] = d[i] - y[i + 1] + y[i];
            m[2*(i - init) + 1] = -d[i+1] + y[i + 1] - y[i];
        }
    }

    template<class T> void get_hermite_coeffs(T *y, T *m, int gamma, int init, int end) {
        int start = init;
        if (init == 0) {
            T d0 = y[1] - y[0], d1 = y[2] - y[1];
            m[0] = checkHermiteSlope(d0, d1);
            start = 1;
        }
        if (end == gamma) {
            T d0 = y[gamma - 1] - y[gamma - 2];
            T d1 = y[gamma - 2] - y[gamma - 3];
            m[gamma - init - 1] = checkHermiteSlope(d0, d1);
            end = gamma - 1;
        }
        for (int i=start; i<end; i++) {
            m[i - init] = (y[i + 1] - y[i - 1]) / T(2);
        }
    }

    template<class T> void get_pchip_coeffs(T *y, T *m, int gamma, int init, int end) {
        int start = init;
        if (init == 0) {
            T d0 = y[1] - y[0];
            T d1 = y[2] - y[1];
            m[0] = checkPchipSlope(d0, d1, 1);
            start = 1;
        }
        if (end == gamma) {
            T d0 = y[gamma - 1] - y[gamma - 2];
            T d1 = y[gamma - 2] - y[gamma - 3];
            m[gamma - init] = checkPchipSlope(d0, d1, 1);
            end = gamma - 1;
        }
        for (int i=start; i<end; i++) {
            T d0 = y[i] - y[i - 1];
            T d1 = y[i + 1] - y[i];
            m[i - init] = checkPchipSlope(d0, d1, 0);
        }
    }

    template<class T> void get_quad_coeffs(T *y, T *m, int init, int end) {
        int start = init;
        if (init == 0) {
            m[0] = y[0];
            m[2] = (y[2] - T(2) * y[1] + y[0]) / T(2);
            m[1] = y[1] - y[0] - m[2];
            start = 1;
        }
        for (int i=start; i<end; i++) {
            m[3*(i - init)] = y[i];
            m[3*(i - init) + 2] = (y[i+1] - T(2) * y[i] + y[i - 1]) / T(2);
            m[3*(i - init) + 1] = y[i] - y[i - 1] + m[3*(i - init) + 2];
        }
    }

    template<class T> T get_linear_score(T y0, T y1, double t) {
        return y0 + (y1 - y0) * t;
    }

    template<class T> T get_spline_score(T y0, T y1, T m0, T m1, double t) {
        return y0 + (y1 - y0 + m0) * t + (m1 - T(2) * m0) * t * t + (m1 - m0) * t * t * t;
    }

    template<class T> T get_hermite_score(T y0, T y1, T m0, T m1, double t) {
        return y0 + m0 * t + (T(3) * (y1 - y0) + m1 - T(2) * m0) * t * t + (m1 + m0 - T(2) * (y1 - y0)) * t * t * t;
    }

    template<class T> T get_pchip_score(T y0, T y1, T m0, T m1, double t) {
        return get_hermite_score(y0, y1, m0, m1, t);
    }

    template<class T> T get_quad_score(T m0, T m1, T m2, double t) {
        return m0 + m1 * t + m2 * t*t;
    }

    template<class T> T get_spline_integral_score(T y0, T y1, T m0, T m1, double t0=0.0, double t1=1.0) {
        return y0 * (t1 - t0) + (y1 - y0 + m0) / T(2) * (t1 * t1 - t0 * t0) + \
               (m1 - T(2) * m0) / T(3) * (t1 * t1 * t1 - t0 * t0 * t0) + \
               (m1 - m0) / T(4) * (t1 * t1 * t1 * t1 - t0 * t0 * t0 * t0);
    }

    template<class T> T get_linear_integral_score(T y0, T y1, double t0=0.0, double t1=1.0) {
        return y0 * (t1 - t0)  + (y1 - y0) / T(2) * (t1 * t1 - t0 * t0);
    }

    template<class T> T get_hermite_integral_score(T y0, T y1, T m0, T m1, double t0=0.0, double t1=1.0) {
        return y0 * (t1 - t0) + m0 / T(2) * (t1 * t1 - t0 * t0) + \
               (T(3) * (y1 - y0) + m1 - T(2) * m0) / T(3) * (t1 * t1 * t1 - t0 * t0 * t0) + \
               (m1 + m0 - T(2) * (y1 - y0)) / T(4) * (t1 * t1 * t1 * t1 - t0 * t0 * t0 * t0);
    }

    template<class T> T get_pchip_integral_score(T y0, T y1, T m0, T m1, double t0=0.0, double t1=1.0) {
        return get_hermite_integral_score(y0, y1, m0, m1, t0, t1);
    }

    template<class T> T get_quad_integral_score(T m0, T m1, T m2, double t0=0.0, double t1=1.0) {
        return m0 * (t1 - t0) + m1 / T(2) * (t1 * t1 - t0 * t0) + m2 / T(3) * (t1 * t1 * t1 - t0 * t0 * t0);
        /*return (t1 - t0) / T(6) * (get_quad_score(m0, m1, m2, t0) +
                                   T(4) * get_quad_score(m0, m1, m2, 0.5 * (t1 - t0)) +
                                   get_quad_score(m0, m1, m2, t1));*/
    }

}
#endif	/* INTERPOLATION_H */

