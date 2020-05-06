#include <stdlib.h>
#include <math.h>
#include "../TYPES/WMMStructs3D.h"
#include <iostream>

#ifndef UTILS_H
#define	UTILS_H

namespace wmm3D {


    void solve_tridiag_equation(double* a, double* b, double* c, double* d, int n) {
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

    void get_spline_coeffs(double *y, double *m, int gamma) {
        double a[gamma+1], b[gamma+1], c[gamma+1], d[gamma+1];
        double h = 1.0 / (double) gamma, h2 = h*h, h_1 = (double) gamma;
        a[0] = 0.0; b[0] = 2.0*h_1; c[0] = h_1; d[0] = 3.0 * (y[1] - y[0]) / h2;
        for (int i=1; i < gamma; i++) {
            a[i] = h_1; c[i] = h_1; b[i] = 4.0*h_1; d[i] = 3.0 * (y[i + 1] - y[i - 1]) / h2;
        }
        a[gamma] = h_1; b[gamma] = 2.0*h_1; c[gamma] = 0.0; d[gamma] = 3.0 * (y[gamma] - y[gamma - 1]) / h2;
        solve_tridiag_equation(a, b, c, d, gamma + 1);
        for (int i=0; i<gamma; i++) {
            m[2*i] = d[i]*h - y[i + 1] + y[i];
            m[2*i + 1] = -d[i+1]*h + y[i + 1] - y[i];
        }
    }

    double checkSlope(double fd0, double fd1, int isborder) {
        double v;
        if (isborder) {
            v = (3.0 * fd0 - fd1) / 2.0;
            if ((v * fd0) <= 0.0)
                v =  0.0;
            else if (((fd1 * fd0) <= 0.0) && (fabs(v) > fabs(3.0 * fd0)))
                v = 3.0 * fd0;
        }
        else {
            v = (fd0 * fd1 <= 0.0) ? 0.0 : 2.0 * fd0 * fd1 / (fd0 + fd1);
        }
        return v;
    }

    double checkHermiteSlope(double d0, double d1) {
        double m = (3.0*d0 - d1)/2.0;
        if (fabs(m) > fabs(3.0*d0))
            m = 3.0*d0;
        return m;
    }

    template <typename T> int sgn(T val) {
        return (T(0) < val) - (val < T(0));
    }

    double GetFValue(Grid3 &image, Node3D &posD, int interp) {

        double value = MAX_VAL;
        Node3 pos(floor(posD.x), floor(posD.y), floor(posD.z)), p(ceil(posD.x), ceil(posD.y), ceil(posD.z)), p2;
        if (pos.x >= image.rows - 2) pos.x = image.rows - 2;
        if (pos.y >= image.cols - 2) pos.y = image.cols - 2;
        if (pos.z >= image.channels - 2) pos.z = image.channels - 2;
        Node3D epsilon = posD - pos;
        if (epsilon.y > 1) epsilon.y = 1;
        if (epsilon.x > 1) epsilon.x = 1;
        if (epsilon.z > 1) epsilon.z = 1;

        double v0 = (image.contains(pos)) ? norm(Node3D(image.at(pos, 0), image.at(pos, 1), image.at(pos, 2))) : MAX_VAL;
        p = pos + Node3(1, 0, 0);
        double v1 = (image.contains(p)) ? norm(Node3D(image.at(p, 0), image.at(p, 1), image.at(p, 2))) : MAX_VAL;
        p = pos + Node3(0, 1, 0);
        double v2 = (image.contains(p)) ? norm(Node3D(image.at(p, 0), image.at(p, 1), image.at(p, 2))) : MAX_VAL;
        p = pos + Node3(0, 0, 1);
        double v3 = (image.contains(p)) ? norm(Node3D(image.at(p, 0), image.at(p, 1), image.at(p, 2))) : MAX_VAL;
        p = pos + Node3(1, 1, 0);
        double v4 = (image.contains(p)) ? norm(Node3D(image.at(p, 0), image.at(p, 1), image.at(p, 2))) : MAX_VAL;
        p = pos + Node3(1, 0, 1);
        double v5 = (image.contains(p)) ? norm(Node3D(image.at(p, 0), image.at(p, 1), image.at(p, 2))) : MAX_VAL;
        p = pos + Node3(0, 1, 1);
        double v6 = (image.contains(p)) ? norm(Node3D(image.at(p, 0), image.at(p, 1), image.at(p, 2))) : MAX_VAL;
        p = pos + Node3(1, 1, 1);
        double v7 = (image.contains(p)) ? norm(Node3D(image.at(p, 0), image.at(p, 1), image.at(p, 2))) : MAX_VAL;

        if (interp == I_PCHIP || interp == I_CUBIC || interp == I_SPLINE) {
            double m[2], v2d[4][4], v[4], v0, v1, dv0, dv1;
            double t = epsilon.z, t_2 = t*t, t_3 = t_2*t;
            int init_x = (pos.x == 0) ? 0 : -1;
            int end_x = (pos.x == image.cols - 2) ? 1 : 2;
            int init_y = (pos.y == 0) ? 0 : -1;
            int end_y = (pos.y == image.rows - 2) ? 1 : 2;
            for (int x=init_x; x<=end_x; x++) {
                for (int y=init_y; y<=end_y; y++) {
                    if (pos.z == 0) {
                        p = pos + Node3(x, y, 1);
                        p2 = pos + Node3(x, y, 0);
                        dv0 = norm(Node3D(image.at(p, 0), image.at(p, 1), image.at(p, 2))) - norm(Node3D(image.at(p2, 0), image.at(p2, 1), image.at(p2, 2)));
                        p = pos + Node3(x, y, 2);
                        p2 = pos + Node3(x, y, 1);
                        dv1 = norm(Node3D(image.at(p, 0), image.at(p, 1), image.at(p, 2))) - norm(Node3D(image.at(p2, 0), image.at(p2, 1), image.at(p2, 2)));
                    }
                    else {
                        p = pos + Node3(x, y, 0);
                        p2 = pos + Node3(x, y, -1);
                        dv0 = norm(Node3D(image.at(p, 0), image.at(p, 1), image.at(p, 2))) - norm(Node3D(image.at(p2, 0), image.at(p2, 1), image.at(p2, 2)));
                        p = pos + Node3(x, y, 1);
                        p2 = pos + Node3(x, y, 0);
                        dv1 = norm(Node3D(image.at(p, 0), image.at(p, 1), image.at(p, 2))) - norm(Node3D(image.at(p2, 0), image.at(p2, 1), image.at(p2, 2)));
                    }
                    m[0] = checkSlope(dv0, dv1, (pos.z == 0));
                    if (pos.z == image.channels - 2) {
                        p = pos + Node3(x, y, 0);
                        p2 = pos + Node3(x, y, -1);
                        dv0 = norm(Node3D(image.at(p, 0), image.at(p, 1), image.at(p, 2))) - norm(Node3D(image.at(p2, 0), image.at(p2, 1), image.at(p2, 2)));
                        p = pos + Node3(x, y, 1);
                        p2 = pos + Node3(x, y, 0);
                        dv1 = norm(Node3D(image.at(p, 0), image.at(p, 1), image.at(p, 2))) - norm(Node3D(image.at(p2, 0), image.at(p2, 1), image.at(p2, 2)));
                    }
                    else {
                        p = pos + Node3(x, y, 1);
                        p2 = pos + Node3(x, y, 0);
                        dv0 = norm(Node3D(image.at(p, 0), image.at(p, 1), image.at(p, 2))) - norm(Node3D(image.at(p2, 0), image.at(p2, 1), image.at(p2, 2)));
                        p = pos + Node3(x, y, 2);
                        p2 = pos + Node3(x, y, 1);
                        dv1 = norm(Node3D(image.at(p, 0), image.at(p, 1), image.at(p, 2))) - norm(Node3D(image.at(p2, 0), image.at(p2, 1), image.at(p2, 2)));
                    }
                    m[1] = checkSlope(dv1, dv0, (pos.z == image.channels - 2));
                    p = pos + Node3(x, y, 0);
                    v0 = norm(Node3D(image.at(p, 0), image.at(p, 1), image.at(p, 2)));
                    p = pos + Node3(x, y, 1);
                    v1 = norm(Node3D(image.at(p, 0), image.at(p, 1), image.at(p, 2)));
                    v2d[x - init_x][y - init_y] = (2.0 * t_3 - 3.0 * t_2 + 1.0) * v0 + (t_3 - 2.0 * t_2 + t) * m[0] +
                           (-2.0 * t_3 + 3.0 * t_2) * v1 + (t_3 - t_2) * m[1];
                }
            }
            t = epsilon.y; t_2 = t*t; t_3 = t_2*t;
            for (int x=0; x<=end_x-init_x; x++) {
                dv0 = v2d[x][1] - v2d[x][0];
                dv1 = v2d[x][2] - v2d[x][1];
                m[0] = checkSlope(dv0, dv1, init_y == 0);
                dv0 = v2d[x][end_y - init_y - 1] - v2d[x][end_y - init_y - 2];
                dv1 = v2d[x][end_y - init_y] - v2d[x][end_y - init_y - 1];
                m[1] = checkSlope(dv1, dv0, end_y == 1);
                v0 = v2d[x][-init_y];
                v1 = v2d[x][1 - init_y];
                v[x] = (2.0 * t_3 - 3.0 * t_2 + 1.0) * v0 + (t_3 - 2.0 * t_2 + t) * m[0] +
                       (-2.0 * t_3 + 3.0 * t_2) * v1 + (t_3 - t_2) * m[1];
            }

            t = epsilon.x; t_2 = t*t; t_3 = t_2*t;
            dv0 = v[1] - v[0];
            dv1 = v[2] - v[1];
            m[0] = checkSlope(dv0, dv1, init_x == 0);
            dv0 = v[end_x - init_x - 1] - v[end_x - init_x - 2];
            dv1 = v[end_x - init_x] - v[end_x - init_x - 1];
            m[1] = checkSlope(dv1, dv0, end_x == 1);
            v0 = v[-init_x];
            v1 = v[1 - init_x];
            value = (2.0 * t_3 - 3.0 * t_2 + 1.0) * v0 + (t_3 - 2.0 * t_2 + t) * m[0] +
                    (-2.0 * t_3 + 3.0 * t_2) * v1 + (t_3 - t_2) * m[1];
        }
        else if (interp == I_HERMITE) {
            double m[2], v2d[4][4], v[4], v0, v1, dv0, dv1;
            double t = epsilon.z, t_2 = t*t, t_3 = t_2*t;
            int init_x = (pos.x == 0) ? 0 : -1;
            int end_x = (pos.x == image.cols - 2) ? 1 : 2;
            int init_y = (pos.y == 0) ? 0 : -1;
            int end_y = (pos.y == image.rows - 2) ? 1 : 2;
            for (int x=init_x; x<=end_x; x++) {
                for (int y=init_y; y<=end_y; y++) {
                    if (pos.z == 0) {
                        p = pos + Node3(x, y, 1);
                        p2 = pos + Node3(x, y, 0);
                        dv0 = norm(Node3D(image.at(p, 0), image.at(p, 1), image.at(p, 2))) - norm(Node3D(image.at(p2, 0), image.at(p2, 1), image.at(p2, 2)));
                        p = pos + Node3(x, y, 2);
                        p2 = pos + Node3(x, y, 1);
                        dv1 = norm(Node3D(image.at(p, 0), image.at(p, 1), image.at(p, 2))) - norm(Node3D(image.at(p2, 0), image.at(p2, 1), image.at(p2, 2)));
                        m[0] = checkHermiteSlope(dv0, dv1);
                    }
                    else {
                        p = pos + Node3(x, y, 1);
                        p2 = pos + Node3(x, y, -1);
                        m[0] = 0.5 * (norm(Node3D(image.at(p, 0), image.at(p, 1), image.at(p, 2))) - norm(Node3D(image.at(p2, 0), image.at(p2, 1), image.at(p2, 2))));
                    }
                    if (pos.z == image.channels - 2) {
                        p = pos + Node3(x, y, 0);
                        p2 = pos + Node3(x, y, -1);
                        dv0 = norm(Node3D(image.at(p, 0), image.at(p, 1), image.at(p, 2))) - norm(Node3D(image.at(p2, 0), image.at(p2, 1), image.at(p2, 2)));
                        p = pos + Node3(x, y, 1);
                        p2 = pos + Node3(x, y, 0);
                        dv1 = norm(Node3D(image.at(p, 0), image.at(p, 1), image.at(p, 2))) - norm(Node3D(image.at(p2, 0), image.at(p2, 1), image.at(p2, 2)));
                        m[1] = checkHermiteSlope(dv1, dv0);
                    }
                    else {
                        p = pos + Node3(x, y, 1);
                        p2 = pos + Node3(x, y, -1);
                        m[1] = 0.5 * (norm(Node3D(image.at(p, 0), image.at(p, 1), image.at(p, 2))) - norm(Node3D(image.at(p2, 0), image.at(p2, 1), image.at(p2, 2))));
                    }
                    p = pos + Node3(x, y, 0);
                    v0 = norm(Node3D(image.at(p, 0), image.at(p, 1), image.at(p, 2)));
                    p = pos + Node3(x, y, 1);
                    v1 = norm(Node3D(image.at(p, 0), image.at(p, 1), image.at(p, 2)));
                    v2d[x - init_x][y - init_y] = (2.0 * t_3 - 3.0 * t_2 + 1.0) * v0 + (t_3 - 2.0 * t_2 + t) * m[0] +
                           (-2.0 * t_3 + 3.0 * t_2) * v1 + (t_3 - t_2) * m[1];
                }
            }
            t = epsilon.y; t_2 = t*t; t_3 = t_2*t;
            for (int x=0; x<=end_x-init_x; x++) {
                if (init_y == 0) {
                    dv0 = v2d[x][1] - v2d[x][0];
                    dv1 = v2d[x][2] - v2d[x][1];
                    m[0] = checkHermiteSlope(dv0, dv1);
                }
                else {
                    m[0] = 0.5 * (v2d[x][2] - v2d[x][0]);
                }
                if (end_y == 1) {
                    dv0 = v2d[x][end_y - init_y - 1] - v2d[x][end_y - init_y - 2];
                    dv1 = v2d[x][end_y - init_y] - v2d[x][end_y - init_y - 1];
                    m[1] = checkHermiteSlope(dv1, dv0);
                }
                else {
                    m[1] = 0.5 * (v2d[x][end_y - init_y] - v2d[x][end_y - init_y - 2]);
                }
                v0 = v2d[x][-init_y];
                v1 = v2d[x][1 - init_y];
                v[x] = (2.0 * t_3 - 3.0 * t_2 + 1.0) * v0 + (t_3 - 2.0 * t_2 + t) * m[0] +
                       (-2.0 * t_3 + 3.0 * t_2) * v1 + (t_3 - t_2) * m[1];
            }

            t = epsilon.x; t_2 = t*t; t_3 = t_2*t;
            if (init_x == 0) {
                dv0 = v[1] - v[0];
                dv1 = v[2] - v[1];
                m[0] = checkHermiteSlope(dv0, dv1);
            }
            else {
                m[0] = 0.5 * (v[2] - v[0]);
            }
            if (end_x == 1) {
                dv0 = v[end_x - init_x - 1] - v[end_x - init_x - 2];
                dv1 = v[end_x - init_x] - v[end_x - init_x - 1];
                m[1] = checkHermiteSlope(dv1, dv0);
            }
            else {
                m[1] = 0.5 * (v[end_x - init_x] - v[end_x - init_x - 2]);
            }
            v0 = v[-init_x];
            v1 = v[1 - init_x];
            value = (2.0 * t_3 - 3.0 * t_2 + 1.0) * v0 + (t_3 - 2.0 * t_2 + t) * m[0] +
                    (-2.0 * t_3 + 3.0 * t_2) * v1 + (t_3 - t_2) * m[1];
        }
        else {
            value = (1.0 - epsilon.x)*(1.0 - epsilon.y)*(1.0 - epsilon.z)*v0 +
                    epsilon.x*(1.0 - epsilon.y)*(1.0 - epsilon.z)*v1 +
                    (1.0 - epsilon.x)*epsilon.y*(1.0 - epsilon.z)*v2 +
                    (1.0 - epsilon.x)*(1.0 - epsilon.y)*epsilon.z*v3 +
                    epsilon.x*epsilon.y*(1.0 - epsilon.z)*v4 +
                    epsilon.x*(1.0 - epsilon.y)*epsilon.z*v5 +
                    (1.0 - epsilon.x)*epsilon.y*epsilon.z*v6 +
                    epsilon.x*epsilon.y*epsilon.z*v7;
        }

        return value;
    }

    double GetIntegralValue(Grid3 &image, Node3D &posD, Node3D &neighD, double fposD, double fneighD,
                            int interp, int inner_nodes=1) {
        double values = fposD + fneighD, step = 1.0 / (inner_nodes + 1.0);
        Node3D p;
        if (inner_nodes > 0) {
            for (int i=1; i <= inner_nodes; i++) {
                p = (1.0 - i * step) * posD + i * step * neighD;
                values += 2.0 * GetFValue(image, p, interp);
            }
        }
        values /= ((inner_nodes + 1.0) * 2.0);
        return values;
    }
}
#endif	/* UTILS_H */

