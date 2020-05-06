#include <stdlib.h>
#include <math.h>
#include "../TYPES/WMMStructs.h"
#include "../TYPES/interpolation.h"
#include <iostream>

#ifndef UTILS_H
#define	UTILS_H

namespace wmm {


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

    void setCoeffs2D(double *y, double *m, int interp, int gamma) {
        switch (interp) {
            case wmm::I_LINEAR:
                break;
            case wmm::I_QUADATRIC:
                m[0] = y[0];
                m[2] = (y[2] + y[0] - 2.0*y[1])/2.0;
                m[1] = y[1] - y[0] - m[2];
                for (int i=1; i < gamma; i++) {
                    m[3*i] = y[i];
                    m[3*i + 2] = (y[i+1] + y[i-1] - 2.0*y[i])/2.0;
                    m[3*i + 1] = y[i] - y[i-1] + m[3*i + 2];
                }
                break;
            case wmm::I_SPLINE:
            {
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
                break;
            }
            case wmm::I_HERMITE: {
                double d0 = (y[1] - y[0]), d1 = (y[2] - y[1]);
                m[0] = checkHermiteSlope(d0, d1);
                d0 = (y[gamma] - y[gamma - 1]);
                d1 = (y[gamma - 1] - y[gamma - 2]);
                m[gamma] = checkHermiteSlope(d0, d1);
                for (int i=1; i<gamma; i++) {
                    m[i] = 0.5*(y[i + 1] - y[i - 1]);
                }
                break;
            }
            default: {//I_PCHIP
                double d0 = (y[1] - y[0]);
                double d1 = (y[2] - y[1]);
                m[0] = checkSlope(d0, d1, 1);
                m[1] = (d0 * d1 <= 0.0) ? 0.0 : 2.0 * d0 * d1 / (d0 + d1);
                for (int i=2; i<gamma; i++) {
                    d0 = d1;
                    d1 = (y[i + 1] - y[i]);
                    m[i] = (d0 * d1 <= 0.0) ? 0.0 : 2.0 * d0 * d1 / (d0 + d1);
                }
                d0 = (y[gamma] - y[gamma - 1]);
                d1 = (y[gamma - 1] - y[gamma - 2]);
                m[gamma] = checkSlope(d0, d1, 1);
            }
        }

    }

    template <typename T> int sgn(T val) {
        return (T(0) < val) - (val < T(0));
    }

    double GetFValue(Grid &image, NodeD &posD, int interp) {

        double value = MAX_VAL;
        Node pos(floor(posD.y), floor(posD.x)), p(ceil(posD.y), ceil(posD.x)), p2;
        if (pos.y >= image.rows - 2) pos.y = image.rows - 2;
        if (pos.x >= image.cols - 2) pos.x = image.cols - 2;
        NodeD epsilon = posD - pos;
        if (epsilon.y > 1) epsilon.y = 1;
        if (epsilon.x > 1) epsilon.x = 1;

        double v0 = (image.contains(pos)) ? norm(NodeD(image.at(pos, 0), image.at(pos, 1))) : MAX_VAL;
        p = pos + Node(1, 0);
        double v1 = (image.contains(p)) ? norm(NodeD(image.at(p, 0), image.at(p, 1))) : MAX_VAL;
        p = pos + Node(0, 1);
        double v2 = (image.contains(p)) ? norm(NodeD(image.at(p, 0), image.at(p, 1))) : MAX_VAL;
        p = pos + Node(1, 1);
        double v3 = (image.contains(p)) ? norm(NodeD(image.at(p, 0), image.at(p, 1))) : MAX_VAL;

        if (interp == I_LINEAR) {
            value = (1.0 - epsilon.y)*(1.0 - epsilon.x)*v0 + epsilon.y*(1.0 - epsilon.x)*v1 +
                    (1.0 - epsilon.y)*epsilon.x*v2 + epsilon.y*epsilon.x*v3;
        }
        /*else if (interp == I_QUADATRIC) {
            double m[4];
            m[0] = v0;
            p = pos + Node(1, 0);
            m[1] = v1 - m[0];
            p = pos + Node(0, 1);
            m[2] = v2 - m[0];
            p = pos + Node(1, 1);
            m[3] = v3 - m[2] - m[1] - m[0];

            value = m[0] + epsilon.y*m[1] + epsilon.x*m[2] + epsilon.y*epsilon.x*m[3];
        }*/
        else if (interp == I_QUADATRIC) {
            double m[9];
            int dy = (pos.y > 0) ? -1 : 2;
            int dx = (pos.x > 0) ? -1 : 2;
            p = pos + Node(1, 0);
            double p1 = norm(NodeD(image.at(p, 0), image.at(p, 1)));
            p = pos + Node(dy, 0);
            double p2 = norm(NodeD(image.at(p, 0), image.at(p, 1)));
            p = pos + Node(0, 1);
            double p3 = norm(NodeD(image.at(p, 0), image.at(p, 1)));
            p = pos + Node(0, dx);
            double p4 = norm(NodeD(image.at(p, 0), image.at(p, 1)));
            m[0] = norm(NodeD(image.at(pos, 0), image.at(pos, 1)));
            m[4] = (p2 - (1.0 - dx) * m[0] - dx * p1) / (dx * (dx - 1.0));
            m[1] = p1 - m[4] - m[0];
            m[5] = (p4 - (1.0 - dy) * m[0] - dy * p3) / (dy * (dy - 1.0));
            m[2] = p3 - m[5] - m[0];
            p = pos + Node(1, 1);
            double a_1 = norm(NodeD(image.at(p, 0), image.at(p, 1))) - m[0] - m[1] - m[2] - m[4] - m[5];
            p = pos + Node(1, dx);
            double a_3 = norm(NodeD(image.at(p, 0), image.at(p, 1))) - m[0] - m[1] - m[4] - dy * m[2] - dy * dy * m[5];
            p = pos + Node(dy, 1);
            double a_2 = norm(NodeD(image.at(p, 0), image.at(p, 1))) - m[0] - dx * m[1] - dx * dx * m[4] - m[2] - m[5];
            p = pos + Node(dy, dx);
            double a_4 = norm(NodeD(image.at(p, 0), image.at(p, 1))) - m[0] - dx * m[1] - dx * dx * m[4] - dy * m[2] - dy * dy * m[5];
            m[8] = (dx * dy * a_1 - dy * a_2 - dx * a_3 + a_4) /
                   (dx * dy - dx * dy * dy - dx * dx * dy + dx * dx * dy * dy);
            m[7] = (dy * a_2 - a_4 - (dx * dx * dy - dx * dx * dy * dy) * m[8]) / (dx * dy - dx * dy * dy);
            m[6] = (dx * a_3 - a_4 - (dx * dy * dy - dx * dx * dy * dy) * m[8]) / (dx * dy - dx * dx * dy);
            m[3] = (a_4 - dx * dx * dy * dy * m[8] - dx * dy * dy * m[7] - dx * dx * dy * m[6]) / (dx * dy);
            value = m[0] + epsilon.y*m[1] + epsilon.x*m[2] + epsilon.y*epsilon.x*m[3] +
                    epsilon.y*epsilon.y*m[4] + epsilon.x*epsilon.x*m[5] +
                    epsilon.y*epsilon.y*epsilon.x*m[6] + epsilon.y*epsilon.x*epsilon.x*m[7] +
                    epsilon.y*epsilon.y*epsilon.x*epsilon.x*m[8];
        }
        else if (interp == I_HERMITE) {
            double m[2], v[4], v0, v1;
            double t = epsilon.y, t_2 = t*t, t_3 = t_2*t;
            int init = (pos.x == 0) ? 0 : -1;
            int end = (pos.x == image.cols - 2) ? 1 : 2;

            for (int i=init; i<=end; i++) {
                if (pos.y == 0) {
                    p = pos + Node(1, i);
                    p2 = pos + Node(0, i);
                    m[0] = norm(NodeD(image.at(p, 0), image.at(p, 1))) - norm(NodeD(image.at(p2, 0), image.at(p2, 1)));
                }
                else {
                    p = pos + Node(1, i);
                    p2 = pos + Node(-1, i);
                    m[0] = 0.5 * (norm(NodeD(image.at(p, 0), image.at(p, 1))) - norm(NodeD(image.at(p2, 0), image.at(p2, 1))));
                }
                if (pos.y == image.rows - 2) {
                    p = pos + Node(1, i);
                    p2 = pos + Node(0, i);
                    m[1] = norm(NodeD(image.at(p, 0), image.at(p, 1))) - norm(NodeD(image.at(p2, 0), image.at(p2, 1)));
                }
                else {
                    p = pos + Node(2, i);
                    p2 = pos + Node(0, i);
                    m[1] = 0.5 * (norm(NodeD(image.at(p, 0), image.at(p, 1))) - norm(NodeD(image.at(p2, 0), image.at(p2, 1))));
                }
                p = pos + Node(0, i);
                v0 = norm(NodeD(image.at(p, 0), image.at(p, 1)));
                p = pos + Node(1, i);
                v1 = norm(NodeD(image.at(p, 0), image.at(p, 1)));
                v[i - init] = (2.0 * t_3 - 3.0 * t_2 + 1.0) * v0 + (t_3 - 2.0 * t_2 + t) * m[0] +
                       (-2.0 * t_3 + 3.0 * t_2) * v1 + (t_3 - t_2) * m[1];
            }
            t = epsilon.x; t_2 = t*t; t_3 = t_2*t;
            m[0] = (init == 0) ? v[1] - v[0] : (v[2] - v[0]) / 2.0;
            m[1] = (end == 1) ? v[2] - v[1] : (v[end - init] - v[end - init - 2]) / 2.0;
            v0 = v[-init];
            v1 = v[1 - init];
            value = (2.0 * t_3 - 3.0 * t_2 + 1.0) * v0 + (t_3 - 2.0 * t_2 + t) * m[0] +
                    (-2.0 * t_3 + 3.0 * t_2) * v1 + (t_3 - t_2) * m[1];

            t = epsilon.x; t_2 = t*t; t_3 = t_2*t;

            init = (pos.y == 0) ? 0 : -1;
            end = (pos.y == image.rows - 2) ? 1 : 2;
            for (int i=init; i<=end; i++) {
                if (pos.x == 0) {
                    p = pos + Node(i, 1);
                    p2 = pos + Node(i, 0);
                    m[0] = norm(NodeD(image.at(p, 0), image.at(p, 1))) - norm(NodeD(image.at(p2, 0), image.at(p2, 1)));
                }
                else {
                    p = pos + Node(i, 1);
                    p2 = pos + Node(i, -1);
                    m[0] = 0.5 * (norm(NodeD(image.at(p, 0), image.at(p, 1))) - norm(NodeD(image.at(p2, 0), image.at(p2, 1))));
                }
                if (pos.x == image.cols - 2) {
                    p = pos + Node(i, 1);
                    p2 = pos + Node(i, 0);
                    m[1] = norm(NodeD(image.at(p, 0), image.at(p, 1))) - norm(NodeD(image.at(p2, 0), image.at(p2, 1)));
                }
                else {
                    p = pos + Node(i, 2);
                    p2 = pos + Node(i, 0);
                    m[1] = 0.5 * (norm(NodeD(image.at(p, 0), image.at(p, 1))) - norm(NodeD(image.at(p2, 0), image.at(p2, 1))));
                }
                p = pos + Node(i, 0);
                v0 = norm(NodeD(image.at(p, 0), image.at(p, 1)));
                p = pos + Node(i, 1);
                v1 = norm(NodeD(image.at(p, 0), image.at(p, 1)));
                v[i - init] = (2.0 * t_3 - 3.0 * t_2 + 1.0) * v0 + (t_3 - 2.0 * t_2 + t) * m[0] +
                       (-2.0 * t_3 + 3.0 * t_2) * v1 + (t_3 - t_2) * m[1];
            }
            t = epsilon.y; t_2 = t*t; t_3 = t_2*t;
            m[0] = (init == 0) ? v[1] - v[0] : (v[2] - v[0]) / 2.0;
            m[1] = (end == 1) ? v[2] - v[1] : (v[end - init] - v[end - init - 2]) / 2.0;
            v0 = v[-init];
            v1 = v[1 - init];
            value = 0.5 * (value + (2.0 * t_3 - 3.0 * t_2 + 1.0) * v0 + (t_3 - 2.0 * t_2 + t) * m[0] +
                    (-2.0 * t_3 + 3.0 * t_2) * v1 + (t_3 - t_2) * m[1]);
        }
        else if (interp == I_PCHIP) {
            double m[2], v[4], v0, v1, dv0, dv1;
            double t = epsilon.y, t_2 = t*t, t_3 = t_2*t;
            int init = (pos.x == 0) ? 0 : -1;
            int end = (pos.x == image.cols - 2) ? 1 : 2;
            for (int i=init; i<=end; i++) {
                if (pos.y == 0) {
                    p = pos + Node(1, i);
                    p2 = pos + Node(0, i);
                    dv0 = norm(NodeD(image.at(p, 0), image.at(p, 1))) - norm(NodeD(image.at(p2, 0), image.at(p2, 1)));
                    p = pos + Node(2, i);
                    p2 = pos + Node(1, i);
                    dv1 = norm(NodeD(image.at(p, 0), image.at(p, 1))) - norm(NodeD(image.at(p2, 0), image.at(p2, 1)));
                }
                else {
                    p = pos + Node(0, i);
                    p2 = pos + Node(-1, i);
                    dv0 = norm(NodeD(image.at(p, 0), image.at(p, 1))) - norm(NodeD(image.at(p2, 0), image.at(p2, 1)));
                    p = pos + Node(1, i);
                    p2 = pos + Node(0, i);
                    dv1 = norm(NodeD(image.at(p, 0), image.at(p, 1))) - norm(NodeD(image.at(p2, 0), image.at(p2, 1)));
                }
                m[0] = checkSlope(dv0, dv1, (pos.y == 0));
                if (pos.y == image.rows - 2) {
                    p = pos + Node(0, i);
                    p2 = pos + Node(-1, i);
                    dv0 = norm(NodeD(image.at(p, 0), image.at(p, 1))) - norm(NodeD(image.at(p2, 0), image.at(p2, 1)));
                    p = pos + Node(1, i);
                    p2 = pos + Node(0, i);
                    dv1 = norm(NodeD(image.at(p, 0), image.at(p, 1))) - norm(NodeD(image.at(p2, 0), image.at(p2, 1)));
                }
                else {
                    p = pos + Node(1, i);
                    p2 = pos + Node(0, i);
                    dv0 = norm(NodeD(image.at(p, 0), image.at(p, 1))) - norm(NodeD(image.at(p2, 0), image.at(p2, 1)));
                    p = pos + Node(2, i);
                    p2 = pos + Node(1, i);
                    dv1 = norm(NodeD(image.at(p, 0), image.at(p, 1))) - norm(NodeD(image.at(p2, 0), image.at(p2, 1)));
                }
                m[1] = checkSlope(dv1, dv0, (pos.y == image.rows - 2));
                p = pos + Node(0, i);
                v0 = norm(NodeD(image.at(p, 0), image.at(p, 1)));
                p = pos + Node(1, i);
                v1 = norm(NodeD(image.at(p, 0), image.at(p, 1)));
                v[i - init] = (2.0 * t_3 - 3.0 * t_2 + 1.0) * v0 + (t_3 - 2.0 * t_2 + t) * m[0] +
                       (-2.0 * t_3 + 3.0 * t_2) * v1 + (t_3 - t_2) * m[1];
            }
            t = epsilon.x; t_2 = t*t; t_3 = t_2*t;
            dv0 = v[1] - v[0];
            dv1 = v[2] - v[1];
            m[0] = checkSlope(dv0, dv1, init == 0);
            dv0 = v[end - init - 1] - v[end - init - 2];
            dv1 = v[end - init] - v[end - init - 1];
            m[1] = checkSlope(dv1, dv0, end == 1);
            v0 = v[-init];
            v1 = v[1 - init];
            value = (2.0 * t_3 - 3.0 * t_2 + 1.0) * v0 + (t_3 - 2.0 * t_2 + t) * m[0] +
                    (-2.0 * t_3 + 3.0 * t_2) * v1 + (t_3 - t_2) * m[1];

            t = epsilon.x; t_2 = t*t; t_3 = t_2*t;

            init = (pos.y == 0) ? 0 : -1;
            end = (pos.y == image.rows - 2) ? 1 : 2;
            for (int i=init; i<=end; i++) {
                if (pos.x == 0) {
                    p = pos + Node(i, 1);
                    p2 = pos + Node(i, 0);
                    dv0 = norm(NodeD(image.at(p, 0), image.at(p, 1))) - norm(NodeD(image.at(p2, 0), image.at(p2, 1)));
                    p = pos + Node(i, 2);
                    p2 = pos + Node(i, 1);
                    dv1 = norm(NodeD(image.at(p, 0), image.at(p, 1))) - norm(NodeD(image.at(p2, 0), image.at(p2, 1)));
                }
                else {
                    p = pos + Node(i, 0);
                    p2 = pos + Node(i, -1);
                    dv0 = norm(NodeD(image.at(p, 0), image.at(p, 1))) - norm(NodeD(image.at(p2, 0), image.at(p2, 1)));
                    p = pos + Node(i, 1);
                    p2 = pos + Node(i, 0);
                    dv1 = norm(NodeD(image.at(p, 0), image.at(p, 1))) - norm(NodeD(image.at(p2, 0), image.at(p2, 1)));
                }
                m[0] = checkSlope(dv0, dv1, pos.y == 0);
                if (pos.x == image.cols - 2) {
                    p = pos + Node(i, 0);
                    p2 = pos + Node(i, -1);
                    dv0 = norm(NodeD(image.at(p, 0), image.at(p, 1))) - norm(NodeD(image.at(p2, 0), image.at(p2, 1)));
                    p = pos + Node(i, 1);
                    p2 = pos + Node(i, 0);
                    dv1 = norm(NodeD(image.at(p, 0), image.at(p, 1))) - norm(NodeD(image.at(p2, 0), image.at(p2, 1)));
                }
                else {
                    p = pos + Node(i, 1);
                    p2 = pos + Node(i, 0);
                    dv0 = norm(NodeD(image.at(p, 0), image.at(p, 1))) - norm(NodeD(image.at(p2, 0), image.at(p2, 1)));
                    p = pos + Node(i, 2);
                    p2 = pos + Node(i, 1);
                    dv1 = norm(NodeD(image.at(p, 0), image.at(p, 1))) - norm(NodeD(image.at(p2, 0), image.at(p2, 1)));
                }
                m[1] = checkSlope(dv1, dv0, pos.x == image.cols - 2);
                p = pos + Node(i, 0);
                v0 = norm(NodeD(image.at(p, 0), image.at(p, 1)));
                p = pos + Node(i, 1);
                v1 = norm(NodeD(image.at(p, 0), image.at(p, 1)));

                v[i - init] = (2.0 * t_3 - 3.0 * t_2 + 1.0) * v0 + (t_3 - 2.0 * t_2 + t) * m[0] +
                       (-2.0 * t_3 + 3.0 * t_2) * v1 + (t_3 - t_2) * m[1];
            }
            t = epsilon.y; t_2 = t*t; t_3 = t_2*t;
            dv0 = v[1] - v[0];
            dv1 = v[2] - v[1];
            m[0] = checkSlope(dv0, dv1, init == 0);
            dv0 = v[end - init - 1] - v[end - init - 2];
            dv1 = v[end - init] - v[end - init - 1];
            m[1] = checkSlope(dv1, dv0, end == 1);
            v0 = v[-init];
            v1 = v[1 - init];
            value = 0.5 * (value + (2.0 * t_3 - 3.0 * t_2 + 1.0) * v0 + (t_3 - 2.0 * t_2 + t) * m[0] +
                    (-2.0 * t_3 + 3.0 * t_2) * v1 + (t_3 - t_2) * m[1]);
        }
        else if (interp == I_SPLINE) {
            double m[8], v[4], v0, v1, y[4];
            double t = epsilon.y;
            int init = (pos.x == 0) ? 0 : -1;
            int end = (pos.x == image.cols - 2) ? 1 : 2;
            int init_2 = (pos.y == 0) ? 0 : -1;
            int end_2 = (pos.y == image.rows - 2) ? 1 : 2;

            for (int p1=init; p1<=end; p1++) {
                for (int p2=init_2; p2<=end_2; p2++) {
                    p = pos + Node(p2, p1);
                    y[p2 - init_2] = norm(NodeD(image.at(p, 0), image.at(p, 1)));
                }
                get_spline_coeffs(y, m, end_2 - init_2);
                p = pos + Node(0, p1);
                v0 = norm(NodeD(image.at(p, 0), image.at(p, 1)));
                p = pos + Node(1, p1);
                v1 = norm(NodeD(image.at(p, 0), image.at(p, 1)));
                v[p1 - init] = (1.0 - t) * v0 + t * v1 + t * (1.0 - t) * ((1.0 - t) * m[-2*init_2] + t * m[-2*init_2 + 1]);
            }
            t = epsilon.x;
            get_spline_coeffs(v, m, end - init);
            v0 = v[-init];
            v1 = v[1 - init];
            value = (1.0 - t) * v0 + t * v1 + t * (1.0 - t) * ((1.0 - t) * m[-2*init] + t * m[-2*init + 1]);

            t = epsilon.x;
            init = (pos.y == 0) ? 0 : -1;
            end = (pos.y == image.rows - 2) ? 1 : 2;
            init_2 = (pos.x == 0) ? 0 : -1;
            end_2 = (pos.x == image.cols - 2) ? 1 : 2;
            for (int p1=init; p1<=end; p1++) {
                for (int p2=init_2; p2<=end_2; p2++) {
                    p = pos + Node(p1, p2);
                    y[p2 - init_2] = norm(NodeD(image.at(p, 0), image.at(p, 1)));
                }
                get_spline_coeffs(y, m, end_2 - init_2);
                p = pos + Node(p1, 0);
                v0 = norm(NodeD(image.at(p, 0), image.at(p, 1)));
                p = pos + Node(p1, 1);
                v1 = norm(NodeD(image.at(p, 0), image.at(p, 1)));
                v[p1 - init] = (1.0 - t) * v0 + t * v1 + t * (1.0 - t) * ((1.0 - t) * m[-2*init_2] + t * m[-2*init_2 + 1]);
            }

            t = epsilon.y;
            get_spline_coeffs(v, m, end - init);
            v0 = v[-init];
            v1 = v[1 - init];
            value = 0.5 * (value + (1.0 - t) * v0 + t * v1 + t * (1.0 - t) * ((1.0 - t) * m[-2*init] + t * m[-2*init + 1]));

        }

        if (((value < v0 && value < v1 && value < v2 && value < v3) || (value > v0 && value > v1 && value > v2 && value > v3)) &&
            (interp != I_LINEAR) && (interp != I_PCHIP)) {
            //value = (1.0 - epsilon.y)*(1.0 - epsilon.x)*v0 + epsilon.y*(1.0 - epsilon.x)*v1 +
            //        (1.0 - epsilon.y)*epsilon.x*v2 + epsilon.y*epsilon.x*v3;
            value = GetFValue(image, posD, I_PCHIP);
        }

        return value;
    }



    double GetGradientFromGrid(Grid &image, NodeD &posD, int interp) {
        Node pos(floor(posD.y), floor(posD.x));
        auto diff = posD - pos;
        double epsilon = fmax(diff.x, diff.y);
        Node dirs(sgn(diff.y), sgn(diff.x));
        Node aux = pos - dirs;
        int init = (image.contains(aux)) ? -1 : 0, end = 4 + init;
        aux = pos + (end - 1) * dirs;
        if (!image.contains(aux)) {
            aux = pos + (init - 1) * dirs;
            if (image.contains(aux)) {
                init -= 1;
                end -= 1;
            }
        }
        double f[4], m[12], value, step = 1.0 / (double) (end - init - 1);
        for (int i=init; i < end; i++) {
            aux = pos + i * dirs;
            f[i - init] = (image.contains(aux)) ? norm(NodeD(image.at(aux, 0), image.at(aux, 1))) : MAX_VAL;
        }
        setCoeffs2D(f, m, interp, end - init - 1);
        switch (interp) {
            case I_LINEAR:
                value = (1.0 - epsilon)*f[-init] + epsilon*f[1-init];
                break;
            case I_QUADATRIC:
                value = epsilon*epsilon*m[-3*init + 2] + epsilon*m[-3*init + 1] + m[-3*init];
                break;
            case I_SPLINE:
                value = (1.0 - epsilon) * f[-init] + epsilon * f[1-init] +
                        epsilon * (1.0 - epsilon) *
                        (m[-2*init] * (1.0 - epsilon) + m[-2*init + 1] * epsilon);
                //value = y0 + epsilon*(-wave.m[2*side + 1]/6.0 - wave.m[2*side]/3.0 + y1 - y0 +
                //        epsilon*(wave.m[2*side]/2.0 + epsilon*(wave.m[2*side + 1] - wave.m[2*side])/6.0));
                break;
            default: //I_HERMITE and I_PCHIP
                double t_2 = epsilon * epsilon;
                double t_3 = epsilon * t_2;
                value = (2.0 * t_3 - 3.0 * t_2 + 1.0) * f[-init] + step * (t_3 - 2.0 * t_2 + epsilon) * m[-init] +
                       (-2.0 * t_3 + 3.0 * t_2) * f[1-init] + step * (t_3 - t_2) * m[1-init];
        }

        /*if (epsilon > 0 && epsilon < 1) {
            std::cout << "posD: " << posD.y << ", " << posD.x << ", f0: " << f[-init] << ", f1: " << f[1-init] << ", f-: " << f[0] << ", f+: " << f[3]  << ", eps: " << epsilon << ", f: " << value << std::endl;
        }*/

        if ((value < f[-init] && value < f[1-init]) || (value > f[-init] && value > f[1-init])) {
            value = (1.0 - epsilon)*f[-init] + epsilon*f[1-init];
        }
        return value;
    }

    double GetIntegralValue(Grid &image, NodeD &posD, NodeD &neighD, double fposD, double fneighD,
                            int interp, int inner_nodes=1) {
        double values = fposD + fneighD, step = 1.0 / (inner_nodes + 1.0);
        NodeD p;
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

