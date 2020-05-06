#include "../TYPES/WMMStructs.h"
#include <math.h>
#include <stdlib.h>
#include <vector>
#include <map>
#include <iostream>


namespace wmm {

    wmm::Node3D GetInterpCpolar(wmm::AugWmm &wave, wmm::Node3D &f0, wmm::Node3D &f1, double epsilon, int interp, int side) {
    
        wmm::Node3D c_polar;
        
        if (wave.dir < 0) {
            return f0;
        }
        
        switch (interp) {
            case wmm::I_LINEAR:
                c_polar = (1.0 - epsilon)*f0 + epsilon*f1;
                break;
            case wmm::I_QUADATRIC:
                c_polar = epsilon*epsilon*wave.fm[2*side + 2] + epsilon*wave.fm[2*side + 1] + wave.fm[0];
                break;
            case wmm::I_SPLINE:
                c_polar = (1.0 - epsilon) * f0 + epsilon * f1 + 
                        epsilon * (1.0 - epsilon) * 
                        ((1.0 - epsilon) * wave.fm[2*side] + epsilon * wave.fm[2*side + 1]);
                //c_polar = f0 + epsilon*(-1.0/6.0*wave.fm[2*side + 1] - 2.0/3.0*wave.fm[2*side] + f1 - f0 + epsilon*(
                //        0.5*wave.fm[2*side] + epsilon/6.0*(wave.fm[2*side + 1] - wave.fm[2*side])));
                break;
            default: //I_HERMITE and I_PCHIP
                double t_2 = epsilon * epsilon;
                double t_3 = epsilon * t_2;
                c_polar = (2.0 * t_3 - 3.0 * t_2 + 1.0) * f0 + (t_3 - 2.0 * t_2 + epsilon) * wave.fm[2*side] +
                     (-2.0 * t_3 + 3.0 * t_2) * f1 + (t_3 - t_2) * wave.fm[2*side + 1];
        }
        if ((c_polar.z < f0.z && c_polar.z < f1.z) || (c_polar.z > f0.z && c_polar.z > f1.z)) {
            c_polar = (1.0 - epsilon)*f0 + epsilon*f1;
        }
        
        return c_polar;
    
    }
    
    double GetInterpV(wmm::AugWmm &wave, double t_epsilon, int interp, int side, int gamma) {
        // f0, f1 and f0 must be in polar coordinates
        
        int segment = (t_epsilon >= 1.0) ? gamma - 1 : (int) floor(gamma * t_epsilon);
        double value, y0 = (segment == 0) ? wave.v[0] : wave.v[side * gamma + segment];
        double y1 = wave.v[side * gamma + segment + 1];
        double step = 1.0/ (double) gamma;
        double epsilon = (t_epsilon - segment * step) / step;
                
        if (wave.dir == -1) {
            return y0;
        }
        
        switch (interp) {
            case wmm::I_LINEAR:
                value = (1.0 - epsilon)*y0 + epsilon*y1;
                break;
            case wmm::I_QUADATRIC:
                value = epsilon*epsilon*wave.m[2*side + 2] + epsilon*wave.m[2*side + 1] + wave.m[0];
                break;
            case wmm::I_SPLINE:
                value = (1.0 - epsilon) * y0 + epsilon * y1 + 
                        epsilon * (1.0 - epsilon) * 
                        (wave.m[2*side*gamma + 2*segment] * (1.0 - epsilon) + wave.m[2*side*gamma + 2*segment + 1] * epsilon);
                //value = y0 + epsilon*(-wave.m[2*side + 1]/6.0 - wave.m[2*side]/3.0 + y1 - y0 + 
                //        epsilon*(wave.m[2*side]/2.0 + epsilon*(wave.m[2*side + 1] - wave.m[2*side])/6.0));
                break;
            default: //I_HERMITE and I_PCHIP
                double t_2 = epsilon * epsilon;
                double t_3 = epsilon * t_2;
                value = (2.0 * t_3 - 3.0 * t_2 + 1.0) * y0 + step * (t_3 - 2.0 * t_2 + epsilon) * wave.m[side*(gamma + 1) + segment] +
                       (-2.0 * t_3 + 3.0 * t_2) * y1 + step * (t_3 - t_2) * wave.m[side*(gamma + 1) + segment + 1];
        }
        
        if ((value < y0 && value < y1) || (value > y0 && value > y1)) {
            value = (1.0 - epsilon)*y0 + epsilon*y1;
        }
        return value;
    }
    
    double GetInterpValue(wmm::AugWmm &wave, wmm::NodeD &dp, wmm::NodeD &dd, wmm::NodeD &dn,
                          wmm::Node3D &f0, wmm::Node3D &f1, wmm::Node3D &fn, double epsilon, int interp, int side, int gamma) {

        // f0, f1 and f0 must be in polar coordinates
        
        double v = GetInterpV(wave, epsilon, interp, side, gamma);
        wmm::NodeD diff = dp + epsilon*dd - dn;
        wmm::NodeD a = 1.0/wmm::norm(diff) * diff, c;
        wmm::Node3D c_polar = GetInterpCpolar(wave, f0, f1, epsilon, interp, side);
        c = wmm::to_cartesian(0.5*(c_polar + fn));
        
        double value = v + wmm::norm(diff) * sqrt(1.0 + (c.x*a.x + c.y*a.y)*(c.x*a.x + c.y*a.y));
        
        return value;
    }


    double GetEpsilonGradient(wmm::NodeD &dd, wmm::NodeD &dp, wmm::NodeD &dn, wmm::NodeD &fn) {

        double epsilon = 0.0;
        double A = -dd.y, B = dd.x, C = dd.y * dp.x - dd.x * dp.y;
        double den = A * fn.y - B * fn.x;
        double t = (A * dn.x + B * dn.y + C) / den;

        wmm::NodeD x(dn.y + t * fn.x, dn.x - t * fn.y);

        if (fabs(dd.x) > 0.0 && fabs(den) > 0.0)
            epsilon = (x.x - dp.x) / dd.x;
        else if (fabs(dd.y) > 0.0 && fabs(den) > 0.0)
            epsilon = (x.y - dp.y) / dd.y;
        else if (fabs(den) == 0.0 && wmm::norm(dd) > 0.0) {
            double dist = fabs(A * dn.x + B * dn.y + C) / sqrt(A * A + B * B);
            epsilon = (wmm::norm(dn - dp) - dist) / (fabs(dd.x) + fabs(dd.y));
        }
        else
            return 0.0;

        if (epsilon < 0.0)
            epsilon = 0.0;
        else if (epsilon > 1.0)
            epsilon = 1.0;
        return epsilon;

    }

    double GetEpsilonHopfLax(wmm::NodeD &dd, wmm::NodeD &dp, wmm::NodeD &dn, double y0, double y1) {

        if (wmm::norm(dn - dp - dd) < TAU) {
            return 0.0;
        }
        else if (wmm::norm(dn - dp) < TAU) {
            return 1.0;
        }
        wmm::NodeD xy = dn - dp;
        double nxy = wmm::norm(xy);
        double nyz = wmm::norm(dd);

        double c_alpha = (xy.x * dd.x + xy.y * dd.y) / (nxy * nyz);
        double c_delta = (y1 - y0) / nyz;

        if (nyz == 0.0 || c_alpha <= c_delta || c_alpha == 1.0) {
            return 0.0;
        }

        wmm::NodeD xz = dn - dp - dd;
        double nxz = wmm::norm(xz);
        double c_beta = (xz.x * dd.x + xz.y * dd.y) / (nxz * nyz);

        if (c_delta <= c_beta) {
            return 1.0;
        }

        double s_delta = sqrt(1.0 - c_delta * c_delta);
        double dist = (c_alpha * c_delta + sqrt(1.0 - c_alpha * c_alpha) * s_delta) * nxy;
        double yzdist = sqrt(nxy * nxy - dist * dist);
        double epsilon = yzdist / (s_delta * nyz);

        return epsilon;

    }

    double GetEpsilonGoldenSearch(wmm::AugWmm &wave, wmm::NodeD &dn, wmm::NodeD &dd, wmm::NodeD &dp,
                                  wmm::Node3D &f0, wmm::Node3D &f1, wmm::Node3D &fn,
                                  double interp, int side, int gamma) {

        double a = 0.0, b = 1.0, x1 = a + (1-RESPHI)*(b - a), x2 = a + RESPHI*(b - a),
                f_x1 = MAX_VAL, f_x2 = MAX_VAL, i1, i2, res;

        wmm::NodeD xtreme = dp + dd;

        wmm::NodeD F_x1, F_x2;

        double epsilon;
        double f_a = GetInterpValue(wave, dp, dd, dn, f0, f1, fn, 0.0, interp, side, gamma);
        double f_b = GetInterpValue(wave, dp, dd, dn, f0, f1, fn, 1.0, interp, side, gamma);
        
        if (f_a < f_b) {
            res = f_a; epsilon = 0.0;
        }
        else {
            res = f_b; epsilon = 1.0;
        }

        f_x1 = GetInterpValue(wave, dp, dd, dn, f0, f1, fn, x1, interp, side, gamma);
        f_x2 = GetInterpValue(wave, dp, dd, dn, f0, f1, fn, x2, interp, side, gamma);

        while (fabs(b - a) > TAU) {
            if(f_x1 < f_x2) {
                b = x2; x2 = x1; f_x2 = f_x1; x1 = a + (1 - RESPHI)*(b - a);
                f_x1 = GetInterpValue(wave, dp, dd, dn, f0, f1, fn, x1, interp, side, gamma);
            }
            else {
                a = x1; x1 = x2; f_x1 = f_x2; x2 = a + RESPHI*(b - a);
                f_x2 = GetInterpValue(wave, dp, dd, dn, f0, f1, fn, x2, interp, side, gamma);
            }
        }

        if (f_x1 < std::min(res, f_x2))
            epsilon = x1;
        else if (f_x2 < std::min(res, f_x1))
            epsilon = x2;

        return epsilon;
        
    }


    double GetVal2D(wmm::Grid &image, wmm::Grid &u_surface, wmm::AugWmm &wave, wmm::NodeD &neigh, wmm::NodeD &fn, wmm::NodeD &h, int interp, int mode, int gamma) {

        wmm::NodeD f0(image.at(wave.p, 0), image.at(wave.p, 1));
        double y0 = wave.v[0], res1, res2, step = 1.0 / (double) gamma;
        int i;

        if (isinf(wmm::norm(f0)) || isnan(wmm::norm(f0)))
            f0 = fn;
        
        wmm::Node3D f0_polar = wmm::to_polar(f0), fn_polar = wmm::to_polar(fn);
        
        double val = wmm::MAX_VAL;
        if (wave.dir < 0) {
            wmm::NodeD diff(h.y * (neigh.y - wave.p.y), h.x * (neigh.x - wave.p.x));
            wmm::NodeD c = wmm::to_cartesian(0.5*(f0_polar + fn_polar));
            wmm::NodeD a = 1.0/wmm::norm(diff) * diff;
            val = y0 + wmm::norm(diff) * sqrt(1.0 + (c.x*a.x + c.y*a.y)*(c.x*a.x + c.y*a.y));
        }
        else {

            wmm::Node p(wave.p.y + wmm::yarray[(wave.dir + 1) % 8] - wmm::yarray[wave.dir],
                    wave.p.x + wmm::xarray[(wave.dir + 1) % 8] - wmm::xarray[wave.dir]);
            res1 = wmm::MAX_VAL;
            wmm::NodeD dp(h.y * wave.p.y, h.x * wave.p.x), dn(h.y * neigh.y, h.x * neigh.x);

            if (u_surface.contains(p)) {
                double y1 = wave.v[gamma];

                wmm::NodeD dd(h.y * (wmm::yarray[(wave.dir + 1) % 8] - wmm::yarray[wave.dir]), h.x * (wmm::xarray[(wave.dir + 1) % 8] - wmm::xarray[wave.dir]));

                wmm::NodeD f1(image.at(p, 0), image.at(p, 1)), dp_aux, dd_aux;
                if (isinf(wmm::norm(f1)) || isnan(wmm::norm(f1)))
                    f1 = fn;
                
                wmm::Node3D f1_polar = wmm::to_polar(f1);

                double epsilon, epsilon_aux, res_aux;
                switch (mode) {
                    case wmm::M_GRADIENT:
                        //epsilon = GetEpsilonGradient(dd, dp, dn, fn);
                        dd_aux = step * dd;
                        res1 = wmm::MAX_VAL;
                        epsilon = 0.0;
                        for (i=0; i<gamma; i++) {
                            dp_aux = dp + i*step*dd;
                            epsilon_aux = (i + GetEpsilonGradient(dd_aux, dp_aux, dn, fn))*step;
                            res_aux = GetInterpValue(wave, dp, dd, dn, f0_polar, f1_polar, fn_polar, epsilon_aux, interp, wmm::S_RIGHT, gamma);
                            if (res_aux < res1) {
                                res1 = res_aux;
                                epsilon = epsilon_aux;
                            }
                        }
                        break;
                    case wmm::M_HOPFLAX:
                        dd_aux = step * dd;
                        res1 = wmm::MAX_VAL;
                        epsilon = 0.0;
                        y0 = wave.v[0];
                        for (i=1; i<=gamma; i++) {
                            y1 = wave.v[i];
                            dp_aux = dp + (i-1)*step*dd;
                            epsilon_aux = (i - 1.0 + GetEpsilonHopfLax(dd_aux, dp_aux, dn, y0, y1))*step;
                            res_aux = GetInterpValue(wave, dp, dd, dn, f0_polar, f1_polar, fn_polar, epsilon_aux, interp, wmm::S_RIGHT, gamma);
                            if (res_aux < res1) {
                                res1 = res_aux;
                                epsilon = epsilon_aux;
                            }
                            y0 = y1;
                        }
                        break;
                    default: //M_GOLDENSEARCH
                        epsilon = GetEpsilonGoldenSearch(wave, dn, dd, dp, f0_polar, f1_polar, fn_polar, interp, wmm::S_RIGHT, gamma);
                }

                res1 = GetInterpValue(wave, dp, dd, dn, f0_polar, f1_polar, fn_polar, epsilon, interp, wmm::S_RIGHT, gamma);
            }

            p = wmm::Node(wave.p.y + wmm::yarray[(wave.dir + 7) % 8] - wmm::yarray[wave.dir],
                    wave.p.x + wmm::xarray[(wave.dir + 7) % 8] - wmm::xarray[wave.dir]);
            res2 = wmm::MAX_VAL;

            if (u_surface.contains(p)) {
                double y1 = wave.v[2*gamma];

                wmm::NodeD dd(h.y * (wmm::yarray[(wave.dir + 7) % 8] - wmm::yarray[wave.dir]), h.x * (wmm::xarray[(wave.dir + 7) % 8] - wmm::xarray[wave.dir]));

                wmm::NodeD f1(image.at(p, 0), image.at(p, 1)), dp_aux, dd_aux;
                if (isinf(wmm::norm(f1)) || isnan(wmm::norm(f1)))
                    f1 = fn;
                
                wmm::Node3D f1_polar = wmm::to_polar(f1);

                double epsilon, epsilon_aux, res_aux;
                switch (mode) {
                    case wmm::M_GRADIENT:
                        //epsilon = GetEpsilonGradient(dd, dp, dn, fn);
                        dd_aux = step * dd;
                        res2 = wmm::MAX_VAL;
                        epsilon = 0.0;
                        for (i=0; i<gamma; i++) {
                            dp_aux = dp + i*step*dd;
                            epsilon_aux = (i + GetEpsilonGradient(dd_aux, dp_aux, dn, fn))*step;
                            res_aux = GetInterpValue(wave, dp, dd, dn, f0_polar, f1_polar, fn_polar, epsilon_aux, interp, wmm::S_LEFT, gamma);
                            if (res_aux < res2) {
                                res2 = res_aux;
                                epsilon = epsilon_aux;
                            }
                        }
                        break;
                    case wmm::M_HOPFLAX:
                        dd_aux = step * dd;
                        res2 = wmm::MAX_VAL;
                        epsilon = 0.0;
                        y0 = wave.v[0];
                        for (i=1; i<=gamma; i++) {
                            y1 = wave.v[gamma + i];
                            dp_aux = dp + (i-1)*step*dd;
                            epsilon_aux = (i - 1.0 + GetEpsilonHopfLax(dd_aux, dp_aux, dn, y0, y1))*step;
                            res_aux = GetInterpValue(wave, dp, dd, dn, f0_polar, f1_polar, fn_polar, epsilon_aux, interp, wmm::S_LEFT, gamma);
                            if (res_aux < res2) {
                                res2 = res_aux;
                                epsilon = epsilon_aux;
                            }
                            y0 = y1;
                        }
                        break;
                    default: //M_GOLDENSEARCH
                        epsilon = GetEpsilonGoldenSearch(wave, dn, dd, dp, f0_polar, f1_polar, fn_polar, interp, wmm::S_LEFT, gamma);
                }

                res2 = GetInterpValue(wave, dp, dd, dn, f0_polar, f1_polar, fn_polar, epsilon, interp, wmm::S_LEFT, gamma);
                

            }

            val = std::min(res1, res2);

        }
        
        return val;

    }
    
    
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

    void setSimpleCoeffs2D(double *y, double *m, int pos, int interp) {

        switch (interp) {
            case wmm::I_LINEAR:
                break;
            case wmm::I_QUADATRIC:
                m[0] = y[pos];
                if (pos%2 == 0) {
                    m[2] = (y[(pos+2)%8] + y[pos] - 2.0*y[(pos+1)%8])/2.0;
                    m[1] = y[(pos+1)%8] - y[pos] - m[2];
                    m[4] = (y[(pos+6)%8] + y[pos] - 2.0*y[(pos+7)%8])/2.0;
                    m[3] = y[(pos+7)%8] - y[pos] - m[4];
                }
                else {
                    m[2] = m[4] = (y[(pos+1)%8] + y[(pos+7)%8] - 2.0*y[pos])/2.0;
                    m[1] = y[pos] - y[(pos+7)%8] + m[2];
                    m[3] = y[pos] - y[(pos+1)%8] + m[4];
                }
                break;
            case wmm::I_SPLINE:
                if (pos%2 == 0) {
                    m[0] = m[2] = 0.0;
                    m[1] = 6.0*(y[pos] - 2.0*y[(pos+1)%8] + y[(pos+2)%8])/4.0;
                    m[3] = 6.0*(y[pos] - 2.0*y[(pos+7)%8] + y[(pos+6)%8])/4.0;
                }
                else {
                    m[0] = m[2] = 6.0*(y[(pos+7)%8] - 2.0*y[pos] + y[(pos+1)%8])/4.0;
                    m[1] = m[3] = 0.0;
                }
                break;
            case wmm::I_HERMITE:
                if (pos%2 == 0) {
                    m[0] = y[(pos + 1) % 8] - y[pos];
                    m[1] = 1.0 / 2.0 * (y[(pos + 2) % 8] - y[pos]);
                    m[2] = y[(pos + 7) % 8] - y[pos];
                    m[3] = 1.0 / 2.0 * (y[(pos + 6) % 8] - y[pos]);
                }
                else {
                    m[0] = 1.0 / 2.0 * (y[(pos + 1) % 8] - y[(pos + 7) % 8]);
                    m[1] = y[(pos + 1) % 8] - y[pos];
                    m[2] = -m[0];
                    m[3] = y[(pos + 7) % 8] - y[pos];
                }
                break;
            default: //I_PCHIP
                if (pos%2 == 0) {
                    double d0 = y[(pos + 1) % 8] - y[pos], d1 = y[(pos + 2) % 8] - y[(pos + 1) % 8];
                    m[0] = (3.0 * d0 - d1) / 2.0;
                    if ((m[0] * d0) <= 0.0)
                        m[0] = 0.0;
                    else if (((d1 * d0) <= 0.0) && (fabs(m[0]) > fabs(3.0 * d0)))
                        m[0] = 3.0 * d0;
                    m[1] = (d0 * d1 <= 0.0) ? 0.0 : 2.0 * d0 * d1 / (d0 + d1);

                    d0 = y[(pos + 7) % 8] - y[pos]; d1 = y[(pos + 6) % 8] - y[(pos + 7) % 8];
                    m[2] = (3.0 * d0 - d1) / 2.0;
                    if ((m[2] * d0) <= 0.0)
                        m[2] = 0.0;
                    else if (((d1 * d0) <= 0.0) && (fabs(m[2]) > fabs(3.0 * d0)))
                        m[2] = 3.0 * d0;
                    m[3] = (d0 * d1 <= 0.0) ? 0.0 : 2.0 * d0 * d1 / (d0 + d1);
                }
                else {
                    double d0 = y[pos] - y[(pos + 7) % 8], d1 = y[(pos + 1) % 8] - y[pos];
                    m[0] = (d0 * d1 <= 0.0) ? 0.0 : 2.0 * d0 * d1 / (d0 + d1);
                    m[1] = (3.0 * d1 - d0) / 2.0;
                    if ((m[1] * d1) <= 0.0)
                        m[1] = 0.0;
                    else if (((d0 * d1) <= 0.0) && (fabs(m[1]) > fabs(3.0 * d1)))
                        m[1] = 3.0 * d1;

                    d0 = y[pos] - y[(pos + 1) % 8]; d1 = y[(pos + 7) % 8] - y[pos];
                    m[2] = (d0 * d1 <= 0.0) ? 0.0 : 2.0 * d0 * d1 / (d0 + d1);
                    m[3] = (3.0 * d1 - d0) / 2.0;
                    if ((m[3] * d1) <= 0.0)
                        m[3] = 0.0;
                    else if (((d0 * d1) <= 0.0) && (fabs(m[3]) > fabs(3.0 * d1)))
                        m[3] = 3.0 * d1;
                }
        }

    }
    
    

    void setCoeffs2D(double *y, double *m, int interp, int gamma) {
        switch (interp) {
            case wmm::I_LINEAR:
                break;
            /*case wmm::I_QUADATRIC:
                m[0] = y[pos];
                if (pos%2 == 0) {
                    m[2] = (y[(pos+2)%8] + y[pos] - 2.0*y[(pos+1)%8])/2.0;
                    m[1] = y[(pos+1)%8] - y[pos] - m[2];
                    m[4] = (y[(pos+6)%8] + y[pos] - 2.0*y[(pos+7)%8])/2.0;
                    m[3] = y[(pos+7)%8] - y[pos] - m[4];
                }
                else {
                    m[2] = m[4] = (y[(pos+1)%8] + y[(pos+7)%8] - 2.0*y[pos])/2.0;
                    m[1] = y[pos] - y[(pos+7)%8] + m[2];
                    m[3] = y[pos] - y[(pos+1)%8] + m[4];
                }
                break;*/
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
                double d0 = gamma*(y[1] - y[0]), d1 = gamma*(y[2] - y[1]);
                m[0] = (3.0*d0 - d1)/2.0;
                if (fabs(m[0]) > fabs(3.0*d0))
                    m[0] = 3.0*d0;
                d0 = gamma*(y[gamma] - y[gamma - 1]);
                d1 = gamma*(y[gamma - 1] - y[gamma - 2]);
                m[gamma] = (3.0*d0 - d1)/2.0;
                if (fabs(m[gamma]) > fabs(3.0*d0))
                    m[gamma] = 3.0*d0;
                for (int i=1; i<gamma; i++) {
                    m[i] = 0.5*gamma*(y[i + 1] - y[i - 1]);
                }                
                break;
            }
            default: {//I_PCHIP
                double d0 = gamma*(y[1] - y[0]);
                double d1 = gamma*(y[2] - y[1]);
                m[0] = (3.0 * d0 - d1) / 2.0;
                if ((m[0] * d0) <= 0.0)
                    m[0] = 0.0;
                else if (((d1 * d0) <= 0.0) && (fabs(m[0]) > fabs(3.0 * d0)))
                    m[0] = 3.0 * d0;
                m[1] = (d0 * d1 <= 0.0) ? 0.0 : 2.0 * d0 * d1 / (d0 + d1);
                for (int i=2; i<gamma; i++) {
                    d0 = d1;
                    d1 = gamma*(y[i + 1] - y[i]);
                    m[i] = (d0 * d1 <= 0.0) ? 0.0 : 2.0 * d0 * d1 / (d0 + d1);
                }
                d0 = gamma*(y[gamma] - y[gamma - 1]);
                d1 = gamma*(y[gamma - 1] - y[gamma - 2]);
                m[gamma] = (3.0 * d0 - d1) / 2.0;
                if ((m[gamma] * d0) <= 0.0)
                    m[gamma] = 0.0;
                else if (((d1 * d0) <= 0.0) && (fabs(m[gamma]) > fabs(3.0 * d0)))
                    m[gamma] = 3.0 * d0;
            }
        }

    }


    wmm::Grid WmmAugAniSurface2D(wmm::Grid &image, std::vector<wmm::Node> &initials, wmm::NodeD &h, int interp, int mode, int N, int M, int gamma) {

        bool isnewpos[8];
        double valcenter[8];
        double imcenter_module[8];
        double imcenter_cos[8];
        double imcenter_sin[8], minimo, step = 1.0/gamma;
        double vinterp[gamma + 1];
        
        int degree = (N == 0) ? 1 : N, range;
                
        double module_values[degree], cos_values[degree], sin_values[degree];

        wmm::Grid u_surface = wmm::Grid(wmm::MAX_VAL, image.rows, image.cols, 1);
        wmm::Grid_<unsigned char> state = wmm::Grid_<unsigned char>(image.rows, image.cols);

        std::multimap<double, wmm::AugWmm > trial_set;
        std::map<int, std::multimap<double, wmm::AugWmm >::iterator> mapa_trial;

        std::multimap<double, wmm::AugWmm >::iterator trial_set_it;
        std::map<int, std::multimap<double, wmm::AugWmm >::iterator>::iterator mapa_trial_it;
        std::pair<double, wmm::AugWmm > pr_trial;
        std::pair<int, std::multimap<double, wmm::AugWmm >::iterator> pr_mapa;

        int key, i, j;
        wmm::AugWmm winner, new_w;
        wmm::Node neigh;
        wmm::NodeD neighD, fn, d1, neigh_aux, c;
        wmm::Node3D polar0, polar1;

        // Initialization
        for (i = 0; i < (int) initials.size(); i++) {
            key = initials[i].y * u_surface.cols + initials[i].x;
            if (mapa_trial.find(key) == mapa_trial.end() && u_surface.contains(initials[i])) {
                winner = wmm::AugWmm(N, M, gamma);
                u_surface.at(initials[i]) = 0.0;
                winner.dir = -1;
                winner.v[0] = 0.0;
                winner.p = initials[i];
                state.at(initials[i]) = wmm::P_TRIAL;
                pr_trial = std::pair<double, wmm::AugWmm >(0.0, winner);
                trial_set_it = trial_set.insert(pr_trial);
                pr_mapa = std::pair<int, std::multimap<double, wmm::AugWmm >::iterator>(key, trial_set_it);
                mapa_trial.insert(pr_mapa);
            }
        }
        
        while (!trial_set.empty()) {

            trial_set_it = trial_set.begin();
            winner = trial_set_it->second;
            key = winner.p.y * u_surface.cols + winner.p.x;
            mapa_trial_it = mapa_trial.find(key);
            //std::cout << "ganador (" << father.p.y << ", " << father.p.x << ", " << father.child_dir << ") key = " << key << std::endl;
            //std::cout << "buscamos (" << winner.p.y << ", " << winner.p.x << ") key = " << key << std::endl;
            if (mapa_trial_it == mapa_trial.end()) {
                printf("ERROR: bad map alloc");
                return u_surface;
            }

            if (mapa_trial_it->second != trial_set_it) {
                printf("ERROR: bad trial/map alloc");
                return u_surface;
            }
            
            trial_set.erase(trial_set_it);
            mapa_trial.erase(mapa_trial_it);

            state.at(winner.p) = wmm::P_ALIVE;
            
            // Neighbour temptative value computation
            for (i = 0; i < 8; i++) {
                neigh = winner.p + wmm::Node(wmm::yarray[i], wmm::xarray[i]);
                isnewpos[i] = false;
                if (u_surface.contains(neigh)) {
                    valcenter[i] = u_surface.at(neigh);
                    fn = wmm::NodeD(image.at(neigh, 0), image.at(neigh, 1));
                }
                else {
                    valcenter[i] = winner.v[0];
                    fn = wmm::NodeD(image.at(winner.p, 0), image.at(winner.p, 1));
                }
                wmm::Node3D polar = wmm::to_polar(fn);
                imcenter_module[i] = polar.z;
                imcenter_cos[i] = polar.x;
                imcenter_sin[i] = polar.y;
                if (u_surface.contains(neigh) && state.at(neigh) != wmm::P_ALIVE) {
                    neighD = wmm::NodeD(winner.p.y + wmm::yarray[i], winner.p.x + wmm::xarray[i]);
                    
                    double val_neigh = GetVal2D(image, u_surface, winner, neighD, fn, h, interp, mode, gamma);
                    if (val_neigh < valcenter[i]) {
                        valcenter[i] = val_neigh;
                        isnewpos[i] = true;
                    }
                }
            }

            // Update
            for (i = 0; i < 8; i++) {
                if (isnewpos[i]) {
                    neigh = winner.p + wmm::Node(wmm::yarray[i], wmm::xarray[i]);
                    key = neigh.y * u_surface.cols + neigh.x;
                    if (state.at(neigh) == wmm::P_TRIAL) {
                        mapa_trial_it = mapa_trial.find(key);
                        del(&(mapa_trial_it->second->second));
                        trial_set.erase(mapa_trial_it->second);
                        mapa_trial.erase(mapa_trial_it);
                    }
                    else {
                        state.at(neigh) = wmm::P_TRIAL;
                    }

                    new_w = wmm::AugWmm(N, M, gamma);
                    new_w.dir = i;
                    new_w.p = neigh;
                    
                    new_w.v[0] = valcenter[i];
                    new_w.v[gamma] = valcenter[(i+1)%8];
                    new_w.v[2*gamma] = valcenter[(i+7)%8];

                    // falta asignar vs y escoger el minimo
                    if (N > 0) {
                        setSimpleCoeffs2D(imcenter_module, module_values, i, interp);
                        setSimpleCoeffs2D(imcenter_cos, cos_values, i, interp);
                        setSimpleCoeffs2D(imcenter_sin, sin_values, i, interp);
                        for (j = 0; j < degree; j++) {
                            new_w.fm[j].z = module_values[j];
                            new_w.fm[j].x = cos_values[j];
                            new_w.fm[j].y = sin_values[j];
                        }
                    }

                    // COMPUTING WAVEFRONT 1st SEGMENT RIGHT
                    neighD = wmm::NodeD(new_w.p.y, new_w.p.x);
                    polar0 = wmm::Node3D(imcenter_cos[i], imcenter_sin[i], imcenter_module[i]);
                    polar1 = wmm::Node3D(imcenter_cos[(i+1)%8], imcenter_sin[(i+1)%8], imcenter_module[(i+1)%8]);
                    vinterp[0] = valcenter[i];
                    vinterp[gamma] = valcenter[(i+1)%8];
                    d1 = wmm::NodeD(wmm::yarray[(i+1)%8], wmm::xarray[(i+1)%8]) - wmm::NodeD(wmm::yarray[i], wmm::xarray[i]);

                    for (j = 1; j < gamma; j++) {
                        c = wmm::to_cartesian(GetInterpCpolar(new_w, polar0, polar1, step*j, interp, wmm::S_RIGHT));
                        neigh_aux = neighD + step*j*d1;
                        if (u_surface.contains(neigh_aux))
                            vinterp[j] = GetVal2D(image, u_surface, winner, neigh_aux, c, h, interp, mode, gamma);
                        else
                            vinterp[j] = valcenter[i];
                        new_w.v[j] = vinterp[j];
                    }
                    setCoeffs2D(vinterp, new_w.m, interp, gamma);


                    // COMPUTING WAVEFRONT 1st SEGMENT LEFT
                    neighD = wmm::NodeD(new_w.p.y, new_w.p.x);
                    polar0 = wmm::Node3D(imcenter_cos[i], imcenter_sin[i], imcenter_module[i]);
                    polar1 = wmm::Node3D(imcenter_cos[(i+7)%8], imcenter_sin[(i+7)%8], imcenter_module[(i+7)%8]);
                    vinterp[0] = valcenter[i];
                    vinterp[gamma] = valcenter[(i+7)%8];
                    d1 = wmm::NodeD(wmm::yarray[(i+7)%8], wmm::xarray[(i+7)%8]) - wmm::NodeD(wmm::yarray[i], wmm::xarray[i]);

                    for (j = 1; j < gamma; j++) {
                        c = wmm::to_cartesian(GetInterpCpolar(new_w, polar0, polar1, step*j, interp, wmm::S_LEFT));
                        neigh_aux = neighD + step*j*d1;
                        if (u_surface.contains(neigh_aux))
                            vinterp[j] = GetVal2D(image, u_surface, winner, neigh_aux, c, h, interp, mode, gamma);
                        else
                            vinterp[j] = valcenter[i];
                        new_w.v[gamma + j] = vinterp[j];
                    }

                    setCoeffs2D(vinterp, &new_w.m[M/2], interp, gamma);
                    //std::cout << "guardamos (" << neigh.y << ", " << neigh.x << ") key = " << key << std::endl;
                    pr_trial = std::pair<double, wmm::AugWmm >(valcenter[i], new_w);
                    trial_set_it = trial_set.insert(pr_trial);
                    pr_mapa = std::pair<int, std::multimap<double, wmm::AugWmm >::iterator>(key, trial_set_it);
                    mapa_trial.insert(pr_mapa);

                    u_surface.at(neigh) = valcenter[i];
                }
            }
            
            del(&winner);

        }
        
        //del(&new_w1);
        //del(&winner);
        free(state.data);

        return u_surface;

    }
    
    
}
