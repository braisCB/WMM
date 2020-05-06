#include "../TYPES/WMMStructs.h"
#include <math.h>
#include <stdlib.h>
#include <vector>
#include <map>
#include <tuple>
#include <iostream>


namespace wmm {

    template<int N> std::pair<double, wmm::NodeD> GetInterpValue(wmm::AniWmm_<double, N> &wave, wmm::NodeD &dp, wmm::NodeD &dd, wmm::NodeD &dn,
            wmm::Node3D &f0, wmm::Node3D &f1, double fn, double epsilon, int interp, int side, double alpha, double momentum) {
        // f0, f1 and f0 must be in polar coordinates
        double value, y0 = wave.v[side], y1 = wave.v[side + 1];
        wmm::NodeD diff = dn - dp - epsilon*dd;
        wmm::NodeD a = 1.0/wmm::norm(diff) * diff, c;
        wmm::Node3D c_polar;
        int s = (side < 2) ? side : side - 1;
        
        switch (interp) {
            case wmm::I_LINEAR:
                c_polar = (1.0 - epsilon)*f0 + epsilon*f1;
                value = (1.0 - epsilon)*y0 + epsilon*y1;
                break;
            case wmm::I_QUADATRIC:
                c_polar = epsilon*epsilon*wave.fm[3*s + 2] + epsilon*wave.fm[3*s + 1] + wave.fm[3*s];
                value = epsilon*epsilon*wave.m[3*s + 2] + epsilon*wave.m[3*s + 1] + wave.m[3*s];
                break;
            case wmm::I_SPLINE:
                c_polar = f0 + epsilon*(-1.0/6.0*wave.fm[2*s + 1] - 2.0/3.0*wave.fm[2*s] + f1 - f0 + epsilon*(
                        0.5*wave.fm[2*s] + epsilon/6.0*(wave.fm[2*s + 1] - wave.fm[2*s])));
                value = y0 + epsilon*(-wave.m[2*s + 1]/6.0 - wave.m[2*s]/3.0 + y1 - y0 + 
                        epsilon*(wave.m[2*s]/2.0 + epsilon*(wave.m[2*s + 1] - wave.m[2*s])/6.0));
                break;
            default: //I_HERMITE and I_PCHIP
                double t_2 = epsilon * epsilon;
                double t_3 = epsilon * t_2;
                c_polar = (2.0 * t_3 - 3.0 * t_2 + 1.0) * f0 + (t_3 - 2.0 * t_2 + epsilon) * wave.fm[2*s] +
                     (-2.0 * t_3 + 3.0 * t_2) * f1 + (t_3 - t_2) * wave.fm[2*s + 1];
                value = (2.0 * t_3 - 3.0 * t_2 + 1.0) * y0 + (t_3 - 2.0 * t_2 + epsilon) * wave.m[2*s] +
                       (-2.0 * t_3 + 3.0 * t_2) * y1 + (t_3 - t_2) * wave.m[2*s + 1];
        }
        
        if (value < y0 && value < y1) {
            value = (1.0 - epsilon)*y0 + epsilon*y1;
        }
        if (c_polar.z < f0.z && c_polar.z < f1.z) {
            c_polar = (1.0 - epsilon)*f0 + epsilon*f1;
        }
        
        c_polar.z = 1.0;
        c = wmm::to_cartesian(c_polar);
        
        //value += wmm::norm(diff) * (1.0 + pow(fabs(c.x*a.y - c.y*a.x), alpha)) * fn;
        //value += wmm::norm(diff) * (1.0 + cos) * fn;
        double cos_ca = c.x*a.x + c.y*a.y;
        double sin_ca = c.x*a.y - c.y*a.x;
        value += wmm::norm(diff) * exp(alpha*(1.0 - cos_ca)) * fn;
        if (fabs(acos(cos_ca)) > momentum) {
            double new_cos = (sin_ca > 0) ? c.x*cos(momentum) - c.y*sin(momentum) : c.x*cos(momentum) + c.y*sin(momentum);
            double new_sin = (sin_ca > 0) ? c.y*cos(momentum) + c.x*sin(momentum) : c.y*cos(momentum) - c.x*sin(momentum);
            //double angle = atan2(c_polar.y, c_polar.x);
            //angle = (side) ? angle + momentum : angle - momentum;
            c = wmm::NodeD(new_sin, new_cos);
            c = 1.0 / wmm::norm(c) * c;
            //c = a;
        }
        else
            c = a;
        
        return std::pair<double, wmm::NodeD>(value, c);
    }


    double GetEpsilonGradient(wmm::NodeD &dp, wmm::NodeD &dd, wmm::NodeD &dn, wmm::NodeD &fn) {

        double epsilon = 0.0;
        double A = -dd.y, B = dd.x, C = dd.y * dp.x - dd.x * dp.y;
        double den = A * fn.x + B * fn.y;
        double t = (A * dn.x + B * dn.y + C) / den;

        wmm::NodeD x(dn.y - t * fn.y, dn.x - t * fn.x);

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

    double GetEpsilonHopfLax(wmm::NodeD &dp, wmm::NodeD &dd, wmm::NodeD &dn, double y0, double y1) {

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

    template<int N> double GetEpsilonGoldenSearch(wmm::AniWmm_<double, N> &wave, wmm::NodeD &dn, wmm::NodeD &dd, wmm::NodeD &dp,
                                                  wmm::Node3D &f0, wmm::Node3D &f1, double fn,
                                                  double interp, int side, double alpha, double momentum) {

        double a = 0.0, b = 1.0, x1 = a + (1-RESPHI)*(b - a), x2 = a + RESPHI*(b - a),
                f_x1 = MAX_VAL, f_x2 = MAX_VAL, i1, i2, res;

        wmm::NodeD xtreme = dp + dd;

        wmm::NodeD F_x1, F_x2;

        double epsilon;
        double f_a = GetInterpValue(wave, dp, dd, dn, f0, f1, fn, 0.0, interp, side, alpha, momentum).first;
        double f_b = GetInterpValue(wave, dp, dd, dn, f0, f1, fn, 1.0, interp, side, alpha, momentum).first;

        if (f_a < f_b) {
            res = f_a; epsilon = 0.0;
        }
        else {
            res = f_b; epsilon = 1.0;
        }

        f_x1 = GetInterpValue(wave, dp, dd, dn, f0, f1, fn, x1, interp, side, alpha, momentum).first;
        f_x2 = GetInterpValue(wave, dp, dd, dn, f0, f1, fn, x2, interp, side, alpha, momentum).first;

        while (fabs(b - a) > TAU) {
            if(f_x1 < f_x2) {
                b = x2; x2 = x1; f_x2 = f_x1; x1 = a + (1 - RESPHI)*(b - a);
                f_x1 = GetInterpValue(wave, dp, dd, dn, f0, f1, fn, x1, interp, side, alpha, momentum).first;
            }
            else {
                a = x1; x1 = x2; f_x1 = f_x2; x2 = a + RESPHI*(b - a);
                f_x2 = GetInterpValue(wave, dp, dd, dn, f0, f1, fn, x2, interp, side, alpha, momentum).first;
            }
        }

        if (f_x1 < std::min(res, f_x2))
            epsilon = x1;
        else if (f_x2 < std::min(res, f_x1))
            epsilon = x2;

        return epsilon;
        
    }


    template<int N> std::tuple<double, double, wmm::NodeD> GetVal2D(wmm::Grid &image, wmm::AniWmm_<double, N> &wave, 
            wmm::Node &neigh, wmm::NodeD &h, int interp, int mode, double alpha, double momentum) {

        wmm::NodeD dn(h.y * neigh.y, h.x * neigh.x);
        double fn(image.at(neigh));
                
        std::tuple<double, double, wmm::NodeD> val(wmm::MAX_VAL, 1.0, wmm::NodeD(0.0, 0.0));
        
        if (wave.dir >=  0) {
            int ncomp = (wave.dir % 2 == 0) ? 2 : 1;
            for (int i=0; i<ncomp; i++) {
                wmm::Node p1(wave.p.y + wmm::yarray[(wave.dir + i + 1) % 8] - wmm::yarray[wave.dir],
                             wave.p.x + wmm::xarray[(wave.dir + i + 1) % 8] - wmm::xarray[wave.dir]);
                wmm::Node diff = neigh - p1;
                if (image.contains(p1) && (i < 2 || (abs(diff.x) <= 1 && abs(diff.y) <= 1))) {
                    wmm::Node p0(wave.p.y + wmm::yarray[(wave.dir + i) % 8] - wmm::yarray[wave.dir],
                                 wave.p.x + wmm::xarray[(wave.dir + i) % 8] - wmm::xarray[wave.dir]);
                    wmm::NodeD dp(h.y * p0.y, h.x * p0.x);
                    double y0 = wave.v[i], y1 = wave.v[i + 1];
                    wmm::NodeD dd(h.y * (wmm::yarray[(wave.dir + i + 1) % 8] - wmm::yarray[(wave.dir + i) % 8]),
                                  h.x * (wmm::xarray[(wave.dir + i + 1) % 8] - wmm::xarray[(wave.dir + i) % 8]));
                     
                    wmm::NodeD f0 = wave.f[i], f1 = wave.f[i+1];
                    
                    wmm::Node3D f0_polar = wmm::to_polar(f0), f1_polar = wmm::to_polar(f1);
                    wmm::NodeD grad = wmm::to_cartesian(0.5*(f0_polar + f1_polar));

                    double epsilon;
                    switch (mode) {
                        case wmm::M_GRADIENT:
                            epsilon = GetEpsilonGradient(dp, dd, dn, grad);
                            break;
                        case wmm::M_HOPFLAX:
                            epsilon = GetEpsilonHopfLax(dp, dd, dn, y0, y1);
                            break;
                        default: //M_GOLDENSEARCH
                            epsilon = GetEpsilonGoldenSearch(wave, dn, dd, dp, f0_polar, f1_polar, fn, interp, i, alpha, momentum);
                    }

                    auto res = GetInterpValue(wave, dp, dd, dn, f0_polar, f1_polar, fn, epsilon, interp, i, alpha, momentum);
                    if (res.first < std::get<0>(val))
                        val = std::make_tuple(res.first, epsilon, res.second);
                }
            }
            
            for (int i=0; i<ncomp; i++) {
                wmm::Node p1(wave.p.y + wmm::yarray[(wave.dir - i + 7) % 8] - wmm::yarray[wave.dir],
                             wave.p.x + wmm::xarray[(wave.dir - i + 7) % 8] - wmm::xarray[wave.dir]);
                wmm::Node diff = neigh - p1;
                if (image.contains(p1) && (i < 2 || (abs(diff.x) <= 1 && abs(diff.y) <= 1))) {
                    wmm::Node p0(wave.p.y + wmm::yarray[(wave.dir - i + 8) % 8] - wmm::yarray[wave.dir],
                                 wave.p.x + wmm::xarray[(wave.dir - i + 8) % 8] - wmm::xarray[wave.dir]);
                    wmm::NodeD dp(h.y * p0.y, h.x * p0.x);
                    double y0 = wave.v[i + 3], y1 = wave.v[i + 4];
                    wmm::NodeD dd(h.y * (wmm::yarray[(wave.dir - i + 7) % 8] - wmm::yarray[(wave.dir - i + 8) % 8]),
                                  h.x * (wmm::xarray[(wave.dir - i + 7) % 8] - wmm::xarray[(wave.dir - i + 8) % 8]));
                     
                    wmm::NodeD f0 = wave.f[i + 3], f1 = wave.f[i + 4];
                    
                    wmm::Node3D f0_polar = wmm::to_polar(f0), f1_polar = wmm::to_polar(f1);
                    wmm::NodeD grad = wmm::to_cartesian(0.5*(f0_polar + f1_polar));

                    double epsilon;
                    switch (mode) {
                        case wmm::M_GRADIENT:
                            epsilon = GetEpsilonGradient(dd, dp, dn, grad);
                            break;
                        case wmm::M_HOPFLAX:
                            epsilon = GetEpsilonHopfLax(dd, dp, dn, y0, y1);
                            break;
                        default: //M_GOLDENSEARCH
                            epsilon = GetEpsilonGoldenSearch(wave, dn, dd, dp, f0_polar, f1_polar, fn, interp, 3 + i, alpha, momentum);
                    }

                    auto res = GetInterpValue(wave, dp, dd, dn, f0_polar, f1_polar, fn, epsilon, interp, 3 + i, alpha, momentum);
                    if (res.first < std::get<0>(val))
                        val = std::make_tuple(res.first, epsilon, res.second);
                }
            }
        }
        else {
            wmm::NodeD diff(h.y * (neigh.y - wave.p.y), h.x * (neigh.x - wave.p.x));
            wmm::Node3D polar_diff = wmm::to_polar(diff);
            polar_diff.z = 1.0;
            wmm::Node3D f0_polar = wmm::to_polar(wave.f[0]);
            double y0 = wave.v[0];
            wmm::NodeD c = wmm::to_cartesian(momentum * f0_polar + (1.0 - momentum) * polar_diff);
            c = 1.0 / wmm::norm(c) * c;
            val = std::make_tuple(y0 + wmm::norm(diff) * fn, 0.0, c);
        }

        return val;

    }

    void setCoeffs2D(double *y, double *m, int pos, int interp, int right) {

        switch (interp) {
            case wmm::I_LINEAR:
                break;
            case wmm::I_QUADATRIC:
                m[0] = y[pos];
                if (pos%2 == 0) {
                    m[2] = (right) ? (y[(pos+2)%8] + y[pos] - 2.0*y[(pos+1)%8])/2.0: (y[(pos+6)%8] + y[pos] - 2.0*y[(pos+7)%8])/2.0;
                    m[1] = (right) ? y[(pos+1)%8] - y[pos] - m[2] : y[(pos+7)%8] - y[pos] - m[2];
                }
                else {
                    m[2] = (y[(pos+1)%8] + y[(pos+7)%8] - 2.0*y[pos])/2.0;
                    m[1] = (right) ? y[pos] - y[(pos+7)%8] + m[2] : y[pos] - y[(pos+1)%8] + m[2];
                }
                break;
            case wmm::I_SPLINE:
                if (pos%2 == 0) {
                    m[0] = 0.0;
                    m[1] = (right) ? 6.0*(y[pos] - 2.0*y[(pos+1)%8] + y[(pos+2)%8])/4.0 : 6.0*(y[pos] - 2.0*y[(pos+7)%8] + y[(pos+6)%8])/4.0;
                }
                else {
                    m[0] = 6.0*(y[(pos+7)%8] - 2.0*y[pos] + y[(pos+1)%8])/4.0;
                    m[1] = 0.0;
                }
                break;
            case wmm::I_HERMITE:
                if (pos%2 == 0) {
                    m[0] = (right) ? y[(pos + 1) % 8] - y[pos] : y[(pos + 7) % 8] - y[pos];
                    m[1] = (right) ? 1.0 / 2.0 * (y[(pos + 2) % 8] - y[pos]) : 1.0 / 2.0 * (y[(pos + 6) % 8] - y[pos]);
                }
                else {
                    m[0] = (right) ? 1.0 / 2.0 * (y[(pos + 1) % 8] - y[(pos + 7) % 8]) : -1.0 / 2.0 * (y[(pos + 1) % 8] - y[(pos + 7) % 8]);
                    m[1] = (right) ? y[(pos + 1) % 8] - y[pos] : y[(pos + 7) % 8] - y[pos];
                }
                break;
            default: //I_PCHIP
                if (pos%2 == 0) {
                    double d0 = (right) ? y[(pos + 1) % 8] - y[pos] : y[(pos + 7) % 8] - y[pos];
                    double d1 = (right) ? y[(pos + 2) % 8] - y[(pos + 1) % 8] : d1 = y[(pos + 6) % 8] - y[(pos + 7) % 8];
                    m[0] = (3.0 * d0 - d1) / 2.0;
                    if ((m[0] * d0) <= 0.0)
                        m[0] = 0.0;
                    else if (((d1 * d0) <= 0.0) && (fabs(m[0]) > fabs(3.0 * d0)))
                        m[0] = 3.0 * d0;
                    m[1] = (d0 * d1 <= 0.0) ? 0.0 : 2.0 * d0 * d1 / (d0 + d1);
                }
                else {
                    double d0 = (right) ? y[pos] - y[(pos + 7) % 8] : y[pos] - y[(pos + 1) % 8];
                    double d1 = (right) ? y[(pos + 1) % 8] - y[pos] : y[(pos + 7) % 8] - y[pos];
                    m[0] = (d0 * d1 <= 0.0) ? 0.0 : 2.0 * d0 * d1 / (d0 + d1);
                    m[1] = (3.0 * d1 - d0) / 2.0;
                    if ((m[1] * d1) <= 0.0)
                        m[1] = 0.0;
                    else if (((d0 * d1) <= 0.0) && (fabs(m[1]) > fabs(3.0 * d1)))
                        m[1] = 3.0 * d1;
                }
        }
        
    }


    template<int wmm_degree> std::pair<wmm::Grid, wmm::Grid_<int> > WmmAniSurfaceMomentum2D(wmm::Grid &image, std::vector<wmm::Node> &initials,
            wmm::NodeD &h, int interp, int mode, double alpha, double momentum, int nangles) {

        bool isnewpos[8];
        double valcenter[8];
        double imcenter_module[8];
        double imcenter_cos[8];
        double imcenter_sin[8];
        int angles[8];
        
        int degree = (wmm_degree == 0) ? 1 : wmm_degree, angle, winner_dir, winner_angle;
                
        double module_values[degree], cos_values[degree], sin_values[degree];

        wmm::Grid u_surface = wmm::Grid(wmm::MAX_VAL, image.rows, image.cols, nangles);
        wmm::Grid u_epsilon = wmm::Grid(wmm::MAX_VAL, image.rows, image.cols, nangles);
        wmm::Grid_<int> u_dirs = wmm::Grid_<int>(-1, image.rows, image.cols, nangles);
        wmm::Grid f_surface = wmm::Grid(image.rows, image.cols, 2*nangles);
        wmm::Grid_<unsigned char> state = wmm::Grid_<unsigned char>(image.rows, image.cols, nangles);

        std::multimap<double, wmm::AniWmm_<double, wmm_degree> > trial_set;
        std::map<int, typename std::multimap<double, wmm::AniWmm_<double, wmm_degree> >::iterator> mapa_trial;

        typename std::multimap<double, wmm::AniWmm_<double, wmm_degree> >::iterator trial_set_it;
        typename std::map<int, typename std::multimap<double, wmm::AniWmm_<double, wmm_degree> >::iterator>::iterator mapa_trial_it;
        std::pair<double, wmm::AniWmm_<double, wmm_degree> > pr_trial;
        std::pair<int, typename std::multimap<double, wmm::AniWmm_<double, wmm_degree> >::iterator> pr_mapa;

        int key, i, j, rango_min, rango_max;
        wmm::AniWmm_<double, wmm_degree> winner, new_w;
        wmm::Node neigh, father;

        // Initialization
        for (i = 0; i < (int) initials.size(); i++) {
            key = initials[i].y + u_surface.rows*initials[i].x;
            if (mapa_trial.find(key) == mapa_trial.end() && u_surface.contains(initials[i])) {
                u_surface.at(initials[i]) = 0.0;
                winner.dir = -1;
                winner.angle = 0;
                winner.v[0] = 0.0;
                winner.p = initials[i];
                for (int f=0; f < nangles; f++) {
                    state.at(initials[i], f) = wmm::P_ALIVE;
                }    
                pr_trial = std::pair<double, wmm::AniWmm_<double, wmm_degree> >(0.0, winner);
                trial_set_it = trial_set.insert(pr_trial);
                pr_mapa = std::pair<int, typename std::multimap<double, wmm::AniWmm_<double, wmm_degree> >::iterator>(key, trial_set_it);
                mapa_trial.insert(pr_mapa);
            }
        }
        
        while (!trial_set.empty()) {
            
            trial_set_it = trial_set.begin();
            key = trial_set_it->second.p.y + u_surface.rows*(trial_set_it->second.p.x + u_surface.cols*trial_set_it->second.angle);
            // key = trial_set_it->second.p.y * u_surface.cols + trial_set_it->second.p.x;
            mapa_trial_it = mapa_trial.find(key);

            if (mapa_trial_it == mapa_trial.end()) {
                printf("ERROR: bad map alloc");
                exit(-1);
            }

            if (mapa_trial_it->second != trial_set_it) {
                printf("ERROR: bad trial/map alloc");
                exit(-1);
            }

            winner = trial_set_it->second;
            
            //std::cout << "WINNER: (" << winner.p.y << ", " << winner.p.x << ") -> dir = " << winner.dir << ", angle " << winner.angle << std::endl;
            //std::cout << "val = " << winner.v[0] << ", f = (" << winner.f[0].y << ", " << winner.f[0].x << ")" << std::endl;
            /*for (int i = 0; i < degree; i++) {
                std::cout << winner.m[i] << ", ";
            }
            std::cout << std::endl;*/
            
            trial_set.erase(trial_set_it);
            mapa_trial.erase(mapa_trial_it);

            state.at(winner.p, winner.angle) = wmm::P_ALIVE;
            
            if (winner.dir < 0) {
                rango_min = 1;
                rango_max = 9;
            }
            else {
                rango_min = - 1 - (winner.dir % 2 == 0);
                rango_max = 2 + (winner.dir % 2 == 0);
                
                valcenter[(winner.dir + 4) % 8] = winner.v[6];
                wmm::Node3D polar = wmm::to_polar(winner.f[6]);
                imcenter_module[(winner.dir + 4) % 8] = polar.z;
                imcenter_cos[(winner.dir + 4) % 8] = polar.x;
                imcenter_sin[(winner.dir + 4) % 8] = polar.y;
                valcenter[(winner.dir + rango_min - 1 + 8) % 8] = winner.v[4];
                polar = wmm::to_polar(winner.f[4]);
                imcenter_module[(winner.dir + rango_min - 1 + 8) % 8] = polar.z;
                imcenter_cos[(winner.dir + rango_min - 1 + 8) % 8] = polar.x;
                imcenter_sin[(winner.dir + rango_min - 1 + 8) % 8] = polar.y;
                valcenter[(winner.dir + rango_max) % 8] = winner.v[1];
                polar = wmm::to_polar(winner.f[1]);
                imcenter_module[(winner.dir + rango_max) % 8] = polar.z;
                imcenter_cos[(winner.dir + rango_max) % 8] = polar.x;
                imcenter_sin[(winner.dir + rango_max) % 8] = polar.y;
                if (winner.dir % 2 != 0) {
                    valcenter[(winner.dir + 3) % 8] = winner.v[2];
                    polar = wmm::to_polar(winner.f[2]);
                    imcenter_module[(winner.dir + 3) % 8] = polar.z;
                    imcenter_cos[(winner.dir + 3) % 8] = polar.x;
                    imcenter_sin[(winner.dir + 3) % 8] = polar.y;
                    valcenter[(winner.dir + 5) % 8] = winner.v[5];
                    polar = wmm::to_polar(winner.f[5]);
                    imcenter_module[(winner.dir + 5) % 8] = polar.z;
                    imcenter_cos[(winner.dir + 5) % 8] = polar.x;
                    imcenter_sin[(winner.dir + 5) % 8] = polar.y;
                }
            }
            
            //std::cout << "rango = ( " << rango_min << ", " << rango_max << ")" << std::endl;
            
            // Neighbour temptative value computation
            for (j = rango_min; j < rango_max; j++) {
                i = (winner.dir + j + 8) % 8;
                neigh = winner.p + wmm::Node(wmm::yarray[i], wmm::xarray[i]);
                isnewpos[i] = false;
                if (u_surface.contains(neigh)) {
                    auto val_neigh = GetVal2D(image, winner, neigh, h, interp, mode, alpha, momentum);
                    if (isnan(std::get<2>(val_neigh).y) || isnan(std::get<2>(val_neigh).x))
                        continue;
                    angle = std::max(0, std::min(nangles - 1, (int) (nangles * (atan2(std::get<2>(val_neigh).y, std::get<2>(val_neigh).x) + M_PI) / (2.0 * M_PI))));
                    valcenter[i] = u_surface.at(neigh, angle);
                    double epsilon = u_epsilon.at(neigh, angle);
                    angles[i] = angle;
                    if (state.at(neigh, angle) != wmm::P_ALIVE && (std::get<0>(val_neigh) < valcenter[i] || 
                            (fabs(std::get<0>(val_neigh) - valcenter[i]) < wmm::TAU && std::get<1>(val_neigh) < epsilon))) {
                        u_epsilon.at(neigh, angle) = std::get<1>(val_neigh);
                        u_dirs.at(neigh, angle) = 8*winner.angle + ((i + 4) % 8);
                    }
                    if (state.at(neigh, angle) != wmm::P_ALIVE && std::get<0>(val_neigh) < valcenter[i]) {
                        valcenter[i] = std::get<0>(val_neigh);
                        f_surface.at(neigh, 2*angle) = std::get<2>(val_neigh).y;
                        f_surface.at(neigh, 2*angle + 1) = std::get<2>(val_neigh).x;
                        u_surface.at(neigh, angle) = valcenter[i];
                        if (winner.dir < 0 || abs(j) < 2)
                            isnewpos[i] = true;
                    }
                    wmm::Node3D polar = wmm::to_polar(wmm::NodeD(f_surface.at(neigh, 2*angle), f_surface.at(neigh, 2*angle + 1)));
                    imcenter_module[i] = polar.z;
                    imcenter_cos[i] = polar.x;
                    imcenter_sin[i] = polar.y;
                }
                else {
                    valcenter[i] = winner.v[0];
                    angles[i] = winner.angle;
                    wmm::Node3D polar = wmm::to_polar(winner.f[0]);
                    imcenter_module[i] = polar.z;
                    imcenter_cos[i] = polar.x;
                    imcenter_sin[i] = polar.y;
                }
                //std::cout << "neigh = ( " << neigh.y << ", " << neigh.x << "), val -> " << valcenter[i] << ", angle = " << 
                //        angles[i] << ", isnew = " << isnewpos[i] << std::endl;
            }

            // Update
            for (i = 0; i < 8; i++) {
                if (isnewpos[i]) {
                    neigh = winner.p + wmm::Node(wmm::yarray[i], wmm::xarray[i]);
                    key = neigh.y + u_surface.rows*(neigh.x + u_surface.cols*angles[i]);
                    if (state.at(neigh, angles[i]) == wmm::P_TRIAL) {
                        mapa_trial_it = mapa_trial.find(key);
                        trial_set.erase(mapa_trial_it->second);
                        mapa_trial.erase(mapa_trial_it);
                    }
                    else {
                        state.at(neigh, angles[i]) = wmm::P_TRIAL;
                    }
                    new_w.p = neigh;
                    new_w.dir = i;
                    new_w.angle = angles[i];
                    
                    new_w.v[0] = new_w.v[3] = valcenter[new_w.dir];
                    new_w.v[1] = valcenter[(new_w.dir+1)%8];
                    new_w.v[2] = valcenter[(new_w.dir+2)%8];
                    new_w.v[4] = valcenter[(new_w.dir+7)%8];
                    new_w.v[5] = valcenter[(new_w.dir+6)%8];
                    new_w.v[6] = winner.v[0];
                    new_w.f[0] = new_w.f[3] = wmm::to_cartesian(wmm::Node3D(imcenter_cos[i], imcenter_sin[i], imcenter_module[i]));
                    new_w.f[1] = wmm::to_cartesian(wmm::Node3D(imcenter_cos[(i+1)%8], imcenter_sin[(i+1)%8], imcenter_module[(i+1)%8]));
                    new_w.f[2] = wmm::to_cartesian(wmm::Node3D(imcenter_cos[(i+2)%8], imcenter_sin[(i+2)%8], imcenter_module[(i+2)%8]));
                    new_w.f[4] = wmm::to_cartesian(wmm::Node3D(imcenter_cos[(i+7)%8], imcenter_sin[(i+7)%8], imcenter_module[(i+7)%8]));
                    new_w.f[5] = wmm::to_cartesian(wmm::Node3D(imcenter_cos[(i+6)%8], imcenter_sin[(i+6)%8], imcenter_module[(i+6)%8]));
                    new_w.f[6] = winner.f[0];
                    
                    setCoeffs2D(valcenter, new_w.m, new_w.dir, interp, 1);
                    setCoeffs2D(valcenter, &(new_w.m[degree/4]), (new_w.dir + 1) % 8, interp, 1);
                    setCoeffs2D(valcenter, &(new_w.m[2*degree/4]), new_w.dir, interp, 0);
                    setCoeffs2D(valcenter, &(new_w.m[3*degree/4]), (new_w.dir + 7) % 8, interp, 0);
                    if (wmm_degree > 0) {
                        setCoeffs2D(imcenter_module, module_values, new_w.dir, interp, 1);
                        setCoeffs2D(imcenter_module, &(module_values[degree/4]), (new_w.dir + 1) % 8, interp, 1);
                        setCoeffs2D(imcenter_module, &(module_values[2*degree/4]), new_w.dir, interp, 0);
                        setCoeffs2D(imcenter_module, &(module_values[3*degree/4]), (new_w.dir + 7) % 8, interp, 0);
                        setCoeffs2D(imcenter_cos, cos_values, new_w.dir, interp, 1);
                        setCoeffs2D(imcenter_cos, &(cos_values[degree/4]), (new_w.dir + 1) % 8, interp, 1);
                        setCoeffs2D(imcenter_cos, &(cos_values[2*degree/4]), new_w.dir, interp, 0);
                        setCoeffs2D(imcenter_cos, &(cos_values[3*degree/4]), (new_w.dir + 7) % 8, interp, 0);
                        setCoeffs2D(imcenter_sin, sin_values, new_w.dir, interp, 1);
                        setCoeffs2D(imcenter_sin, &(sin_values[degree/4]), (new_w.dir + 1) % 8, interp, 1);
                        setCoeffs2D(imcenter_sin, &(sin_values[2*degree/4]), new_w.dir, interp, 0);
                        setCoeffs2D(imcenter_sin, &(sin_values[3*degree/4]), (new_w.dir + 7) % 8, interp, 0);
                        for (j = 0; j < wmm_degree; j++) {
                            new_w.fm[j].z = module_values[j];
                            new_w.fm[j].x = cos_values[j];
                            new_w.fm[j].y = sin_values[j];
                        }
                    }
                    
                    //std::cout << "NEW: (" << new_w.p.y << ", " << new_w.p.x << ") -> dir = " << new_w.dir << ", angle " << new_w.angle << std::endl;
                    //std::cout << "val = " << new_w.v[0] << ", f = (" << new_w.f[0].y << ", " << new_w.f[0].x << ")" << std::endl;

                    pr_trial = std::pair<double, wmm::AniWmm_<double, wmm_degree> >(valcenter[i], new_w);
                    trial_set_it = trial_set.insert(pr_trial);
                    pr_mapa = std::pair<int, typename std::multimap<double, wmm::AniWmm_<double, wmm_degree> >::iterator>(key, trial_set_it);
                    mapa_trial.insert(pr_mapa);
                }
                isnewpos[i] = false;
            }

        }
        
        free(state.data);
        free(f_surface.data);
        free(u_epsilon.data);

        return std::pair<wmm::Grid, wmm::Grid_<int> >(u_surface, u_dirs);

    }
    
    
}
