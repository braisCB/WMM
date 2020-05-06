#include "../TYPES/WMMStructs.h"
#include "../TYPES/utils.h"
#include <math.h>
#include <stdlib.h>
#include <vector>
#include <map>
#include <iostream>


namespace wmm {

    template<int N> double GetInterpValue(Grid &image, Wmm_<double, N> &wave, NodeD &dp, NodeD &dd, NodeD &dn, NodeD &h,
                                          NodeD &f0, NodeD &f1, NodeD &fn, double epsilon, int interp,
                                          int side) {

        NodeD d = dp + epsilon * dd, f_aux;
        double value, y0 = wave.v[0], y1 = wave.v[side + 1], ft;

        switch (interp) {
            default:
                f_aux = (1.0 - epsilon)*f0 + epsilon*f1;
                ft = norm(f_aux);
                value = (1.0 - epsilon)*y0 + epsilon*y1;
                break;
        }

        double ft2 = GetIntegralValue(image, d, dn, ft, norm(fn), interp);
        value += norm(h * (dn - d)) * ft2; //(ft + norm(fn))/2.0;

        return value;
    }

    double GetEpsilonGradient(NodeD dd, NodeD dp, NodeD dn, NodeD &fn) {

        double epsilon = 0.0;
        double A = -dd.y, B = dd.x, C = dd.y * dp.x - dd.x * dp.y;
        double den = A * fn.x + B * fn.y;
        double t = (A * dn.x + B * dn.y + C) / den;

        NodeD x(dn.y - t * fn.y, dn.x - t * fn.x);

        if (fabs(dd.x) > 0.0 && fabs(den) > 0.0)
            epsilon = (x.x - dp.x) / dd.x;
        else if (fabs(dd.y) > 0.0 && fabs(den) > 0.0)
            epsilon = (x.y - dp.y) / dd.y;
        else if (fabs(den) == 0.0 && norm(dd) > 0.0) {
            double dist = fabs(A * dn.x + B * dn.y + C) / sqrt(A * A + B * B);
            epsilon = (norm(dn - dp) - dist) / (fabs(dd.x) + fabs(dd.y));
        }
        else
            return 0.0;

        if (epsilon < 0.0)
            epsilon = 0.0;
        else if (epsilon > 1.0)
            epsilon = 1.0;
        return epsilon;

    }

    /*template<int N> double GetEpsilonHopfLax(Grid &image, Wmm_<double, N> &wave, NodeD dn, NodeD dd, NodeD dp,
                                             NodeD &f0, NodeD &f1, int interp, int side) {

        NodeD wave_d = dp + dd;
        NodeD ori = dn - dp, w0, w1;
        double grad_0 = -GetGradientValue(wave, 0.5, interp, side);
        double grad_1 = -GetGradientValue(wave, 1., interp, side);

        w0 = (dd.x == 0) ? NodeD(-sgn(ori.x)*grad_0, dd.y) : NodeD(-dd.x, sgn(ori.y)*grad_0);
        w1 = (dd.x == 0) ? NodeD(-sgn(ori.x)*grad_1, dd.y) : NodeD(-dd.x, sgn(ori.y)*grad_1);

        double A = -w0.y, B = w0.x, C = w0.y * dp.x - w0.x * dp.y;
        double den = A * w1.x + B * w1.y;
        NodeD fn;
        if (den == 0)
            fn = w1;
        else {
            double t = (A * wave_d.x + B * wave_d.y + C) / den;
            NodeD x(wave_d.y - t * w1.y, wave_d.x - t * w1.x);
            fn = x - dn;
        }
        return GetEpsilonGradient(dd, dp, dn, w0);

    }*/

    double GetEpsilonHopfLax(NodeD dd, NodeD dp, NodeD dn, double y0, double y1) {

        if (norm(dn - dp - dd) < TAU) {
            return 0.0;
        }
        else if (norm(dn - dp) < TAU) {
            return 1.0;
        }
        NodeD xy = dn - dp;
        double nxy = norm(xy);
        double nyz = norm(dd);

        double c_alpha = (xy.x * dd.x + xy.y * dd.y) / (nxy * nyz);
        double c_delta = (y1 - y0) / nyz;

        if (nyz == 0.0 || c_alpha <= c_delta || c_alpha == 1.0) {
            return 0.0;
        }

        NodeD xz = dn - dp - dd;
        double nxz = norm(xz);
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

    template<int N> double GetEpsilonGoldenSearch(Grid &image, Wmm_<double, N> &wave, NodeD &dn, NodeD &dd, NodeD &dp,
                                                  NodeD &h, NodeD &f0, NodeD &f1, NodeD &fn, int interp, int side) {

        double a = 0.0, b = 1.0, x1 = a + (1-RESPHI)*(b - a), x2 = a + RESPHI*(b - a),
                f_x1 = MAX_VAL, f_x2 = MAX_VAL, res;

        NodeD F_x1, F_x2;

        double epsilon;
        double f_a = GetInterpValue(image, wave, dp, dd, dn, h, f0, f1, fn, 0.0, interp, side);
        double f_b = GetInterpValue(image, wave, dp, dd, dn, h, f0, f1, fn, 1.0, interp, side);

        if (f_a < f_b) {
            res = f_a; epsilon = 0.0;
        }
        else {
            res = f_b; epsilon = 1.0;
        }

        f_x1 = GetInterpValue(image, wave, dp, dd, dn, h, f0, f1, fn, x1, interp, side);
        f_x2 = GetInterpValue(image, wave, dp, dd, dn, h, f0, f1, fn, x2, interp, side);

        while (fabs(b - a) > TAU) {
            if(f_x1 < f_x2) {
                b = x2; x2 = x1; f_x2 = f_x1; x1 = a + (1 - RESPHI)*(b - a);
                f_x1 = GetInterpValue(image, wave, dp, dd, dn, h, f0, f1, fn, x1, interp, side);
            }
            else {
                a = x1; x1 = x2; f_x1 = f_x2; x2 = a + RESPHI*(b - a);
                f_x2 = GetInterpValue(image, wave, dp, dd, dn, h, f0, f1, fn, x2, interp, side);
            }
        }

        if (f_x1 < res) {
            epsilon = x1;
            res = f_x1;
        }
        if (f_x2 < res) {
            epsilon = x2;
        }

        return epsilon;
        
    }


    template<int N> double GetVal2D(Grid &image, Grid &u_surface, Wmm_<double, N> &wave, Node &neigh, NodeD &h, int interp, int mode) {

        NodeD f0(image.at(wave.p, 0), image.at(wave.p, 1)), fn(image.at(neigh, 0), image.at(neigh, 1));
        double y0 = wave.v[0];

        if (isinf(norm(f0)) || isnan(norm(f0)))
            f0 = fn;

        NodeD dp(wave.p.y, wave.p.x), dn(neigh.y, neigh.x);

        double val = MAX_VAL;
        if (wave.dir < 0) {
            NodeD diff(h.y * (neigh.y - wave.p.y), h.x * (neigh.x - wave.p.x));
            NodeD dirs(0.0, 0.0);
            double ft = GetIntegralValue(image, dp, dn, norm(f0), norm(fn), interp);
            val = y0 + norm(diff) * ft; //(norm(f0) + norm(fn)) / 2.0;
        }
        else {

            Node p(wave.p.y + yarray[(wave.dir + 1) % 8] - yarray[wave.dir],
                    wave.p.x + xarray[(wave.dir + 1) % 8] - xarray[wave.dir]);
            double res1 = MAX_VAL;

            if (u_surface.contains(p)) {
                double y1 = wave.v[1];

                NodeD dd((yarray[(wave.dir + 1) % 8] - yarray[wave.dir]), (xarray[(wave.dir + 1) % 8] - xarray[wave.dir]));

                NodeD f1(image.at(p, 0), image.at(p, 1));
                if (isinf(norm(f1)) || isnan(norm(f1)))
                    f1 = fn;

                double epsilon;
                switch (mode) {
                    case M_GRADIENT:
                        epsilon = GetEpsilonGradient(h*dd, h*dp, h*dn, fn);
                        break;
                    case M_HOPFLAX:
                        epsilon = GetEpsilonHopfLax(h*dd, h*dp, h*dn, y0, y1);
                        //epsilon = GetEpsilonHopfLax(image, wave, h*dn, h*dd, h*dp, f0, f1, interp, S_RIGHT);
                        break;
                    default: //M_GOLDENSEARCH
                        epsilon = GetEpsilonGoldenSearch(image, wave, dn, dd, dp, h, f0, f1, fn, interp, S_RIGHT);
                }
                res1 = GetInterpValue(image, wave, dp, dd, dn, h, f0, f1, fn, epsilon, interp, S_RIGHT);
            }

            p = Node(wave.p.y + yarray[(wave.dir + 7) % 8] - yarray[wave.dir],
                    wave.p.x + xarray[(wave.dir + 7) % 8] - xarray[wave.dir]);
            double res2 = MAX_VAL;

            if (u_surface.contains(p)) {
                double y1 = wave.v[2];

                NodeD dd((yarray[(wave.dir + 7) % 8] - yarray[wave.dir]), (xarray[(wave.dir + 7) % 8] - xarray[wave.dir]));

                NodeD f1(image.at(p, 0), image.at(p, 1));
                if (isinf(norm(f1)) || isnan(norm(f1)))
                    f1 = fn;

                double epsilon;
                switch (mode) {
                    case M_GRADIENT:
                        epsilon = GetEpsilonGradient(h*dd, h*dp, h*dn, fn);
                        break;
                    case M_HOPFLAX:
                        epsilon = GetEpsilonHopfLax(h*dd, h*dp, h*dn, y0, y1);
                        //epsilon = GetEpsilonHopfLax(image, wave, h*dn, h*dd, h*dp, f0, f1, interp, S_LEFT);
                        break;
                    default: //M_GOLDENSEARCH
                        epsilon = GetEpsilonGoldenSearch(image, wave, dn, dd, dp, h, f0, f1, fn, interp, S_LEFT);
                }
                res2 = GetInterpValue(image, wave, dp, dd, dn, h, f0, f1, fn, epsilon, interp, S_LEFT);

            }

            val = std::min(res1, res2);

        }

        return val;

    }


    template<int wmm_degree> void WmmIsoSurface2D(Grid &image, std::vector<Node> &initials, NodeD &h,
                                                       int interp, int mode, Grid &u_surface) {

        bool isnewpos[8];
        double valcenter[8];
        double imcenter[8];

        Grid_<unsigned char> state = Grid_<unsigned char>(image.rows, image.cols);

        std::multimap<double, Wmm_<double, wmm_degree> > trial_set;
        std::map<int, typename std::multimap<double, Wmm_<double, wmm_degree> >::iterator> mapa_trial;

        typename std::multimap<double, Wmm_<double, wmm_degree> >::iterator trial_set_it;
        typename std::map<int, typename std::multimap<double, Wmm_<double, wmm_degree> >::iterator>::iterator mapa_trial_it;
        std::pair<double, Wmm_<double, wmm_degree> > pr_trial;
        std::pair<int, typename std::multimap<double, Wmm_<double, wmm_degree> >::iterator> pr_mapa;

        int key, i;
        Wmm_<double, wmm_degree> winner, new_w;
        Node neigh, father;

        // Initialization
        for (i = 0; i < (int) initials.size(); i++) {
            key = initials[i].y * u_surface.cols + initials[i].x;
            if (mapa_trial.find(key) == mapa_trial.end() && u_surface.contains(initials[i])) {
                u_surface.at(initials[i]) = 0.0;
                winner.dir = -1;
                winner.v[0] = 0.0;
                winner.p = initials[i];
                state.at(initials[i]) = P_TRIAL;
                pr_trial = std::pair<double, Wmm_<double, wmm_degree> >(0.0, winner);
                trial_set_it = trial_set.insert(pr_trial);
                pr_mapa = std::pair<int, typename std::multimap<double, Wmm_<double, wmm_degree> >::iterator>(key, trial_set_it);
                mapa_trial.insert(pr_mapa);
            }
        }

        while (!trial_set.empty()) {

            trial_set_it = trial_set.begin();
            key = trial_set_it->second.p.y * u_surface.cols + trial_set_it->second.p.x;
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

            trial_set.erase(trial_set_it);
            mapa_trial.erase(mapa_trial_it);

            state.at(winner.p) = P_ALIVE;

            if (winner.dir >= 0) {
                father = winner.p + Node(yarray[(winner.dir + 4) % 8], xarray[(winner.dir + 4) % 8]);
                for (i = 0; i < 8; i++) {
                    neigh = father + Node(yarray[i], xarray[i]);
                    isnewpos[i] = false;
                    valcenter[i] = u_surface.contains(neigh) ? u_surface.at(neigh) : MAX_VAL;
                    imcenter[i] = u_surface.contains(neigh) ? norm(NodeD(image.at(neigh, 0), image.at(neigh, 1))) : MAX_VAL;
                }
                /*for (i = 0; i < 8; i++) {
                    if (valcenter[i] == MAX_VAL) {
                        if (valcenter[(i + 2) % 8] == MAX_VAL && valcenter[(i + 6) % 8] == MAX_VAL)
                            valcenter[i] = winner.v[0];
                        else if (valcenter[(i + 2) % 8] == MAX_VAL)
                            valcenter[i] = 1.5*valcenter[(i + 7) % 8] - valcenter[(i + 6) % 8];
                        else
                            valcenter[i] = 1.5*valcenter[(i + 1) % 8] - valcenter[(i + 2) % 8];
                    }
                    if (imcenter[i] == MAX_VAL) {
                        if (imcenter[(i + 2) % 8] == MAX_VAL && imcenter[(i + 6) % 8] == MAX_VAL)
                            imcenter[i] = norm(NodeD(image.at(winner.p, 0), image.at(winner.p, 1)));
                        else if (imcenter[(i + 2) % 8] == MAX_VAL)
                            imcenter[i] = 1.5*imcenter[(i + 7) % 8] - imcenter[(i + 6) % 8];
                        else
                            imcenter[i] = 1.5*imcenter[(i + 1) % 8] - imcenter[(i + 2) % 8];
                    }
                }*/
                i = winner.dir;
                winner.v[0] = valcenter[i];
                winner.v[1] = valcenter[(i+1)%8];
                winner.v[2] = valcenter[(i+7)%8];

                //setCoeffs(valcenter, winner.m, i, interp);
                //setCoeffs(imcenter, winner.fm, i, interp);
            }

            // Neighbour temptative value computation
            for (i = 0; i < 8; i++) {
                neigh = winner.p + Node(yarray[i], xarray[i]);
                isnewpos[i] = false;
                valcenter[i] = u_surface.contains(neigh) ? u_surface.at(neigh) : MAX_VAL;
                imcenter[i] = u_surface.contains(neigh) ? norm(NodeD(image.at(neigh, 0), image.at(neigh, 1))) : MAX_VAL;
                if (u_surface.contains(neigh) && state.at(neigh) != P_ALIVE) {
                    double val_neigh = GetVal2D(image, u_surface, winner, neigh, h, interp, mode);
                    // if ((val_neigh - valcenter[i]) < -TAU) {
                    if (val_neigh < valcenter[i]) {
                        valcenter[i] = val_neigh;
                        isnewpos[i] = true;
                    }
                }
            }

            // Update
            for (i = 0; i < 8; i++) {
                if (isnewpos[i]) {
                    neigh = winner.p + Node(yarray[i], xarray[i]);
                    u_surface.at(neigh) = valcenter[i];
                    if (valcenter[(i+1) % 8] == MAX_VAL || valcenter[(i+7) % 8] == MAX_VAL)
                        continue;
                    key = neigh.y * u_surface.cols + neigh.x;
                    if (state.at(neigh) == P_TRIAL) {
                        mapa_trial_it = mapa_trial.find(key);
                        trial_set.erase(mapa_trial_it->second);
                        mapa_trial.erase(mapa_trial_it);
                    }
                    else {
                        state.at(neigh) = P_TRIAL;
                    }
                    new_w.p = neigh;
                    new_w.dir = i;
                    new_w.v[0] = valcenter[i];

                    pr_trial = std::pair<double, Wmm_<double, wmm_degree> >(valcenter[i], new_w);
                    trial_set_it = trial_set.insert(pr_trial);
                    pr_mapa = std::pair<int, typename std::multimap<double, Wmm_<double, wmm_degree> >::iterator>(key, trial_set_it);
                    mapa_trial.insert(pr_mapa);

                    /*if ((neigh.x == 7) && (neigh.y == 0 || neigh.y == 100)) {
                        std::cout << "winner: " << winner.p.y << ", " << winner.p.x << " : " << winner.dir << " - " << winner.v[0] << std::endl;
                        std::cout << "vs: " << winner.v[1] << ", " << winner.v[2] << std::endl;
                        std::cout << "fs0: " << winner.fm[0] << ", " << winner.v[1] << std::endl;
                        std::cout << "fs1: " << winner.fm[2] << ", " << winner.v[3] << std::endl;
                        std::cout << "update: " << new_w.p.y << ", " << new_w.p.x << " - " << new_w.v[0] << std::endl;
                    }*/
                }
            }

        }

        free(state.data);

        return;

    }


    template<int wmm_degree> Grid WmmIsoSurface2D(Grid &image, std::vector<Node> &initials, NodeD &h,
                                                       int interp, int mode) {
        Grid u_surface = Grid(MAX_VAL, image.rows, image.cols, 1);
        WmmIsoSurface2D<wmm_degree>(image, initials, h, interp, mode, u_surface);
        return u_surface;
    }
    
    
}
