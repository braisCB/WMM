#include "../TYPES/WMMStructs.h"
#include <math.h>
#include <stdlib.h>
#include <vector>
#include <map>
#include <iostream>


namespace wmm {

    double GetInterpValue(wmm::NodeD &p0, wmm::NodeD &p1, wmm::NodeD &pn,
                          double y0, double y1, wmm::NodeD &c,
                          double epsilon) {
        
        wmm::NodeD diff = (1.0 - epsilon) * p0 + epsilon * p1 - pn;
        wmm::NodeD a = 1.0/wmm::norm(diff) * diff;
        double value = (1.0 - epsilon)*y0 + epsilon*y1 + wmm::norm(diff) * sqrt(1.0 + (c.x*a.x + c.y*a.y)*(c.x*a.x + c.y*a.y));
        
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
    
    
    double GetEpsilonGoldenSearch(wmm::NodeD &p0, wmm::NodeD &p1, wmm::NodeD &pn,
                                       double y0, double y1, wmm::NodeD &c) {
        
        double a = 0.0, b = 1.0, x1 = a + (1-RESPHI)*(b - a), x2 = a + RESPHI*(b - a),
                f_x1 = MAX_VAL, f_x2 = MAX_VAL, res, epsilon;
    

        wmm::NodeD F_x1, F_x2;

        double f_a = GetInterpValue(p0, p1, pn, y0, y1, c, 0.0);
        double f_b = GetInterpValue(p0, p1, pn, y0, y1, c, 1.0);

        if (f_a < f_b) {
            res = f_a;
            epsilon = 0.0;
        }
        else {
            res = f_b;
            epsilon = 1.0;
        }

        f_x1 = GetInterpValue(p0, p1, pn, y0, y1, c, x1);
        f_x2 = GetInterpValue(p0, p1, pn, y0, y1, c, x2);

        while (fabs(b - a) > TAU) {
            if(f_x1 < f_x2) {
                b = x2; x2 = x1; f_x2 = f_x1; x1 = a + (1 - RESPHI)*(b - a);
                f_x1 = GetInterpValue(p0, p1, pn, y0, y1, c, x1);
            }
            else {
                a = x1; x1 = x2; f_x1 = f_x2; x2 = a + RESPHI*(b - a);
                f_x2 = GetInterpValue(p0, p1, pn, y0, y1, c, x2);
            }
        }

        if (f_x1 < std::min(res, f_x2))
            epsilon = x1;
        else if (f_x2 < std::min(res, f_x1))
            epsilon = x2;

        return epsilon;
        
    }

    bool is_considered(wmm::Grid_<unsigned char> &state, wmm::Node &p) {
    
        if (!state.contains(p))
            return false;
        
        if (state.at(p) != wmm::P_ALIVE)
            return false;
        
        wmm::Node neigh;
        bool considered = false;
        
        for (int i = 0; i < 8; i++) {
            neigh = p + wmm::Node(wmm::yarray[i], wmm::xarray[i]);
            if (state.contains(neigh) && state.at(neigh) != wmm::P_ALIVE) {
                considered = true;
                break;
            }
        }
        
        return considered;
    
    }

    double GetVal2D(wmm::Grid &image, wmm::Grid &u_surface, wmm::Grid_<unsigned char> &state,
            wmm::Node &neigh, wmm::NodeD &h, double radius, int mode) {

        wmm::NodeD c(image.at(neigh, 0), image.at(neigh, 1)), p0, p1, pn = wmm::NodeD(h.y*neigh.y, h.x*neigh.x), dd;
        wmm::Node n0, n1;
        double y0, y1;
                
        double val = wmm::MAX_VAL, candidate_val;
        
        int i_radius = (int) radius;
        double h_radius = radius * std::max(h.x, h.y) + 1e-8, epsilon;
        
        for (int i=-i_radius; i<=i_radius; i++) {
            for (int j=-i_radius; j<=i_radius; j++) {
                if (i == 0 && j == 0)
                    continue;
                n0 = wmm::Node(neigh.y + j, neigh.x + i);
                p0 = wmm::NodeD(h.y*n0.y, h.x*n0.x);
                if (wmm::norm(p0 - pn) <= h_radius && is_considered(state, n0)) {
                    y0 = u_surface.at(n0);
                    candidate_val = GetInterpValue(p0, p0, pn, y0, y0, c, 0.0);
                    val = std::min(val, candidate_val);
                    for (int k = 0; k < 8; k++) {
                        n1 = n0 + wmm::Node(wmm::yarray[k], wmm::xarray[k]);
                        if (is_considered(state, n1)) {
                            p1 = wmm::NodeD(h.y*n1.y, h.x*n1.x);
                            y1 = u_surface.at(n1);
                            dd = p1 - p0;
                            switch (mode) {
                                case wmm::M_GRADIENT:
                                    epsilon = GetEpsilonGradient(dd, p0, pn, c);
                                    break;
                                case wmm::M_HOPFLAX:
                                    epsilon = GetEpsilonHopfLax(dd, p0, pn, y0, y1);
                                    break;
                                default: //M_GOLDENSEARCH
                                    epsilon = GetEpsilonGoldenSearch(p0, p1, pn, y0, y1, c);
                            }
                            candidate_val = GetInterpValue(p0, p1, pn, y0, y1, c, epsilon);
                            val = std::min(val, candidate_val);
                        }
                    }
                }
            }
        }

        return val;

    }


    void OUMSurface2D(wmm::Grid &image, std::vector<wmm::Node> &initials, wmm::NodeD &h, double radius, int mode,
                      wmm::Grid &u_surface) {

        bool isnewpos[8];
        double valcenter[8];

        wmm::Grid_<unsigned char> state = wmm::Grid_<unsigned char>(image.rows, image.cols);

        std::multimap<double, wmm::Node > trial_set;
        std::map<int, std::multimap<double, wmm::Node >::iterator> mapa_trial;        

        int key, i, j, i_radius = (int) radius;
        wmm::Node winner, neigh;

        // Initialization
        for (i = 0; i < (int) initials.size(); i++) {
            key = initials[i].y * u_surface.cols + initials[i].x;
            if (mapa_trial.find(key) == mapa_trial.end() && u_surface.contains(initials[i])) {
                u_surface.at(initials[i]) = 0.0;
                state.at(initials[i]) = wmm::P_TRIAL;
                auto pr_trial = std::pair<double, wmm::Node >(0.0, initials[i]);
                auto trial_set_it = trial_set.insert(pr_trial);
                auto pr_mapa = std::pair<int, std::multimap<double, wmm::Node >::iterator>(key, trial_set_it);
                mapa_trial.insert(pr_mapa);
            }
        }

        while (!trial_set.empty()) {

            auto trial_set_it = trial_set.begin();
            key = trial_set_it->second.y * u_surface.cols + trial_set_it->second.x;
            auto mapa_trial_it = mapa_trial.find(key);

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

            state.at(winner) = wmm::P_ALIVE;
            
            // Neighbour temptative value computation
            for (i = 0; i < 8; i++) {
                neigh = winner + wmm::Node(wmm::yarray[i], wmm::xarray[i]);
                isnewpos[i] = false;
                if (u_surface.contains(neigh)) {
                    valcenter[i] = u_surface.at(neigh);
                }
                else {
                    valcenter[i] = wmm::MAX_VAL;
                }
                if (u_surface.contains(neigh) && state.at(neigh) != wmm::P_ALIVE) {
                    double val_neigh = GetVal2D(image, u_surface, state, neigh, h, radius, mode);
                    if (val_neigh < valcenter[i]) {
                        valcenter[i] = val_neigh;
                        isnewpos[i] = true;
                    }
                }
            }

            // Update
            for (i = 0; i < 8; i++) {
                if (isnewpos[i]) {
                    neigh = winner + wmm::Node(wmm::yarray[i], wmm::xarray[i]);
                    key = neigh.y * u_surface.cols + neigh.x;
                    if (state.at(neigh) == wmm::P_TRIAL) {
                        mapa_trial_it = mapa_trial.find(key);
                        trial_set.erase(mapa_trial_it->second);
                        mapa_trial.erase(mapa_trial_it);
                    }
                    else {
                        state.at(neigh) = wmm::P_TRIAL;
                    }

                    auto pr_trial = std::pair<double, wmm::Node >(valcenter[i], neigh);
                    auto trial_set_it = trial_set.insert(pr_trial);
                    auto pr_mapa = std::pair<int, std::multimap<double, wmm::Node >::iterator>(key, trial_set_it);
                    mapa_trial.insert(pr_mapa);

                    u_surface.at(neigh) = valcenter[i];
                }
            }
            
            for (i = -i_radius; i <= i_radius; i++) {
                for (j = -i_radius; j <= i_radius; j++) {
                    if (i == 0 && j == 0)
                        continue;
                    neigh = wmm::Node(winner.y + j, winner.x + i);
                    if (state.contains(neigh) && state.at(neigh) == wmm::P_TRIAL) {
                        double old_val = u_surface.at(neigh);
                        double new_val = GetVal2D(image, u_surface, state, neigh, h, radius, mode);
                        if (new_val < old_val) {
                            key = neigh.y * u_surface.cols + neigh.x;
                            mapa_trial_it = mapa_trial.find(key);
                            trial_set.erase(mapa_trial_it->second);
                            mapa_trial.erase(mapa_trial_it);

                            auto pr_trial = std::pair<double, wmm::Node >(new_val, neigh);
                            auto trial_set_it = trial_set.insert(pr_trial);
                            auto pr_mapa = std::pair<int, std::multimap<double, wmm::Node >::iterator>(key, trial_set_it);
                            mapa_trial.insert(pr_mapa);

                            u_surface.at(neigh) = new_val;
                        }
                    }
                }
            }

        }

        free(state.data);

        return;

    }


    wmm::Grid OUMSurface2D(wmm::Grid &image, std::vector<wmm::Node> &initials, wmm::NodeD &h, double radius, int mode) {
        wmm::Grid u_surface = wmm::Grid(wmm::MAX_VAL, image.rows, image.cols, 1);
        OUMSurface2D(image, initials, h, radius, mode, u_surface);
        return u_surface;
    }
    
    
}
