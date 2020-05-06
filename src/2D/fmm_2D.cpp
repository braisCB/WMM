#include "../TYPES/WMMStructs.h"
#include <math.h>
#include <stdlib.h>
#include <vector>
#include <map>
#include <iostream>


namespace wmm {

    double StencilS1(wmm::Grid_<unsigned char> &state, wmm::Grid &u_surface, wmm::Node &neigh, wmm::NodeD &h, double F, int order) {

        double T_2a = wmm::MAX_VAL, T_2b = wmm::MAX_VAL, T_1a = wmm::MAX_VAL, T_1b = wmm::MAX_VAL, gh1, gh2;
        double T_21a = wmm::MAX_VAL, T_21b = wmm::MAX_VAL, T_11a = wmm::MAX_VAL, T_11b = wmm::MAX_VAL;
        bool T_1n = false, T_2n = false;
        wmm::Node p;

        gh1 = 1.0/(h.y*h.y);
        p = wmm::Node(neigh.y - 1, neigh.x);
        if (state.contains(p)) {
            if (state.at(p) == wmm::P_ALIVE) {
                T_1a = u_surface.at(p);
                if (order == 2) {
                    p = wmm::Node(neigh.y - 2, neigh.x);
                    if (state.contains(p) && state.at(p) == wmm::P_ALIVE) {
                        T_11a = u_surface.at(p);
                        if (T_1a > T_11a)
                            T_1n = true;
                    }
                }
            }
        }
        p = wmm::Node(neigh.y + 1, neigh.x);
        if (state.contains(p)) {
            if (state.at(p) == wmm::P_ALIVE) {
                T_1b = u_surface.at(p);
                if (order == 2) {
                    p = wmm::Node(neigh.y + 2, neigh.x);
                    if (state.contains(p) && state.at(p) == wmm::P_ALIVE) {
                        T_11b = u_surface.at(p);
                        if (T_1b > T_11b) {
                            T_1n = true;
                        }
                    }
                }
            }
        }

        gh2 = 1.0/(h.x*h.x);
        p = wmm::Node(neigh.y, neigh.x + 1);
        if (state.contains(p)) {
            if (state.at(p) == wmm::P_ALIVE) {
                T_2a = u_surface.at(p);
                if (order == 2) {
                    p = wmm::Node(neigh.y, neigh.x + 2);
                    if (state.contains(p) && state.at(p) == wmm::P_ALIVE) {
                        T_21a = u_surface.at(p);
                        if (T_2a > T_21a)
                            T_2n = true;
                    }
                }
            }
        }
        p = wmm::Node(neigh.y, neigh.x - 1);
        if (state.contains(p)) {
            if (state.at(p) == wmm::P_ALIVE) {
                T_2b = u_surface.at(p);
                if (order == 2) {
                    p = wmm::Node(neigh.y, neigh.x - 2);
                    if (state.contains(p) && state.at(p) == wmm::P_ALIVE) {
                        T_21b = u_surface.at(p);
                        if (T_2b > T_21b) {
                            T_2n = true;
                        }
                    }
                }
            }
        }

        int real_order = (T_1n && T_2n) ? 2 : 1;
        double T_1, T_2;
        switch (real_order) {
            case 1:
                T_1 = std::min(T_1a, T_1b);
                T_2 = std::min(T_2a, T_2b);
                break;
            case 2:
                if (T_1a >= T_11a && T_1b >= T_11b)
                    T_1 = std::min((4.0*T_1a - T_11a)/3.0, (4.0*T_1b - T_11b)/3.0);
                else if (T_1a >= T_11a)
                    T_1 = (4.0*T_1a - T_11a)/3.0;
                else
                    T_1 = (4.0*T_1b - T_11b)/3.0;
                if (T_2a >= T_21a && T_2b >= T_21b)
                    T_2 = std::min((4.0*T_2a - T_21a)/3.0, (4.0*T_2b - T_21b)/3.0);
                else if (T_2a >= T_21a)
                    T_2 = (4.0*T_2a - T_21a)/3.0;
                else
                    T_2 = (4.0*T_2b - T_21b)/3.0;
                gh1 *= 9.0/4.0;
                gh2 *= 9.0/4.0;
                break;
        }

        if (T_1 == wmm::MAX_VAL && T_2 == wmm::MAX_VAL)
            return wmm::MAX_VAL;

        double new_val;

        double a_1 = gh1 + gh2;
        double b_1 = -2.0*(T_1*gh1 + T_2*gh2);
        double c_1 = T_1*T_1*gh1 + T_2*T_2*gh2 - 1.0/(F*F);
        double disc = b_1*b_1 - 4.0*a_1*c_1;
        new_val = -1.0;
        if (disc >= 0.0) {
            new_val = (-b_1 + sqrt(disc))/(2.0*a_1);
            if (new_val < std::min(T_1, T_2))
                new_val = wmm::MAX_VAL;
        } 
        if ((new_val <= T_1 && T_1 < wmm::MAX_VAL) || (new_val <= T_2 && T_2 < wmm::MAX_VAL) || new_val >= wmm::MAX_VAL) {
            double s1 = T_1 + 1.0/(F*sqrt(gh1)), s2 = T_2 + 1.0/(F*sqrt(gh2));
            new_val = std::min(s1, s2);
        }

        return new_val;

    }

    double GetVal2D(wmm::Grid_<unsigned char> &state, wmm::Grid &u_surface, wmm::Node &neigh, wmm::NodeD &h, double F, int order) {

        double stencilS1_val = wmm::StencilS1(state, u_surface, neigh, h, F, order);
        return stencilS1_val;

    }

    void FMMIsoSurface2D(wmm::Grid &image, std::vector<wmm::Node> &initials, wmm::NodeD &h, int order, wmm::Grid &u_surface) {

        wmm::Grid_<unsigned char> state = wmm::Grid_<unsigned char>(image.rows, image.cols);

        std::multimap<double, wmm::Node > trial_set;
        std::map<int, std::multimap<double, wmm::Node >::iterator> mapa_trial;

        std::multimap<double, wmm::Node >::iterator trial_set_it;
        std::map<int, std::multimap<double, wmm::Node >::iterator>::iterator mapa_trial_it;
        std::pair<double, wmm::Node > pr_trial;
        std::pair<int, std::multimap<double, wmm::Node >::iterator> pr_mapa;

        int key, i;
        wmm::Node winner, neigh;
        
        double valcenter[8];
        bool isnewpos[8];

        // Initialization
        for (i = 0; i < (int) initials.size(); i++) {
            key = initials[i].y * u_surface.cols + initials[i].x;
            if (mapa_trial.find(key) == mapa_trial.end() && u_surface.contains(initials[i])) {
                u_surface.at(initials[i]) = 0.0;
                state.at(initials[i]) = wmm::P_TRIAL;
                pr_trial = std::pair<double, wmm::Node>(0.0, initials[i]);
                trial_set_it = trial_set.insert(pr_trial);
                pr_mapa = std::pair<int, std::multimap<double, wmm::Node >::iterator>(key, trial_set_it);
                mapa_trial.insert(pr_mapa);
            }
        }

        while (!trial_set.empty()) {

            trial_set_it = trial_set.begin();
            key = trial_set_it->second.y * u_surface.cols + trial_set_it->second.x;
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

            state.at(winner) = wmm::P_ALIVE;

            // Neighbour temptative value computation
            for (i = 0; i < 8; i++) {
                neigh = winner + wmm::Node(wmm::yarray[i], wmm::xarray[i]);
                isnewpos[i] = false;
                valcenter[i] = u_surface.contains(neigh) ? u_surface.at(neigh) : wmm::MAX_VAL;
                if (u_surface.contains(neigh) && state.at(neigh) != wmm::P_ALIVE) {
                    double imcenter = u_surface.contains(neigh) ? 1.0/wmm::norm(wmm::NodeD(image.at(neigh, 0), image.at(neigh, 1))) : wmm::MAX_VAL;
                    double val_neigh = GetVal2D(state, u_surface, neigh, h, imcenter, order);
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

                    pr_trial = std::pair<double, wmm::Node >(valcenter[i], neigh);
                    trial_set_it = trial_set.insert(pr_trial);
                    pr_mapa = std::pair<int, std::multimap<double, wmm::Node >::iterator>(key, trial_set_it);
                    mapa_trial.insert(pr_mapa);

                    u_surface.at(neigh) = valcenter[i];
                }
            }

        }

        free(state.data);

        return;

    }


    wmm::Grid FMMIsoSurface2D(wmm::Grid &image, std::vector<wmm::Node> &initials, wmm::NodeD &h, int order) {
        wmm::Grid u_surface = wmm::Grid(wmm::MAX_VAL, image.rows, image.cols, 1);
        FMMIsoSurface2D(image, initials, h, order, u_surface);
        return u_surface;
    }
    
}
