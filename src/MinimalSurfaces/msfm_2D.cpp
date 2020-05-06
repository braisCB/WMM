#include "../TYPES/WMMStructs.h"
#include <math.h>
#include <stdlib.h>
#include <vector>
#include <map>
#include <iostream>


namespace wmm {

    const wmm::Node stencils[2][2] = { {wmm::Node(1, 0), wmm::Node(0, 1)}, {wmm::Node(1, 1), wmm::Node(-1, 1)}};

    wmm::Node3D GetVal2D(wmm::Grid_<unsigned char> &state, wmm::Grid &u_surface, wmm::Grid &d_surface, wmm::Grid &t_surface,
                         wmm::Node &node, wmm::NodeD &h, double freq, double vel, int order, double f_factor) {

        wmm::Node3D scores = wmm::Node3D(wmm::MAX_VAL, wmm::MAX_VAL, wmm::MAX_VAL);
        wmm::Node neigh;
        double a;
        wmm::Node3D stencil_scores, b, c, disc, F = wmm::Node3D((1. - f_factor) * freq + f_factor * vel, 1., vel);
        wmm::NodeD stencil_h, stencil_hs;
        int factor;

        for (int i=0; i < 2; i++) {
            wmm::Node3D scores_H = wmm::Node3D(wmm::MAX_VAL, wmm::MAX_VAL, wmm::MAX_VAL);
            wmm::Node3D scores_V = wmm::Node3D(wmm::MAX_VAL, wmm::MAX_VAL, wmm::MAX_VAL);
            wmm::Node3D scores_H2 = wmm::Node3D(wmm::MAX_VAL, wmm::MAX_VAL, wmm::MAX_VAL);
            wmm::Node3D scores_V2 = wmm::Node3D(wmm::MAX_VAL, wmm::MAX_VAL, wmm::MAX_VAL);
            bool V_pos = true, H_pos = true;

            neigh = node + stencils[i][0];
            if (state.contains(neigh) && state.at(neigh) == wmm::P_ALIVE) {
                double uval = u_surface.at(neigh);
                if (uval < scores_V.x)
                    scores_V = wmm::Node3D(uval, d_surface.at(neigh), t_surface.at(neigh));
            }
            neigh = node - stencils[i][0];
            if (state.contains(neigh) && state.at(neigh) == wmm::P_ALIVE) {
                double uval = u_surface.at(neigh);
                if (uval < scores_V.x) {
                    scores_V = wmm::Node3D(uval, d_surface.at(neigh), t_surface.at(neigh));
                    V_pos = false;
                }
            }
            neigh = node + stencils[i][1];
            if (state.contains(neigh) && state.at(neigh) == wmm::P_ALIVE) {
                double uval = u_surface.at(neigh);
                if (uval < scores_H.x)
                    scores_H = wmm::Node3D(uval, d_surface.at(neigh), t_surface.at(neigh));
            }
            neigh = node - stencils[i][1];
            if (state.contains(neigh) && state.at(neigh) == wmm::P_ALIVE) {
                double uval = u_surface.at(neigh);
                if (uval < scores_H.x) {
                    scores_H = wmm::Node3D(uval, d_surface.at(neigh), t_surface.at(neigh));
                    H_pos = false;
                }
            }

            // check order 2
            int real_order = ((scores_H.x == wmm::MAX_VAL) || (scores_V.x == wmm::MAX_VAL)) ? 1 : order;
            if (real_order == 2) {
                factor = (V_pos) ? 2 : -2;
                neigh = node + factor * stencils[i][0];
                if (state.contains(neigh) && state.at(neigh) == wmm::P_ALIVE) {
                    double uval = u_surface.at(neigh), dval = d_surface.at(neigh), tval = t_surface.at(neigh);
                    if (uval < scores_V.x && dval < scores_V.y && tval < scores_V.z)
                        scores_V2 = wmm::Node3D(uval, dval, tval);
                }
                factor = (H_pos) ? 2 : -2;
                neigh = node + factor * stencils[i][1];
                if (state.contains(neigh) && state.at(neigh) == wmm::P_ALIVE) {
                    double uval = u_surface.at(neigh), dval = d_surface.at(neigh), tval = t_surface.at(neigh);
                    if (uval < scores_H.x && dval < scores_H.y && tval < scores_H.z)
                        scores_H2 = wmm::Node3D(uval, dval, tval);
                }
                real_order = ((scores_H2.x == wmm::MAX_VAL) || (scores_V2.x == wmm::MAX_VAL)) ? 1 : order;
            }

            stencil_h = wmm::NodeD(
                wmm::norm(wmm::NodeD(stencils[i][0].y * h.y, stencils[i][0].x * h.x)), 
                wmm::norm(wmm::NodeD(stencils[i][1].y * h.y, stencils[i][1].x * h.x))
            );
            stencil_hs = wmm::NodeD(stencil_h.y * stencil_h.y, stencil_h.x * stencil_h.x);

            if (scores_V.x + stencil_h.y / F.x < scores_H.x + stencil_h.x / F.x)
                stencil_scores = scores_V + stencil_h.y / F;
            else
                stencil_scores = scores_H + stencil_h.x / F;

            if (scores_V.x < wmm::MAX_VAL && scores_H.x < wmm::MAX_VAL) {
                a = 1. / stencil_hs.y + 1. / stencil_hs.x;
                b = -2. * (scores_V / stencil_hs.y + scores_H / stencil_hs.x);
                c = scores_H * scores_H / stencil_hs.x + scores_V * scores_V / stencil_hs.y - 1. / (F * F);
                disc = b * b - 4. * a * c;
                if (disc.x > 0. and disc.y > 0. and disc.z > 0)
                    stencil_scores = (wmm::node_sqrt(disc) - b) / (2. * a);
                if (real_order == 2) {
                    scores_H = (4. * scores_H - scores_H2) / 3.;
                    scores_V = (4. * scores_V - scores_V2) / 3.;
                    a = 9. / 4. * (1. / stencil_hs.y + 1. / stencil_hs.x);
                    b = -9. / 2. * (scores_V / stencil_hs.y + scores_H / stencil_hs.x);
                    c = 9. / 4. * (scores_H * scores_H / stencil_hs.x + scores_V * scores_V / stencil_hs.y) - 1. / (F * F);
                    disc = b * b - 4. * a * c;
                    if (disc.x > 0. and disc.y > 0. and disc.z > 0)
                        stencil_scores = (wmm::node_sqrt(disc) - b) / (2. * a);
                }
            }

            if (stencil_scores.x < scores.x)
                scores = stencil_scores;

        }
        
        return scores;

    }


    void MSFMTimeSurface2D(wmm::Grid &frequency, wmm::Grid &velocity, std::vector<wmm::Node> &initials, wmm::NodeD &h,
                           int order, int backwards, int min_conn, double factor, wmm::Grid &u_surface,
                           wmm::Grid &d_surface, wmm::Grid &t_surface) {

        wmm::Grid_<unsigned char> state = wmm::Grid_<unsigned char>(frequency.rows, frequency.cols);

        std::multimap<double, wmm::Node > trial_set;
        std::map<int, std::multimap<double, wmm::Node >::iterator> mapa_trial;

        std::multimap<double, wmm::Node >::iterator trial_set_it;
        std::map<int, std::multimap<double, wmm::Node >::iterator>::iterator mapa_trial_it;
        std::pair<double, wmm::Node > pr_trial;
        std::pair<int, std::multimap<double, wmm::Node >::iterator> pr_mapa;

        int key, i;
        wmm::Node winner, neigh;
        
        wmm::Node3D valcenter[8];
        bool isnewpos[8];

        // Initialization
        for (i = 0; i < (int) initials.size(); i++) {
            key = initials[i].y * u_surface.cols + initials[i].x;
            if (mapa_trial.find(key) == mapa_trial.end() && u_surface.contains(initials[i])) {
                u_surface.at(initials[i]) = 0.0;
                d_surface.at(initials[i]) = 0.0;
                t_surface.at(initials[i]) = 0.0;
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
                if (!u_surface.contains(neigh) || state.at(neigh) == wmm::P_ALIVE)
                    continue;
                double freq = wmm::MAX_VAL, vel = wmm::MAX_VAL;
                if (u_surface.contains(neigh)) {
                    if (frequency.channels == 1) {
                        freq = frequency.at(neigh);
                        vel = velocity.at(neigh);
                    }
                    else {
                        freq = (backwards) ? frequency.at(neigh, (i+4) % 8) : frequency.at(winner, i);
                        vel = (backwards) ? velocity.at(neigh, (i+4) % 8) : velocity.at(winner, i);
                    }
                }
                if (freq <= min_conn)
                    continue;
                wmm::Node3D val_neigh = GetVal2D(state, u_surface, d_surface, t_surface, neigh, h, freq, vel, order, factor);
                if (val_neigh.x < u_surface.at(neigh)) {
                    valcenter[i] = val_neigh;
                    isnewpos[i] = true;
                }
            }

            // Update
            for (i = 0; i < 8; i++) {
                if (!isnewpos[i])
                    continue;
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

                pr_trial = std::pair<double, wmm::Node >(valcenter[i].x, neigh);
                trial_set_it = trial_set.insert(pr_trial);
                pr_mapa = std::pair<int, std::multimap<double, wmm::Node >::iterator>(key, trial_set_it);
                mapa_trial.insert(pr_mapa);

                u_surface.at(neigh) = valcenter[i].x;
                d_surface.at(neigh) = valcenter[i].y;
                t_surface.at(neigh) = valcenter[i].z;
            }

        }

        free(state.data);

        return;

    }
    
}
