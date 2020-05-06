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
    
    

    double GetGoldenSearchValue(wmm::NodeD &p0, wmm::NodeD &p1, wmm::NodeD &pn,
                              double y0, double y1, wmm::NodeD &c) {
        
        double a = 0.0, b = 1.0, x1 = a + (1-RESPHI)*(b - a), x2 = a + RESPHI*(b - a),
                f_x1 = MAX_VAL, f_x2 = MAX_VAL, i1, i2, res;


        wmm::NodeD F_x1, F_x2;

        double f_a = GetInterpValue(p0, p1, pn, y0, y1, c, 0.0);
        double f_b = GetInterpValue(p0, p1, pn, y0, y1, c, 1.0);

        if (f_a < f_b) {
            res = f_a;
        }
        else {
            res = f_b;
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
            res = f_x1;
        else if (f_x2 < std::min(res, f_x1))
            res = f_x2;

        return res;
        
    }

    double GetVal2D(wmm::Grid_<unsigned char> &state, wmm::Grid &u_surface, wmm::Node &neigh, wmm::NodeD &h, wmm::NodeD &c, int mode) {

        wmm::Node node0, node1;
        wmm::NodeD p0, p1, pn = wmm::NodeD(h.y * neigh.y, h.x * neigh.x);
        double y0, y1;
        double res = wmm::MAX_VAL;
        int i;
        for (i = 0; i < 8; i++) {
            node0 = neigh + wmm::Node(wmm::yarray[i], wmm::xarray[i]);
            if (u_surface.contains(node0) && state.at(node0) == wmm::P_ALIVE) {
                p0 = wmm::NodeD(h.y * node0.y, h.x * node0.x);
                y0 = u_surface.at(node0);
                node1 = neigh + wmm::Node(wmm::yarray[(i + 1) % 8], wmm::xarray[(i + 1) % 8]);
                double value;
                if (u_surface.contains(node1) && state.at(node1) == wmm::P_ALIVE) {
                    p1 = wmm::NodeD(h.y * node1.y, h.x * node1.x);
                    y1 = u_surface.at(node1);
                    switch (mode) {
                        default: //M_GOLDENSEARCH
                            value = GetGoldenSearchValue(p0, p1, pn, y0, y1, c);
                    }
                }
                else {
                    value = GetInterpValue(p0, p0, pn, y0, y0, c, 0.0);
                }
                if (value < res) res = value;
            }
        }
        
        for (i = 1; i < 8; i+=2) {
            node0 = neigh + wmm::Node(wmm::yarray[i], wmm::xarray[i]);
            node1 = neigh + wmm::Node(wmm::yarray[(i + 2) % 8], wmm::xarray[(i + 2) % 8]);
            if (u_surface.contains(node0) && state.at(node0) == wmm::P_ALIVE &&
                   u_surface.contains(node1) && state.at(node1) == wmm::P_ALIVE) {
                p0 = wmm::NodeD(h.y * node0.y, h.x * node0.x);
                p1 = wmm::NodeD(h.y * node1.y, h.x * node1.x);
                y0 = u_surface.at(node0);
                y1 = u_surface.at(node1);
                switch (mode) {
                    default: //M_GOLDENSEARCH
                        double value = GetGoldenSearchValue(p0, p1, pn, y0, y1, c);
                        if (value < res) res = value;
                }
            }
        }
        
        return res;

    }

    void FMMAniSurface2D(wmm::Grid &image, std::vector<wmm::Node> &initials, wmm::NodeD &h, int mode, wmm::Grid &u_surface) {

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
                    wmm::NodeD imcenter = u_surface.contains(neigh) ? wmm::NodeD(image.at(neigh, 0), image.at(neigh, 1)) : wmm::NodeD(wmm::MAX_VAL, wmm::MAX_VAL);
                    double val_neigh = GetVal2D(state, u_surface, neigh, h, imcenter, mode);
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


    wmm::Grid FMMAniSurface2D(wmm::Grid &image, std::vector<wmm::Node> &initials, wmm::NodeD &h, int mode) {
        wmm::Grid u_surface = wmm::Grid(wmm::MAX_VAL, image.rows, image.cols, 1);
        FMMAniSurface2D(&image, &initials, &h, mode, &u_surface);
        return u_surface;
    }

}
