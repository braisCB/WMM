#include "../TYPES/WMMStructs.h"
#include <math.h>
#include <stdlib.h>
#include <vector>
#include <map>
#include <iostream>


namespace wmm {

    template<typename _Tp> struct TimeWmm0_ {
        Node p;
        Node3_<_Tp> v[3];
        int dir;
        int angle;
    };

    template<typename _Tp, int N> struct TimeWmm_ : public TimeWmm0_<_Tp> {
        Node3_<_Tp> m[N];
    };

    template<int N> Node3D GetInterpValue(TimeWmm_<double, N> &wave, NodeD &dp, NodeD &dd,
                                               NodeD &dn, Node3D &F,
                                               double epsilon, int interp, int side) {

        Node3D value, y0 = wave.v[0], y1 = wave.v[side + 1];

        switch (interp) {
            case I_LINEAR:
                value = (1.0 - epsilon)*y0 + epsilon*y1;
                break;
            case I_QUADATRIC:
                value = wave.m[2*side + 2]*epsilon*epsilon + wave.m[2*side + 1]*epsilon + wave.m[0];
                break;
            case I_SPLINE:  
                value = y0 + epsilon*(-1. * wave.m[2*side + 1]/6.0 - wave.m[2*side]/3.0 + y1 - y0 +
                        epsilon*(wave.m[2*side]/2.0 + epsilon*(wave.m[2*side + 1] - wave.m[2*side])/6.0));
                break;
            default: //I_HERMITE and I_PCHIP
                double t_2 = epsilon * epsilon;
                double t_3 = epsilon * t_2;
                value = (2.0 * t_3 - 3.0 * t_2 + 1.0) * y0 + (t_3 - 2.0 * t_2 + epsilon) * wave.m[2*side] +
                        (-2.0 * t_3 + 3.0 * t_2) * y1 + (t_3 - t_2) * wave.m[2*side + 1];

        }

        if ((value.x < y0.x) && (value.x < y1.x))
            value = (1.0 - epsilon)*y0 + epsilon*y1;
        value = value + norm(dn - dp - epsilon*dd) / F;

        return value;
    }

    double GetEpsilonHopfLax(NodeD &dd, NodeD &dp, NodeD &dn, double y0, double y1) {

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

    template<int N> double GetEpsilonGoldenSearch(TimeWmm_<double, N> &wave, NodeD &dn, NodeD &dd,
                                                  NodeD &dp, Node3D &F, double interp, int side) {

        double a = 0.0, b = 1.0, x1 = a + (1-RESPHI)*(b - a), x2 = a + RESPHI*(b - a);

        NodeD xtreme = dp + dd;

        double epsilon;
        Node3D f_a = wave.v[0] + norm(dn - dp) / F, f_x1, f_x2, res;
        Node3D f_b = wave.v[1 + side] + norm(dn - xtreme) / F;

        if (f_a.x < f_b.x) {
            res = f_a; epsilon = 0.0;
        }
        else {
            res = f_b; epsilon = 1.0;
        }

        f_x1 = GetInterpValue(wave, dp, dd, dn, F, x1, interp, side);
        f_x2 = GetInterpValue(wave, dp, dd, dn, F, x2, interp, side);

        while (fabs(b - a) > TAU) {
            if(f_x1.x < f_x2.x) {
                b = x2; x2 = x1; f_x2 = f_x1; x1 = a + (1 - RESPHI)*(b - a);
                f_x1 = GetInterpValue(wave, dp, dd, dn, F, x1, interp, side);
            }
            else {
                a = x1; x1 = x2; f_x1 = f_x2; x2 = a + RESPHI*(b - a);
                f_x2 = GetInterpValue(wave, dp, dd, dn, F, x2, interp, side);
            }
        }

        if (f_x1.x < std::min(res.x, f_x2.x))
            epsilon = x1;
        else if (f_x2.x < std::min(res.x, f_x1.x))
            epsilon = x2;

        return epsilon;
        
    }


    template<int N> Node3D GetVal2D(Grid_<unsigned char> &state, Grid &u_surface, Grid &d_surface,
                                    Grid &t_surface, TimeWmm_<double, N> &wave, Node &neigh,
                                    NodeD &h, double freq, double vel, int interp, int mode, double factor) {

        Node3D y0 = wave.v[0];

        Node3D F = Node3D((1. - factor) * freq + factor * vel, 1., vel);

        Node3D val = Node3D(MAX_VAL, MAX_VAL, MAX_VAL);
        if (wave.dir < 0) {
            NodeD diff(h.y * (neigh.y - wave.p.y), h.x * (neigh.x - wave.p.x));
            val = y0 + norm(diff) / F;
        }
        else {
            Node p(wave.p.y + yarray[(wave.dir + 1) % 8] - yarray[wave.dir],
                    wave.p.x + xarray[(wave.dir + 1) % 8] - xarray[wave.dir]);
            Node3D res1 = Node3D(MAX_VAL, MAX_VAL, MAX_VAL);
            NodeD dp(h.y * wave.p.y, h.x * wave.p.x), dn(h.y * neigh.y, h.x * neigh.x);

            if (u_surface.contains(p)) {
                Node3D y1 = wave.v[1];

                NodeD dd(h.y * (yarray[(wave.dir + 1) % 8] - yarray[wave.dir]), h.x * (xarray[(wave.dir + 1) % 8] - xarray[wave.dir]));

                double epsilon;
                switch (mode) {
                    case M_HOPFLAX:
                        epsilon = GetEpsilonHopfLax(dd, dp, dn, y0.x, y1.x);
                        break;
                    default: //M_GOLDENSEARCH
                        epsilon = GetEpsilonGoldenSearch(wave, dn, dd, dp, F, interp, S_RIGHT);
                }

                res1 = GetInterpValue(wave, dp, dd, dn, F, epsilon, interp, S_RIGHT);
                //std::cout << "RIGHT: " << wave.p.y << ", " << wave.p.x << " - " << epsilon << " . " << neigh.y << ", " << neigh.x << " -> " << res1.x << std::endl;
            }

            p = Node(wave.p.y + yarray[(wave.dir + 7) % 8] - yarray[wave.dir],
                     wave.p.x + xarray[(wave.dir + 7) % 8] - xarray[wave.dir]);
            Node3D res2 = Node3D(MAX_VAL, MAX_VAL, MAX_VAL);

            if (u_surface.contains(p)) {
                Node3D y2 = wave.v[2];

                NodeD dd(h.y * (yarray[(wave.dir + 7) % 8] - yarray[wave.dir]), h.x * (xarray[(wave.dir + 7) % 8] - xarray[wave.dir]));

                double epsilon;
                switch (mode) {
                    case M_HOPFLAX:
                        epsilon = GetEpsilonHopfLax(dd, dp, dn, y0.x, y2.x);
                        break;
                    default: //M_GOLDENSEARCH
                        epsilon = GetEpsilonGoldenSearch(wave, dn, dd, dp, F, interp, S_LEFT);
                }

                res2 = GetInterpValue(wave, dp, dd, dn, F, epsilon, interp, S_LEFT);
                //std::cout << "LEFT: " << wave.p.y << ", " << wave.p.x << " - " << epsilon << " . " << neigh.y << ", " << neigh.x << " -> " << res2.x << std::endl;

            }

            val = (res1.x < res2.x) ? res1 : res2;

        }

        return val;

    }



    void setCoeffs2D(Node3D *y, Node3D *m, int pos, int interp) {

        switch (interp) {
            case I_LINEAR:
                break;
            case I_QUADATRIC:
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
            case I_SPLINE:
                if (pos%2 == 0) {
                    m[0] = m[2] = Node3D(0., 0., 0.);
                    m[1] = 6.0*(y[pos] - 2.0*y[(pos+1)%8] + y[(pos+2)%8])/4.0;
                    m[3] = 6.0*(y[pos] - 2.0*y[(pos+7)%8] + y[(pos+6)%8])/4.0;
                }
                else {
                    m[0] = m[2] = 6.0*(y[(pos+7)%8] - 2.0*y[pos] + y[(pos+1)%8])/4.0;
                    m[1] = m[3] = Node3D(0., 0., 0.);
                }
                break;
            case I_HERMITE:
                if (pos%2 == 0) {
                    m[0] = y[(pos + 1) % 8] - y[pos];
                    m[1] = 1.0 / 2.0 * (y[(pos + 2) % 8] - y[pos]);
                    m[2] = y[(pos + 7) % 8] - y[pos];
                    m[3] = 1.0 / 2.0 * (y[(pos + 6) % 8] - y[pos]);
                }
                else {
                    m[0] = 1.0 / 2.0 * (y[(pos + 1) % 8] - y[(pos + 7) % 8]);
                    m[1] = y[(pos + 1) % 8] - y[pos];
                    m[2] = -1. * m[0];
                    m[3] = y[(pos + 7) % 8] - y[pos];
                }
                break;
            default: //I_PCHIP
                if (pos%2 == 0) {
                    Node3D d0 = y[(pos + 1) % 8] - y[pos], d1 = y[(pos + 2) % 8] - y[(pos + 1) % 8];
                    m[0] = (3.0 * d0 - d1) / 2.0;
                    if ((m[0].x * d0.x) <= 0.0)
                        m[0].x = 0.;
                    else if (((d1.x * d0.x) <= 0.0) && (fabs(m[0].x) > fabs(3.0 * d0.x)))
                        m[0].x = 3.0 * d0.x;
                    if ((m[0].y * d0.y) <= 0.0)
                        m[0].y = 0.;
                    else if (((d1.y * d0.y) <= 0.0) && (fabs(m[0].y) > fabs(3.0 * d0.y)))
                        m[0].y = 3.0 * d0.y;
                    if ((m[0].z * d0.z) <= 0.0)
                        m[0].z = 0.;
                    else if (((d1.z * d0.z) <= 0.0) && (fabs(m[0].z) > fabs(3.0 * d0.z)))
                        m[0].z = 3.0 * d0.z;
                    m[1].x = (d0.x * d1.x <= 0.0) ? 0. : 2.0 * d0.x * d1.x / (d0.x + d1.x);
                    m[1].y = (d0.y * d1.y <= 0.0) ? 0. : 2.0 * d0.y * d1.y / (d0.y + d1.y);
                    m[1].z = (d0.z * d1.z <= 0.0) ? 0. : 2.0 * d0.z * d1.z / (d0.z + d1.z);

                    d0 = y[(pos + 7) % 8] - y[pos]; d1 = y[(pos + 6) % 8] - y[(pos + 7) % 8];
                    m[2] = (3.0 * d0 - d1) / 2.0;
                    if ((m[2].x * d0.x) <= 0.0)
                        m[2].x = 0.;
                    else if (((d1.x * d0.x) <= 0.0) && (fabs(m[2].x) > fabs(3.0 * d0.x)))
                        m[2].x = 3.0 * d0.x;
                    if ((m[2].y * d0.y) <= 0.0)
                        m[2].y = 0.;
                    else if (((d1.y * d0.y) <= 0.0) && (fabs(m[2].y) > fabs(3.0 * d0.y)))
                        m[2].y = 3.0 * d0.y;
                    if ((m[2].z * d0.z) <= 0.0)
                        m[2].z = 0.;
                    else if (((d1.z * d0.z) <= 0.0) && (fabs(m[2].z) > fabs(3.0 * d0.z)))
                        m[2].z = 3.0 * d0.z;
                    m[3].x = (d0.x * d1.x <= 0.0) ? 0. : 2.0 * d0.x * d1.x / (d0.x + d1.x);
                    m[3].y = (d0.y * d1.y <= 0.0) ? 0. : 2.0 * d0.y * d1.y / (d0.y + d1.y);
                    m[3].z = (d0.z * d1.z <= 0.0) ? 0. : 2.0 * d0.z * d1.z / (d0.z + d1.z);
                }
                else {
                    Node3D d0 = y[pos] - y[(pos + 7) % 8], d1 = y[(pos + 1) % 8] - y[pos];
                    m[0].x = (d0.x * d1.x <= 0.0) ? 0.0 : 2.0 * d0.x * d1.x / (d0.x + d1.x);
                    m[0].y = (d0.y * d1.y <= 0.0) ? 0.0 : 2.0 * d0.y * d1.y / (d0.y + d1.y);
                    m[0].z = (d0.z * d1.z <= 0.0) ? 0.0 : 2.0 * d0.z * d1.z / (d0.z + d1.z);
                    m[1] = (3.0 * d1 - d0) / 2.0;
                    if ((m[1].x * d1.x) <= 0.0)
                        m[1].x = 0.;
                    else if (((d1.x * d0.x) <= 0.0) && (fabs(m[1].x) > fabs(3.0 * d1.x)))
                        m[1].x = 3.0 * d1.x;
                    if ((m[1].y * d1.y) <= 0.0)
                        m[1].y = 0.;
                    else if (((d1.y * d0.y) <= 0.0) && (fabs(m[1].y) > fabs(3.0 * d1.y)))
                        m[1].y = 3.0 * d1.y;
                    if ((m[1].z * d1.z) <= 0.0)
                        m[1].z = 0.;
                    else if (((d1.z * d0.z) <= 0.0) && (fabs(m[1].z) > fabs(3.0 * d1.z)))
                        m[1].z = 3.0 * d1.z;

                    d0 = y[pos] - y[(pos + 1) % 8]; d1 = y[(pos + 7) % 8] - y[pos];
                    m[2].x = (d0.x * d1.x <= 0.0) ? 0.0 : 2.0 * d0.x * d1.x / (d0.x + d1.x);
                    m[2].y = (d0.y * d1.y <= 0.0) ? 0.0 : 2.0 * d0.y * d1.y / (d0.y + d1.y);
                    m[2].z = (d0.z * d1.z <= 0.0) ? 0.0 : 2.0 * d0.z * d1.z / (d0.z + d1.z);
                    m[3] = (3.0 * d1 - d0) / 2.0;
                    if ((m[3].x * d1.x) <= 0.0)
                        m[3].x = 0.;
                    else if (((d1.x * d0.x) <= 0.0) && (fabs(m[3].x) > fabs(3.0 * d1.x)))
                        m[3].x = 3.0 * d1.x;
                    if ((m[3].y * d1.y) <= 0.0)
                        m[3].y = 0.;
                    else if (((d1.y * d0.y) <= 0.0) && (fabs(m[3].y) > fabs(3.0 * d1.y)))
                        m[3].y = 3.0 * d1.y;
                    if ((m[3].z * d1.z) <= 0.0)
                        m[3].z = 0.;
                    else if (((d1.z * d0.z) <= 0.0) && (fabs(m[3].z) > fabs(3.0 * d1.z)))
                        m[3].z = 3.0 * d1.z;
                }
        }

    }


    template<int wmm_degree> void WMMTimeSurface2D(Grid &frequency, Grid &velocity,
                                                   std::vector<Node> &initials, NodeD &h, int interp,
                                                   int mode, int backwards, int min_conn, double factor,
                                                   Grid &u_surface, Grid &d_surface, Grid &t_surface) {

        bool isnewpos[8];
        Node3D valcenter[8];

        Grid_<unsigned char> state = Grid_<unsigned char>(frequency.rows, frequency.cols);

        std::multimap<double, TimeWmm_<double, wmm_degree> > trial_set;
        std::map<int, typename std::multimap<double, TimeWmm_<double, wmm_degree> >::iterator> mapa_trial;

        typename std::multimap<double, TimeWmm_<double, wmm_degree> >::iterator trial_set_it;
        typename std::map<int, typename std::multimap<double, TimeWmm_<double, wmm_degree> >::iterator>::iterator mapa_trial_it;
        std::pair<double, TimeWmm_<double, wmm_degree> > pr_trial;
        std::pair<int, typename std::multimap<double, TimeWmm_<double, wmm_degree> >::iterator> pr_mapa;

        int key, i;
        TimeWmm_<double, wmm_degree> winner, new_w;
        Node neigh;

        // Initialization
        for (i = 0; i < (int) initials.size(); i++) {
            key = initials[i].y * u_surface.cols + initials[i].x;
            if (mapa_trial.find(key) == mapa_trial.end() && u_surface.contains(initials[i])) {
                u_surface.at(initials[i]) = 0.0;
                d_surface.at(initials[i]) = 0.0;
                t_surface.at(initials[i]) = 0.0;
                winner.dir = -1;
                winner.v[0] = Node3D(0., 0., 0.);
                winner.p = initials[i];
                state.at(initials[i]) = P_TRIAL;
                pr_trial = std::pair<double, TimeWmm_<double, wmm_degree> >(0.0, winner);
                trial_set_it = trial_set.insert(pr_trial);
                pr_mapa = std::pair<int, typename std::multimap<double, TimeWmm_<double, wmm_degree> >::iterator>(key, trial_set_it);
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

            /*std::cout << "WINNER: " << winner.p.y << ", " << winner.p.x << " - " << u_surface.at(winner.p) << " -.- "
                      << winner.v[0].y << ", "
                      << winner.v[1].y << ", " << winner.v[2].y << " -.- "
                      << winner.m[0].y << ", " << winner.m[1].y << ", " << winner.m[2].y << ", " << winner.m[3].y
                      << std::endl;*/

            // Neighbour temptative value computation
            for (i = 0; i < 8; i++) {
                neigh = winner.p + Node(yarray[i], xarray[i]);
                isnewpos[i] = false;
                if (!u_surface.contains(neigh)) {
                    valcenter[i] = Node3D(MAX_VAL, MAX_VAL, MAX_VAL);
                    continue;
                }
                valcenter[i] = Node3D(u_surface.at(neigh), d_surface.at(neigh), t_surface.at(neigh));
                if (state.at(neigh) == P_ALIVE)
                    continue;
                double freq = MAX_VAL, vel = MAX_VAL;
                if (frequency.channels == 1) {
                    freq = frequency.at(neigh);
                    vel = velocity.at(neigh);
                }
                else {
                    freq = (backwards) ? frequency.at(neigh, (i+4) % 8) : frequency.at(winner.p, i);
                    vel = (backwards) ? velocity.at(neigh, (i+4) % 8) : velocity.at(winner.p, i);
                }
                if (freq <= min_conn)
                    continue;
                Node3D val_neigh = GetVal2D(state, u_surface, d_surface, t_surface, winner, neigh, h, freq, vel, interp,
                                            mode, factor);
                if ((val_neigh.x - u_surface.at(neigh)) < -TAU) {
                    valcenter[i] = val_neigh;
                    isnewpos[i] = true;
                }
            }

            for (i = 0; i < 8; i+=2) {
                if (valcenter[i].x == MAX_VAL) {
                    Node n = winner.p + Node(yarray[(i + 2) % 8], xarray[(i + 2) % 8]);
                    if (frequency.contains(n)) {
                        valcenter[i] = valcenter[(i + 2) % 8];
                    }
                    else {
                        valcenter[i] = valcenter[(i + 6) % 8];
                    }
                    /*if (valcenter[(i + 2) % 8].x == MAX_VAL && valcenter[(i + 6) % 8].x == MAX_VAL)
                        valcenter[i] = winner.v[0];
                    else if (valcenter[(i + 2) % 8].x == MAX_VAL)
                        valcenter[i] = 1.5*valcenter[(i + 7) % 8] - valcenter[(i + 6) % 8];
                    else
                        valcenter[i] = 1.5*valcenter[(i + 1) % 8] - valcenter[(i + 2) % 8];*/
                }
            }

            // Update
            for (i = 0; i < 8; i++) {
                //std::cout << winner.p.y + yarray[i] << ", " << winner.p.x + xarray[i] << ", valcenter[" << i << "] = " << valcenter[i].x << std::endl;
                if (isnewpos[i]) {
                    neigh = winner.p + Node(yarray[i], xarray[i]);
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
                    new_w.v[1] = valcenter[(i+1)%8];
                    new_w.v[2] = valcenter[(i+7)%8];

                    setCoeffs2D(valcenter, new_w.m, i, interp);

                    /*std::cout << "UPDATE: " << new_w.p.y << ", " << new_w.p.x << " - " << new_w.v[0].x << ", "
                      << new_w.v[1].x << ", " << new_w.v[2].x << " -.- "
                      << new_w.m[0].x << ", " << new_w.m[1].x << ", " << new_w.m[2].x << ", " << new_w.m[3].x
                      << std::endl;*/

                    pr_trial = std::pair<double, TimeWmm_<double, wmm_degree> >(new_w.v[0].x, new_w);
                    trial_set_it = trial_set.insert(pr_trial);
                    pr_mapa = std::pair<int, typename std::multimap<double, TimeWmm_<double, wmm_degree> >::iterator>(key, trial_set_it);
                    mapa_trial.insert(pr_mapa);

                    u_surface.at(new_w.p) = valcenter[i].x;
                    d_surface.at(new_w.p) = valcenter[i].y;
                    t_surface.at(new_w.p) = valcenter[i].z;
                }
            }

        }

        free(state.data);

        return;

    }
    
    
}
