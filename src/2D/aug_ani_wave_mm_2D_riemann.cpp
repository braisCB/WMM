#include "../TYPES/WMMStructs.h"
#include "../TYPES/utils.h"
#include <math.h>
#include <stdlib.h>
#include <vector>
#include <map>
#include <unordered_map>
#include <iostream>


namespace wmm {


    Node3D GetFValue(Grid &image, AugWmm &wave, int ndir, int interp, double epsilon, NodeD &h) {

        Node3D ft;

        Node p1 = wave.p + wave.planes[ndir].dir;

        Node3D f0 = (image.channels == 3) ? Node3D(image.at(wave.p, 0), image.at(wave.p, 1), image.at(wave.p, 2)) :
                                            Node3D(image.at(wave.p, 0), image.at(wave.p, 1), 1.0), f1;
        if (image.contains(p1)) {
            f1 = (image.channels == 3) ? Node3D(image.at(p1, 0), image.at(p1, 1), image.at(p1, 2)) :
                                         Node3D(image.at(p1, 0), image.at(p1, 1), 1.0);
        }
        else {
            f1 = f0;
        }

        if (interp == I_LINEAR) {
            ft = (1.0 - epsilon)*f0 + epsilon*f1;
        }
        else if (interp == I_QUADATRIC) {
            ft = wave.planes[ndir].fm[0] + epsilon*wave.planes[ndir].fm[1] + epsilon*epsilon*wave.planes[ndir].fm[2];
        }
        else if (interp == I_SPLINE) {
            ft = (1.0 - epsilon) * f0 + epsilon * f1 +  epsilon * (1.0 - epsilon) *
                 ((1.0 - epsilon) * wave.planes[ndir].fm[0] + epsilon * wave.planes[ndir].fm[1]);
        }
        else {
            double t_2 = epsilon * epsilon;
            double t_3 = epsilon * t_2;
            ft = (2.0 * t_3 - 3.0 * t_2 + 1.0) * f0 + (t_3 - 2.0 * t_2 + epsilon) * wave.planes[ndir].fm[0] +
                 (-2.0 * t_3 + 3.0 * t_2) * f1 + (t_3 - t_2) * wave.planes[ndir].fm[1];
        }
        if ((ft.z < f0.z && ft.z < f1.z) || (ft.z > f0.z && ft.z > f1.z) || (ft.y < f0.y && ft.y < f1.y) ||
                (ft.y > f0.y && ft.y > f1.y) || (ft.x > f0.x && ft.x > f1.x) || (ft.x < f0.x && ft.x < f1.x)) {
            ft = (1.0 - epsilon)*f0 + epsilon*f1;
        }

        return ft;
    }

    double GetUValue(Grid &image, AugWmm &wave, int ndir, int interp, double epsilon_t, NodeD &h) {

        double value = MAX_VAL;
        double gamma = wave.gamma;

        int segment_1 = (epsilon_t >= 1.0) ? gamma - 1 : (int) floor(gamma * epsilon_t);

        double step = 1.0/ gamma;
        double epsilon = (epsilon_t - segment_1 * step) / step;

        if (interp == I_LINEAR) {
            value = (1.0 - epsilon)*wave.planes[ndir].v[segment_1] + epsilon*wave.planes[ndir].v[segment_1 + 1];
        }
        else if (interp == I_QUADATRIC) {
            double m[3];
            m[0] = wave.planes[ndir].v[segment_1];
            int dx = (segment_1 > 0) ? -1 : 2;
            m[2] = (wave.planes[ndir].v[segment_1 + dx] + (dx - 1.) * m[0] - dx * wave.planes[ndir].v[segment_1 + 1]) / (dx * (dx - 1.));
            m[1] = wave.planes[ndir].v[segment_1 + 1] - m[0] - m[2];
            value = m[0] + epsilon * m[1] + epsilon * epsilon * m[2];
        }
        else if (interp == I_HERMITE) {
            double m[2];
            double t = epsilon, t_2 = t*t, t_3 = t_2*t;
            if (segment_1 == 0) {
                double d0 = wave.planes[ndir].v[segment_1 + 1] - wave.planes[ndir].v[segment_1];
                double d1 = wave.planes[ndir].v[segment_1 + 2] - wave.planes[ndir].v[segment_1 + 1];
                m[0] = checkHermiteSlope(d0, d1);
            }
            else {
                m[0] = (wave.planes[ndir].v[segment_1 + 1] - wave.planes[ndir].v[segment_1 - 1]) / 2.;
            }
            if (segment_1 == (gamma - 1)) {
                double d1 = wave.planes[ndir].v[segment_1 + 1] - wave.planes[ndir].v[segment_1];
                double d0 = wave.planes[ndir].v[segment_1] - wave.planes[ndir].v[segment_1 - 1];
                m[1] = checkHermiteSlope(d1, d0);
            }
            else {
                m[1] = (wave.planes[ndir].v[segment_1 + 2] - wave.planes[ndir].v[segment_1]) / 2.;
            }
            /*m[0] = (segment_1 == 0) ? wave.planes[ndir].v[segment_1 + 1] - wave.planes[ndir].v[segment_1] :
                                      (wave.planes[ndir].v[segment_1 + 1] - wave.planes[ndir].v[segment_1 - 1]) / 2.;
            m[1] = (segment_1 == (gamma - 1)) ? wave.planes[ndir].v[segment_1 + 1] - wave.planes[ndir].v[segment_1] :
                                                (wave.planes[ndir].v[segment_1 + 2] - wave.planes[ndir].v[segment_1]) / 2.;*/
            value = (2.0 * t_3 - 3.0 * t_2 + 1.0) * wave.planes[ndir].v[segment_1] + (t_3 - 2.0 * t_2 + t) * m[0] +
                    (-2.0 * t_3 + 3.0 * t_2) * wave.planes[ndir].v[segment_1 + 1] + (t_3 - t_2) * m[1];
        }
        else if (interp == I_PCHIP) {
            double m[2], dv0, dv1;
            double t = epsilon, t_2 = t*t, t_3 = t_2*t;
            if (segment_1 == 0) {
                dv0 = wave.planes[ndir].v[segment_1 + 1] - wave.planes[ndir].v[segment_1];
                dv1 = wave.planes[ndir].v[segment_1 + 2] - wave.planes[ndir].v[segment_1 + 1];
            }
            else {
                dv0 = wave.planes[ndir].v[segment_1] - wave.planes[ndir].v[segment_1 - 1];
                dv1 = wave.planes[ndir].v[segment_1 + 1] - wave.planes[ndir].v[segment_1];
            }
            m[0] = checkSlope(dv0, dv1, (segment_1 == 0));
            if (segment_1 == (gamma - 1)) {
                dv0 = wave.planes[ndir].v[segment_1] - wave.planes[ndir].v[segment_1 - 1];
                dv1 = wave.planes[ndir].v[segment_1 + 1] - wave.planes[ndir].v[segment_1];
            }
            else {
                dv0 = wave.planes[ndir].v[segment_1 + 1] - wave.planes[ndir].v[segment_1];
                dv1 = wave.planes[ndir].v[segment_1 + 2] - wave.planes[ndir].v[segment_1 + 1];
            }
            m[1] = checkSlope(dv1, dv0, (segment_1 == (gamma - 1)));
            value = (2.0 * t_3 - 3.0 * t_2 + 1.0) * wave.planes[ndir].v[segment_1] + (t_3 - 2.0 * t_2 + t) * m[0] +
                    (-2.0 * t_3 + 3.0 * t_2) * wave.planes[ndir].v[segment_1 + 1] + (t_3 - t_2) * m[1];
        }
        else if (interp == I_SPLINE) {
            double m[8], v0, v1, y[4];
            double t = epsilon;
            int init_2 = (segment_1 == 0) ? 0 : -1;
            int end_2 = (segment_1 == (gamma - 1)) ? 1 : 2;

            for (int p2=init_2; p2<=end_2; p2++) {
                y[p2 - init_2] = wave.planes[ndir].v[segment_1 + p2];
            }
            get_spline_coeffs(y, m, end_2 - init_2);
            v0 = wave.planes[ndir].v[segment_1];
            v1 = wave.planes[ndir].v[segment_1 + 1];
            value = (1.0 - t) * v0 + t * v1 + t * (1.0 - t) * ((1.0 - t) * m[-2*init_2] + t * m[-2*init_2 + 1]);
        }
        if (((value < wave.planes[ndir].v[segment_1] && value < wave.planes[ndir].v[segment_1 + 1]) ||
            (value > wave.planes[ndir].v[segment_1] && value > wave.planes[ndir].v[segment_1 + 1])) &&
            interp != I_PCHIP && interp != I_LINEAR) {
                value = GetUValue(image, wave, ndir, I_PCHIP, epsilon_t, h);
        }

        return value;
    }


    double GetInterpValue(Grid &image, AugWmm &wave, int ndir,
                          int interp, double epsilon, NodeD &neighD, Node3D &fn,
                          NodeD &h) {

        NodeD plane_posD(wave.p.y + epsilon*wave.planes[ndir].dir.y, wave.p.x + epsilon*wave.planes[ndir].dir.x);
        Node neigh_f((int) floor(neighD.y), (int) floor(neighD.x));
        Node neigh_c((int) ceil(neighD.y), (int) ceil(neighD.x));
        double eps = neighD.y - neigh_f.y + neighD.x - neigh_f.x;


        double a_00 = (1. - eps) * image.at(neigh_f, 0) + eps * image.at(neigh_c, 0);
        double a_01 = (1. - eps) * image.at(neigh_f, 1) + eps * image.at(neigh_c, 1);
        double a_10 = (1. - eps) * image.at(neigh_f, 2) + eps * image.at(neigh_c, 2);
        double a_11 = (1. - eps) * image.at(neigh_f, 3) + eps * image.at(neigh_c, 3);

        double value = MAX_VAL;
        // std::cout << "iniciando F " << interp << std::endl;
        // NodeD f0 = GetFValue(image, wave, ndir, interp, epsilon, h);
        // std::cout << "iniciando U: " << ft.x << ", " << ft.y << ", " << ft.z << std::endl;
        double vt = GetUValue(image, wave, ndir, interp, epsilon, h);
        // std::cout << "finish him!" << std::endl;

        if (vt < wave.v && vt < wave.planes[ndir].v[wave.gamma] && interp != I_PCHIP && interp != I_LINEAR) {
            if (interp == I_SPLINE || interp == I_HERMITE) {
                vt = GetUValue(image, wave, ndir, I_PCHIP, epsilon, h);
            }
            else {
                vt = GetUValue(image, wave, ndir, I_LINEAR, epsilon, h);
            }
        }
        // double ft = GetIntegralValue(image, plane_posD, neighD, norm(f0), norm(fn), interp, wave.gamma - 1);
        NodeD diff = h * (neighD - plane_posD);
        NodeD a_diff(a_00 * diff.y + a_01 * diff.x, a_10 * diff.y + a_11 * diff.x);
        double dist = sqrt(diff.y * a_diff.y + diff.x * a_diff.x);
        value = vt + dist;
        //value = vt + norm(h * (neighD - plane_posD)) * (norm(f0) + norm(fn))/2.0;
        /*if (fabs(neighD.y - 17.25) < 0.01 && (fabs(neighD.x - 4) < 0.01 || fabs(neighD.x - 96) < 0.01)) {
            std::cout << "value: " << neighD.y << ", " << neighD.x << " : " << vt << ", eps : " << epsilon <<
            ", d = " << norm(neighD - plane_posD) << ", v = " << value << std::endl;
        }*/

        return value;
    }


    double GetEpsilonGradient(AugWmm &wave, int ndir, NodeD &neigh, Node3D &c, NodeD &h) {

        NodeD fn = NodeD(-c.x, c.y);
        NodeD ph = h * wave.p, nh = h * neigh;
        NodeD d0h = h*wave.planes[ndir].dir;

        double A_1 = d0h.y;
        double B_1 = -d0h.x;
        double C_1 = d0h.x * ph.y - d0h.y * ph.x;

        double A_2 = fn.y;
        double B_2 = -fn.x;
        double C_2 = fn.x * nh.y - fn.y * nh.x;

        double den = A_1*B_2 - A_2*B_1;
        NodeD x = NodeD(A_2 * C_1 - A_1 * C_2, B_1 * C_2 - B_2 * C_1) / den;

        double epsilon;

        if (fabs(den) == 0.0)
            epsilon = 0.0;
        else {
            NodeD diff = (x - ph) / h;
            epsilon = diff.x * wave.planes[ndir].dir.x + diff.y * wave.planes[ndir].dir.y;
        }

        if (epsilon < 0.0)
            epsilon = 0.0;
        else if (epsilon > 1.0)
            epsilon = 1.0;

        return epsilon;

    }

    /*double GetEpsilonHopfLax(AugWmm &wave, int ndir, NodeD &neigh, NodeD &h) {

        NodeD ph = h * wave.p, nh = h * neigh;
        NodeD d0h = h * wave.planes[ndir].dir;

        NodeD inc0, wave_p(wave.p.y, wave.p.x);
        NodeD ori = neigh - wave_p;
        int half = std::max(0, wave.gamma / 2 - 1);
        double inc_v0 = (wave.planes[ndir].v[half] - wave.planes[ndir].v[half + 2]) * wave.gamma / 2.0;

        if (d0h.y == 0) {
            inc0 = NodeD(ori.y*inc_v0, d0h.x);
        }
        else {
            inc0 = NodeD(d0h.y, ori.x*inc_v0);
        }
        NodeD fn = NodeD(-inc0.x, inc0.y);

        double A_1 = d0h.y;
        double B_1 = -d0h.x;
        double C_1 = d0h.x * ph.y - d0h.y * ph.x;

        double A_2 = fn.y;
        double B_2 = -fn.x;
        double C_2 = fn.x * nh.y - fn.y * nh.x;

        double den = A_1*B_2 - A_2*B_1;
        NodeD x = NodeD(A_2 * C_1 - A_1 * C_2, B_1 * C_2 - B_2 * C_1) / den;

        double epsilon;

        if (fabs(den) == 0.0)
            epsilon = 0.0;
        else {
            NodeD diff = (x - ph) / h;
            epsilon = diff.x * wave.planes[ndir].dir.x + diff.y * wave.planes[ndir].dir.y;
        }

        if (epsilon < 0.0)
            epsilon = 0.0;
        else if (epsilon > 1.0)
            epsilon = 1.0;

        return epsilon;

    }*/

    /*double GetEpsilonHopfLax(NodeD &wave_p, NodeD &pd, NodeD dir, double v[2], NodeD &neigh, NodeD &h) {

        NodeD ph = h * pd, nh = h * neigh;
        NodeD d0h = h * dir;

        NodeD inc0, inc1;
        NodeD ori = neigh - wave_p;
        double inc_v0 = v[0] - v[1];
        if (d0h.x == 0.0) {
            inc0 = NodeD(d0h.y, ori.x*inc_v0);
        }
        else {
            inc0 = NodeD(ori.y*inc_v0, d0h.x);
        }
        NodeD fn = NodeD(-inc0.x, inc0.y);

        double A_1 = d0h.y;
        double B_1 = -d0h.x;
        double C_1 = d0h.x * ph.y - d0h.y * ph.x;

        double A_2 = fn.y;
        double B_2 = -fn.x;
        double C_2 = fn.x * nh.y - fn.y * nh.x;

        double den = A_1*B_2 - A_2*B_1;
        NodeD x = NodeD(A_2 * C_1 - A_1 * C_2, B_1 * C_2 - B_2 * C_1) / den;

        double epsilon;

        if (fabs(den) == 0.0)
            epsilon = 0.0;
        else {
            NodeD diff = (x - ph) / (d0h.x + d0h.y);
            epsilon = diff.x * dir.x + diff.y * dir.y;
        }

        if (epsilon < 0.0)
            epsilon = 0.0;
        else if (epsilon > 1.0)
            epsilon = 1.0;

        return epsilon;

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

    double GetEpsilonGoldenSearch(Grid &image, AugWmm &wave, int ndir, int interp, NodeD &neigh, Node3D &fn, NodeD &h) {

        double br = (sqrt(5.0) + 1.0)/2.0;
        double a_1 = 0.0, b_1 = 1.0, x1_1 = b_1 - (b_1 - a_1) / br, x2_1 = a_1 + (b_1 - a_1) / br,
                f_x1 = MAX_VAL, f_x2 = MAX_VAL;

        double epsilon = 0.0;

        double res = GetInterpValue(image, wave, ndir, interp, 0.0, neigh, fn, h);
        double f_b = GetInterpValue(image, wave, ndir, interp, 1.0, neigh, fn, h);
        if (f_b < res) {
            res = f_b; epsilon = 1.0;
        }

        while (fabs(b_1 - a_1) > TAU) {
            f_x1 = GetInterpValue(image, wave, ndir, interp, x1_1, neigh, fn, h);
            f_x2 = GetInterpValue(image, wave, ndir, interp, x2_1, neigh, fn, h);
            if(f_x1 < f_x2) {
                b_1 = x2_1;
            }
            else {
                a_1 = x1_1;
            }
            x1_1 = b_1 - (b_1 - a_1) / br, x2_1 = a_1 + (b_1 - a_1) / br;
        }

        /*if (fabs(neigh.y - 17.25) < 0.01 && (fabs(neigh.x - 4) < 0.01 || fabs(neigh.x - 96) < 0.01)) {
            std::cout << "wvalue: " << neigh.y << ", " << neigh.x << " : res = " << res << ", f_x1 : " << f_x1 <<
            ", f_x2 = " << f_x2 << std::endl;
        }*/

        if (f_x1 < res) {
            epsilon = x1_1;
            res = f_x1;
        }
        if (f_x2 < res) {
            epsilon = x2_1;
        }

        /*if (fabs(neigh.y - 17.25) < 0.01 && (fabs(neigh.x - 4) < 0.01 || fabs(neigh.x - 96) < 0.01)) {
            std::cout << "epsilon: " << epsilon << std::endl;
        }*/

        return epsilon;

    }

    double GetVal2D(Grid &image, Grid &u_surface, AugWmm &wave, NodeD &neigh, Node3D &fn, NodeD &h, int interp, int mode) {
        double y0 = wave.v;

        Node neigh_f((int) floor(neigh.y), (int) floor(neigh.x));
        Node neigh_c((int) ceil(neigh.y), (int) ceil(neigh.x));
        double epsilon = neigh.y - neigh_f.y + neigh.x - neigh_f.x;


        double a_00 = (1. - epsilon) * image.at(neigh_f, 0) + epsilon * image.at(neigh_c, 0);
        double a_01 = (1. - epsilon) * image.at(neigh_f, 1) + epsilon * image.at(neigh_c, 1);
        double a_10 = (1. - epsilon) * image.at(neigh_f, 2) + epsilon * image.at(neigh_c, 2);
        double a_11 = (1. - epsilon) * image.at(neigh_f, 3) + epsilon * image.at(neigh_c, 3);

        NodeD wave_p(wave.p.y, wave.p.x);
        NodeD diff = h * (neigh - wave_p);
        NodeD a_diff(a_00 * diff.y + a_01 * diff.x, a_10 * diff.y + a_11 * diff.x);
        double dist = sqrt(diff.y * a_diff.y + diff.x * a_diff.x);
        double val = y0 + dist;

        // double ft = GetIntegralValue(image, wave_p, neigh, norm(f0), norm(fn), interp, wave.gamma - 1);
        // double val = y0 + norm(diff) * ft;
        //double val = y0 + norm(diff) * (norm(f0) + norm(fn)) / 2.0;

        if (wave.ndirs > 0) {

            for (int ndir = 0; ndir < wave.ndirs; ndir++) {

                Node p0 = wave.p + wave.planes[ndir].dir;

                if (!image.contains(p0))
                    continue;

                double epsilon, epsilon_aux;
                if (mode == M_GRADIENT)
                    epsilon = GetEpsilonGradient(wave, ndir, neigh, fn, h);
                else if (mode == M_HOPFLAX) {
                    //epsilon = GetEpsilonHopfLax(wave, ndir, neigh, h);
                    double res, step = 1.0 / (double) wave.gamma, res_aux;
                    // double v[2];
                    //NodeD dd_aux = step * wave.planes[ndir].dir, dp_aux;
                    res = MAX_VAL;
                    epsilon = 0.0;
                    for (int i=0; i<wave.gamma; i++) {
                        double y0 = wave.planes[ndir].v[i];
                        double y1 = wave.planes[ndir].v[i+1];
                        NodeD dd = step * h * wave.planes[ndir].dir;
                        NodeD dp = h * wave_p + i*dd;
                        NodeD dn = h * neigh;
                        epsilon_aux = GetEpsilonHopfLax(dd, dp, dn, y0, y1);
                        /*v[0] = wave.planes[ndir].v[i];
                        v[1] = wave.planes[ndir].v[i+1];
                        dp_aux = wave_p + i*dd_aux;
                        epsilon_aux = GetEpsilonHopfLax(wave_p, dp_aux, dd_aux, v, neigh, h);*/
                        epsilon_aux = (i + epsilon_aux)*step;
                        res_aux = GetInterpValue(image, wave, ndir, interp, epsilon_aux, neigh, fn, h);
                        if (res_aux < res) {
                            res = res_aux;
                            epsilon = epsilon_aux;
                        }
                    }
                }
                else
                    epsilon = GetEpsilonGoldenSearch(image, wave, ndir, interp, neigh, fn, h);
                // std::cout << "test " << neigh.y << ", " << neigh.x << " - " << ndir << " : " << epsilon << " - "<< GetInterpValue(image, wave, ndir, interp, epsilon, neigh, fn, h) << std::endl;
                val = std::min(val, GetInterpValue(image, wave, ndir, interp, epsilon, neigh, fn, h));
            }
        }
        return val;

    }


    void WmmAugAniSurface2D(Grid &image, std::vector<Node> &initials, NodeD &h,
                            int interp, int mode, int N, int gamma, Grid &u_surface) {
        bool isnewpos[8], wascomputed[8];
        double valcenter[8], value;
        Node3D imcenter[8];
        Node dirs[2], neigh, neigh_aux;
        int ndirs = 2;
        double step = 1.0 / (double) gamma, epsilon1;

        Grid_<unsigned char> state = Grid_<unsigned char>(image.rows, image.cols);

        std::multimap<double, AugWmm > trial_set;
        std::unordered_map<int, std::multimap<double, AugWmm >::iterator> mapa_trial;
        std::unordered_map<int, double> augmented_values;

        std::multimap<double, AugWmm >::iterator trial_set_it;
        std::unordered_map<int, std::multimap<double, AugWmm >::iterator>::iterator mapa_trial_it;
        std::pair<double, AugWmm > pr_trial;
        std::pair<int, std::multimap<double, AugWmm >::iterator> pr_mapa;

        int key, i, aug_key;
        AugWmm winner;
        NodeD neighD, d1, neighD_aux;
        Node3D fi;

        // Initialization
        for (i = 0; i < (int) initials.size(); i++) {
            key = initials[i].y + u_surface.rows*initials[i].x;
            if (mapa_trial.find(key) == mapa_trial.end() && u_surface.contains(initials[i])) {
                u_surface.at(initials[i]) = 0.0;
                winner = AugWmm(N, gamma, 0);
                winner.v = 0.0;
                winner.p = initials[i];
                state.at(initials[i]) = P_TRIAL;
                pr_trial = std::pair<double, AugWmm >(0.0, winner);
                trial_set_it = trial_set.insert(pr_trial);
                pr_mapa = std::pair<int, typename std::multimap<double, AugWmm >::iterator>(key, trial_set_it);
                mapa_trial.insert(pr_mapa);
            }
        }
        while (!trial_set.empty()) {

            trial_set_it = trial_set.begin();
            key = trial_set_it->second.p.y + u_surface.rows*trial_set_it->second.p.x;
            mapa_trial_it = mapa_trial.find(key);

            if (mapa_trial_it == mapa_trial.end()) {
                printf("ERROR: bad map alloc");
                return;
            }

            if (mapa_trial_it->second != trial_set_it) {
                printf("ERROR: bad trial/map alloc");
                return;
            }

            winner = trial_set_it->second;
            double score = trial_set_it->first;

            trial_set.erase(trial_set_it);
            mapa_trial.erase(mapa_trial_it);

            state.at(winner.p) = (score < winner.v) ? P_INTERMEDIATE : P_ALIVE;

            // std::cout << "winner: " << winner.p.y << ", " << winner.p.x << " - " << winner.v << std::endl;
            if (winner.ndirs > 0) {
                ndirs = 2;
                neigh = winner.p;
                neighD = NodeD(neigh.y, neigh.x);
                for (int m=0; m < ndirs; m++) {
                    d1 = NodeD(winner.planes[m].dir.y, winner.planes[m].dir.x);
                    // std::cout << "m: " << m << " : " << d1.y << ", " << d1.x << std::endl;
                    for (int v=0; v <= gamma; v++) {
                        epsilon1 = v*step;
                        neighD_aux = neighD + epsilon1*d1;
                        if (!image.contains(neighD_aux)) {
                            winner.planes[m].v[v] = MAX_VAL;
                        }
                        else {
                            neigh_aux = Node(round(gamma * neighD_aux.y), round(gamma * neighD_aux.x));
                            aug_key = neigh_aux.y + gamma * u_surface.rows*neigh_aux.x;
                            winner.planes[m].v[v] = (augmented_values.find(aug_key) == augmented_values.end()) ? MAX_VAL : augmented_values[aug_key];
                            // std::cout << aug_key << " : " << m << " - " << v << " - " << neigh_aux.y << ", " << neigh_aux.x << " : " << winner.planes[m].v[v] << std::endl;
                        }
                    }
                }
            }


            // Neighbour temptative value computation
            for (int t=0; t < 8; t++) {
                isnewpos[t] = false;
                wascomputed[t] = false;
                neigh = winner.p + Node(yarray[t], xarray[t]);
                if (!u_surface.contains(neigh)) {
                    fi = Node3D(MAX_VAL, MAX_VAL, MAX_VAL);
                    // fi = Node3D(MAX_VAL, MAX_VAL, MAX_VAL);
                    valcenter[t] = MAX_VAL;
                    continue;
                }
                fi = (image.channels == 3) ? Node3D(image.at(neigh, 0), image.at(neigh, 1), image.at(neigh, 2)) :
                                             Node3D(image.at(neigh, 0), image.at(neigh, 1), 1.0);
                valcenter[t] = u_surface.at(neigh);
                if (isinf(norm(fi)) || isnan(norm(fi)))
                    fi = Node3D(MAX_VAL, MAX_VAL, MAX_VAL);
                imcenter[t] = fi;
                if (u_surface.contains(neigh) && state.at(neigh) != P_ALIVE) {
                    neighD = NodeD(neigh.y, neigh.x);
                    double val_neigh = GetVal2D(image, u_surface, winner, neighD, fi, h, interp, mode);
                    // std::cout << neighD.y << ", " << neighD.x << " : " << val_neigh << " - " << valcenter[t] << std::endl;
                    if (val_neigh < valcenter[t]) {
                        valcenter[t] = val_neigh;
                        isnewpos[t] = true;
                    }
                }
            }

            // Update
            for (int t=0; t < 8; t++) {
                int i = yarray[t] + 1, j = xarray[t] + 1;
                if (!isnewpos[t])
                    continue;

                Node d(i - 1, j - 1);
                neigh = winner.p + d;

                key = neigh.y + u_surface.rows*neigh.x;
                mapa_trial_it = mapa_trial.find(key);
                if (mapa_trial_it != mapa_trial.end()) {
                    trial_set.erase(mapa_trial_it->second);
                    mapa_trial.erase(mapa_trial_it);
                }
                state.at(neigh) = (state.at(neigh) == P_INTERMEDIATE) ? P_INTERMEDIATE : P_TRIAL;

                AugWmm new_w(N, gamma, 2);
                new_w.p = neigh;
                dirs[0] = Node(yarray[(t+1)%8], xarray[(t+1)%8]) - d;
                dirs[1] = Node(yarray[(t+7)%8], xarray[(t+7)%8]) - d;
                new_w.v = valcenter[t];
                double v_border[2] = {valcenter[(t+1)%8], valcenter[(t+7)%8]}, score = new_w.v;
                for (int m=0; m < 2; m++) {
                    new_w.planes[m].dir = dirs[m];
                    int t_aux = (t + (m == 0) - (m == 1) + 8) % 8;
                    if (wascomputed[t_aux])
                        continue;
                    neighD = NodeD(neigh.y, neigh.x);
                    d1 = NodeD(dirs[m].y, dirs[m].x);
                    for (int v=0; v <= gamma; v++) {
                        epsilon1 = v*step;
                        neighD_aux = neighD + epsilon1*d1;
                        if (!image.contains(neighD_aux)) {
                            continue;
                        }
                        neigh_aux = Node(round(gamma * neighD_aux.y), round(gamma * neighD_aux.x));
                        aug_key = neigh_aux.y + gamma * u_surface.rows*neigh_aux.x;
                        if (v == 0)
                            value = valcenter[t];
                        else if (v == gamma)
                            value = v_border[m];
                        else {
                            if (u_surface.contains(neighD_aux)) {
                                fi = GetFValue(image, new_w, m, interp, epsilon1, h);
                                value = GetVal2D(image, u_surface, winner, neighD_aux, fi, h, interp, mode);
                            }
                            else {
                                value = MAX_VAL;
                            }
                        }
                        if (augmented_values.find(aug_key) == augmented_values.end() || (value < augmented_values[aug_key])) {
                            augmented_values[aug_key] = value;
                        }
                        if (value < score) score = value;
                    }
                }
                wascomputed[t] = false;
                score = (state.at(neigh) == P_INTERMEDIATE) ? new_w.v : score;
                pr_trial = std::pair<double, AugWmm >(score, new_w);
                trial_set_it = trial_set.insert(pr_trial);
                pr_mapa = std::pair<int, std::multimap<double, AugWmm >::iterator>(key, trial_set_it);
                mapa_trial.insert(pr_mapa);

                u_surface.at(new_w.p) = valcenter[t];
            }

            del(&winner);
        }

        free(state.data);
        return;

    }

    Grid WmmAugAniSurface2D(Grid &image, std::vector<Node> &initials, NodeD &h, int interp, int mode,
                              int N, int gamma) {
        Grid u_surface = Grid(MAX_VAL, image.rows, image.cols, 1);
        WmmAugAniSurface2D(image, initials, h, interp, mode, N, gamma, u_surface);
        return u_surface;
    }

    
}
