#include "../TYPES/DistanceWMMStructs.h"
#include "../TYPES/utils.h"
#include <math.h>
#include <stdlib.h>
#include <vector>
#include <map>
#include <set>
#include <unordered_map>
#include <iostream>


namespace wmm {

    Node3D checkNodeSlope(Node3D &fd0, Node3D &fd1, int isborder) {
        Node3D node = Node3D(
            checkSlope(fd0.x, fd1.x, isborder),
            checkSlope(fd0.y, fd1.y, isborder),
            checkSlope(fd0.z, fd1.z, isborder)
        );
        return node;
    }

    Node3D checkHermiteSlope(Node3D &fd0, Node3D &fd1) {
        Node3D node = Node3D(
            checkHermiteSlope(fd0.x, fd1.x),
            checkHermiteSlope(fd0.y, fd1.y),
            checkHermiteSlope(fd0.z, fd1.z)
        );
        return node;
    }

    NodeD GetFValue(Grid &image, AugWmm &wave, int ndir, int interp, double epsilon, NodeD &h) {

        NodeD value;
        Node3D ft;

        Node p1 = wave.p + wave.planes[ndir].dir;

        NodeD f0_3 = NodeD(image.at(wave.p, 0), image.at(wave.p, 1));
        NodeD f1_3 = (image.contains(p1)) ? NodeD(image.at(p1, 0), image.at(p1, 1)) : f0_3;

        Node3D f0 = to_polar(f0_3);
        Node3D f1 = to_polar(f1_3);

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
        if ((ft.z < f0.z && ft.z < f1.z) || (ft.z > f0.z && ft.z > f1.z)) {
            ft = (1.0 - epsilon)*f0 + epsilon*f1;
        }
        value = to_cartesian(ft);

        return value;
    }

    std::pair<double, double> GetUValue(Grid &image, AugWmm &wave, int ndir, int interp, double epsilon_t, NodeD &h) {

        double value = MAX_VAL, dist = MAX_VAL;
        double gamma = wave.gamma;

        int segment_1 = (epsilon_t >= 1.0) ? gamma - 1 : (int) floor(gamma * epsilon_t);

        double step = 1.0/ gamma;
        double epsilon = (epsilon_t - segment_1 * step) / step;

        if (interp == I_LINEAR) {
            value = (1.0 - epsilon)*wave.planes[ndir].v[segment_1] + epsilon*wave.planes[ndir].v[segment_1 + 1];
            dist = (1.0 - epsilon)*wave.planes[ndir].d[segment_1] + epsilon*wave.planes[ndir].d[segment_1 + 1];
        }
        else if (interp == I_QUADATRIC) {
            double m[3];
            int dx = (segment_1 > 0) ? -1 : 2;
            m[0] = wave.planes[ndir].v[segment_1];
            m[2] = (wave.planes[ndir].v[segment_1 + dx] + (dx - 1.) * m[0] - dx * wave.planes[ndir].v[segment_1 + 1]) / (dx * (dx - 1.));
            m[1] = wave.planes[ndir].v[segment_1 + 1] - m[0] - m[2];
            value = m[0] + epsilon * m[1] + epsilon * epsilon * m[2];
            m[0] = wave.planes[ndir].d[segment_1];
            m[2] = (wave.planes[ndir].d[segment_1 + dx] + (dx - 1.) * m[0] - dx * wave.planes[ndir].d[segment_1 + 1]) / (dx * (dx - 1.));
            m[1] = wave.planes[ndir].d[segment_1 + 1] - m[0] - m[2];
            dist = m[0] + epsilon * m[1] + epsilon * epsilon * m[2];
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
            value = (2.0 * t_3 - 3.0 * t_2 + 1.0) * wave.planes[ndir].v[segment_1] + (t_3 - 2.0 * t_2 + t) * m[0] +
                    (-2.0 * t_3 + 3.0 * t_2) * wave.planes[ndir].v[segment_1 + 1] + (t_3 - t_2) * m[1];
            if (segment_1 == 0) {
                double d0 = wave.planes[ndir].d[segment_1 + 1] - wave.planes[ndir].d[segment_1];
                double d1 = wave.planes[ndir].d[segment_1 + 2] - wave.planes[ndir].d[segment_1 + 1];
                m[0] = checkHermiteSlope(d0, d1);
            }
            else {
                m[0] = (wave.planes[ndir].d[segment_1 + 1] - wave.planes[ndir].d[segment_1 - 1]) / 2.;
            }
            if (segment_1 == (gamma - 1)) {
                double d1 = wave.planes[ndir].d[segment_1 + 1] - wave.planes[ndir].d[segment_1];
                double d0 = wave.planes[ndir].d[segment_1] - wave.planes[ndir].d[segment_1 - 1];
                m[1] = checkHermiteSlope(d1, d0);
            }
            else {
                m[1] = (wave.planes[ndir].d[segment_1 + 2] - wave.planes[ndir].d[segment_1]) / 2.;
            }
            dist = (2.0 * t_3 - 3.0 * t_2 + 1.0) * wave.planes[ndir].d[segment_1] + (t_3 - 2.0 * t_2 + t) * m[0] +
                   (-2.0 * t_3 + 3.0 * t_2) * wave.planes[ndir].d[segment_1 + 1] + (t_3 - t_2) * m[1];
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
            if (segment_1 == 0) {
                dv0 = wave.planes[ndir].d[segment_1 + 1] - wave.planes[ndir].d[segment_1];
                dv1 = wave.planes[ndir].d[segment_1 + 2] - wave.planes[ndir].d[segment_1 + 1];
            }
            else {
                dv0 = wave.planes[ndir].d[segment_1] - wave.planes[ndir].d[segment_1 - 1];
                dv1 = wave.planes[ndir].d[segment_1 + 1] - wave.planes[ndir].d[segment_1];
            }
            m[0] = checkSlope(dv0, dv1, (segment_1 == 0));
            if (segment_1 == (gamma - 1)) {
                dv0 = wave.planes[ndir].d[segment_1] - wave.planes[ndir].d[segment_1 - 1];
                dv1 = wave.planes[ndir].d[segment_1 + 1] - wave.planes[ndir].d[segment_1];
            }
            else {
                dv0 = wave.planes[ndir].d[segment_1 + 1] - wave.planes[ndir].d[segment_1];
                dv1 = wave.planes[ndir].d[segment_1 + 2] - wave.planes[ndir].d[segment_1 + 1];
            }
            m[1] = checkSlope(dv1, dv0, (segment_1 == (gamma - 1)));
            dist = (2.0 * t_3 - 3.0 * t_2 + 1.0) * wave.planes[ndir].d[segment_1] + (t_3 - 2.0 * t_2 + t) * m[0] +
                   (-2.0 * t_3 + 3.0 * t_2) * wave.planes[ndir].d[segment_1 + 1] + (t_3 - t_2) * m[1];
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
            for (int p2=init_2; p2<=end_2; p2++) {
                y[p2 - init_2] = wave.planes[ndir].d[segment_1 + p2];
            }
            get_spline_coeffs(y, m, end_2 - init_2);
            v0 = wave.planes[ndir].d[segment_1];
            v1 = wave.planes[ndir].d[segment_1 + 1];
            dist = (1.0 - t) * v0 + t * v1 + t * (1.0 - t) * ((1.0 - t) * m[-2*init_2] + t * m[-2*init_2 + 1]);
        }
        if (((value < wave.planes[ndir].v[segment_1] && value < wave.planes[ndir].v[segment_1 + 1]) ||
            (value > wave.planes[ndir].v[segment_1] && value > wave.planes[ndir].v[segment_1 + 1])) &&
            ((dist < wave.planes[ndir].d[segment_1] && dist < wave.planes[ndir].d[segment_1 + 1]) ||
            (dist > wave.planes[ndir].d[segment_1] && dist > wave.planes[ndir].d[segment_1 + 1])) &&
            interp != I_PCHIP && interp != I_LINEAR) {
                return GetUValue(image, wave, ndir, I_PCHIP, epsilon_t, h);
        }

        return std::pair<double, double>(value, dist);
    }


    void setCoeffs2D(Node3D *y, Node3D *m, int pos, int interp, int right) {

        int ori = -2 * right + 1;
        if (interp == I_LINEAR)
            return;
        if (interp == I_QUADATRIC) {
            m[0] = y[pos];
            int dx = (pos % 2 == 1) ? -1 : 2;
            m[2] = (y[(pos + dx*ori + 8) % 8] + (dx - 1.) * m[0] - dx * y[(pos + ori + 8) % 8]) / (dx * (dx - 1.));
            m[1] = y[(pos + ori + 8) % 8] - m[0] - m[2];
        }
        else if (interp == I_HERMITE) {
            Node3D d0, d1;
            if (pos % 2 == 0) {
                d0 = y[(pos + ori + 8) % 8] - y[pos];
                d1 = y[(pos + 2 * ori + 8) % 8] - y[(pos + ori + 8) % 8];
                m[0] = checkHermiteSlope(d0, d1);
                m[1] = (y[(pos + ori + 8) % 8] - y[(pos - ori + 8) % 8]) / 2.;
            }
            else {
                d0 = y[pos] - y[(pos - ori + 8) % 8];
                d1 = y[(pos + ori + 8) % 8] - y[pos];
                m[0] = (y[(pos + ori + 8) % 8] - y[(pos - ori + 8) % 8]) / 2.;
                m[1] = checkHermiteSlope(d1, d0);
            }
            /*m[0] = (pos % 2 == 0) ? y[(pos + ori + 8) % 8] - y[pos] :
                                    (y[(pos + ori + 8) % 8] - y[(pos - ori + 8) % 8]) / 2.;
            m[1] = (pos % 2 == 1) ? y[(pos + ori + 8) % 8] - y[pos] :
                                    (y[(pos + ori + 8) % 8] - y[(pos - ori + 8) % 8]) / 2.;*/
        }
        else if (interp == I_PCHIP) {
            Node3D dv0, dv1;
            if (pos % 2 == 0) {
                dv0 = y[(pos + ori + 8) % 8] - y[pos];
                dv1 = y[(pos + 2 * ori + 8) % 8] - y[(pos + ori + 8) % 8];
            }
            else {
                dv0 = y[pos] - y[(pos - ori + 8) % 8];
                dv1 = y[(pos + ori + 8) % 8] - y[pos];
            }
            m[0] = checkNodeSlope(dv0, dv1, (pos % 2 == 0));
            m[1] = checkNodeSlope(dv1, dv0, (pos % 2 == 1));
        }
        else if (interp == I_SPLINE) {
            if (pos%2 == 0) {
                m[0] = 0.0;
                m[1] = 6.0*(y[pos] - 2.0*y[(pos + ori + 8)%8] + y[(pos + 2*ori + 8)%8])/4.0;
            }
            else {
                m[0] = 6.0*(y[(pos + ori + 8)%8] - 2.0*y[pos] + y[(pos - ori + 8)%8])/4.0;
                m[1] = 0.0;
            }
        }
    }


    std::pair<double, double> GetInterpValue(Grid &image, AugWmm &wave, int ndir,
                          int interp, double epsilon, NodeD &neighD, NodeD &fn,
                          NodeD &h) {

        NodeD plane_posD(wave.p.y + epsilon*wave.planes[ndir].dir.y, wave.p.x + epsilon*wave.planes[ndir].dir.x);

        double value = MAX_VAL, dist = MAX_VAL;
        // std::cout << "iniciando F " << interp << std::endl;
        NodeD f0 = GetFValue(image, wave, ndir, interp, epsilon, h);
        // std::cout << "iniciando U: " << ft.x << ", " << ft.y << ", " << ft.z << std::endl;
        std::pair<double, double> vt = GetUValue(image, wave, ndir, interp, epsilon, h);
        // std::cout << "finish him!" << std::endl;

        if (((vt.first < wave.v && vt.first < wave.planes[ndir].v[wave.gamma]) || (vt.second < wave.d && vt.second < wave.planes[ndir].d[wave.gamma])) && interp != I_PCHIP && interp != I_LINEAR) {
            if (interp == I_SPLINE || interp == I_HERMITE) {
                vt = GetUValue(image, wave, ndir, I_PCHIP, epsilon, h);
            }
            else {
                vt = GetUValue(image, wave, ndir, I_LINEAR, epsilon, h);
            }
        }
        // double ft = GetIntegralValue(image, plane_posD, neighD, norm(f0), norm(fn), interp, wave.gamma - 1);
        value = vt.first + norm(h * (neighD - plane_posD)) * (norm(f0) + norm(fn))/2.0;
        dist = vt.second + norm(h * (neighD - plane_posD));
        //value = vt + norm(h * (neighD - plane_posD)) * (norm(f0) + norm(fn))/2.0;
        /*if (fabs(neighD.y - 17.25) < 0.01 && (fabs(neighD.x - 4) < 0.01 || fabs(neighD.x - 96) < 0.01)) {
            std::cout << "value: " << neighD.y << ", " << neighD.x << " : " << vt << ", eps : " << epsilon <<
            ", d = " << norm(neighD - plane_posD) << ", v = " << value << std::endl;
        }*/

        return std::pair<double, double>(value, dist);
    }


    double GetEpsilonGradient(AugWmm &wave, int ndir, NodeD &neigh, NodeD &fn, NodeD &h) {

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

    double GetEpsilonGoldenSearch(Grid &image, AugWmm &wave, int ndir, int interp, NodeD &neigh, NodeD &fn, NodeD &h) {

        double br = (sqrt(5.0) + 1.0)/2.0;
        double a_1 = 0.0, b_1 = 1.0, x1_1 = b_1 - (b_1 - a_1) / br, x2_1 = a_1 + (b_1 - a_1) / br,
                f_x1 = MAX_VAL, f_x2 = MAX_VAL;

        double epsilon = 0.0;

        double res = GetInterpValue(image, wave, ndir, interp, 0.0, neigh, fn, h).first;
        double f_b = GetInterpValue(image, wave, ndir, interp, 1.0, neigh, fn, h).first;
        if (f_b < res) {
            res = f_b; epsilon = 1.0;
        }

        while (fabs(b_1 - a_1) > TAU) {
            f_x1 = GetInterpValue(image, wave, ndir, interp, x1_1, neigh, fn, h).first;
            f_x2 = GetInterpValue(image, wave, ndir, interp, x2_1, neigh, fn, h).first;
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

    std::pair<double, double> GetVal2D(Grid &image, Grid &u_surface, AugWmm &wave, NodeD &neigh, NodeD &fn, NodeD &h, int interp, int mode) {
        NodeD f0(image.at(wave.p, 0), image.at(wave.p, 1));
        double y0 = wave.v, d0 = wave.d;

        if (isinf(norm(f0)) || isnan(norm(f0)))
            f0 = fn;

        NodeD wave_p(wave.p.y, wave.p.x);
        NodeD diff = h * (neigh - wave_p);
        // double ft = GetIntegralValue(image, wave_p, neigh, norm(f0), norm(fn), interp, wave.gamma - 1);
        std::pair<double, double> val(y0 + norm(diff) * (norm(f0) + norm(fn)) / 2.0, d0 + norm(diff)), val_aux;
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
                    std::pair<double, double> res(MAX_VAL, MAX_VAL), res_aux;
                    double step = 1.0 / (double) wave.gamma;
                    // double v[2];
                    //NodeD dd_aux = step * wave.planes[ndir].dir, dp_aux;
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
                        if (res_aux.first < res.first) {
                            res = res_aux;
                            epsilon = epsilon_aux;
                        }
                    }
                }
                else
                    epsilon = GetEpsilonGoldenSearch(image, wave, ndir, interp, neigh, fn, h);
                // std::cout << "test " << neigh.y << ", " << neigh.x << " - " << ndir << " : " << epsilon << " - "<< GetInterpValue(image, wave, ndir, interp, epsilon, neigh, fn, h) << std::endl;
                val_aux = GetInterpValue(image, wave, ndir, interp, epsilon, neigh, fn, h);
                if (val_aux.first < val.first)
                    val = val_aux;
            }
        }
        return val;

    }


    void WmmAugSurface2D(Grid &image, std::vector<Node> &initials, std::vector<Node> &finals, NodeD &h,
                         int interp, int mode, int N, int gamma, Grid &u_surface) {
        bool isnewpos[8], wascomputed[8];
        double valcenter[8], distance[8];
        std::pair<double, double> value;
        Node3D imcenter[8];
        Node dirs[2], neigh, neigh_aux;
        int ndirs = 2;
        double step = 1.0 / (double) gamma, epsilon1;

        Grid_<unsigned char> state = Grid_<unsigned char>(image.rows, image.cols);

        std::multimap<double, AugWmm > trial_set;
        std::unordered_map<int, std::multimap<double, AugWmm >::iterator> mapa_trial;
        std::unordered_map<int, double> augmented_values;
        std::unordered_map<int, double> distance_augmented_values;

        std::multimap<double, AugWmm >::iterator trial_set_it;
        std::unordered_map<int, std::multimap<double, AugWmm >::iterator>::iterator mapa_trial_it;
        std::set<int> final_keys;
        std::pair<double, AugWmm > pr_trial;
        std::pair<int, std::multimap<double, AugWmm >::iterator> pr_mapa;

        int key, i, aug_key;
        AugWmm winner;
        NodeD fi, neighD, d1, neighD_aux;

        // Initialization
        for (i = 0; i < (int) initials.size(); i++) {
            key = initials[i].y + u_surface.rows*initials[i].x;
            if (mapa_trial.find(key) == mapa_trial.end() && u_surface.contains(initials[i])) {
                u_surface.at(initials[i], 0) = 0.0;
                u_surface.at(initials[i], 1) = 0.0;
                winner = AugWmm(N, gamma, 0);
                winner.v = 0.0;
                winner.d = 0.0;
                winner.p = initials[i];
                state.at(initials[i]) = P_TRIAL;
                pr_trial = std::pair<double, AugWmm >(0.0, winner);
                trial_set_it = trial_set.insert(pr_trial);
                pr_mapa = std::pair<int, typename std::multimap<double, AugWmm >::iterator>(key, trial_set_it);
                mapa_trial.insert(pr_mapa);
            }
        }

        for (i = 0; i < (int) finals.size(); i++) {
            key = finals[i].y * u_surface.cols + finals[i].x;
            if (final_keys.find(key) == final_keys.end() && u_surface.contains(finals[i])) {
                final_keys.insert(key);
            }
        }

        if (final_keys.empty())
            final_keys.insert(-1);

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

            if (final_keys.find(key) != final_keys.end()) {
                // std::cout << "Reaching point " << trial_set_it->second.p.y << ", " << trial_set_it->second.p.x << std::endl;
                final_keys.erase(key);
            }

            if (final_keys.empty())
                break;

            winner = trial_set_it->second;

            trial_set.erase(trial_set_it);
            mapa_trial.erase(mapa_trial_it);

            state.at(winner.p) = P_ALIVE;

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
                            winner.planes[m].d[v] = MAX_VAL;
                        }
                        else {
                            neigh_aux = Node(round(gamma * neighD_aux.y), round(gamma * neighD_aux.x));
                            aug_key = neigh_aux.y + gamma * u_surface.rows*neigh_aux.x;
                            winner.planes[m].v[v] = (augmented_values.find(aug_key) == augmented_values.end()) ? MAX_VAL : augmented_values[aug_key];
                            winner.planes[m].d[v] = (distance_augmented_values.find(aug_key) == distance_augmented_values.end()) ? MAX_VAL : distance_augmented_values[aug_key];
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
                    fi = NodeD(MAX_VAL, MAX_VAL);
                    // fi = Node3D(MAX_VAL, MAX_VAL, MAX_VAL);
                    valcenter[t] = MAX_VAL;
                    distance[t] = MAX_VAL;
                    continue;
                }
                fi = NodeD(image.at(neigh, 0), image.at(neigh, 1));
                valcenter[t] = u_surface.at(neigh, 0);
                distance[t] = u_surface.at(neigh, 1);
                if (isinf(norm(fi)) || isnan(norm(fi)))
                    fi = NodeD(MAX_VAL, MAX_VAL);
                imcenter[t] = to_polar(fi);
                if (u_surface.contains(neigh) && state.at(neigh) != P_ALIVE) {
                    neighD = NodeD(neigh.y, neigh.x);
                    std::pair<double, double> val_neigh = GetVal2D(image, u_surface, winner, neighD, fi, h, interp, mode);
                    // std::cout << neighD.y << ", " << neighD.x << " : " << val_neigh << " - " << valcenter[t] << std::endl;
                    if (val_neigh.first < valcenter[t]) {
                        valcenter[t] = val_neigh.first;
                        distance[t] = val_neigh.second;
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
                if (state.at(neigh) == P_TRIAL) {
                    mapa_trial_it = mapa_trial.find(key);
                    trial_set.erase(mapa_trial_it->second);
                    mapa_trial.erase(mapa_trial_it);
                }
                else {
                    state.at(neigh) = P_TRIAL;
                }

                AugWmm new_w(N, gamma, 2);
                new_w.p = neigh;
                dirs[0] = Node(yarray[(t+1)%8], xarray[(t+1)%8]) - d;
                dirs[1] = Node(yarray[(t+7)%8], xarray[(t+7)%8]) - d;
                new_w.v = valcenter[t];
                new_w.d = distance[t];
                double v_border[2] = {valcenter[(t+1)%8], valcenter[(t+7)%8]};
                double d_border[2] = {distance[(t+1)%8], distance[(t+7)%8]};
                /*if ((neigh.y == 17) && (neigh.x == 4 || neigh.x == 96)) {
                    std::cout << "winner: " << winner.p.y << ", " << winner.p.x << " - " << winner.v << std::endl;
                    if (winner.ndirs > 0) {
                        for (int m=0; m < ndirs; m++) {
                            std::cout << "m: " << m << " : " << winner.planes[m].dir.y << ", " << winner.planes[m].dir.x << std::endl;
                            for (int v=0; v <= gamma; v++) {
                                std::cout << "v : " << m << " - " << v << " - " << " : " << winner.planes[m].v[v] << std::endl;
                            }
                            for (int v=0; v <= 1; v++) {
                                std::cout << "fm : " << m << " - " << v << " - " << " : " <<
                                winner.planes[m].fm[v].z << ", " << winner.planes[m].fm[v].y << ", " << winner.planes[m].fm[v].x << std::endl;
                            }
                        }
                    }
                    std::cout << "update: " << new_w.p.y << ", " << new_w.p.x << " - " << new_w.v << std::endl;
                }*/
                // std::cout << "update: " << new_w.p.x << ", " << new_w.p.y << ", " << new_w.p.z << " - " << new_w.v << std::endl;
                for (int m=0; m < 2; m++) {
                    new_w.planes[m].dir = dirs[m];
                    setCoeffs2D(imcenter, new_w.planes[m].fm, t, interp, m);
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
                            value = std::pair<double, double>(valcenter[t], distance[t]);
                        else if (v == gamma)
                            value = std::pair<double, double>(v_border[m], d_border[m]);
                        else {
                            if (u_surface.contains(neighD_aux)) {
                                fi = GetFValue(image, new_w, m, interp, epsilon1, h);
                                value = GetVal2D(image, u_surface, winner, neighD_aux, fi, h, interp, mode);
                            }
                            else {
                                value = std::pair<double, double>(MAX_VAL, MAX_VAL);
                            }
                        }
                        /*if ((neigh.y == 17) && (neigh.x == 4 || neigh.x == 96)) {
                            std::cout << aug_key << ", " << neighD_aux.y << ", " << neighD_aux.x << " <-> " << aug_key <<
                                         " : " << v << " - " << value << std::endl;
                        }*/
                        if (augmented_values.find(aug_key) == augmented_values.end() || (value.first < augmented_values[aug_key])) {
                            augmented_values[aug_key] = value.first;
                        }
                        if (distance_augmented_values.find(aug_key) == distance_augmented_values.end() || (value.second < distance_augmented_values[aug_key])) {
                            distance_augmented_values[aug_key] = value.second;
                        }
                    }
                }
                wascomputed[t] = false;
                pr_trial = std::pair<double, AugWmm >(valcenter[t], new_w);
                trial_set_it = trial_set.insert(pr_trial);
                pr_mapa = std::pair<int, std::multimap<double, AugWmm >::iterator>(key, trial_set_it);
                mapa_trial.insert(pr_mapa);

                u_surface.at(new_w.p, 0) = valcenter[t];
                u_surface.at(new_w.p, 1) = distance[t];
            }

            // del(&winner);
        }

        free(state.data);
        return;

    }

    
}
