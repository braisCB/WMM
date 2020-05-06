#include "../TYPES/WMMStructs3D.h"
#include "../TYPES/utils3D.h"
#include <math.h>
#include <stdlib.h>
#include <vector>
#include <map>
#include <iostream>


namespace wmm3D {

    template<int N> double GetInterpValue(Grid3 &image, wmm3_<double, N> &wave, int ndir,
                                          int interp, std::pair<double, double> epsilon, Node3 &neigh,
                                          Node3D &h) {

        Node3D fn(image.at(neigh, 0), image.at(neigh, 1), image.at(neigh, 2));

        Node3 p1 = wave.p + wave.planes[ndir].dirs[0];
        Node3 p2 = wave.p + wave.planes[ndir].dirs[1];
        Node3 p3 = wave.p + wave.planes[ndir].dirs[0] + wave.planes[ndir].dirs[1];

        double value = MAX_VAL, ft, vt;

        if (interp == I_BILINEAR) {
            vt = (1.0 - epsilon.first)*(1.0 - epsilon.second)*wave.v +
                 epsilon.first*(1.0 - epsilon.second)*wave.planes[ndir].v[0] +
                 (1.0 - epsilon.first)*epsilon.second*wave.planes[ndir].v[1] +
                 epsilon.first*epsilon.second*wave.planes[ndir].v[2];
        }
        else if (interp == I_QUADATRIC) {
            vt = wave.planes[ndir].m[0] + epsilon.first*wave.planes[ndir].m[1] +
                 epsilon.second*wave.planes[ndir].m[2] + epsilon.first*epsilon.second*wave.planes[ndir].m[3];
        }
        else if (interp == I_CUBIC) {
            vt = wave.planes[ndir].m[0] + epsilon.first*wave.planes[ndir].m[1] +
                 epsilon.second*wave.planes[ndir].m[2] + epsilon.first*epsilon.second*wave.planes[ndir].m[3] +
                 epsilon.first*epsilon.first*wave.planes[ndir].m[4] + epsilon.second*epsilon.second*wave.planes[ndir].m[5] +
                 epsilon.first*epsilon.first*epsilon.second*wave.planes[ndir].m[6] +
                 epsilon.first*epsilon.second*epsilon.second*wave.planes[ndir].m[7] +
                 epsilon.first*epsilon.first*epsilon.second*epsilon.second*wave.planes[ndir].m[8];
        }
        else if (interp == I_HERMITE) {
            double m[2], v[3], y0, y1;
            double t = epsilon.first, t_2 = t*t, t_3 = t_2*t;
            int pos = (int) wave.planes[ndir].m[0];
            for (int p=0; p<3; p++) {
                if (pos == 0) {
                    m[0] = wave.planes[ndir].m[2 + 3*p + 1] - wave.planes[ndir].m[2 + 3*p];
                    m[1] = 1.0 / 2.0 * (wave.planes[ndir].m[2 + 3*p + 2] - wave.planes[ndir].m[2 + 3*p]);
                }
                else {
                    m[0] = 1.0 / 2.0 * (wave.planes[ndir].m[2 + 3*p + 2] - wave.planes[ndir].m[2 + 3*p]);
                    m[1] = wave.planes[ndir].m[2 + 3*p + 2] - wave.planes[ndir].m[2 + 3*p + 1];
                }
                y0 = wave.planes[ndir].m[2 + 3*p + pos];
                y1 = wave.planes[ndir].m[2 + 3*p + 1 + pos];
                v[p] = (2.0 * t_3 - 3.0 * t_2 + 1.0) * y0 + (t_3 - 2.0 * t_2 + t) * m[0] +
                       (-2.0 * t_3 + 3.0 * t_2) * y1 + (t_3 - t_2) * m[1];
            }
            t = epsilon.second; t_2 = t*t; t_3 = t_2*t;
            pos = wave.planes[ndir].m[1];
            if (pos == 0) {
                m[0] = v[1] - v[0];
                m[1] = 1.0 / 2.0 * (v[2] - v[0]);
            }
            else {
                m[0] = 1.0 / 2.0 * (v[2] - v[0]);
                m[1] = v[2] - v[1];
            }
            y0 = v[pos];
            y1 = v[1 + pos];
            vt = (2.0 * t_3 - 3.0 * t_2 + 1.0) * y0 + (t_3 - 2.0 * t_2 + t) * m[0] +
                   (-2.0 * t_3 + 3.0 * t_2) * y1 + (t_3 - t_2) * m[1];

            
            t = epsilon.second; t_2 = t*t; t_3 = t_2*t;
            pos = (int) wave.planes[ndir].m[1];
            for (int p=0; p<3; p++) {
                if (pos == 0) {
                    m[0] = wave.planes[ndir].m[2 + 3 + p] - wave.planes[ndir].m[2 + p];
                    m[1] = 1.0 / 2.0 * (wave.planes[ndir].m[2 + 6 + p] - wave.planes[ndir].m[2 + p]);
                }
                else {
                    m[0] = 1.0 / 2.0 * (wave.planes[ndir].m[2 + 6 + p] - wave.planes[ndir].m[2 + p]);
                    m[1] = wave.planes[ndir].m[2 + 6 + p] - wave.planes[ndir].m[2 + 3 + p];
                }
                y0 = wave.planes[ndir].m[2 + 3*pos + p];
                y1 = wave.planes[ndir].m[2 + 3 + 3*pos + p];
                v[p] = (2.0 * t_3 - 3.0 * t_2 + 1.0) * y0 + (t_3 - 2.0 * t_2 + t) * m[0] +
                       (-2.0 * t_3 + 3.0 * t_2) * y1 + (t_3 - t_2) * m[1];
            }
            t = epsilon.first; t_2 = t*t; t_3 = t_2*t;
            pos = wave.planes[ndir].m[0];
            if (pos == 0) {
                m[0] = v[1] - v[0];
                m[1] = 1.0 / 2.0 * (v[2] - v[0]);
            }
            else {
                m[0] = 1.0 / 2.0 * (v[2] - v[0]);
                m[1] = v[2] - v[1];
            }
            y0 = v[pos];
            y1 = v[1 + pos];
            vt = 0.5 * (vt + (2.0 * t_3 - 3.0 * t_2 + 1.0) * y0 + (t_3 - 2.0 * t_2 + t) * m[0] +
                   (-2.0 * t_3 + 3.0 * t_2) * y1 + (t_3 - t_2) * m[1]);
        }
        else if (interp == I_PCHIP) {
            double m[2], v[3], y0, y1, d0, d1;
            double t = epsilon.first, t_2 = t*t, t_3 = t_2*t;
            int pos = (int) wave.planes[ndir].m[0];
            for (int p=0; p<3; p++) {
                d0 = wave.planes[ndir].m[2 + 3*p + 1] - wave.planes[ndir].m[2 + 3*p];
                d1 = wave.planes[ndir].m[2 + 3*p + 2] - wave.planes[ndir].m[2 + 3*p + 1];
                if (pos == 0) {
                    m[0] = (3.0 * d0 - d1) / 2.0;
                    if ((m[0] * d0) <= 0.0)
                        m[0] = 0.0;
                    else if (((d1 * d0) <= 0.0) && (fabs(m[0]) > fabs(3.0 * d0)))
                        m[0] = 3.0 * d0;
                    m[1] = (d0 * d1 <= 0.0) ? 0.0 : 2.0 * d0 * d1 / (d0 + d1);
                }
                else {
                    m[0] = (d0 * d1 <= 0.0) ? 0.0 : 2.0 * d0 * d1 / (d0 + d1);
                    m[1] = (3.0 * d1 - d0) / 2.0;
                    if ((m[1] * d1) <= 0.0)
                        m[1] = 0.0;
                    else if (((d0 * d1) <= 0.0) && (fabs(m[1]) > fabs(3.0 * d1)))
                        m[1] = 3.0 * d1;
                }
                y0 = wave.planes[ndir].m[2 + 3*p + pos];
                y1 = wave.planes[ndir].m[2 + 3*p + 1 + pos];
                v[p] = (2.0 * t_3 - 3.0 * t_2 + 1.0) * y0 + (t_3 - 2.0 * t_2 + t) * m[0] +
                       (-2.0 * t_3 + 3.0 * t_2) * y1 + (t_3 - t_2) * m[1];
            }
            t = epsilon.second; t_2 = t*t; t_3 = t_2*t;
            pos = wave.planes[ndir].m[1];
            d0 = v[1] - v[0];
            d1 = v[2] - v[1];
            if (pos == 0) {
                m[0] = (3.0 * d0 - d1) / 2.0;
                if ((m[0] * d0) <= 0.0)
                    m[0] = 0.0;
                else if (((d1 * d0) <= 0.0) && (fabs(m[0]) > fabs(3.0 * d0)))
                    m[0] = 3.0 * d0;
                m[1] = (d0 * d1 <= 0.0) ? 0.0 : 2.0 * d0 * d1 / (d0 + d1);
            }
            else {
                m[0] = (d0 * d1 <= 0.0) ? 0.0 : 2.0 * d0 * d1 / (d0 + d1);
                m[1] = (3.0 * d1 - d0) / 2.0;
                if ((m[1] * d1) <= 0.0)
                    m[1] = 0.0;
                else if (((d0 * d1) <= 0.0) && (fabs(m[1]) > fabs(3.0 * d1)))
                    m[1] = 3.0 * d1;
            }
            y0 = v[pos];
            y1 = v[1 + pos];
            vt = (2.0 * t_3 - 3.0 * t_2 + 1.0) * y0 + (t_3 - 2.0 * t_2 + t) * m[0] +
                   (-2.0 * t_3 + 3.0 * t_2) * y1 + (t_3 - t_2) * m[1];
            
            t = epsilon.second; t_2 = t*t; t_3 = t_2*t;
            pos = (int) wave.planes[ndir].m[1];
            
            for (int p=0; p<3; p++) {
                d0 = wave.planes[ndir].m[2 + 3 + p] - wave.planes[ndir].m[2 + p];
                d1 = wave.planes[ndir].m[2 + 6 + p] - wave.planes[ndir].m[2 + 3 + p];
                if (pos == 0) {
                    m[0] = (3.0 * d0 - d1) / 2.0;
                    if ((m[0] * d0) <= 0.0)
                        m[0] = 0.0;
                    else if (((d1 * d0) <= 0.0) && (fabs(m[0]) > fabs(3.0 * d0)))
                        m[0] = 3.0 * d0;
                    m[1] = (d0 * d1 <= 0.0) ? 0.0 : 2.0 * d0 * d1 / (d0 + d1);
                }
                else {
                    m[0] = (d0 * d1 <= 0.0) ? 0.0 : 2.0 * d0 * d1 / (d0 + d1);
                    m[1] = (3.0 * d1 - d0) / 2.0;
                    if ((m[1] * d1) <= 0.0)
                        m[1] = 0.0;
                    else if (((d0 * d1) <= 0.0) && (fabs(m[1]) > fabs(3.0 * d1)))
                        m[1] = 3.0 * d1;
                }
                y0 = wave.planes[ndir].m[2 + 3*pos + p];
                y1 = wave.planes[ndir].m[2 + 3 + 3*pos + p];
                v[p] = (2.0 * t_3 - 3.0 * t_2 + 1.0) * y0 + (t_3 - 2.0 * t_2 + t) * m[0] +
                       (-2.0 * t_3 + 3.0 * t_2) * y1 + (t_3 - t_2) * m[1];
            }
            t = epsilon.first; t_2 = t*t; t_3 = t_2*t;
            pos = wave.planes[ndir].m[0];
            d0 = v[1] - v[0];
            d1 = v[2] - v[1];
            if (pos == 0) {
                m[0] = (3.0 * d0 - d1) / 2.0;
                if ((m[0] * d0) <= 0.0)
                    m[0] = 0.0;
                else if (((d1 * d0) <= 0.0) && (fabs(m[0]) > fabs(3.0 * d0)))
                    m[0] = 3.0 * d0;
                m[1] = (d0 * d1 <= 0.0) ? 0.0 : 2.0 * d0 * d1 / (d0 + d1);
            }
            else {
                m[0] = (d0 * d1 <= 0.0) ? 0.0 : 2.0 * d0 * d1 / (d0 + d1);
                m[1] = (3.0 * d1 - d0) / 2.0;
                if ((m[1] * d1) <= 0.0)
                    m[1] = 0.0;
                else if (((d0 * d1) <= 0.0) && (fabs(m[1]) > fabs(3.0 * d1)))
                    m[1] = 3.0 * d1;
            }
            y0 = v[pos];
            y1 = v[1 + pos];
            vt = 0.5 * (vt + (2.0 * t_3 - 3.0 * t_2 + 1.0) * y0 + (t_3 - 2.0 * t_2 + t) * m[0] +
                   (-2.0 * t_3 + 3.0 * t_2) * y1 + (t_3 - t_2) * m[1]);
        }
        else {
            double m[2], v[3], y0, y1;
            double t = epsilon.first;
            int pos = (int) wave.planes[ndir].m[0];
            for (int p=0; p<3; p++) {
                if (pos == 0) {
                    m[0] = 0.0;
                    m[1] = 6.0*(wave.planes[ndir].m[2 + 3*p] - 2.0*wave.planes[ndir].m[2 + 3*p + 1] + wave.planes[ndir].m[2 + 3*p + 2])/4.0;
                }
                else {
                    m[0] = 6.0*(wave.planes[ndir].m[2 + 3*p] - 2.0*wave.planes[ndir].m[2 + 3*p + 1] + wave.planes[ndir].m[2 + 3*p + 2])/4.0;
                    m[1] = 0.0;
                }
                y0 = wave.planes[ndir].m[2 + 3*p + pos];
                y1 = wave.planes[ndir].m[2 + 3*p + 1 + pos];
                v[p] = y0 + t*(-m[1]/6.0 - m[0]/3.0 + y1 - y0 + t*(m[0]/2.0 + t*(m[1] - m[0])/6.0));
            }
            t = epsilon.second;
            pos = wave.planes[ndir].m[1];
            if (pos == 0) {
                m[0] = 0.0;
                m[1] = 6.0*(v[0] - 2.0*v[1] + v[2])/4.0;
            }
            else {
                m[0] = 6.0*(v[0] - 2.0*v[1] + v[2])/4.0;
                m[1] = 0.0;
            }
            y0 = v[pos];
            y1 = v[1 + pos];
            vt = y0 + t*(-m[1]/6.0 - m[0]/3.0 + y1 - y0 + t*(m[0]/2.0 + t*(m[1] - m[0])/6.0));
            
            t = epsilon.second;
            pos = (int) wave.planes[ndir].m[1];
            
            for (int p=0; p<3; p++) {
                if (pos == 0) {
                    m[0] = 0.0;
                    m[1] = 6.0*(wave.planes[ndir].m[2 + p] - 2.0*wave.planes[ndir].m[2 + 3 + p] + wave.planes[ndir].m[2 + 6 + p])/4.0;
                }
                else {
                    m[0] = 6.0*(wave.planes[ndir].m[2 + p] - 2.0*wave.planes[ndir].m[2 + 3 + p] + wave.planes[ndir].m[2 + 6 + p])/4.0;
                    m[1] = 0.0;
                }
                y0 = wave.planes[ndir].m[2 + 3*pos + p];
                y1 = wave.planes[ndir].m[2 + 3 + 3*pos + p];
                v[p] = y0 + t*(-m[1]/6.0 - m[0]/3.0 + y1 - y0 + t*(m[0]/2.0 + t*(m[1] - m[0])/6.0));
            }
            t = epsilon.first;
            pos = wave.planes[ndir].m[0];
            if (pos == 0) {
                m[0] = 0.0;
                m[1] = 6.0*(v[0] - 2.0*v[1] + v[2])/4.0;
            }
            else {
                m[0] = 6.0*(v[0] - 2.0*v[1] + v[2])/4.0;
                m[1] = 0.0;
            }
            y0 = v[pos];
            y1 = v[1 + pos];
            vt = y0 + t*(-m[1]/6.0 - m[0]/3.0 + y1 - y0 + t*(m[0]/2.0 + t*(m[1] - m[0])/6.0));
        }
        
        if ((vt < wave.v && vt < wave.planes[ndir].v[0] && vt < wave.planes[ndir].v[1] && vt < wave.planes[ndir].v[2] &&
            interp != I_PCHIP && interp != I_BILINEAR) ||
            ((vt > wave.v && vt > wave.planes[ndir].v[0] && vt > wave.planes[ndir].v[1] && vt > wave.planes[ndir].v[2] &&
            interp != I_PCHIP && interp != I_BILINEAR))) {
            if (interp == I_HERMITE)
                value = GetInterpValue(image, wave, ndir, I_PCHIP, epsilon, neigh, h);
            else
                value = GetInterpValue(image, wave, ndir, I_BILINEAR, epsilon, neigh, h);
        }
        else {
            Node3D plane_posD(wave.p.x + epsilon.first*wave.planes[ndir].dirs[0].x + epsilon.second*wave.planes[ndir].dirs[1].x,
                                     wave.p.y + epsilon.first*wave.planes[ndir].dirs[0].y + epsilon.second*wave.planes[ndir].dirs[1].y,
                                     wave.p.z + epsilon.first*wave.planes[ndir].dirs[0].z + epsilon.second*wave.planes[ndir].dirs[1].z);
            Node3D neighD(neigh.x, neigh.y, neigh.z);
            Node3D diff = h * (neighD - plane_posD);
            Node3D a = diff / norm(diff);
            value = vt + norm(diff) / sqrt(1.0 + (fn.x*a.x + fn.y*a.y + fn.z*a.z)*(fn.x*a.x + fn.y*a.y + fn.z*a.z));
        }

        return value;
    }


    template<int N> std::pair<double, double> GetEpsilonGradient(wmm3_<double, N> &wave, int ndir,
                                                                 Node3 &neigh, Node3D &fn, Node3D &h) {

        Node3D ph(h.x*wave.p.x, h.y*wave.p.y, h.z*wave.p.z), nh(h.x*neigh.x, h.y*neigh.y, h.z*neigh.z), diff, x;
        Node3D d0h(h.x*wave.planes[ndir].dirs[0].x, h.y*wave.planes[ndir].dirs[0].y, h.z*wave.planes[ndir].dirs[0].z);
        Node3D d1h(h.x*wave.planes[ndir].dirs[1].x, h.y*wave.planes[ndir].dirs[1].y, h.z*wave.planes[ndir].dirs[1].z);

        double A = d0h.y*d1h.z - d0h.z*d1h.y;
        double B = d0h.z*d1h.x - d0h.x*d1h.z;
        double C = d0h.x*d1h.y - d0h.y*d1h.x;
        double D = -A*ph.x - B*ph.y - C*ph.z;

        double den = A*fn.x + B*fn.y + C*fn.z, t;

        std::pair<double, double> epsilon(0.0, 0.0);

        if (fabs(den) > 0.0) {
            t = -(A*nh.x + B*nh.y + C*nh.z + D)/den;
            x = Node3D(nh.x + t*fn.x, nh.y + t*fn.y, nh.z + t*fn.z), diff;
            diff = (x - ph) / h;
        }
        else {
            t = (A*nh.x + B*nh.y + C*nh.z + D)/(A*A + B*B + C*C);
            x = Node3D(nh.x - t*A, nh.y - t*B, nh.z - t*C);
            diff = (x - ph) / h;
        }

        epsilon.first = diff.x*wave.planes[ndir].dirs[0].x + diff.y*wave.planes[ndir].dirs[0].y + diff.z*wave.planes[ndir].dirs[0].z;
        if (epsilon.first < 0.0)
            epsilon.first = 0.0;
        else if (epsilon.first > 1.0)
            epsilon.first = 1.0;

        epsilon.second = diff.x*wave.planes[ndir].dirs[1].x + diff.y*wave.planes[ndir].dirs[1].y + diff.z*wave.planes[ndir].dirs[1].z;
        if (epsilon.second < 0.0)
            epsilon.second = 0.0;
        else if (epsilon.second > 1.0)
            epsilon.second = 1.0;

        //std::cout << neigh.x << ", " << neigh.y << ", " << neigh.z << " : " << epsilon.first << ", " << epsilon.second << std::endl;
        return epsilon;

    }

    template<int N> std::pair<double, double> GetEpsilonHopfLax(wmm3_<double, N> &wave, int ndir,
                                                                    Node3 &neigh, Node3D &h) {

        Node3D ph(h.x*wave.p.x, h.y*wave.p.y, h.z*wave.p.z), nh(h.x*neigh.x, h.y*neigh.y, h.z*neigh.z);
        Node3D d0h(h.x*wave.planes[ndir].dirs[0].x, h.y*wave.planes[ndir].dirs[0].y, h.z*wave.planes[ndir].dirs[0].z);
        Node3D d1h(h.x*wave.planes[ndir].dirs[1].x, h.y*wave.planes[ndir].dirs[1].y, h.z*wave.planes[ndir].dirs[1].z);

        Node3D inc0, inc1;
        Node3 d01h = wave.planes[ndir].dirs[0] + wave.planes[ndir].dirs[1], ori = neigh - wave.p;
        double cos_delta_0 = (wave.planes[ndir].v[0] - wave.v) / norm(d0h);
        double cos_delta_1 = (wave.planes[ndir].v[1] - wave.v) / norm(d1h);
        if (fabs(cos_delta_0) > 1) cos_delta_0 = sgn(cos_delta_0);
        if (fabs(cos_delta_1) > 1) cos_delta_1 = sgn(cos_delta_1);
        double sin_delta_0 = sqrt(1. - cos_delta_0*cos_delta_0);
        double sin_delta_1 = sqrt(1. - cos_delta_1*cos_delta_1);
        if (d01h.x == 0) {
            // inc0 = (wave.planes[ndir].dirs[0].y == 0) ? Node3D(sgn(ori.x)*sin_delta_0, 0., sgn(ori.z)*cos_delta_0) : Node3D(sgn(ori.x)*sin_delta_0, sgn(ori.y)*cos_delta_0, 0.);
            inc0 = (wave.planes[ndir].dirs[0].y == 0) ? Node3D(-wave.planes[ndir].dirs[0].z*cos_delta_0, 0., sgn(ori.x)*sin_delta_0) :
                                                        Node3D(-wave.planes[ndir].dirs[0].y*cos_delta_0, sgn(ori.x)*sin_delta_0, 0.);
            inc1 = (wave.planes[ndir].dirs[1].y == 0) ? Node3D(-wave.planes[ndir].dirs[1].z*cos_delta_1, 0., sgn(ori.x)*sin_delta_1) :
                                                        Node3D(-wave.planes[ndir].dirs[1].y*cos_delta_1, sgn(ori.x)*sin_delta_1, 0.);
        }
        else if (d01h.y == 0) {
            inc0 = (wave.planes[ndir].dirs[0].x == 0) ? Node3D(0., -wave.planes[ndir].dirs[0].z*cos_delta_0, sgn(ori.y)*sin_delta_0) :
                                                        Node3D(sgn(ori.y)*sin_delta_0, -wave.planes[ndir].dirs[0].x*cos_delta_0, 0.);
            inc1 = (wave.planes[ndir].dirs[1].x == 0) ? Node3D(0., -wave.planes[ndir].dirs[1].z*cos_delta_1, sgn(ori.y)*sin_delta_1) :
                                                        Node3D(sgn(ori.y)*sin_delta_1, -wave.planes[ndir].dirs[1].x*cos_delta_1, 0.);
        }
        else {
            inc0 = (wave.planes[ndir].dirs[0].x == 0) ? Node3D(0., sgn(ori.z)*sin_delta_0, -wave.planes[ndir].dirs[0].y*cos_delta_0) :
                                                        Node3D(sgn(ori.z)*sin_delta_0,  0., -wave.planes[ndir].dirs[0].x*cos_delta_0);
            inc1 = (wave.planes[ndir].dirs[1].x == 0) ? Node3D(0., sgn(ori.z)*sin_delta_1, -wave.planes[ndir].dirs[1].y*cos_delta_1) :
                                                        Node3D(sgn(ori.z)*sin_delta_1,  0., -wave.planes[ndir].dirs[1].x*cos_delta_1);
        }
        Node3D fn = Node3D(inc0.y*inc1.z - inc0.z*inc1.y, inc0.z*inc1.x - inc0.x*inc1.z, inc0.x*inc1.y - inc0.y*inc1.x);

        double A = d0h.y*d1h.z - d0h.z*d1h.y;
        double B = d0h.z*d1h.x - d0h.x*d1h.z;
        double C = d0h.x*d1h.y - d0h.y*d1h.x;
        double D = -A*ph.x - B*ph.y - C*ph.z;

        double den = A*fn.x + B*fn.y + C*fn.z;

        double t = -(A*nh.x + B*nh.y + C*nh.z + D)/den;
        Node3D x(nh.x + t*fn.x, nh.y + t*fn.y, nh.z + t*fn.z), diff;

        std::pair<double, double> epsilon(0.0, 0.0);

        if (fabs(den) > 0.0)
            diff = (x - ph) / h;
        else {
            t = (A*nh.x + B*nh.y + C*nh.z + D)/(A*A + B*B + C*C);
            x = Node3D(nh.x - t*A, nh.y - t*B, nh.z - t*C);
            diff = (x - ph) / h;
        }

        epsilon.first = diff.x*wave.planes[ndir].dirs[0].x + diff.y*wave.planes[ndir].dirs[0].y + diff.z*wave.planes[ndir].dirs[0].z;
        if (epsilon.first < 0.0)
            epsilon.first = 0.0;
        else if (epsilon.first > 1.0)
            epsilon.first = 1.0;

        epsilon.second = diff.x*wave.planes[ndir].dirs[1].x + diff.y*wave.planes[ndir].dirs[1].y + diff.z*wave.planes[ndir].dirs[1].z;
        if (epsilon.second < 0.0)
            epsilon.second = 0.0;
        else if (epsilon.second > 1.0)
            epsilon.second = 1.0;

        return epsilon;

    }

    template<int N> std::pair<double, double> GetEpsilonGoldenSearch(Grid3 &image, wmm3_<double, N> &wave,
                                                                     int ndir, int interp, Node3 &neigh,
                                                                     Node3D &h) {

        double br = (sqrt(5.0) + 1.0)/2.0;
        double a_1 = 0.0, b_1 = 1.0, x1_1 = b_1 - (b_1 - a_1) / br, x2_1 = a_1 + (b_1 - a_1) / br,
               a_2 = 0.0, b_2 = 1.0, x1_2 = b_2 - (b_2 - a_2) / br, x2_2 = a_2 + (b_2 - a_2) / br,
               f_x1_1 = MAX_VAL, f_x2_1 = MAX_VAL, f_x1_2 = MAX_VAL, f_x2_2 = MAX_VAL;

        std::pair<double, double> epsilon(0.0, 0.0);

        double res = GetInterpValue(image, wave, ndir, interp, std::pair<double, double>(0.0, 0.0), neigh, h);
        double f_b = GetInterpValue(image, wave, ndir, interp, std::pair<double, double>(1.0, 0.0), neigh, h);
        if (f_b < res) {
            res = f_b; epsilon = std::pair<double, double>(1.0, 0.0);
        }
        f_b = GetInterpValue(image, wave, ndir, interp, std::pair<double, double>(0.0, 1.0), neigh, h);
        if (f_b < res) {
            res = f_b; epsilon = std::pair<double, double>(0.0, 1.0);
        }
        f_b = GetInterpValue(image, wave, ndir, interp, std::pair<double, double>(1.0, 1.0), neigh, h);
        if (f_b < res) {
            res = f_b; epsilon = std::pair<double, double>(1.0, 1.0);
        }

        while (fabs(b_1 - a_1) > TAU && fabs(b_2 - a_2) > TAU) {
            f_x1_1 = GetInterpValue(image, wave, ndir, interp, std::pair<double, double>(x1_1, x1_2), neigh, h);
            f_x2_1 = GetInterpValue(image, wave, ndir, interp, std::pair<double, double>(x2_1, x1_2), neigh, h);
            f_x1_2 = GetInterpValue(image, wave, ndir, interp, std::pair<double, double>(x1_1, x2_2), neigh, h);
            f_x2_2 = GetInterpValue(image, wave, ndir, interp, std::pair<double, double>(x2_1, x2_2), neigh, h);
            if (f_x1_1 <= std::min(f_x2_1, f_x1_2) &&  f_x1_1 <=  f_x2_2) {
                b_1 = x2_1;
                b_2 = x2_2;
            }
            else if (f_x2_1 <= std::min(f_x1_1, f_x1_2) && f_x2_1 <=  f_x2_2) {
                a_1 = x1_1;
                b_2 = x2_2;
            }
            else if (f_x1_2 <= std::min(f_x1_1, f_x2_2) && f_x2_2) {
                b_1 = x2_1;
                a_2 = x1_2;
            }
            else {
                a_1 = x1_1;
                a_2 = x1_2;
            }

            x1_1 = b_1 - (b_1 - a_1) / br, x2_1 = a_1 + (b_1 - a_1) / br;
            x1_2 = b_2 - (b_2 - a_2) / br, x2_2 = a_2 + (b_2 - a_2) / br;
        }

        if (f_x1_1 < res) {
            epsilon = std::pair<double, double>(x1_1, x1_2);
            res = f_x1_1;
        }
        if (f_x2_1 <= res) {
            epsilon = std::pair<double, double>(x2_1, x1_2);
            res = f_x2_1;
        }
        if (f_x1_2 <= res) {
            epsilon = std::pair<double, double>(x1_1, x2_2);
            res = f_x1_2;
        }
        if (f_x2_2 <= res) {
            epsilon = std::pair<double, double>(x2_1, x2_2);
            res = f_x2_2;
        }

        return epsilon;

    }



    /*template<int N> std::pair<double, double> GetEpsilonGoldenSearch(Grid3 &image, wmm3_<double, N> &wave,
                                                                     int ndir, int interp, Node3 &neigh,
                                                                     Node3D &h) {

        double a_1 = 0.0, b_1 = 1.0, x1_1 = a_1 + (1.0-RESPHI)*(b_1 - a_1), x2_1 = a_1 + RESPHI*(b_1 - a_1),
                f_x1 = MAX_VAL, f_x2 = MAX_VAL;

        double a_2 = 0.0, b_2 = 1.0, x1_2 = a_2 + (1.0-RESPHI)*(b_2 - a_2), x2_2 = a_2 + RESPHI*(b_2 - a_2);

        std::pair<double, double> epsilon(0.0, 0.0);

        double res = GetInterpValue(image, wave, ndir, interp, std::pair<double, double>(0.0, 0.0), neigh, h);
        double f_b = GetInterpValue(image, wave, ndir, interp, std::pair<double, double>(1.0, 0.0), neigh, h);
        if (f_b < res) {
            res = f_b; epsilon = std::pair<double, double>(1.0, 0.0);
        }
        f_b = GetInterpValue(image, wave, ndir, interp, std::pair<double, double>(0.0, 1.0), neigh, h);
        if (f_b < res) {
            res = f_b; epsilon = std::pair<double, double>(0.0, 1.0);
        }
        f_b = GetInterpValue(image, wave, ndir, interp, std::pair<double, double>(1.0, 1.0), neigh, h);
        if (f_b < res) {
            res = f_b; epsilon = std::pair<double, double>(1.0, 1.0);
        }

        f_x1 = GetInterpValue(image, wave, ndir, interp, std::pair<double, double>(0.0, x1_2), neigh, h);
        f_x2 = GetInterpValue(image, wave, ndir, interp, std::pair<double, double>(0.0, x2_2), neigh, h);

        while (fabs(b_1 - a_1) > TAU && fabs(b_2 - a_2) > TAU) {
            if(f_x1 < f_x2) {
                f_x1 = GetInterpValue(image, wave, ndir, interp, std::pair<double, double>(x1_1, x1_2), neigh, h);
                f_x2 = GetInterpValue(image, wave, ndir, interp, std::pair<double, double>(x2_1, x1_2), neigh, h);
                b_2 = x2_2; x2_2 = x1_2; x1_2 = a_2 + (1.0 - RESPHI)*(b_2 - a_2);
            }
            else {
                f_x1 = GetInterpValue(image, wave, ndir, interp, std::pair<double, double>(x1_1, x2_2), neigh, h);
                f_x2 = GetInterpValue(image, wave, ndir, interp, std::pair<double, double>(x2_1, x2_2), neigh, h);
                a_2 = x1_2; x1_2 = x2_2; x2_2 = a_2 + RESPHI*(b_2 - a_2);
            }
            if(f_x1 < f_x2) {
                f_x1 = GetInterpValue(image, wave, ndir, interp, std::pair<double, double>(x1_1, x1_2), neigh, h);
                f_x2 = GetInterpValue(image, wave, ndir, interp, std::pair<double, double>(x1_1, x2_2), neigh, h);
                b_1 = x2_1; x2_1 = x1_1; x1_1 = a_1 + (1.0 - RESPHI)*(b_1 - a_1);
            }
            else {
                f_x1 = GetInterpValue(image, wave, ndir, interp, std::pair<double, double>(x2_1, x1_2), neigh, h);
                f_x2 = GetInterpValue(image, wave, ndir, interp, std::pair<double, double>(x2_1, x2_2), neigh, h);
                a_1 = x1_1; x1_1 = x2_1; x2_1 = a_1 + RESPHI*(b_1 - a_1);
            }
        }

        if (f_x1 < res) {
            epsilon = std::pair<double, double>(x1_1, x1_2);
            f_x1 = GetInterpValue(image, wave, ndir, interp, std::pair<double, double>(x1_1, 1.0), neigh, h);
            if (f_x1 < res)
                epsilon = std::pair<double, double>(x1_1, 1.0);
            f_x1 = GetInterpValue(image, wave, ndir, interp, std::pair<double, double>(x1_1, 0.0), neigh, h);
            if (f_x1 < res)
                epsilon = std::pair<double, double>(x1_1, 0.0);
            f_x1 = GetInterpValue(image, wave, ndir, interp, std::pair<double, double>(1.0, x1_2), neigh, h);
            if (f_x1 < res)
                epsilon = std::pair<double, double>(1.0, x1_2);
            f_x1 = GetInterpValue(image, wave, ndir, interp, std::pair<double, double>(0.0, x1_2), neigh, h);
            if (f_x1 < res)
                epsilon = std::pair<double, double>(0.0, x1_2);
        }
        if (f_x2 < std::min(res, f_x1)) {
            epsilon = std::pair<double, double>(x2_1, x2_2);
            f_x1 = GetInterpValue(image, wave, ndir, interp, std::pair<double, double>(x2_1, 1.0), neigh, h);
            if (f_x1 < res)
                epsilon = std::pair<double, double>(x2_1, 1.0);
            f_x1 = GetInterpValue(image, wave, ndir, interp, std::pair<double, double>(x2_1, 0.0), neigh, h);
            if (f_x1 < res)
                epsilon = std::pair<double, double>(x2_1, 0.0);
            f_x1 = GetInterpValue(image, wave, ndir, interp, std::pair<double, double>(1.0, x2_2), neigh, h);
            if (f_x1 < res)
                epsilon = std::pair<double, double>(1.0, x2_2);
            f_x1 = GetInterpValue(image, wave, ndir, interp, std::pair<double, double>(0.0, x2_2), neigh, h);
            if (f_x1 < res)
                epsilon = std::pair<double, double>(0.0, x2_2);
        }


        return epsilon;
        
    }*/


    template<int N> double GetVal3D(Grid3 &image, Grid3 &u_surface, wmm3_<double, N> &wave,
                                    Node3 &neigh, Node3D &h, int interp, int mode) {

        Node3D f0(image.at(wave.p, 0), image.at(wave.p, 1), image.at(wave.p, 2)), fn(image.at(neigh, 0), image.at(neigh, 1), image.at(neigh, 2));
        double y0 = wave.v;

        if (isinf(norm(f0)) || isnan(norm(f0)))
            f0 = fn;

        Node3D diff = h * (neigh - wave.p);
        Node3D a = diff / norm(diff);
        double val = y0 + norm(diff) / sqrt(1.0 + (fn.x*a.x + fn.y*a.y + fn.z*a.z)*(fn.x*a.x + fn.y*a.y + fn.z*a.z));

        if (wave.ndirs > 0) {

            for (int ndir = 0; ndir < wave.ndirs; ndir++) {

                //Node3 p0 = wave.p + wave.planes[ndir].dirs[0];
                //Node3 p1 = wave.p + wave.planes[ndir].dirs[1];
                //Node3 p2 = wave.p - wave.planes[ndir].dirs[0];
                //Node3 p3 = wave.p - wave.planes[ndir].dirs[1];

                //if (!image.contains(p0) || !image.contains(p1) || !image.contains(p2) || !image.contains(p3))
                //    continue;

                std::pair<double, double> epsilon;
                if (mode == M_GRADIENT)
                    epsilon = GetEpsilonGradient(wave, ndir, neigh, fn, h);
                else if (mode == M_HOPFLAX)
                    epsilon = GetEpsilonHopfLax(wave, ndir, neigh, h);
                else
                    epsilon = GetEpsilonGoldenSearch(image, wave, ndir, interp, neigh, h);
                val = std::min(val, GetInterpValue(image, wave, ndir, interp, epsilon, neigh, h));
            }
        }
        return val;

    }

    void setCoeffs3D(double y[3][3][3], double *m, Node3 d, Node3 *dirs, int interp) {

        if (interp == I_BILINEAR)
            return;
        else if (interp == I_QUADATRIC) {
            Node3 pos = d + Node3(1, 1, 1);
            Node3 pos1 = d + dirs[0] + Node3(1, 1, 1);
            Node3 pos2 = d + dirs[1] + Node3(1, 1, 1);
            Node3 pos3 = d + dirs[0] + dirs[1] + Node3(1, 1, 1);
            m[0] = y[pos.x][pos.y][pos.z];
            m[1] = y[pos1.x][pos1.y][pos1.z] - m[0];
            m[2] = y[pos2.x][pos2.y][pos2.z] - m[0];
            m[3] = y[pos3.x][pos3.y][pos3.z] - m[2] - m[1] - m[0];
        }
        else if (interp == I_CUBIC) {
            Node3 pos = d + Node3(1, 1, 1);
            int dx = 2;
            if ((dirs[0].x != 0 && pos.x == 1) || (dirs[0].y != 0 && pos.y == 1) || (dirs[0].z != 0 && pos.z == 1))
                dx = -1;
            int dy = 2;
            if ((dirs[1].x != 0 && pos.x == 1) || (dirs[1].y != 0 && pos.y == 1) || (dirs[1].z != 0 && pos.z == 1))
                dy = -1;
            double p1 = y[pos.x + dirs[0].x][pos.y + dirs[0].y][pos.z + dirs[0].z];
            double p2 = y[pos.x + dx * dirs[0].x][pos.y + dx * dirs[0].y][pos.z + dx * dirs[0].z];
            double p3 = y[pos.x + dirs[1].x][pos.y + dirs[1].y][pos.z + dirs[1].z];
            double p4 = y[pos.x + dy * dirs[1].x][pos.y + dy * dirs[1].y][pos.z + dy * dirs[1].z];
            m[0] = y[pos.x][pos.y][pos.z];
            m[4] = (p2 - (1.0 - dx) * m[0] - dx * p1) / (dx * (dx - 1.0));
            m[1] = p1 - m[4] - m[0];
            m[5] = (p4 - (1.0 - dy) * m[0] - dy * p3) / (dy * (dy - 1.0));
            m[2] = p3 - m[5] - m[0];
            double a_1 = y[pos.x + dirs[0].x + dirs[1].x][pos.y + dirs[0].y + dirs[1].y][pos.z + dirs[0].z + dirs[1].z] -
                    m[0] - m[1] - m[2] - m[4] - m[5];
            double a_3 = y[pos.x + dirs[0].x + dy * dirs[1].x][pos.y + dirs[0].y +
                    dy * dirs[1].y][pos.z + dirs[0].z + dy * dirs[1].z] - m[0] - m[1] - m[4] - dy * m[2] - dy * dy * m[5];
            double a_2 = y[pos.x + dx * dirs[0].x + dirs[1].x][pos.y + dx * dirs[0].y +
                    dirs[1].y][pos.z + dx * dirs[0].z + dirs[1].z] - m[0] - dx * m[1] - dx * dx * m[4] - m[2] - m[5];
            double a_4 = y[pos.x + dx * dirs[0].x + dy * dirs[1].x][pos.y + dx * dirs[0].y +
                    dy * dirs[1].y][pos.z + dx * dirs[0].z + dy * dirs[1].z] -
                    m[0] - dx * m[1] - dx * dx * m[4] - dy * m[2] - dy * dy * m[5];
            m[8] = (dx * dy * a_1 - dy * a_2 - dx * a_3 + a_4) /
                   (dx * dy - dx * dy * dy - dx * dx * dy + dx * dx * dy * dy);
            m[7] = (dy * a_2 - a_4 - (dx * dx * dy - dx * dx * dy * dy) * m[8]) / (dx * dy - dx * dy * dy);
            m[6] = (dx * a_3 - a_4 - (dx * dy * dy - dx * dx * dy * dy) * m[8]) / (dx * dy - dx * dx * dy);
            m[3] = (a_4 - dx * dx * dy * dy * m[8] - dx * dy * dy * m[7] - dx * dx * dy * m[6]) / (dx * dy);
        }
        else {
            Node3 pos = d + Node3(1, 1, 1);
            //int dx = 2;
            m[0] = 0;
            if ((dirs[0].x != 0 && pos.x == 1) || (dirs[0].y != 0 && pos.y == 1) || (dirs[0].z != 0 && pos.z == 1)) {
                //dx = -1;
                m[0] = 1;
            }
            //int dy = 2;
            m[1] = 0;
            if ((dirs[1].x != 0 && pos.x == 1) || (dirs[1].y != 0 && pos.y == 1) || (dirs[1].z != 0 && pos.z == 1)) {
                //dy = -1;
                m[1] = 1;
            }
            int px = pos.x;
            if (dirs[0].x != 0 || dirs[1].x != 0)
                px = (dirs[0].x < 0 || dirs[1].x < 0) ? 2 : 0;
            int py = pos.y;
            if (dirs[0].y != 0 || dirs[1].y != 0)
                py = (dirs[0].y < 0 || dirs[1].y < 0) ? 2 : 0;
            int pz = pos.z;
            if (dirs[0].z != 0 || dirs[1].z != 0)
                pz = (dirs[0].z < 0 || dirs[1].z < 0) ? 2 : 0;
            int cont = 2;
            for (int i=0; i<3; i++) {
                for (int j=0; j<3; j++) {
                    m[cont] = y[px + i*dirs[1].x + j*dirs[0].x][py + i*dirs[1].y + j*dirs[0].y][pz + i*dirs[1].z + j*dirs[0].z];
                    cont++;
                }
            }
        }
    }


    template<int wmm_degree> void wmm3DAniSurface3D(Grid3 &image, std::vector<Node3> &initials,
                                                            Node3D &h, int interp, int mode, Grid3 &u_surface) {

        bool isnewpos[3][3][3], is_valid;
        double valcenter[3][3][3], imcenter[3][3][3];
        Node3 dirs[4], neigh, father;
        int ndirs = 3;

        Grid3_<unsigned char> state = Grid3_<unsigned char>(image.rows, image.cols, image.channels);

        std::multimap<double, wmm3_<double, wmm_degree> > trial_set;
        std::map<int, typename std::multimap<double, wmm3_<double, wmm_degree> >::iterator> mapa_trial;

        typename std::multimap<double, wmm3_<double, wmm_degree> >::iterator trial_set_it;
        typename std::map<int, typename std::multimap<double, wmm3_<double, wmm_degree> >::iterator>::iterator mapa_trial_it;
        std::pair<double, wmm3_<double, wmm_degree> > pr_trial;
        std::pair<int, typename std::multimap<double, wmm3_<double, wmm_degree> >::iterator> pr_mapa;

        int key, i;
        wmm3_<double, wmm_degree> winner;
        Node3D fi;

        // Initialization
        for (i = 0; i < (int) initials.size(); i++) {
            key = initials[i].x + u_surface.rows*(initials[i].y + u_surface.cols*initials[i].z);
            if (mapa_trial.find(key) == mapa_trial.end() && u_surface.contains(initials[i])) {
                u_surface.at(initials[i]) = 0.0;
                winner.v = 0.0;
                winner.p = initials[i];
                winner.ndirs = 0;
                winner.father = Node3(-1., -1., -1);
                state.at(initials[i]) = P_TRIAL;
                pr_trial = std::pair<double, wmm3_<double, wmm_degree> >(0.0, winner);
                trial_set_it = trial_set.insert(pr_trial);
                pr_mapa = std::pair<int, typename std::multimap<double, wmm3_<double, wmm_degree> >::iterator>(key, trial_set_it);
                mapa_trial.insert(pr_mapa);
            }
        }

        while (!trial_set.empty()) {

            trial_set_it = trial_set.begin();
            key = trial_set_it->second.p.x + u_surface.rows*(trial_set_it->second.p.y + u_surface.cols*trial_set_it->second.p.z);
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

            trial_set.erase(trial_set_it);
            mapa_trial.erase(mapa_trial_it);

            state.at(winner.p) = P_ALIVE;

            if (image.is_border(winner.p))
                continue;

            if (winner.father.x >= 0) {
                // Father computation
                is_valid = true;
                for (int i=0; i < 3; i++) {
                    for (int j=0; j < 3; j++) {
                        for (int k=0; k < 3; k++) {
                            if ((i == 1) && (j == 1) && (k == 1)) {
                                fi = Node3D(image.at(winner.father, 0), image.at(winner.father, 1),image.at(winner.father, 2));
                                imcenter[1][1][1] = norm(fi);
                                valcenter[1][1][1] = u_surface.at(winner.father);
                                continue;
                            }
                            Node3 d(i - 1, j - 1, k - 1);
                            neigh = winner.father + d;
                            if (!u_surface.contains(neigh)) {
                                valcenter[i][j][k] = winner.v;
                                fi = Node3D(image.at(winner.p, 0), image.at(winner.p, 1),image.at(winner.p, 2));
                            }
                            else {
                                fi = Node3D(image.at(neigh, 0), image.at(neigh, 1),image.at(neigh, 2));
                                valcenter[i][j][k] = u_surface.at(neigh);
                            }
                            imcenter[i][j][k] = norm(fi);
                        }
                    }
                }

                Node3 d = winner.p - winner.father;
                int cont = 0, cont2 = 0;
                int sum = abs(d.x) + abs(d.y) + abs(d.z);
                ndirs = (sum == 3) ? 3 : 4;
                if (ndirs == 3) {
                    dirs[0] = Node3(-d.x, 0, 0);
                    dirs[1] = Node3(0, -d.y, 0);
                    dirs[2] = Node3(0, 0, -d.z);
                }
                else {
                    if ((abs(d.x) == 1) && (sum == 2)) {
                        dirs[0] = Node3(-d.x, 0, 0);
                        cont++; cont2 = 2;
                    }
                    else if (abs(d.x) == 0) {
                        dirs[0] = Node3(1, 0, 0);
                        dirs[2] = Node3(-1, 0, 0);
                        cont++; cont2 = 3;
                    }
                    if ((abs(d.y) == 1) && (sum == 2)) {
                        dirs[cont] = Node3(0, -d.y, 0);
                        cont++;
                    }
                    else if (abs(d.y) == 0) {
                        dirs[cont] = Node3(0, 1, 0);
                        dirs[cont+2] = Node3(0, -1, 0);
                        cont++;
                    }
                    if ((abs(d.z) == 1) && (sum == 2)) {
                        dirs[cont2] = Node3(0, 0, -d.z);
                    }
                    else if (abs(d.z) == 0) {
                        if (cont > 1)
                            dirs[2] = dirs[1];
                        dirs[1] = Node3(0, 0, 1);
                        dirs[3] = Node3(0, 0, -1);
                    }
                }

                winner.ndirs = ndirs;
                winner.v = valcenter[d.x + 1][d.y + 1][d.z + 1];
                for (int m=0; m < ndirs; m++) {
                    winner.planes[m].dirs[0] = dirs[m];
                    winner.planes[m].dirs[1] = dirs[(m+1)%ndirs];
                    Node3 pos1 = d + dirs[m] + Node3(1, 1, 1);
                    Node3 pos2 = d + dirs[(m+1)%ndirs] + Node3(1, 1, 1);
                    Node3 esquina = d + dirs[m] + dirs[(m+1)%ndirs] + Node3(1, 1, 1);
                    winner.planes[m].v[0] = valcenter[pos1.x][pos1.y][pos1.z];
                    winner.planes[m].v[1] = valcenter[pos2.x][pos2.y][pos2.z];
                    winner.planes[m].v[2] = valcenter[esquina.x][esquina.y][esquina.z];
                    setCoeffs3D(valcenter, winner.planes[m].m, d, winner.planes[m].dirs, interp);
                    setCoeffs3D(imcenter, winner.planes[m].fm, d, winner.planes[m].dirs, interp);
                }
            }
            
            // Neighbour temptative value computation
            for (int i=0; i < 3; i++) {
                for (int j=0; j < 3; j++) {
                    for (int k=0; k < 3; k++) {
                        isnewpos[i][j][k] = false;
                        if ((i == 1) && (j == 1) && (k == 1)) {
                            fi = Node3D(image.at(winner.p, 0), image.at(winner.p, 1),image.at(winner.p, 2));
                            imcenter[1][1][1] = norm(fi);
                            valcenter[1][1][1] = u_surface.at(winner.p);
                            continue;
                        }
                        Node3 d(i - 1, j - 1, k - 1);
                        neigh = winner.p + d;
                        if (!u_surface.contains(neigh)) {
                            fi = Node3D(image.at(winner.p, 0), image.at(winner.p, 1),image.at(winner.p, 2));
                            // fi = Node3D(MAX_VAL, MAX_VAL, MAX_VAL);
                            valcenter[i][j][k] = winner.v;
                        }
                        else {
                            fi = Node3D(image.at(neigh, 0), image.at(neigh, 1),image.at(neigh, 2));
                            valcenter[i][j][k] = u_surface.at(neigh);
                        }
                        imcenter[i][j][k] = norm(fi);
                        if (u_surface.contains(neigh) && state.at(neigh) != P_ALIVE) {
                            double val_neigh = GetVal3D(image, u_surface, winner, neigh, h, interp, mode);
                            if (val_neigh < valcenter[i][j][k]) {
                                valcenter[i][j][k] = val_neigh;
                                isnewpos[i][j][k] = true;
                            }
                        }
                    }
                }
            }

            // Update
            for (int i=0; i < 3; i++) {
                for (int j=0; j < 3; j++) {
                    for (int k=0; k < 3; k++) {
                        if (!isnewpos[i][j][k])
                            continue;

                        Node3 d(i - 1, j - 1, k - 1);
                        neigh = winner.p + d;

                        key = neigh.x + u_surface.rows*(neigh.y + u_surface.cols*neigh.z);
                        if (state.at(neigh) == P_TRIAL) {
                            mapa_trial_it = mapa_trial.find(key);
                            trial_set.erase(mapa_trial_it->second);
                            mapa_trial.erase(mapa_trial_it);
                        }
                        else {
                            state.at(neigh) = P_TRIAL;
                        }

                        wmm3_<double, wmm_degree> new_w;
                        new_w.p = neigh;
                        new_w.father = winner.p;

                        pr_trial = std::pair<double, wmm3_<double, wmm_degree> >(valcenter[i][j][k], new_w);
                        trial_set_it = trial_set.insert(pr_trial);
                        pr_mapa = std::pair<int, typename std::multimap<double, wmm3_<double, wmm_degree> >::iterator>(key, trial_set_it);
                        mapa_trial.insert(pr_mapa);

                        u_surface.at(new_w.p) = valcenter[i][j][k];
                    }
                }
            }

        }

        free(state.data);
        return;

    }

    template<int wmm_degree> Grid3 wmm3DAniSurface3D(Grid3 &image, std::vector<Node3> &initials,
                                                            Node3D &h, int interp, int mode) {
        Grid3 u_surface = Grid3(MAX_VAL, image.rows, image.cols, image.channels, 1);
        wmm3DAniSurface3D<wmm_degree>(image, initials, h, interp, mode, u_surface);
        return u_surface;
    }
    
}
