#include "../TYPES/WMMStructs3D.h"
#include <math.h>
#include <stdlib.h>
#include <vector>
#include <map>
#include <iostream>


namespace wmm3D {



    double GetInterpValue(wmm3D::Grid3 &image, wmm3D::AugWmm3 &wave, int ndir,
                          int interp, std::pair<double, double> epsilon, wmm3D::Node3D &neigh, wmm3D::Node3D fn, 
                          wmm3D::Node3D &h, int gamma) {
        int segment1 = (epsilon.first >= 1.0) ? gamma - 1 : (int) floor(gamma * epsilon.first);
        int segment2 = (epsilon.second >= 1.0) ? gamma - 1 : (int) floor(gamma * epsilon.second);
        int r1 = segment2*(gamma + 1) + segment1;
        double step = 1.0/ (double) gamma;
        std::pair<double, double> t_epsilon = std::pair<double, double>((epsilon.first - segment1 * step) / step,
                                                                        (epsilon.second - segment2 * step) / step);
        
        wmm3D::Node3D neighD(h.x*neigh.x, h.y*neigh.y, h.z*neigh.z);
        wmm3D::Node3D plane_posD(h.x*(wave.p.x + epsilon.first*wave.planes[ndir].dirs[0].x + epsilon.second*wave.planes[ndir].dirs[1].x),
                                 h.y*(wave.p.y + epsilon.first*wave.planes[ndir].dirs[0].y + epsilon.second*wave.planes[ndir].dirs[1].y),
                                 h.z*(wave.p.z + epsilon.first*wave.planes[ndir].dirs[0].z + epsilon.second*wave.planes[ndir].dirs[1].z));
        double value = wmm3D::MAX_VAL, ft, vt;
        if (interp == wmm3D::I_BILINEAR) {
            wmm3D::Node3 p1 = wave.p + wave.planes[ndir].dirs[0];
            wmm3D::Node3 p2 = wave.p + wave.planes[ndir].dirs[1];
            wmm3D::Node3 p3 = wave.p + wave.planes[ndir].dirs[0] + wave.planes[ndir].dirs[1];
            wmm3D::Node3D f0 = wmm3D::Node3D(image.at(wave.p, 0), image.at(wave.p, 1), image.at(wave.p, 2)), f1, f2, f3;
            if (image.contains(p1))
                f1 = wmm3D::Node3D(image.at(p1, 0), image.at(p1, 1), image.at(p1, 2));
            else
                f1 = f0;
            if (image.contains(p2))
                f2 = wmm3D::Node3D(image.at(p2, 0), image.at(p2, 1), image.at(p2, 2));
            else
                f2 = f0;
            if (image.contains(p3))
                f3 = wmm3D::Node3D(image.at(p3, 0), image.at(p3, 1), image.at(p3, 2));
            else
                f3 = f0;
            ft = (1.0 - epsilon.first)*(1.0 - epsilon.second)*wmm3D::norm(f0) +
                 epsilon.first*(1.0 - epsilon.second)*wmm3D::norm(f1) +
                 (1.0 - epsilon.first)*epsilon.second*wmm3D::norm(f2) +
                 epsilon.first*epsilon.second*wmm3D::norm(f3);
            vt = (1.0 - t_epsilon.first)*(1.0 - t_epsilon.second)*wave.planes[ndir].v[r1] +
                 t_epsilon.first*(1.0 - t_epsilon.second)*wave.planes[ndir].v[r1 + 1] +
                 (1.0 - t_epsilon.first)*t_epsilon.second*wave.planes[ndir].v[r1 + gamma + 1] +
                 t_epsilon.first*t_epsilon.second*wave.planes[ndir].v[r1 + gamma + 2];
        }
        
        value = vt + wmm3D::norm(neighD - plane_posD) * (ft + wmm3D::norm(fn))/2.0;
        /*if (value < wave.v && interp != wmm3D::I_BILINEAR) {
            value = GetInterpValue(image, wave, ndir, wmm3D::I_BILINEAR, epsilon, neigh, h);
        }*/
        return value;
    }
    
    
    std::pair<double, double> GetEpsilonHopfLax(wmm3D::AugWmm3 &wave, int ndir,
                                                wmm3D::Node3D &neigh, wmm3D::Node3D &fn, wmm3D::Node3D &h) {

        wmm3D::Node3D ph(h.x*wave.p.x, h.y*wave.p.y, h.z*wave.p.z), nh(h.x*neigh.x, h.y*neigh.y, h.z*neigh.z);
        wmm3D::Node3D d0h(h.x*wave.planes[ndir].dirs[0].x, h.y*wave.planes[ndir].dirs[0].y, h.z*wave.planes[ndir].dirs[0].z);
        wmm3D::Node3D d1h(h.x*wave.planes[ndir].dirs[1].x, h.y*wave.planes[ndir].dirs[1].y, h.z*wave.planes[ndir].dirs[1].z);

        double A = d0h.y*d1h.z - d0h.z*d1h.y;
        double B = d0h.z*d1h.x - d0h.x*d1h.z;
        double C = d0h.x*d1h.y - d0h.y*d1h.x;
        double D = -A*ph.x - B*ph.y - C*ph.z;

        double den = A*fn.x + B*fn.y + C*fn.z;

        double t = -(A*nh.x + B*nh.y + C*nh.z + D)/den;
        wmm3D::Node3D x(nh.x + t*fn.x, nh.y + t*fn.y, nh.z + t*fn.z), diff;

        std::pair<double, double> epsilon(0.0, 0.0);

        if (fabs(den) > 0.0)
            diff = x - ph;
        else if (fabs(den) == 0.0 && wmm3D::norm(d0h + d1h) > 0.0) {
            t = (A*nh.x + B*nh.y + C*nh.z + D)/(A*A + B*B + C*C);
            x = wmm3D::Node3D(nh.x - t*A, nh.y - t*B, nh.z - t*C);
            diff = x - ph;
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


    std::pair<double, double> GetEpsilonGradient(wmm3D::AugWmm3 &wave, int ndir, 
                                                 wmm3D::Node3D &neigh, wmm3D::Node3D &fn, wmm3D::Node3D &h) {

        double lambda = 0.0;
        std::pair<double, double> epsilon;
        if (wave.planes[ndir].dirs[0].x == 0 && wave.planes[ndir].dirs[1].x == 0 && fn.x != 0.0) {
            lambda = (neigh.x - wave.p.x)/fn.x;
        }
        else if (wave.planes[ndir].dirs[0].y == 0 && wave.planes[ndir].dirs[1].y == 0 && fn.y != 0.0) {
            lambda = (neigh.y - wave.p.y)/fn.y;
        }
        else if (wave.planes[ndir].dirs[0].z == 0 && wave.planes[ndir].dirs[1].z == 0 && fn.z != 0.0) {
            lambda = (neigh.z - wave.p.z)/fn.z;
        }
        else {
            epsilon = GetEpsilonHopfLax(wave, ndir, neigh, fn, h);
            return epsilon;
        }
        
        wmm3D::Node3D diff(neigh.x - lambda*fn.x - wave.p.x, neigh.y - lambda*fn.y - wave.p.y, neigh.z - lambda*fn.z - wave.p.z);
        
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


    std::pair<double, double> GetEpsilonGoldenSearch(wmm3D::Grid3 &image, wmm3D::AugWmm3 &wave,
                                                     int ndir, int interp, wmm3D::Node3D &neigh, wmm3D::Node3D fn,
                                                     wmm3D::Node3D &h, int gamma) {

        double a_1 = 0.0, b_1 = 1.0, x1_1 = a_1 + (1.0-RESPHI)*(b_1 - a_1), x2_1 = a_1 + RESPHI*(b_1 - a_1),
                f_x1 = MAX_VAL, f_x2 = MAX_VAL;

        double a_2 = 0.0, b_2 = 1.0, x1_2 = a_2 + (1.0-RESPHI)*(b_2 - a_2), x2_2 = a_2 + RESPHI*(b_2 - a_2);

        std::pair<double, double> epsilon(0.0, 0.0);

        double res = GetInterpValue(image, wave, ndir, interp, std::pair<double, double>(0.0, 0.0), neigh, fn, h, gamma);
        double f_b = GetInterpValue(image, wave, ndir, interp, std::pair<double, double>(1.0, 0.0), neigh, fn, h, gamma);
        if (f_b < res) {
            res = f_b; epsilon = std::pair<double, double>(1.0, 0.0);
        }
        f_b = GetInterpValue(image, wave, ndir, interp, std::pair<double, double>(0.0, 1.0), neigh, fn, h, gamma);
        if (f_b < res) {
            res = f_b; epsilon = std::pair<double, double>(0.0, 1.0);
        }
        f_b = GetInterpValue(image, wave, ndir, interp, std::pair<double, double>(1.0, 1.0), neigh, fn, h, gamma);
        if (f_b < res) {
            res = f_b; epsilon = std::pair<double, double>(1.0, 1.0);
        }

        f_x1 = GetInterpValue(image, wave, ndir, interp, std::pair<double, double>(0.0, x1_2), neigh, fn, h, gamma);
        f_x2 = GetInterpValue(image, wave, ndir, interp, std::pair<double, double>(0.0, x2_2), neigh, fn, h, gamma);

        while (fabs(b_1 - a_1) > TAU && fabs(b_2 - a_2) > TAU) {
            if(f_x1 < f_x2) {
                f_x1 = GetInterpValue(image, wave, ndir, interp, std::pair<double, double>(x1_1, x1_2), neigh, fn, h, gamma);
                f_x2 = GetInterpValue(image, wave, ndir, interp, std::pair<double, double>(x2_1, x1_2), neigh, fn, h, gamma);
                b_2 = x2_2; x2_2 = x1_2; x1_2 = a_2 + (1.0 - RESPHI)*(b_2 - a_2);
            }
            else {
                f_x1 = GetInterpValue(image, wave, ndir, interp, std::pair<double, double>(x1_1, x2_2), neigh, fn, h, gamma);
                f_x2 = GetInterpValue(image, wave, ndir, interp, std::pair<double, double>(x2_1, x2_2), neigh, fn, h, gamma);
                a_2 = x1_2; x1_2 = x2_2; x2_2 = a_2 + RESPHI*(b_2 - a_2);
            }
            if(f_x1 < f_x2) {
                f_x1 = GetInterpValue(image, wave, ndir, interp, std::pair<double, double>(x1_1, x1_2), neigh, fn, h, gamma);
                f_x2 = GetInterpValue(image, wave, ndir, interp, std::pair<double, double>(x1_1, x2_2), neigh, fn, h, gamma);
                b_1 = x2_1; x2_1 = x1_1; x1_1 = a_1 + (1.0 - RESPHI)*(b_1 - a_1);
            }
            else {
                f_x1 = GetInterpValue(image, wave, ndir, interp, std::pair<double, double>(x2_1, x1_2), neigh, fn, h, gamma);
                f_x2 = GetInterpValue(image, wave, ndir, interp, std::pair<double, double>(x2_1, x2_2), neigh, fn, h, gamma);
                a_1 = x1_1; x1_1 = x2_1; x2_1 = a_1 + RESPHI*(b_1 - a_1);
            }
        }

        if (f_x1 < std::min(res, f_x2))
            epsilon = std::pair<double, double>(x1_1, x1_2);
        else if (f_x2 < std::min(res, f_x1))
            epsilon = std::pair<double, double>(x2_1, x2_2);

        return epsilon;
        
    }


    double GetVal3D(wmm3D::Grid3 &image, wmm3D::Grid3 &u_surface, wmm3D::AugWmm3 &wave,
                    wmm3D::Node3D &neigh, wmm3D::Node3D &fn, wmm3D::Node3D &h, int interp, int mode, int gamma) {

        wmm3D::Node3D f0(image.at(wave.p, 0), image.at(wave.p, 1), image.at(wave.p, 2));
        double y0 = wave.planes[0].v[0];

        if (isinf(wmm3D::norm(f0)) || isnan(wmm3D::norm(f0)))
            f0 = fn;

        wmm3D::Node3D diff(h.x * (neigh.x - wave.p.x), h.y * (neigh.y - wave.p.y), h.z * (neigh.z - wave.p.z));
        double val = y0 + wmm3D::norm(diff) * (wmm3D::norm(f0) + wmm3D::norm(fn)) / 2.0;
        if (wave.ndirs > 0.0) {

            for (int ndir = 0; ndir < wave.ndirs; ndir++) {
                wmm3D::Node3 p0 = wave.p + wave.planes[ndir].dirs[0];
                wmm3D::Node3 p1 = wave.p + wave.planes[ndir].dirs[1];

                if (!image.contains(p0) || !image.contains(p1))
                    continue;
                std::pair<double, double> epsilon;
                if (mode == wmm3D::M_GRADIENT)
                    epsilon = GetEpsilonGradient(wave, ndir, neigh, fn, h);
                else if (mode == wmm3D::M_HOPFLAX)
                    epsilon = GetEpsilonHopfLax(wave, ndir, neigh, fn, h);
                else
                    epsilon = GetEpsilonGoldenSearch(image, wave, ndir, interp, neigh, fn, h, gamma);
                val = std::min(val, GetInterpValue(image, wave, ndir, interp, epsilon, neigh, fn, h, gamma));
            }
        }
        return val;

    }


    void setCoeffs3D(double y[3][3][3], double *m, wmm3D::Node3 d, wmm3D::Node3 *dirs, int interp) {

        if (interp == wmm3D::I_QUADATRIC) {
            wmm3D::Node3 pos = d + wmm3D::Node3(1, 1, 1);
            wmm3D::Node3 pos1 = d + dirs[0] + wmm3D::Node3(1, 1, 1);
            wmm3D::Node3 pos2 = d + dirs[1] + wmm3D::Node3(1, 1, 1);
            wmm3D::Node3 pos3 = d + dirs[0] + dirs[1] + wmm3D::Node3(1, 1, 1);
            m[0] = y[pos.x][pos.y][pos.z];
            m[1] = y[pos1.x][pos1.y][pos1.z] - m[0];
            m[2] = y[pos2.x][pos2.y][pos2.z] - m[0];
            m[3] = y[pos3.x][pos3.y][pos3.z] - m[2] - m[1] - m[0];
        }
        else if (interp == wmm3D::I_CUBIC) {
            wmm3D::Node3 pos = d + wmm3D::Node3(1, 1, 1);
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
        else if (interp == wmm3D::I_BICUBIC) {
            wmm3D::Node3 pos = d + wmm3D::Node3(1, 1, 1);
            int dx = 2;
            if ((dirs[0].x != 0 && pos.x == 1) || (dirs[0].y != 0 && pos.y == 1) || (dirs[0].z != 0 && pos.z == 1))
                dx = -1;
            int dy = 2;
            if ((dirs[1].x != 0 && pos.x == 1) || (dirs[1].y != 0 && pos.y == 1) || (dirs[1].z != 0 && pos.z == 1))
                dy = -1;
            double f00 = y[pos.x][pos.y][pos.z];
            double f10 = y[pos.x + dirs[0].x][pos.y + dirs[0].y][pos.z + dirs[0].z];
            double f01 = y[pos.x + dirs[1].x][pos.y + dirs[1].y][pos.z + dirs[1].z];
            double f11 = y[pos.x + dirs[0].x + dirs[1].x][pos.y + dirs[0].y + dirs[1].y][pos.z + dirs[0].z + dirs[1].z];
            double fx00 = (dx == 2) ? f10 - f00 : (f10 - y[pos.x - dirs[0].x][pos.y - dirs[0].y][pos.z - dirs[0].z])/2.0;
            double fx10 = (dx == 2) ? (y[pos.x + dx * dirs[0].x][pos.y + dx * dirs[0].y][pos.z + dx * dirs[0].z] - f00)/2.0 : f10 - f00;
            double fx01 = (dx == 2) ? f11 - f01 : (f11 - y[pos.x - dirs[0].x + dirs[1].x][pos.y - dirs[0].y + dirs[1].y][pos.z - dirs[0].z + dirs[1].z])/2.0;
            double fx11 = (dx == 2) ? (y[pos.x + dx * dirs[0].x + dirs[1].x][pos.y + dx * dirs[0].y + dirs[1].y][pos.z + dx * dirs[0].z + dirs[1].z] - f01)/2.0 : f11 - f01;
            
            double fy00 = (dy == 2) ? f01 - f00 : (f01 - y[pos.x - dirs[1].x][pos.y - dirs[1].y][pos.z - dirs[1].z])/2.0;
            double fy10 = (dy == 2) ? (y[pos.x + dy * dirs[1].x][pos.y + dy * dirs[1].y][pos.z + dy * dirs[1].z] - f00)/2.0 : f01 - f00;
            double fy01 = (dy == 2) ? f11 - f10 : (f11 - y[pos.x + dirs[0].x - dirs[1].x][pos.y + dirs[0].y - dirs[1].y][pos.z + dirs[0].z - dirs[1].z])/2.0;
            double fy11 = (dy == 2) ? (y[pos.x + dirs[0].x + dy * dirs[1].x][pos.y + dirs[0].y + dy * dirs[1].y][pos.z + dirs[0].z + dy * dirs[1].z] - f10)/2.0 : f11 - f10;
            
            double fxy00 = 0.5*(fy10 - fy00 + fx01 - fx00);
            double fxy10 = 0.5*(fy10 - fy00 + fx11 - fx10);
            double fxy01 = 0.5*(fy11 - fy01 + fx01 - fx00);
            double fxy11 = 0.5*(fy11 - fy01 + fx11 - fx10);
            
            double a[4][4] = {{f00, f01, fy00, fy01}, {f10, f11, fy10, fy11}, {fx00, fx01, fxy00, fxy01}, {fx10, fx11, fxy10, fxy11}};
            
            double b[4][4] = {{a[0][0], a[0][2], -3*a[0][0] + 3*a[0][1] - 2*a[0][2] - a[0][3], 2*a[0][0] - 2*a[0][1] + a[0][2] + a[0][3]},
                              {a[1][0], a[1][2], -3*a[1][0] + 3*a[1][1] - 2*a[1][2] - a[1][3], 2*a[1][0] - 2*a[1][1] + a[1][2] + a[1][3]},
                              {a[2][0], a[2][2], -3*a[2][0] + 3*a[2][1] - 2*a[2][2] - a[2][3], 2*a[2][0] - 2*a[2][1] + a[2][2] + a[2][3]},
                              {a[3][0], a[3][2], -3*a[3][0] + 3*a[3][1] - 2*a[3][2] - a[3][3], 2*a[3][0] - 2*a[3][1] + a[3][2] + a[3][3]}};
            
            double m_aux[16] = {b[0][0], b[0][1], b[0][2], b[0][3],
                 b[2][0], b[2][1], b[2][2], b[2][3],
                 -3*b[0][0] + 3*b[1][0] - 2*b[2][0] - b[3][0], -3*b[0][1] + 3*b[1][1] - 2*b[2][1] - b[3][1], 
                 -3*b[0][2] + 3*b[1][2] - 2*b[2][2] - b[3][2], -3*b[0][3] + 3*b[1][3] - 2*b[2][3] - b[3][3],
                 2*b[0][0] - 2*b[1][0] + b[2][0] + b[3][0], 2*b[0][1] - 2*b[1][1] + b[2][1] + b[3][1],
                 2*b[0][2] - 2*b[1][2] + b[2][2] + b[3][2], 2*b[0][3] - 2*b[1][3] + b[2][3] + b[3][3]};
            
            std::copy(m_aux, m_aux + 16, m);
        }
    }


    void augWmm3DIsoSurface3D(wmm3D::Grid3 &image, std::vector<wmm3D::Node3> &initials, wmm3D::Node3D &h,
                                      int interp, int mode, int N, int M, int gamma, wmm3D::Grid3 &u_surface) {

        bool isnewpos[3][3][3];
        double valcenter[3][3][3], imcenter[3][3][3];
        wmm3D::Node5D imnodes[3][3][3];
        wmm3D::Node3 dirs[4], neigh;
        wmm3D::Node3D neighD, neighD_aux, d1, d2;
        int ndirs = 3, cont = 0;
        double step = 1.0 / (double) gamma, epsilon1, epsilon2;

        wmm3D::Grid3_<unsigned char> state = wmm3D::Grid3_<unsigned char>(image.rows, image.cols, image.channels);

        std::multimap<double, wmm3D::AugWmm3 > trial_set;
        std::map<int, typename std::multimap<double, wmm3D::AugWmm3 >::iterator> mapa_trial;

        typename std::multimap<double, wmm3D::AugWmm3 >::iterator trial_set_it;
        typename std::map<int, typename std::multimap<double, wmm3D::AugWmm3 >::iterator>::iterator mapa_trial_it;
        std::pair<double, wmm3D::AugWmm3> pr_trial;
        std::pair<int, typename std::multimap<double, wmm3D::AugWmm3 >::iterator> pr_mapa;

        int key, i;
        wmm3D::AugWmm3 winner, new_w;
        wmm3D::Node3D fi;


        // Initialization
        for (i = 0; i < (int) initials.size(); i++) {
            key = initials[i].y + u_surface.rows*(initials[i].x + u_surface.cols*initials[i].z);
            if (mapa_trial.find(key) == mapa_trial.end() && u_surface.contains(initials[i])) {
                u_surface.at(initials[i]) = 0.0;
                winner = wmm3D::AugWmm3(N, M, gamma, 1);
                winner.planes[0].v[0] = 0.0;
                winner.p = initials[i];
                winner.ndirs = 0;
                state.at(initials[i]) = wmm3D::P_TRIAL;
                pr_trial = std::pair<double, wmm3D::AugWmm3 >(0.0, winner);
                trial_set_it = trial_set.insert(pr_trial);
                pr_mapa = std::pair<int, typename std::multimap<double, wmm3D::AugWmm3 >::iterator>(key, trial_set_it);
                mapa_trial.insert(pr_mapa);
            }
        }

        while (!trial_set.empty()) {

            trial_set_it = trial_set.begin();
            key = trial_set_it->second.p.y + u_surface.rows*(trial_set_it->second.p.x + u_surface.cols*trial_set_it->second.p.z);
            mapa_trial_it = mapa_trial.find(key);

            if (mapa_trial_it == mapa_trial.end()) {
                printf("ERROR: bad map alloc\n");
                return;
            }

            if (mapa_trial_it->second != trial_set_it) {
                printf("ERROR: bad trial/map alloc\n");
                return;
            }

            winner = trial_set_it->second;
            /*std::cout << "WINNER = ( " << winner.p.x << ", " << winner.p.y << ", " << winner.p.z << " ) v = " << winner.planes[0].v[0] << std::endl;
            for (int f=0; f<winner.ndirs; f++) {
                std::cout << "PLANE " << f << " -> ";
                for (int t=0; t<(gamma + 1)*(gamma + 1); t++) {
                    std::cout << winner.planes[f].v[t] << " ";
                }
                std::cout << std::endl;
            }*/
            /*std::cout << "VAL = " << winner.v << std::endl;
            std::cout << "NDIRS = " << winner.ndirs << std::endl;
            for (int f=0; f<winner.ndirs; f++) {
                std::cout << "PLANE " << f << std::endl;
                std::cout <<"\t DIR 0 = ( " << winner.planes[f].dirs[0].x << ", " << winner.planes[f].dirs[0].y << ", " << winner.planes[f].dirs[0].z << " )" << std::endl;
                std::cout <<"\t DIR 1 = ( " << winner.planes[f].dirs[1].x << ", " << winner.planes[f].dirs[1].y << ", " << winner.planes[f].dirs[1].z << " )" << std::endl;
                std::cout << "\t VALS = " << winner.planes[f].v[0] << ", " << winner.planes[f].v[1] << ", " << winner.planes[f].v[2] << std::endl;
            }
            std::cout << std::endl << std::endl;*/

            trial_set.erase(trial_set_it);
            mapa_trial.erase(mapa_trial_it);

            state.at(winner.p) = wmm3D::P_ALIVE;
            // Neighbour temptative value computation
            for (int i=0; i < 3; i++) {
                for (int j=0; j < 3; j++) {
                    for (int k=0; k < 3; k++) {
                        isnewpos[i][j][k] = false;
                        if ((i == 1) && (j == 1) && (k == 1)) {
                            imnodes[1][1][1] = wmm3D::to_spherical(wmm3D::Node3D(image.at(winner.p, 0), image.at(winner.p, 1),image.at(winner.p, 2)));
                            imcenter[1][1][1] = wmm3D::norm(fi);
                            valcenter[1][1][1] = u_surface.at(winner.p);
                            continue;
                        }
                        wmm3D::Node3 d(i - 1, j - 1, k - 1);
                        neigh = winner.p + d;
                        isnewpos[i][j][k] = false;
                        if (!u_surface.contains(neigh)) {
                            fi = wmm3D::Node3D(image.at(winner.p, 0), image.at(winner.p, 1),image.at(winner.p, 2));
                            // fi = wmm3D::Node3D(MAX_VAL, MAX_VAL, MAX_VAL);
                            valcenter[i][j][k] = winner.planes[0].v[0];
                        }
                        else {
                            fi = wmm3D::Node3D(image.at(neigh, 0), image.at(neigh, 1),image.at(neigh, 2));
                            valcenter[i][j][k] = u_surface.at(neigh);
                        }
                        imcenter[i][j][k] = wmm3D::norm(fi);
                        imnodes[i][j][k] = wmm3D::to_spherical(fi);
                        if (u_surface.contains(neigh) && state.at(neigh) != P_ALIVE) {
                            neighD = wmm3D::Node3D(neigh.x, neigh.y, neigh.z);
                            double val_neigh = GetVal3D(image, u_surface, winner, neighD, fi, h, interp, mode, gamma);
                            if (val_neigh < valcenter[i][j][k]) {
                                valcenter[i][j][k] = val_neigh;
                                isnewpos[i][j][k] = true;
                            }
                        }
                        //std::cout << winner.p << neigh << imcenter[i][j][k] << valcenter[i][j][k] << h << std::endl;
                    }
                }
            }

            // Update
            for (int i=0; i < 3; i++) {
                for (int j=0; j < 3; j++) {
                    for (int k=0; k < 3; k++) {
                        if (!isnewpos[i][j][k])
                            continue;

                        wmm3D::Node3 d(i - 1, j - 1, k - 1);
                        neigh = winner.p + d;

                        key = neigh.y + u_surface.rows*(neigh.x + u_surface.cols*neigh.z);
                        
                        if (state.at(neigh) == wmm3D::P_TRIAL) {
                            mapa_trial_it = mapa_trial.find(key);
                            trial_set.erase(mapa_trial_it->second);
                            mapa_trial.erase(mapa_trial_it);
                        }
                        else {
                            state.at(neigh) = wmm3D::P_TRIAL;
                        }

                        int cont = 0, cont2 = 0;
                        int sum = abs(d.x) + abs(d.y) + abs(d.z);
                        ndirs = (sum == 3) ? 3 : 4;
                        if (ndirs == 3) {
                            dirs[0] = wmm3D::Node3(-d.x, 0, 0);
                            dirs[1] = wmm3D::Node3(0, -d.y, 0);
                            dirs[2] = wmm3D::Node3(0, 0, -d.z);
                        }
                        else {
                            if ((abs(d.x) == 1) && (sum == 2)) {
                                dirs[0] = wmm3D::Node3(-d.x, 0, 0);
                                cont++; cont2 = 2;
                            }
                            else if (abs(d.x) == 0) {
                                dirs[0] = wmm3D::Node3(1, 0, 0);
                                dirs[2] = wmm3D::Node3(-1, 0, 0);
                                cont++; cont2 = 3;
                            }
                            if ((abs(d.y) == 1) && (sum == 2)) {
                                dirs[cont] = wmm3D::Node3(0, -d.y, 0);
                                cont++;
                            }
                            else if (abs(d.y) == 0) {
                                dirs[cont] = wmm3D::Node3(0, 1, 0);
                                dirs[cont+2] = wmm3D::Node3(0, -1, 0);
                                cont++;
                            }
                            if ((abs(d.z) == 1) && (sum == 2)) {
                                dirs[cont2] = wmm3D::Node3(0, 0, -d.z);
                            }
                            else if (abs(d.z) == 0) {
                                if (cont > 1)
                                    dirs[2] = dirs[1];
                                dirs[1] = wmm3D::Node3(0, 0, 1);
                                dirs[3] = wmm3D::Node3(0, 0, -1);
                            }
                        }
                        
                        new_w = wmm3D::AugWmm3(N, M, gamma, ndirs);
                        new_w.p = neigh;
                        
                        for (int m=0; m < ndirs; m++) {
                            new_w.planes[m].dirs[0] = dirs[m];
                            new_w.planes[m].dirs[1] = dirs[(m+1)%ndirs];
                            
                            wmm3D::Node3 pos1 = neigh + dirs[m] - winner.p + wmm3D::Node3(1, 1, 1);
                            wmm3D::Node3 pos2 = neigh + dirs[(m+1)%ndirs] - winner.p + wmm3D::Node3(1, 1, 1);
                            wmm3D::Node3 esquina = neigh + dirs[m] + dirs[(m+1)%ndirs] - winner.p + wmm3D::Node3(1, 1, 1);
                            //std::cout << neigh << pos1 << pos2 << esquina << std::endl;
                            cont = 0;
                            neighD = wmm3D::Node3D(neigh.x, neigh.y, neigh.z);
                            d1 = wmm3D::Node3D(dirs[m].x, dirs[m].y, dirs[m].z);
                            d2 = wmm3D::Node3D(dirs[(m+1)%ndirs].x, dirs[(m+1)%ndirs].y, dirs[(m+1)%ndirs].z);
                            for (int v=0; v <= gamma; v++) {
                                for (int w=0; w <= gamma; w++) {
                                    if (v == 0 && w == 0)
                                        new_w.planes[m].v[cont] = valcenter[i][j][k];
                                    else if (v == 0 && w == gamma)
                                        new_w.planes[m].v[cont] = valcenter[pos1.x][pos1.y][pos1.z];
                                    else if (v == gamma && w == 0)
                                        new_w.planes[m].v[cont] = valcenter[pos2.x][pos2.y][pos2.z];
                                    else if (v == gamma && w == gamma)
                                        new_w.planes[m].v[cont] = valcenter[esquina.x][esquina.y][esquina.z];
                                    else {
                                        epsilon1 = w*step;
                                        epsilon2 = v*step;
                                        neighD_aux = neighD + epsilon1*d1 + epsilon2*d2;
                                        if (u_surface.contains(neighD_aux)) {
                                            fi = wmm3D::to_cartesian((1.0 - epsilon1)*(1.0 - epsilon2)*imnodes[i][j][k] +
                                                 epsilon1*(1.0 - epsilon2)*imnodes[pos1.x][pos1.y][pos1.z] +
                                                 (1.0 - epsilon1)*epsilon2*imnodes[pos2.x][pos2.y][pos2.z] +
                                                 epsilon1*epsilon2*imnodes[esquina.x][esquina.y][esquina.z]);
                                            new_w.planes[m].v[cont] = GetVal3D(image, u_surface, winner, neighD_aux, fi, h, interp, mode, gamma);
                                        }
                                        else {
                                            new_w.planes[m].v[cont] = valcenter[1][1][1];
                                        }
                                        
                                    }
                                    cont++;
                                }
                            }
                            /*if (interp != wmm3D::I_BILINEAR) {
                                setCoeffs3D(valcenter, new_w.planes[m].m, d, new_w.planes[m].dirs, interp);
                                setCoeffs3D(imcenter, new_w.planes[m].fm, d, new_w.planes[m].dirs, interp);
                            }*/
                                
                        }
                        pr_trial = std::pair<double, wmm3D::AugWmm3 >(valcenter[i][j][k], new_w);
                        trial_set_it = trial_set.insert(pr_trial);
                        pr_mapa = std::pair<int, typename std::multimap<double, wmm3D::AugWmm3 >::iterator>(key, trial_set_it);
                        mapa_trial.insert(pr_mapa);
                        u_surface.at(new_w.p) = valcenter[i][j][k];
                        /*if (new_w.p.x == 1 && new_w.p.y == 0 && new_w.p.z == 2) {
                            std::cout << "UPDATE = ( " << new_w.p.x << ", " << new_w.p.y << ", " << new_w.p.z << " )" <<
                            std::endl;
                            std::cout << "VAL = " << new_w.v << std::endl;
                            std::cout << "NDIRS = " << new_w.ndirs << std::endl;
                            for (int f = 0; f < new_w.ndirs; f++) {
                                std::cout << "PLANE " << f << std::endl;
                                std::cout << "\t DIR 0 = ( " << new_w.planes[f].dirs[0].x << ", " <<
                                new_w.planes[f].dirs[0].y << ", " << new_w.planes[f].dirs[0].z << " )" << std::endl;
                                std::cout << "\t DIR 1 = ( " << new_w.planes[f].dirs[1].x << ", " <<
                                new_w.planes[f].dirs[1].y << ", " << new_w.planes[f].dirs[1].z << " )" << std::endl;
                                std::cout << "\t VALS = " << new_w.planes[f].v[0] << ", " << new_w.planes[f].v[1] <<
                                ", " << new_w.planes[f].v[2] << std::endl;
                            }
                            std::cout << std::endl << std::endl;
                        }*/
                    }
                }
            }
            //del(&winner);

        }

        free(state.data);
        return;

    }

    wmm3D::Grid3 augWmm3DIsoSurface3D(wmm3D::Grid3 &image, std::vector<wmm3D::Node3> &initials,
                                   wmm3D::Node3D &h, int interp, int mode, int N, int M, int gamma) {
        wmm3D::Grid3 u_surface = wmm3D::Grid3(wmm3D::MAX_VAL, image.rows, image.cols, image.channels, 1);
        augWmm3DIsoSurface3D(image, initials, h, interp, mode, N, M, gamma, u_surface);
        return u_surface;
    }
    
}
