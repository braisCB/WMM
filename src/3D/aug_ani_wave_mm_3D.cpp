#include "../TYPES/WMMStructs3D.h"
#include "../TYPES/utils3D.h"
#include <math.h>
#include <stdlib.h>
#include <vector>
#include <map>
#include <unordered_map>
#include <iostream>


namespace wmm3D {

    Node5D checkNodeSlope(Node5D &fd0, Node5D &fd1, int isborder) {
        Node5D node = Node5D(
            checkSlope(fd0.v, fd1.v, isborder),
            checkSlope(fd0.w, fd1.w, isborder),
            checkSlope(fd0.x, fd1.x, isborder),
            checkSlope(fd0.y, fd1.y, isborder),
            checkSlope(fd0.z, fd1.z, isborder)
        );
        return node;
    }

    Node3D GetFValue(Grid3 &image, AugWmm3 &wave, int ndir, int interp, std::pair<double, double> epsilon, Node3D &h) {

        Node3D value;
        Node5D ft;

        Node3 p1 = wave.p + wave.planes[ndir].dirs[0];
        Node3 p2 = wave.p + wave.planes[ndir].dirs[1];
        Node3 p3 = wave.p + wave.planes[ndir].dirs[0] + wave.planes[ndir].dirs[1];
        Node3D f0_3, f1_3, f2_3, f3_3;

        f0_3 = Node3D(image.at(wave.p, 0), image.at(wave.p, 1), image.at(wave.p, 2));
        f1_3 = (image.contains(p1)) ? Node3D(image.at(p1, 0), image.at(p1, 1), image.at(p1, 2)) : f0_3;
        f2_3 = (image.contains(p2)) ? Node3D(image.at(p2, 0), image.at(p2, 1), image.at(p2, 2)) : f0_3;
        f3_3 = (image.contains(p3)) ? Node3D(image.at(p3, 0), image.at(p3, 1), image.at(p3, 2)) : f0_3;

        Node5D f0_5 = to_spherical(f0_3);
        Node5D f1_5 = to_spherical(f1_3);
        Node5D f2_5 = to_spherical(f2_3);
        Node5D f3_5 = to_spherical(f3_3);

        if (interp == I_BILINEAR) {
            ft = (1.0 - epsilon.first)*(1.0 - epsilon.second)*f0_5 +
                 epsilon.first*(1.0 - epsilon.second)*f1_5 +
                 (1.0 - epsilon.first)*epsilon.second*f2_5 +
                 epsilon.first*epsilon.second*f3_5;
        }
        else if (interp == I_QUADATRIC) {
            ft = wave.planes[ndir].fm[0] + epsilon.first*wave.planes[ndir].fm[1] +
                 epsilon.second*wave.planes[ndir].fm[2] + epsilon.first*epsilon.second*wave.planes[ndir].fm[3];
        }
        else if (interp == I_CUBIC) {
            ft = wave.planes[ndir].fm[0] + epsilon.first*wave.planes[ndir].fm[1] +
                 epsilon.second*wave.planes[ndir].fm[2] + epsilon.first*epsilon.second*wave.planes[ndir].fm[3] +
                 epsilon.first*epsilon.first*wave.planes[ndir].fm[4] + epsilon.second*epsilon.second*wave.planes[ndir].fm[5] +
                 epsilon.first*epsilon.first*epsilon.second*wave.planes[ndir].fm[6] +
                 epsilon.first*epsilon.second*epsilon.second*wave.planes[ndir].fm[7] +
                 epsilon.first*epsilon.first*epsilon.second*epsilon.second*wave.planes[ndir].fm[8];
        }
        else if (interp == I_HERMITE) {
            Node5D fm[2], f[3], f0, f1;
            double t = epsilon.first, t_2 = t*t, t_3 = t_2*t;
            int pos = (int) wave.planes[ndir].fm[0].v;
            for (int p=0; p<3; p++) {
                if (pos == 0) {
                    fm[0] = wave.planes[ndir].fm[2 + 3*p + 1] - wave.planes[ndir].fm[2 + 3*p];
                    fm[1] = 1.0 / 2.0 * (wave.planes[ndir].fm[2 + 3*p + 2] - wave.planes[ndir].fm[2 + 3*p]);
                }
                else {
                    fm[0] = 1.0 / 2.0 * (wave.planes[ndir].fm[2 + 3*p + 2] - wave.planes[ndir].fm[2 + 3*p]);
                    fm[1] = wave.planes[ndir].fm[2 + 3*p + 2] - wave.planes[ndir].fm[2 + 3*p + 1];
                }
                f0 = wave.planes[ndir].fm[2 + 3*p + pos];
                f1 = wave.planes[ndir].fm[2 + 3*p + 1 + pos];
                f[p] = (2.0 * t_3 - 3.0 * t_2 + 1.0) * f0 + (t_3 - 2.0 * t_2 + t) * fm[0] +
                       (-2.0 * t_3 + 3.0 * t_2) * f1 + (t_3 - t_2) * fm[1];
            }
            t = epsilon.second; t_2 = t*t; t_3 = t_2*t;
            pos = (int) wave.planes[ndir].fm[1].v;
            if (pos == 0) {
                fm[0] = f[1] - f[0];
                fm[1] = 1.0 / 2.0 * (f[2] - f[0]);
            }
            else {
                fm[0] = 1.0 / 2.0 * (f[2] - f[0]);
                fm[1] = f[2] - f[1];
            }
            f0 = f[pos];
            f1 = f[1 + pos];
            ft = (2.0 * t_3 - 3.0 * t_2 + 1.0) * f0 + (t_3 - 2.0 * t_2 + t) * fm[0] +
                   (-2.0 * t_3 + 3.0 * t_2) * f1 + (t_3 - t_2) * fm[1];

            t = epsilon.second; t_2 = t*t; t_3 = t_2*t;
            pos = (int) wave.planes[ndir].fm[1].v;
            for (int p=0; p<3; p++) {
                if (pos == 0) {
                    fm[0] = wave.planes[ndir].fm[2 + 3 + p] - wave.planes[ndir].fm[2 + p];
                    fm[1] = 1.0 / 2.0 * (wave.planes[ndir].fm[2 + 6 + p] - wave.planes[ndir].fm[2 + p]);
                }
                else {
                    fm[0] = 1.0 / 2.0 * (wave.planes[ndir].fm[2 + 6 + p] - wave.planes[ndir].fm[2 + p]);
                    fm[1] = wave.planes[ndir].fm[2 + 6 + p] - wave.planes[ndir].fm[2 + 3 + p];
                }
                f0 = wave.planes[ndir].fm[2 + 3*pos + p];
                f1 = wave.planes[ndir].fm[2 + 3 + 3*pos + p];
                f[p] = (2.0 * t_3 - 3.0 * t_2 + 1.0) * f0 + (t_3 - 2.0 * t_2 + t) * fm[0] +
                       (-2.0 * t_3 + 3.0 * t_2) * f1 + (t_3 - t_2) * fm[1];
            }
            t = epsilon.first; t_2 = t*t; t_3 = t_2*t;
            pos = (int) wave.planes[ndir].fm[0].v;
            if (pos == 0) {
                fm[0] = f[1] - f[0];
                fm[1] = 1.0 / 2.0 * (f[2] - f[0]);
            }
            else {
                fm[0] = 1.0 / 2.0 * (f[2] - f[0]);
                fm[1] = f[2] - f[1];
            }
            f0 = f[pos];
            f1 = f[1 + pos];
            ft = 0.5 * (ft + (2.0 * t_3 - 3.0 * t_2 + 1.0) * f0 + (t_3 - 2.0 * t_2 + t) * fm[0] +
                   (-2.0 * t_3 + 3.0 * t_2) * f1 + (t_3 - t_2) * fm[1]);
        }
        else if (interp == I_PCHIP) {
            Node5D fm[2], f[3], f0, f1, fd0, fd1;
            double t = epsilon.first, t_2 = t*t, t_3 = t_2*t;
            int pos = (int) wave.planes[ndir].fm[0].v;
            for (int p=0; p<3; p++) {
                fd0 = wave.planes[ndir].fm[2 + 3*p + 1] - wave.planes[ndir].fm[2 + 3*p];
                fd1 = wave.planes[ndir].fm[2 + 3*p + 2] - wave.planes[ndir].fm[2 + 3*p + 1];
                //std::cout << "fd0: " << fd0.v << ", " << fd0.w << ", " << fd0.x << ", " << fd0.y << ", " << fd0.z << std::endl;
                //std::cout << "fd1: " << fd1.v << ", " << fd1.w << ", " << fd1.x << ", " << fd1.y << ", " << fd1.z << std::endl;
                if (pos == 0) {
                    fm[0] = checkNodeSlope(fd0, fd1, 1);
                    fm[1] = checkNodeSlope(fd0, fd1, 0);
                }
                else {
                    fm[0] = checkNodeSlope(fd1, fd0, 0);
                    fm[1] = checkNodeSlope(fd1, fd0, 1);
                }
                //std::cout << "fm[0]: " << fm[0].v << ", " << fm[0].w << ", " << fm[0].x << ", " << fm[0].y << ", " << fm[0].z << std::endl;
                //std::cout << "fm[1]: " << fm[1].v << ", " << fm[1].w << ", " << fm[1].x << ", " << fm[1].y << ", " << fm[1].z << std::endl;
                f0 = wave.planes[ndir].fm[2 + 3*p + pos];
                f1 = wave.planes[ndir].fm[2 + 3*p + 1 + pos];
                //std::cout << "f0: " << f0.v << ", " << f0.w << ", " << f0.x << ", " << f0.y << ", " << f0.z << std::endl;
                //std::cout << "f1: " << f1.v << ", " << f1.w << ", " << f1.x << ", " << f1.y << ", " << f1.z << std::endl;
                f[p] = (2.0 * t_3 - 3.0 * t_2 + 1.0) * f0 + (t_3 - 2.0 * t_2 + t) * fm[0] +
                       (-2.0 * t_3 + 3.0 * t_2) * f1 + (t_3 - t_2) * fm[1];
                /*std::cout << "N: " << wave.N << std::endl;
                for (int i = 0; i < wave.N; i++)
                    std::cout << wave.planes[ndir].fm[i].v << ", " << wave.planes[ndir].fm[i].w << ", " << wave.planes[ndir].fm[i].x << ", " << wave.planes[ndir].fm[i].y << ", " << wave.planes[ndir].fm[i].z << std::endl;
                std::cout << "f[p]: " << f[p].v << ", " << f[p].w << ", " << f[p].x << ", " << f[p].y << ", " << f[p].z << std::endl;*/
            }
            t = epsilon.second; t_2 = t*t; t_3 = t_2*t;
            pos = (int) wave.planes[ndir].fm[1].v;
            fd0 = f[1] - f[0];
            fd1 = f[2] - f[1];
            if (pos == 0) {
                fm[0] = checkNodeSlope(fd0, fd1, 1);
                fm[1] = checkNodeSlope(fd0, fd1, 0);
            }
            else {
                fm[0] = checkNodeSlope(fd1, fd0, 0);
                fm[1] = checkNodeSlope(fd1, fd0, 1);
            }
            f0 = f[pos];
            f1 = f[1 + pos];
            ft = (2.0 * t_3 - 3.0 * t_2 + 1.0) * f0 + (t_3 - 2.0 * t_2 + t) * fm[0] +
                   (-2.0 * t_3 + 3.0 * t_2) * f1 + (t_3 - t_2) * fm[1];

            //std::cout << "ft: " << ft.v << ", " << ft.w << ", " << ft.x << ", " << ft.y << ", " << ft.z << std::endl;

            t = epsilon.second; t_2 = t*t; t_3 = t_2*t;
            pos = (int) wave.planes[ndir].fm[1].v;

            for (int p=0; p<3; p++) {
                fd0 = wave.planes[ndir].fm[2 + 3 + p] - wave.planes[ndir].fm[2 + p];
                fd1 = wave.planes[ndir].fm[2 + 6 + p] - wave.planes[ndir].fm[2 + 3 + p];
                if (pos == 0) {
                    fm[0] = checkNodeSlope(fd0, fd1, 1);
                    fm[1] = checkNodeSlope(fd0, fd1, 0);
                }
                else {
                    fm[0] = checkNodeSlope(fd1, fd0, 0);
                    fm[1] = checkNodeSlope(fd1, fd0, 1);
                }
                f0 = wave.planes[ndir].fm[2 + 3*pos + p];
                f1 = wave.planes[ndir].fm[2 + 3 + 3*pos + p];
                f[p] = (2.0 * t_3 - 3.0 * t_2 + 1.0) * f0 + (t_3 - 2.0 * t_2 + t) * fm[0] +
                       (-2.0 * t_3 + 3.0 * t_2) * f1 + (t_3 - t_2) * fm[1];
            }
            t = epsilon.first; t_2 = t*t; t_3 = t_2*t;
            pos = (int) wave.planes[ndir].fm[0].v;
            fd0 = f[1] - f[0];
            fd1 = f[2] - f[1];
            if (pos == 0) {
                fm[0] = checkNodeSlope(fd0, fd1, 1);
                fm[1] = checkNodeSlope(fd0, fd1, 0);
            }
            else {
                fm[0] = checkNodeSlope(fd1, fd0, 0);
                fm[1] = checkNodeSlope(fd1, fd0, 1);
            }
            f0 = f[pos];
            f1 = f[1 + pos];
            ft = 0.5 * (ft + (2.0 * t_3 - 3.0 * t_2 + 1.0) * f0 + (t_3 - 2.0 * t_2 + t) * fm[0] +
                   (-2.0 * t_3 + 3.0 * t_2) * f1 + (t_3 - t_2) * fm[1]);
        }
        else if (interp == I_SPLINE) {
            Node5D fm[2], f[3], f0, f1;
            double t = epsilon.first;
            int pos = (int) wave.planes[ndir].fm[0].v;
            for (int p=0; p<3; p++) {
                if (pos == 0) {
                    fm[0] = 0.0;
                    fm[1] = 6.0/4.0*(wave.planes[ndir].fm[2 + 3*p] - 2.0*wave.planes[ndir].fm[2 + 3*p + 1] + wave.planes[ndir].fm[2 + 3*p + 2]);
                }
                else {
                    fm[0] = 6.0/4.0*(wave.planes[ndir].fm[2 + 3*p] - 2.0*wave.planes[ndir].fm[2 + 3*p + 1] + wave.planes[ndir].fm[2 + 3*p + 2]);
                    fm[1] = 0.0;
                }
                f0 = wave.planes[ndir].fm[2 + 3*p + pos];
                f1 = wave.planes[ndir].fm[2 + 3*p + 1 + pos];
                f[p] = f0 + t*(-1.0 * fm[1]/6.0 - fm[0]/3.0 + f1 - f0 + t*(fm[0]/2.0 + t*(fm[1] - fm[0])/6.0));
            }
            t = epsilon.second;
            pos = (int) wave.planes[ndir].fm[1].v;
            if (pos == 0) {
                fm[0] = 0.0;
                fm[1] = 6.0*(f[0] - 2.0*f[1] + f[2])/4.0;
            }
            else {
                fm[0] = 6.0*(f[0] - 2.0*f[1] + f[2])/4.0;
                fm[1] = 0.0;
            }
            f0 = f[pos];
            f1 = f[1 + pos];
            ft = f0 + t*(-1.0 * fm[1]/6.0 - fm[0]/3.0 + f1 - f0 + t*(fm[0]/2.0 + t*(fm[1] - fm[0])/6.0));

            t = epsilon.second;
            pos = (int) wave.planes[ndir].fm[1].v;

            for (int p=0; p<3; p++) {
                if (pos == 0) {
                    fm[0] = 0.0;
                    fm[1] = 6.0*(wave.planes[ndir].fm[2 + p] - 2.0*wave.planes[ndir].fm[2 + 3 + p] + wave.planes[ndir].fm[2 + 6 + p])/4.0;
                }
                else {
                    fm[0] = 6.0*(wave.planes[ndir].fm[2 + p] - 2.0*wave.planes[ndir].fm[2 + 3 + p] + wave.planes[ndir].fm[2 + 6 + p])/4.0;
                    fm[1] = 0.0;
                }
                f0 = wave.planes[ndir].fm[2 + 3*pos + p];
                f1 = wave.planes[ndir].fm[2 + 3 + 3*pos + p];
                f[p] = f0 + t*(-1.0 * fm[1]/6.0 - fm[0]/3.0 + f1 - f0 + t*(fm[0]/2.0 + t*(fm[1] - fm[0])/6.0));
            }
            t = epsilon.first;
            pos = (int) wave.planes[ndir].fm[0].v;
            if (pos == 0) {
                fm[0] = 0.0;
                fm[1] = 6.0*(f[0] - 2.0*f[1] + f[2])/4.0;
            }
            else {
                fm[0] = 6.0*(f[0] - 2.0*f[1] + f[2])/4.0;
                fm[1] = 0.0;
            }
            f0 = f[pos];
            f1 = f[1 + pos];
            ft = f0 + t*(-1.0 * fm[1]/6.0 - fm[0]/3.0 + f1 - f0 + t*(fm[0]/2.0 + t*(fm[1] - fm[0])/6.0));
        }
        value = to_cartesian(ft);
        if (((ft.z < f0_3.z && ft.z < f1_3.z && ft.z < f2_3.z && ft.z < f3_3.z) ||
              (ft.z > f0_3.z && ft.z > f1_3.z && ft.z > f2_3.z && ft.z > f3_3.z))  && interp != I_PCHIP && interp != I_BILINEAR) {
            value = GetFValue(image, wave, ndir, I_BILINEAR, epsilon, h);
        }

        return value;
    }

    double GetUValue(Grid3 &image, AugWmm3 &wave, int ndir, int interp, std::pair<double, double> epsilon_t, Node3D &h) {

        double value = MAX_VAL;
        double gamma = wave.gamma;

        int segment_1 = (epsilon_t.first >= 1.0) ? gamma - 1 : (int) floor(gamma * epsilon_t.first);
        int segment_2 = (epsilon_t.second >= 1.0) ? gamma - 1 : (int) floor(gamma * epsilon_t.second);

        double step = 1.0/ gamma;
        std::pair<double, double> epsilon = std::pair<double, double>(
            (epsilon_t.first - segment_1 * step) / step,
            (epsilon_t.second - segment_2 * step) / step
        );

        if (interp == I_BILINEAR) {
            double v0 = wave.planes[ndir].v[segment_1][segment_2];
            double v1 = wave.planes[ndir].v[segment_1 + 1][segment_2];
            double v2 = wave.planes[ndir].v[segment_1][segment_2 + 1];
            double v3 = wave.planes[ndir].v[segment_1 + 1][segment_2 + 1];

            value = (1.0 - epsilon.first)*(1.0 - epsilon.second)*v0 +
                 epsilon.first*(1.0 - epsilon.second)*v1 +
                 (1.0 - epsilon.first)*epsilon.second*v2 +
                 epsilon.first*epsilon.second*v3;
            /*if (wave.p.x == 0 && wave.p.y == 0 && wave.p.z == 1) {
                std::cout << "values" << std::endl;
                for (int i=0; i <= (int) gamma; i++) {
                    for (int j=0; j <= (int) gamma; j++) {
                        std::cout << wave.planes[ndir].v[i][j] << ", ";
                    }
                    std::cout << std::endl;
                }
                std::cout << "vs: " << v0 << ", " << v1 << ", " << v2 << ", " << v3 << " - " << value << std::endl;
            }*/
            //std::cout << "vs: " << epsilon.first << ", " << epsilon.second << ", " << v0 << ", " << v1 << ", " << v2 << ", " << v3 << " - " << value << std::endl;
        }
        else if (interp == I_QUADATRIC) {
            double m[4];
            m[0] = wave.planes[ndir].v[segment_1][segment_2];
            m[1] = wave.planes[ndir].v[segment_1 + 1][segment_2] - m[0];
            m[2] = wave.planes[ndir].v[segment_1][segment_2 + 1] - m[0];
            m[3] = wave.planes[ndir].v[segment_1 + 1][segment_2 + 1] - m[2] - m[1] - m[0];
            value = m[0] + epsilon.first*m[1] + epsilon.second*m[2] + epsilon.first*epsilon.second*m[3];
        }
        else if (interp == I_CUBIC) {
            double m[9];
            int dx = (segment_1 > 0) ? -1 : 2;
            int dy = (segment_2 > 0) ? -1 : 2;
            double p1 = wave.planes[ndir].v[segment_1 + 1][segment_2];
            double p2 = wave.planes[ndir].v[segment_1 + dx][segment_2];
            double p3 = wave.planes[ndir].v[segment_1][segment_2 + 1];
            double p4 = wave.planes[ndir].v[segment_1][segment_2 + dy];
            m[0] = wave.planes[ndir].v[segment_1][segment_2];
            m[4] = (p2 - (1.0 - dx) * m[0] - dx * p1) / (dx * (dx - 1.0));
            m[1] = p1 - m[4] - m[0];
            m[5] = (p4 - (1.0 - dy) * m[0] - dy * p3) / (dy * (dy - 1.0));
            m[2] = p3 - m[5] - m[0];
            double a_1 = wave.planes[ndir].v[segment_1 + 1][segment_2 + 1] -
                    m[0] - m[1] - m[2] - m[4] - m[5];
            double a_3 = wave.planes[ndir].v[segment_1 + 1][segment_2 + dy] - m[0] - m[1] - m[4] - dy * m[2] - dy * dy * m[5];
            double a_2 = wave.planes[ndir].v[segment_1 + dx][segment_2 + 1] - m[0] - dx * m[1] - dx * dx * m[4] - m[2] - m[5];
            double a_4 = wave.planes[ndir].v[segment_1 + dx][segment_2 + dy] - m[0] - dx * m[1] - dx * dx * m[4] - dy * m[2] - dy * dy * m[5];
            m[8] = (dx * dy * a_1 - dy * a_2 - dx * a_3 + a_4) /
                   (dx * dy - dx * dy * dy - dx * dx * dy + dx * dx * dy * dy);
            m[7] = (dy * a_2 - a_4 - (dx * dx * dy - dx * dx * dy * dy) * m[8]) / (dx * dy - dx * dy * dy);
            m[6] = (dx * a_3 - a_4 - (dx * dy * dy - dx * dx * dy * dy) * m[8]) / (dx * dy - dx * dx * dy);
            m[3] = (a_4 - dx * dx * dy * dy * m[8] - dx * dy * dy * m[7] - dx * dx * dy * m[6]) / (dx * dy);
            value = m[0] + epsilon.first*m[1] + epsilon.second*m[2] + epsilon.first*epsilon.second*m[3] +
                    epsilon.first*epsilon.first*m[4] + epsilon.second*epsilon.second*m[5] +
                    epsilon.first*epsilon.first*epsilon.second*m[6] + epsilon.first*epsilon.second*epsilon.second*m[7] +
                    epsilon.first*epsilon.first*epsilon.second*epsilon.second*m[8];
        }
        else if (interp == I_HERMITE) {
            double m[2], v[4], v0, v1;
            double t = epsilon.first, t_2 = t*t, t_3 = t_2*t;
            int init = (segment_2 == 0) ? 0 : -1;
            int end = (segment_2 == (gamma - 1)) ? 1 : 2;
            for (int p=init; p<=end; p++) {
                m[0] = (segment_1 == 0) ? wave.planes[ndir].v[segment_1 + 1][segment_2 + p] - wave.planes[ndir].v[segment_1][segment_2 + p] :
                                          (wave.planes[ndir].v[segment_1 + 1][segment_2 + p] - wave.planes[ndir].v[segment_1 - 1][segment_2 + p]) / 2.;
                m[1] = (segment_1 == (gamma - 1)) ? wave.planes[ndir].v[segment_1 + 1][segment_2 + p] - wave.planes[ndir].v[segment_1][segment_2 + p] :
                                                    (wave.planes[ndir].v[segment_1 + 2][segment_2 + p] - wave.planes[ndir].v[segment_1][segment_2 + p]) / 2.;
                v0 = wave.planes[ndir].v[segment_1][segment_2 + p];
                v1 = wave.planes[ndir].v[segment_1 + 1][segment_2 + p];
                v[p - init] = (2.0 * t_3 - 3.0 * t_2 + 1.0) * v0 + (t_3 - 2.0 * t_2 + t) * m[0] +
                       (-2.0 * t_3 + 3.0 * t_2) * v1 + (t_3 - t_2) * m[1];
            }
            t = epsilon.second; t_2 = t*t; t_3 = t_2*t;
            m[0] = (init == 0) ? v[1] - v[0] : (v[2] - v[0]) / 2.0;
            m[1] = (end == 1) ? v[2] - v[1] : (v[end - init] - v[end - init - 2]) / 2.0;
            v0 = v[-init];
            v1 = v[1 - init];
            value = (2.0 * t_3 - 3.0 * t_2 + 1.0) * v0 + (t_3 - 2.0 * t_2 + t) * m[0] +
                    (-2.0 * t_3 + 3.0 * t_2) * v1 + (t_3 - t_2) * m[1];

            t = epsilon.second; t_2 = t*t; t_3 = t_2*t;

            init = (segment_1 == 0) ? 0 : -1;
            end = (segment_1 == (gamma - 1)) ? 1 : 2;
            for (int p=init; p<=end; p++) {
                m[0] = (segment_2 == 0) ? wave.planes[ndir].v[segment_1 + p][segment_2 + 1] - wave.planes[ndir].v[segment_1 + p][segment_2] :
                                          (wave.planes[ndir].v[segment_1 + p][segment_2 + 1] - wave.planes[ndir].v[segment_1 + p][segment_2 - 1]) / 2.;
                m[1] = (segment_2 == (gamma - 1)) ? wave.planes[ndir].v[segment_1 + p][segment_2 + 1] - wave.planes[ndir].v[segment_1 + p][segment_2] :
                                                    (wave.planes[ndir].v[segment_1 + p][segment_2 + 2] - wave.planes[ndir].v[segment_1 + p][segment_2]) / 2.;
                v0 = wave.planes[ndir].v[segment_1 + p][segment_2];
                v1 = wave.planes[ndir].v[segment_1 + p][segment_2 + 1];
                v[p - init] = (2.0 * t_3 - 3.0 * t_2 + 1.0) * v0 + (t_3 - 2.0 * t_2 + t) * m[0] +
                       (-2.0 * t_3 + 3.0 * t_2) * v1 + (t_3 - t_2) * m[1];
            }
            t = epsilon.first; t_2 = t*t; t_3 = t_2*t;
            m[0] = (init == 0) ? v[1] - v[0] : (v[2] - v[0]) / 2.0;
            m[1] = (end == 1) ? v[2] - v[1] : (v[end - init] - v[end - init - 2]) / 2.0;
            v0 = v[-init];
            v1 = v[1 - init];
            value = 0.5 * (value + (2.0 * t_3 - 3.0 * t_2 + 1.0) * v0 + (t_3 - 2.0 * t_2 + t) * m[0] +
                    (-2.0 * t_3 + 3.0 * t_2) * v1 + (t_3 - t_2) * m[1]);
        }
        else if (interp == I_PCHIP) {
            double m[2], v[4], v0, v1, dv0, dv1;
            double t = epsilon.first, t_2 = t*t, t_3 = t_2*t;
            int init = (segment_2 == 0) ? 0 : -1;
            int end = (segment_2 == (gamma - 1)) ? 1 : 2;
            for (int p=init; p<=end; p++) {
                if (segment_1 == 0) {
                    dv0 = wave.planes[ndir].v[segment_1 + 1][segment_2 + p] - wave.planes[ndir].v[segment_1][segment_2 + p];
                    dv1 = wave.planes[ndir].v[segment_1 + 2][segment_2 + p] - wave.planes[ndir].v[segment_1 + 1][segment_2 + p];
                }
                else {
                    dv0 = wave.planes[ndir].v[segment_1][segment_2 + p] - wave.planes[ndir].v[segment_1 - 1][segment_2 + p];
                    dv1 = wave.planes[ndir].v[segment_1 + 1][segment_2 + p] - wave.planes[ndir].v[segment_1][segment_2 + p];
                }
                m[0] = checkSlope(dv0, dv1, (segment_1 == 0));
                if (segment_1 == (gamma - 1)) {
                    dv0 = wave.planes[ndir].v[segment_1][segment_2 + p] - wave.planes[ndir].v[segment_1 - 1][segment_2 + p];
                    dv1 = wave.planes[ndir].v[segment_1 + 1][segment_2 + p] - wave.planes[ndir].v[segment_1][segment_2 + p];
                }
                else {
                    dv0 = wave.planes[ndir].v[segment_1 + 1][segment_2 + p] - wave.planes[ndir].v[segment_1][segment_2 + p];
                    dv1 = wave.planes[ndir].v[segment_1 + 2][segment_2 + p] - wave.planes[ndir].v[segment_1 + 1][segment_2 + p];
                }
                m[1] = checkSlope(dv1, dv0, (segment_1 == (gamma - 1)));
                v0 = wave.planes[ndir].v[segment_1][segment_2 + p];
                v1 = wave.planes[ndir].v[segment_1 + 1][segment_2 + p];
                v[p - init] = (2.0 * t_3 - 3.0 * t_2 + 1.0) * v0 + (t_3 - 2.0 * t_2 + t) * m[0] +
                       (-2.0 * t_3 + 3.0 * t_2) * v1 + (t_3 - t_2) * m[1];
            }
            t = epsilon.second; t_2 = t*t; t_3 = t_2*t;
            dv0 = v[1] - v[0];
            dv1 = v[2] - v[1];
            m[0] = checkSlope(dv0, dv1, (init == 0));
            dv0 = v[end - init - 1] - v[end - init - 2];
            dv1 = v[end - init] - v[end - init - 1];
            m[1] = checkSlope(dv1, dv0, (end == 1));
            v0 = v[-init];
            v1 = v[1 - init];
            value = (2.0 * t_3 - 3.0 * t_2 + 1.0) * v0 + (t_3 - 2.0 * t_2 + t) * m[0] +
                    (-2.0 * t_3 + 3.0 * t_2) * v1 + (t_3 - t_2) * m[1];

            t = epsilon.second; t_2 = t*t; t_3 = t_2*t;

            init = (segment_1 == 0) ? 0 : -1;
            end = (segment_1 == (gamma - 1)) ? 1 : 2;
            for (int p=init; p<=end; p++) {
                if (segment_2 == 0) {
                    dv0 = wave.planes[ndir].v[segment_1 + p][segment_2 + 1] - wave.planes[ndir].v[segment_1 + p][segment_2];
                    dv1 = wave.planes[ndir].v[segment_1 + p][segment_2 + 2] - wave.planes[ndir].v[segment_1 + p][segment_2 + 1];
                }
                else {
                    dv0 = wave.planes[ndir].v[segment_1 + p][segment_2] - wave.planes[ndir].v[segment_1 + p][segment_2 - 1];
                    dv1 = wave.planes[ndir].v[segment_1 + p][segment_2 + 1] - wave.planes[ndir].v[segment_1 + p][segment_2];
                }
                m[0] = checkSlope(dv0, dv1, (segment_2 == 0));
                if (segment_2 == (gamma - 1)) {
                    dv0 = wave.planes[ndir].v[segment_1 + p][segment_2] - wave.planes[ndir].v[segment_1 + p][segment_2 - 1];
                    dv1 = wave.planes[ndir].v[segment_1 + p][segment_2 + 1] - wave.planes[ndir].v[segment_1 + p][segment_2];
                }
                else {
                    dv0 = wave.planes[ndir].v[segment_1 + p][segment_2 + 1] - wave.planes[ndir].v[segment_1 + p][segment_2];
                    dv1 = wave.planes[ndir].v[segment_1 + p][segment_2 + 2] - wave.planes[ndir].v[segment_1 + p][segment_2 + 1];
                }
                m[1] = checkSlope(dv1, dv0, (segment_2 == (gamma - 1)));
                v0 = wave.planes[ndir].v[segment_1 + p][segment_2];
                v1 = wave.planes[ndir].v[segment_1 + p][segment_2 + 1];
                v[p - init] = (2.0 * t_3 - 3.0 * t_2 + 1.0) * v0 + (t_3 - 2.0 * t_2 + t) * m[0] +
                       (-2.0 * t_3 + 3.0 * t_2) * v1 + (t_3 - t_2) * m[1];
            }
            t = epsilon.first; t_2 = t*t; t_3 = t_2*t;
            dv0 = v[1] - v[0];
            dv1 = v[2] - v[1];
            m[0] = checkSlope(dv0, dv1, (init == 0));
            dv0 = v[end - init - 1] - v[end - init - 2];
            dv1 = v[end - init] - v[end - init - 1];
            m[1] = checkSlope(dv1, dv0, (end == 1));
            v0 = v[-init];
            v1 = v[1 - init];
            value = 0.5 * (value + (2.0 * t_3 - 3.0 * t_2 + 1.0) * v0 + (t_3 - 2.0 * t_2 + t) * m[0] +
                    (-2.0 * t_3 + 3.0 * t_2) * v1 + (t_3 - t_2) * m[1]);
        }
        else if (interp == I_SPLINE) {

            double m[8], v[4], v0, v1, y[4];
            double t = epsilon.first;
            int init = (segment_2 == 0) ? 0 : -1;
            int end = (segment_2 == (gamma - 1)) ? 1 : 2;
            int init_2 = (segment_1 == 0) ? 0 : -1;
            int end_2 = (segment_1 == (gamma - 1)) ? 1 : 2;

            for (int p=init; p<=end; p++) {
                for (int p2=init_2; p2<=end_2; p2++) {
                    y[p2 - init_2] = wave.planes[ndir].v[segment_1 + p2][segment_2 + p];
                    /*if (wave.p.x == 0 && wave.p.y == 1 && wave.p.z == 1) {
                        std::cout << "y: " << y[p2 - init_2] << std::endl;
                    }*/
                }
                get_spline_coeffs(y, m, end_2 - init_2);
                v0 = wave.planes[ndir].v[segment_1][segment_2 + p];
                v1 = wave.planes[ndir].v[segment_1 + 1][segment_2 + p];
                v[p - init] = (1.0 - t) * v0 + t * v1 + t * (1.0 - t) * ((1.0 - t) * m[-2*init_2] + t * m[-2*init_2 + 1]);

            }
            t = epsilon.second;
            get_spline_coeffs(v, m, end - init);
            v0 = v[-init];
            v1 = v[1 - init];
            value = (1.0 - t) * v0 + t * v1 + t * (1.0 - t) * ((1.0 - t) * m[-2*init] + t * m[-2*init + 1]);

            t = epsilon.second;
            init = (segment_1 == 0) ? 0 : -1;
            end = (segment_1 == (gamma - 1)) ? 1 : 2;
            init_2 = (segment_2 == 0) ? 0 : -1;
            end_2 = (segment_2 == (gamma - 1)) ? 1 : 2;
            for (int p=init; p<=end; p++) {
                for (int p2=init_2; p2<=end_2; p2++) {
                    y[p2 - init_2] = wave.planes[ndir].v[segment_1 + p][segment_2 + p2];
                }
                get_spline_coeffs(y, m, end_2 - init_2);
                v0 = wave.planes[ndir].v[segment_1 + p][segment_2];
                v1 = wave.planes[ndir].v[segment_1 + p][segment_2 + 1];
                v[p - init] = (1.0 - t) * v0 + t * v1 + t * (1.0 - t) * ((1.0 - t) * m[-2*init_2] + t * m[-2*init_2 + 1]);
            }

            t = epsilon.first;
            get_spline_coeffs(v, m, end - init);
            v0 = v[-init];
            v1 = v[1 - init];
            value = 0.5 * (value + (1.0 - t) * v0 + t * v1 + t * (1.0 - t) * ((1.0 - t) * m[-2*init] + t * m[-2*init + 1]));

        }
        return value;
    }


    double GetInterpValue(Grid3 &image, AugWmm3 &wave, int ndir,
                          int interp, std::pair<double, double> epsilon, Node3D &neighD, Node3D &fn,
                          Node3D &h, int inner_nodes=1) {

        Node3D plane_posD((wave.p.x + epsilon.first*wave.planes[ndir].dirs[0].x + epsilon.second*wave.planes[ndir].dirs[1].x),
                                 (wave.p.y + epsilon.first*wave.planes[ndir].dirs[0].y + epsilon.second*wave.planes[ndir].dirs[1].y),
                                 (wave.p.z + epsilon.first*wave.planes[ndir].dirs[0].z + epsilon.second*wave.planes[ndir].dirs[1].z));

        double value = MAX_VAL;
        // std::cout << "iniciando F " << interp << std::endl;
        //Node3D ft = GetFValue(image, wave, ndir, interp, epsilon, h);
        // std::cout << "iniciando U: " << ft.x << ", " << ft.y << ", " << ft.z << std::endl;
        double vt = GetUValue(image, wave, ndir, interp, epsilon, h);
        // std::cout << "finish him!" << std::endl;

        if (((vt < wave.v && vt < wave.planes[ndir].v[0][0] && vt < wave.planes[ndir].v[wave.gamma][0] && vt < wave.planes[ndir].v[0][wave.gamma]) ||
                (vt > wave.v && vt > wave.planes[ndir].v[0][0] && vt > wave.planes[ndir].v[wave.gamma][0] && vt > wave.planes[ndir].v[0][wave.gamma])) &&
                interp != I_PCHIP && interp != I_BILINEAR) {
            vt = GetUValue(image, wave, ndir, I_PCHIP, epsilon, h);
        }
        //double fti = GetIntegralValue(image, plane_posD, neighD, norm(ft), norm(fn), interp, inner_nodes);
        Node3D diff = h * (neighD - plane_posD);
        Node3D a = diff / norm(diff);
        value = vt + norm(diff) / sqrt(1.0 + (fn.x*a.x + fn.y*a.y + fn.z*a.z)*(fn.x*a.x + fn.y*a.y + fn.z*a.z));
        // value = vt + norm(h * (neighD - plane_posD)) * fti;
            /*if (wave.p.x == 0 && wave.p.y == 0 && wave.p.z == 1) {
                std::cout << "value: " << vt << ", " << norm(neighD - plane_posD) << ", " << value << std::endl;
            }*/

        return value;
    }


    std::pair<double, double> GetEpsilonGradient(AugWmm3 &wave, int ndir, Node3D &neigh, Node3D &fn, Node3D &h) {

        Node3D ph(h.x*wave.p.x, h.y*wave.p.y, h.z*wave.p.z), nh(h.x*neigh.x, h.y*neigh.y, h.z*neigh.z);
        Node3D d0h(h.x*wave.planes[ndir].dirs[0].x, h.y*wave.planes[ndir].dirs[0].y, h.z*wave.planes[ndir].dirs[0].z);
        Node3D d1h(h.x*wave.planes[ndir].dirs[1].x, h.y*wave.planes[ndir].dirs[1].y, h.z*wave.planes[ndir].dirs[1].z);

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
            x = wmm3D::Node3D(nh.x - t*A, nh.y - t*B, nh.z - t*C);
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


    /*std::pair<double, double> GetEpsilonHopfLax(AugWmm3 &wave, int ndir, Node3D &neigh, Node3D &h) {

        Node3D ph(h.x*wave.p.x, h.y*wave.p.y, h.z*wave.p.z), nh(h.x*neigh.x, h.y*neigh.y, h.z*neigh.z);
        Node3D d0h(h.x*wave.planes[ndir].dirs[0].x, h.y*wave.planes[ndir].dirs[0].y, h.z*wave.planes[ndir].dirs[0].z);
        Node3D d1h(h.x*wave.planes[ndir].dirs[1].x, h.y*wave.planes[ndir].dirs[1].y, h.z*wave.planes[ndir].dirs[1].z);

        Node3D inc0, inc1, wave_p(wave.p.x, wave.p.y, wave.p.z);
        Node3 d01h = wave.planes[ndir].dirs[0] + wave.planes[ndir].dirs[1];
        Node3D ori = neigh - wave_p;
        int half = std::max(0, wave.gamma / 2 - 1);
        double inc_v0 = (wave.planes[ndir].v[half][half] - wave.planes[ndir].v[half + 2][half]) * wave.gamma / 2.0;
        double inc_v1 = (wave.planes[ndir].v[half][half] - wave.planes[ndir].v[half][half + 2]) * wave.gamma / 2.0;
        if (d01h.x == 0) {
            inc0 = (wave.planes[ndir].dirs[0].y == 0) ? Node3D(ori.x*inc_v0, 0., d0h.z) : Node3D(ori.x*inc_v0, d0h.y, 0.);
            inc1 = (wave.planes[ndir].dirs[1].y == 0) ? Node3D(ori.x*inc_v1, 0., d1h.z) : Node3D(ori.x*inc_v1, d1h.y, 0.);
        }
        else if (d01h.y == 0) {
            inc0 = (wave.planes[ndir].dirs[0].x == 0) ? Node3D(0., ori.y*inc_v0, d0h.z) : Node3D(d0h.x, ori.y*inc_v0, 0.);
            inc1 = (wave.planes[ndir].dirs[1].x == 0) ? Node3D(0., ori.y*inc_v1, d1h.z) : Node3D(d1h.x, ori.y*inc_v1, 0.);
        }
        else {
            inc0 = (wave.planes[ndir].dirs[0].x == 0) ? Node3D(0., d0h.y, ori.z*inc_v0) : Node3D(d0h.x, 0., ori.z*inc_v0);
            inc1 = (wave.planes[ndir].dirs[1].x == 0) ? Node3D(0., d1h.y, ori.z*inc_v1) : Node3D(d1h.x, 0., ori.z*inc_v1);
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
            diff = x - ph;
        else if (fabs(den) == 0.0 && norm(d0h + d1h) > 0.0) {
            t = (A*nh.x + B*nh.y + C*nh.z + D)/(A*A + B*B + C*C);
            x = Node3D(nh.x - t*A, nh.y - t*B, nh.z - t*C);
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

    }*/

    std::pair<double, double> GetEpsilonHopfLax(Node3D &wave_p, Node3D dirs[2], double v[4], Node3D &neigh, Node3D &h, double step) {

        Node3D ph(h.x*wave_p.x, h.y*wave_p.y, h.z*wave_p.z), nh(h.x*neigh.x, h.y*neigh.y, h.z*neigh.z);
        Node3D d0h = h * dirs[0];
        Node3D d1h = h * dirs[1];

        Node3D inc0, inc1;
        Node3D d01h = dirs[0] + dirs[1];
        Node3D ori = neigh - wave_p;
        double cos_delta_0 = (v[1] - v[0]) / norm(step * d0h);
        double cos_delta_1 = (v[2] - v[0]) / norm(step * d1h);
        if (fabs(cos_delta_0) > 1) cos_delta_0 = sgn(cos_delta_0);
        if (fabs(cos_delta_1) > 1) cos_delta_1 = sgn(cos_delta_1);
        double sin_delta_0 = sqrt(1. - cos_delta_0*cos_delta_0);
        double sin_delta_1 = sqrt(1. - cos_delta_1*cos_delta_1);
        if (d01h.x == 0) {
            // inc0 = (wave.planes[ndir].dirs[0].y == 0) ? Node3D(sgn(ori.x)*sin_delta_0, 0., sgn(ori.z)*cos_delta_0) : Node3D(sgn(ori.x)*sin_delta_0, sgn(ori.y)*cos_delta_0, 0.);
            inc0 = (dirs[0].y == 0) ? Node3D(-sgn(dirs[0].z)*cos_delta_0, 0., sgn(ori.x)*sin_delta_0) :
                                      Node3D(-sgn(dirs[0].y)*cos_delta_0, sgn(ori.x)*sin_delta_0, 0.);
            inc1 = (dirs[1].y == 0) ? Node3D(-sgn(dirs[1].z)*cos_delta_1, 0., sgn(ori.x)*sin_delta_1) :
                                      Node3D(-sgn(dirs[1].y)*cos_delta_1, sgn(ori.x)*sin_delta_1, 0.);
        }
        else if (d01h.y == 0) {
            inc0 = (dirs[0].x == 0) ? Node3D(0., -sgn(dirs[0].z)*cos_delta_0, sgn(ori.y)*sin_delta_0) :
                                      Node3D(sgn(ori.y)*sin_delta_0, -sgn(dirs[0].x)*cos_delta_0, 0.);
            inc1 = (dirs[1].x == 0) ? Node3D(0., -sgn(dirs[1].z)*cos_delta_1, sgn(ori.y)*sin_delta_1) :
                                      Node3D(sgn(ori.y)*sin_delta_1, -sgn(dirs[1].x)*cos_delta_1, 0.);
        }
        else {
            inc0 = (dirs[0].x == 0) ? Node3D(0., sgn(ori.z)*sin_delta_0, -sgn(dirs[0].y)*cos_delta_0) :
                                      Node3D(sgn(ori.z)*sin_delta_0,  0., -sgn(dirs[0].x)*cos_delta_0);
            inc1 = (dirs[1].x == 0) ? Node3D(0., sgn(ori.z)*sin_delta_1, -sgn(dirs[1].y)*cos_delta_1) :
                                      Node3D(sgn(ori.z)*sin_delta_1,  0., -sgn(dirs[1].x)*cos_delta_1);
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
            diff = (x - ph) / (step * h);
        else {
            t = (A*nh.x + B*nh.y + C*nh.z + D)/(A*A + B*B + C*C);
            x = Node3D(nh.x - t*A, nh.y - t*B, nh.z - t*C);
            diff = (x - ph) / (step * h);
        }

        epsilon.first = diff.x*dirs[0].x + diff.y*dirs[0].y + diff.z*dirs[0].z;
        if (epsilon.first < 0.0)
            epsilon.first = 0.0;
        else if (epsilon.first > 1.0)
            epsilon.first = 1.0;

        epsilon.second = diff.x*dirs[1].x + diff.y*dirs[1].y + diff.z*dirs[1].z;
        if (epsilon.second < 0.0)
            epsilon.second = 0.0;
        else if (epsilon.second > 1.0)
            epsilon.second = 1.0;

        return epsilon;

    }

    std::pair<double, double> GetEpsilonGoldenSearch(Grid3 &image, AugWmm3 &wave, int ndir, int interp, Node3D &neigh, Node3D &fn,
                                                     Node3D &h, int gamma) {

        double a_1 = 0.0, b_1 = 1.0, x1_1 = a_1 + (1.0-RESPHI)*(b_1 - a_1), x2_1 = a_1 + RESPHI*(b_1 - a_1),
                f_x1 = MAX_VAL, f_x2 = MAX_VAL;

        double a_2 = 0.0, b_2 = 1.0, x1_2 = a_2 + (1.0-RESPHI)*(b_2 - a_2), x2_2 = a_2 + RESPHI*(b_2 - a_2);

        std::pair<double, double> epsilon(0.0, 0.0);

        double res = GetInterpValue(image, wave, ndir, interp, std::pair<double, double>(0.0, 0.0), neigh, fn, h, gamma - 1);
        double f_b = GetInterpValue(image, wave, ndir, interp, std::pair<double, double>(1.0, 0.0), neigh, fn, h, gamma - 1);
        if (f_b < res) {
            res = f_b; epsilon = std::pair<double, double>(1.0, 0.0);
        }
        f_b = GetInterpValue(image, wave, ndir, interp, std::pair<double, double>(0.0, 1.0), neigh, fn, h, gamma - 1);
        if (f_b < res) {
            res = f_b; epsilon = std::pair<double, double>(0.0, 1.0);
        }
        f_b = GetInterpValue(image, wave, ndir, interp, std::pair<double, double>(1.0, 1.0), neigh, fn, h, gamma - 1);
        if (f_b < res) {
            res = f_b; epsilon = std::pair<double, double>(1.0, 1.0);
        }

        f_x1 = GetInterpValue(image, wave, ndir, interp, std::pair<double, double>(0.0, x1_2), neigh, fn, h, gamma - 1);
        f_x2 = GetInterpValue(image, wave, ndir, interp, std::pair<double, double>(0.0, x2_2), neigh, fn, h, gamma - 1);

        while (fabs(b_1 - a_1) > TAU && fabs(b_2 - a_2) > TAU) {
            if(f_x1 < f_x2) {
                f_x1 = GetInterpValue(image, wave, ndir, interp, std::pair<double, double>(x1_1, x1_2), neigh, fn, h, gamma - 1);
                f_x2 = GetInterpValue(image, wave, ndir, interp, std::pair<double, double>(x2_1, x1_2), neigh, fn, h, gamma - 1);
                b_2 = x2_2; x2_2 = x1_2; x1_2 = a_2 + (1.0 - RESPHI)*(b_2 - a_2);
            }
            else {
                f_x1 = GetInterpValue(image, wave, ndir, interp, std::pair<double, double>(x1_1, x2_2), neigh, fn, h, gamma - 1);
                f_x2 = GetInterpValue(image, wave, ndir, interp, std::pair<double, double>(x2_1, x2_2), neigh, fn, h, gamma - 1);
                a_2 = x1_2; x1_2 = x2_2; x2_2 = a_2 + RESPHI*(b_2 - a_2);
            }
            if(f_x1 < f_x2) {
                f_x1 = GetInterpValue(image, wave, ndir, interp, std::pair<double, double>(x1_1, x1_2), neigh, fn, h, gamma - 1);
                f_x2 = GetInterpValue(image, wave, ndir, interp, std::pair<double, double>(x1_1, x2_2), neigh, fn, h, gamma - 1);
                b_1 = x2_1; x2_1 = x1_1; x1_1 = a_1 + (1.0 - RESPHI)*(b_1 - a_1);
            }
            else {
                f_x1 = GetInterpValue(image, wave, ndir, interp, std::pair<double, double>(x2_1, x1_2), neigh, fn, h, gamma - 1);
                f_x2 = GetInterpValue(image, wave, ndir, interp, std::pair<double, double>(x2_1, x2_2), neigh, fn, h, gamma - 1);
                a_1 = x1_1; x1_1 = x2_1; x2_1 = a_1 + RESPHI*(b_1 - a_1);
            }
        }

        if (f_x1 < res)
            epsilon = std::pair<double, double>(x1_1, x1_2);
        if (f_x2 < std::min(res, f_x1))
            epsilon = std::pair<double, double>(x2_1, x2_2);

        return epsilon;
        
    }

    double GetVal3D(Grid3 &image, Grid3 &u_surface, AugWmm3 &wave, Node3D &neigh, Node3D &fn, Node3D &h, int interp, int mode) {

        Node3D f0(image.at(wave.p, 0), image.at(wave.p, 1), image.at(wave.p, 2));
        double y0 = wave.v;

        if (isinf(norm(f0)) || isnan(norm(f0)))
            f0 = fn;

        Node3D wave_p(wave.p.x, wave.p.y, wave.p.z);
        Node3D diff = h * (neigh - wave_p);
        Node3D a = diff / norm(diff);
        double val = y0 + norm(diff) / sqrt(1.0 + (fn.x*a.x + fn.y*a.y + fn.z*a.z)*(fn.x*a.x + fn.y*a.y + fn.z*a.z));

        if (wave.ndirs > 0) {

            for (int ndir = 0; ndir < wave.ndirs; ndir++) {

                Node3 p0 = wave.p + wave.planes[ndir].dirs[0];
                Node3 p1 = wave.p + wave.planes[ndir].dirs[1];

                if (!image.contains(p0) || !image.contains(p1))
                    continue;

                std::pair<double, double> epsilon, epsilon_aux;
                if (mode == M_GRADIENT)
                    epsilon = GetEpsilonGradient(wave, ndir, neigh, fn, h);
                else if (mode == M_HOPFLAX) {
                    //epsilon = GetEpsilonHopfLax(wave, ndir, neigh, h);
                    double res, step = 1.0 / (double) wave.gamma, res_aux, v[4];
                    Node3D dd_aux[2], dp_aux;
                    for (int i=0; i<2; i++) {
                        dd_aux[i] = Node3D(wave.planes[ndir].dirs[i].x,
                                           wave.planes[ndir].dirs[i].y,
                                           wave.planes[ndir].dirs[i].z);
                    }
                    res = MAX_VAL;
                    epsilon = std::pair<double, double>(0.0, 0.0);
                    for (int i=0; i<wave.gamma; i++) {
                        for (int j=0; j<wave.gamma; j++) {
                            v[0] = wave.planes[ndir].v[i][j];
                            v[1] = wave.planes[ndir].v[i+1][j];
                            v[2] = wave.planes[ndir].v[i][j+1];
                            v[3] = wave.planes[ndir].v[i+1][j+1];
                            dp_aux = wave_p + i*dd_aux[0] + j*dd_aux[1];
                            epsilon_aux = GetEpsilonHopfLax(dp_aux, dd_aux, v, neigh, h, step);
                            epsilon_aux.first = (i + epsilon_aux.first)*step;
                            epsilon_aux.second = (j + epsilon_aux.second)*step;
                            res_aux = GetInterpValue(image, wave, ndir, interp, epsilon_aux, neigh, fn, h, wave.gamma - 1);
                            if (res_aux < res) {
                                res = res_aux;
                                epsilon = epsilon_aux;
                            }
                        }
                    }
                }
                else
                    epsilon = GetEpsilonGoldenSearch(image, wave, ndir, interp, neigh, fn, h, wave.gamma);
                // std::cout << ndir << " : " << epsilon.first << ", " << epsilon.second << " - "<< GetInterpValue(image, wave, ndir, interp, epsilon, neigh, fn, h) << std::endl;
                val = std::min(val, GetInterpValue(image, wave, ndir, interp, epsilon, neigh, fn, h, wave.gamma - 1));
            }
        }
        return val;

    }

    void setCoeffs3D(Node5D y[3][3][3], Node5D *m, Node3 d, Node3 *dirs, int interp) {

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
            Node5D p1 = y[pos.x + dirs[0].x][pos.y + dirs[0].y][pos.z + dirs[0].z];
            Node5D p2 = y[pos.x + dx * dirs[0].x][pos.y + dx * dirs[0].y][pos.z + dx * dirs[0].z];
            Node5D p3 = y[pos.x + dirs[1].x][pos.y + dirs[1].y][pos.z + dirs[1].z];
            Node5D p4 = y[pos.x + dy * dirs[1].x][pos.y + dy * dirs[1].y][pos.z + dy * dirs[1].z];
            m[0] = y[pos.x][pos.y][pos.z];
            m[4] = (p2 - (1.0 - dx) * m[0] - dx * p1) / (dx * (dx - 1.0));
            m[1] = p1 - m[4] - m[0];
            m[5] = (p4 - (1.0 - dy) * m[0] - dy * p3) / (dy * (dy - 1.0));
            m[2] = p3 - m[5] - m[0];
            Node5D a_1 = y[pos.x + dirs[0].x + dirs[1].x][pos.y + dirs[0].y + dirs[1].y][pos.z + dirs[0].z + dirs[1].z] -
                    m[0] - m[1] - m[2] - m[4] - m[5];
            Node5D a_3 = y[pos.x + dirs[0].x + dy * dirs[1].x][pos.y + dirs[0].y +
                    dy * dirs[1].y][pos.z + dirs[0].z + dy * dirs[1].z] - m[0] - m[1] - m[4] - dy * m[2] - dy * dy * m[5];
            Node5D a_2 = y[pos.x + dx * dirs[0].x + dirs[1].x][pos.y + dx * dirs[0].y +
                    dirs[1].y][pos.z + dx * dirs[0].z + dirs[1].z] - m[0] - dx * m[1] - dx * dx * m[4] - m[2] - m[5];
            Node5D a_4 = y[pos.x + dx * dirs[0].x + dy * dirs[1].x][pos.y + dx * dirs[0].y +
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
            m[0] = Node5D(0, 0, 0, 0, 0);
            if ((dirs[0].x != 0 && pos.x == 1) || (dirs[0].y != 0 && pos.y == 1) || (dirs[0].z != 0 && pos.z == 1)) {
                m[0] = Node5D(1, 1, 1, 1, 1);
            }
            m[1] = Node5D(0, 0, 0, 0, 0);
            if ((dirs[1].x != 0 && pos.x == 1) || (dirs[1].y != 0 && pos.y == 1) || (dirs[1].z != 0 && pos.z == 1)) {
                m[1] = Node5D(1, 1, 1, 1, 1);
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


    void augAniWmm3DIsoSurface3D(Grid3 &image, std::vector<Node3> &initials, Node3D &h,
                                      int interp, int mode, int N, int gamma, Grid3 &u_surface) {

        bool isnewpos[3][3][3], wascomputed[3][3][3];
        double valcenter[3][3][3], value;
        Node5D imcenter[3][3][3];
        Node3 dirs[4], neigh, neigh_aux;
        int ndirs = 3;
        double step = 1.0 / (double) gamma, epsilon1, epsilon2;

        Grid3_<unsigned char> state = Grid3_<unsigned char>(image.rows, image.cols, image.channels);

        std::multimap<double, AugWmm3 > trial_set;
        std::unordered_map<int, std::multimap<double, AugWmm3 >::iterator> mapa_trial;
        std::unordered_map<int, double> augmented_values;

        std::multimap<double, AugWmm3 >::iterator trial_set_it;
        std::unordered_map<int, std::multimap<double, AugWmm3 >::iterator>::iterator mapa_trial_it;
        std::pair<double, AugWmm3 > pr_trial;
        std::pair<int, std::multimap<double, AugWmm3 >::iterator> pr_mapa;

        int key, i, aug_key;
        AugWmm3 winner;
        Node3D fi, neighD, d1, d2, neighD_aux;

        // Initialization
        for (i = 0; i < (int) initials.size(); i++) {
            key = initials[i].x + u_surface.rows*(initials[i].y + u_surface.cols*initials[i].z);
            if (mapa_trial.find(key) == mapa_trial.end() && u_surface.contains(initials[i])) {
                u_surface.at(initials[i]) = 0.0;
                winner = AugWmm3(N, gamma, 0);
                winner.v = 0.0;
                winner.p = initials[i];
                state.at(initials[i]) = P_TRIAL;
                pr_trial = std::pair<double, AugWmm3 >(0.0, winner);
                trial_set_it = trial_set.insert(pr_trial);
                pr_mapa = std::pair<int, typename std::multimap<double, AugWmm3 >::iterator>(key, trial_set_it);
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

            // std::cout << "winner: " << winner.p.x << ", " << winner.p.y << ", " << winner.p.z << " - " << winner.v << std::endl;
            if (winner.ndirs > 0) {
                ndirs = winner.ndirs;
                neigh = winner.p;
                neighD = wmm3D::Node3D(neigh.x, neigh.y, neigh.z);
                for (int m=0; m < ndirs; m++) {
                    d1 = wmm3D::Node3D(winner.planes[m].dirs[0].x, winner.planes[m].dirs[0].y, winner.planes[m].dirs[0].z);
                    d2 = wmm3D::Node3D(winner.planes[m].dirs[1].x, winner.planes[m].dirs[1].y, winner.planes[m].dirs[1].z);
                    // std::cout << "m: " << m << ", 1: " << d1.x << ", " << d1.y << ", " << d1.z << std::endl;
                    // std::cout << "m: " << m << ", 2: " << d2.x << ", " << d2.y << ", " << d2.z << std::endl;
                    for (int v=0; v <= gamma; v++) {
                        for (int w=0; w <= gamma; w++) {
                            epsilon1 = v*step;
                            epsilon2 = w*step;
                            neighD_aux = neighD + epsilon1*d1 + epsilon2*d2;
                            if (!image.contains(neighD_aux)) {
                                winner.planes[m].v[v][w] = MAX_VAL;
                            }
                            else {
                                neigh_aux = Node3(round(gamma * neighD_aux.x), round(gamma * neighD_aux.y), round(gamma * neighD_aux.z));
                                aug_key = neigh_aux.x + gamma * u_surface.rows*(neigh_aux.y + gamma * u_surface.cols*neigh_aux.z);
                                winner.planes[m].v[v][w] = augmented_values[aug_key];
                                // std::cout << aug_key << " : " << m << " - " << v << ", " << w << " - " << neigh_aux.x << ", " << neigh_aux.y << ", " << neigh_aux.z << " : " << winner.planes[m].v[v][w] << std::endl;
                            }
                        }
                    }
                }
            }


            // Neighbour temptative value computation
            for (int i=0; i < 3; i++) {
                for (int j=0; j < 3; j++) {
                    for (int k=0; k < 3; k++) {
                        isnewpos[i][j][k] = false;
                        wascomputed[i][j][k] = false;
                        if ((i == 1) && (j == 1) && (k == 1)) {
                            fi = Node3D(image.at(winner.p, 0), image.at(winner.p, 1),image.at(winner.p, 2));
                            imcenter[1][1][1] = to_spherical(fi);
                            valcenter[1][1][1] = u_surface.at(winner.p);
                            continue;
                        }
                        Node3 d(i - 1, j - 1, k - 1);
                        neigh = winner.p + d;
                        isnewpos[i][j][k] = false;
                        if (!u_surface.contains(neigh)) {
                            fi = Node3D(image.at(winner.p, 0), image.at(winner.p, 1),image.at(winner.p, 2));
                            // fi = Node3D(MAX_VAL, MAX_VAL, MAX_VAL);
                            valcenter[i][j][k] = MAX_VAL;
                        }
                        else {
                            fi = Node3D(image.at(neigh, 0), image.at(neigh, 1),image.at(neigh, 2));
                            valcenter[i][j][k] = u_surface.at(neigh);
                        }
                        if (isinf(norm(fi)) || isnan(norm(fi)))
                            fi = Node3D(MAX_VAL, MAX_VAL, MAX_VAL);
                        imcenter[i][j][k] = to_spherical(fi);
                        if (u_surface.contains(neigh) && state.at(neigh) != P_ALIVE) {
                            neighD = Node3D(neigh.x, neigh.y, neigh.z);
                            double val_neigh = GetVal3D(image, u_surface, winner, neighD, fi, h, interp, mode);
                            // std::cout << neighD.x << ", " << neighD.y << ", " << neighD.z << " : " << val_neigh << std::endl;
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
                        wascomputed[i][j][k] = true;

                        key = neigh.x + u_surface.rows*(neigh.y + u_surface.cols*neigh.z);
                        if (state.at(neigh) == P_TRIAL) {
                            mapa_trial_it = mapa_trial.find(key);
                            trial_set.erase(mapa_trial_it->second);
                            mapa_trial.erase(mapa_trial_it);
                        }
                        else {
                            state.at(neigh) = P_TRIAL;
                        }

                        int cont = 0, cont2 = 0;
                        int sum = abs(d.x) + abs(d.y) + abs(d.z);
                        ndirs = (sum == 3) ? 3 : 4;

                        AugWmm3 new_w(N, gamma, ndirs);
                        new_w.p = neigh;
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
                        new_w.ndirs = ndirs;
                        new_w.v = valcenter[d.x + 1][d.y + 1][d.z + 1];
                        // std::cout << "update: " << new_w.p.x << ", " << new_w.p.y << ", " << new_w.p.z << " - " << new_w.v << std::endl;
                        for (int m=0; m < ndirs; m++) {
                            new_w.planes[m].dirs[0] = dirs[m];
                            new_w.planes[m].dirs[1] = dirs[(m+1)%ndirs];
                            Node3 pos1 = d + dirs[m] + Node3(1, 1, 1);
                            Node3 pos2 = d + dirs[(m+1)%ndirs] + Node3(1, 1, 1);
                            Node3 esquina = d + dirs[m] + dirs[(m+1)%ndirs] + Node3(1, 1, 1);
                            setCoeffs3D(imcenter, new_w.planes[m].fm, d, new_w.planes[m].dirs, interp);
                            if (wascomputed[pos1.x][pos1.y][pos1.z] || wascomputed[pos2.x][pos2.y][pos2.z] ||
                                    wascomputed[esquina.x][esquina.y][esquina.z])
                                continue;
                            neighD = wmm3D::Node3D(neigh.x, neigh.y, neigh.z);
                            d1 = wmm3D::Node3D(dirs[m].x, dirs[m].y, dirs[m].z);
                            d2 = wmm3D::Node3D(dirs[(m+1)%ndirs].x, dirs[(m+1)%ndirs].y, dirs[(m+1)%ndirs].z);
                            for (int v=0; v <= gamma; v++) {
                                for (int w=0; w <= gamma; w++) {
                                    epsilon1 = v*step;
                                    epsilon2 = w*step;
                                    neighD_aux = neighD + epsilon1*d1 + epsilon2*d2;
                                    if (!image.contains(neighD_aux)) {
                                        continue;
                                    }
                                    neigh_aux = Node3(round(gamma * neighD_aux.x), round(gamma * neighD_aux.y), round(gamma * neighD_aux.z));
                                    aug_key = neigh_aux.x + gamma * u_surface.rows*(neigh_aux.y + gamma * u_surface.cols*neigh_aux.z);
                                    if (v == 0 && w == 0)
                                        value = valcenter[i][j][k];
                                    else if (v == gamma && w == 0)
                                        value = valcenter[pos1.x][pos1.y][pos1.z];
                                    else if (v == 0 && w == gamma)
                                        value = valcenter[pos2.x][pos2.y][pos2.z];
                                    else if (v == gamma && w == gamma)
                                        value = valcenter[esquina.x][esquina.y][esquina.z];
                                    else {
                                        std::pair<double, double> epsilon = std::pair<double, double>(epsilon1, epsilon2);
                                        if (u_surface.contains(neighD_aux)) {
                                            fi = GetFValue(image, new_w, m, interp, epsilon, h);
                                            value = GetVal3D(image, u_surface, winner, neighD_aux, fi, h, interp, mode);
                                        }
                                        else {
                                            value = MAX_VAL;
                                        }
                                    }

                                    if (augmented_values.find(aug_key) == augmented_values.end() || (value < augmented_values[aug_key])) {
                                        // std::cout << aug_key << ", " << neighD_aux.x << ", " << neighD_aux.y << ", " << neighD_aux.z << " <-> " << aug_key << " : " << v << ", " << w << " - " << value << std::endl;
                                        augmented_values[aug_key] = value;
                                    }
                                }
                            }
                            cont++;
                        }
                        pr_trial = std::pair<double, AugWmm3 >(valcenter[i][j][k], new_w);
                        trial_set_it = trial_set.insert(pr_trial);
                        pr_mapa = std::pair<int, std::multimap<double, AugWmm3 >::iterator>(key, trial_set_it);
                        mapa_trial.insert(pr_mapa);

                        u_surface.at(new_w.p) = valcenter[i][j][k];
                    }
                }
            }

            del(&winner);
        }

        free(state.data);
        return;

    }

    Grid3 augAniWmm3DIsoSurface3D(Grid3 &image, std::vector<Node3> &initials,
                                      Node3D &h, int interp, int mode, int N, int gamma) {
        Grid3 u_surface = Grid3(MAX_VAL, image.rows, image.cols, image.channels, 1);
        augAniWmm3DIsoSurface3D(image, initials, h, interp, mode, N, gamma, u_surface);
        return u_surface;
    }
    
}
