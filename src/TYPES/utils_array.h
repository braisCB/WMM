#include <stdlib.h>
#include <math.h>
#include "../TYPES/structs.h"
#include "../TYPES/interpolation.h"

#ifndef UTILS_H
#define	UTILS_H

namespace interp {

    const double MAX_VAL = 100000000;

    template<typename _Tp> Array_<_Tp> GetGradientFromGrid(const Array_<_Tp> &gradient, const Array_<int> &pos) {
        if (pos.size != 2)
            throw "GetGradientFromGrid: pos should contain only 2 elements";
        int output_size = gradient.shape[gradient.ndim - 1];
        Array_<int> pos_aux(gradient.ndim);
        Array_<_Tp> result(output_size);
        for (int i=0; i<pos.size; i++)
            pos_aux.data[i] = pos.data[i];
        for (int i=0; i<output_size; i++) {
            pos_aux.data[gradient.ndim - 1] = i;
            result.data[i] = gradient[pos_aux];
        }
        return result;
    }

    template<typename _Tp> Array_<_Tp> GetInterpolatedGradient(
            const Array_<_Tp> &gradient, const Array_<double> &posD, const Array_<int> &dirs,
            int interp
    ) {
        Array_<int> pos;
        posD.copyTo(pos);
        auto diff = posD - pos;
        double epsilon = (double) diff.sum();
        Array_<_Tp> result = GetInterpolatedGradient(gradient, pos, dirs, epsilon, interp);
        return
    }

    template<typename _Tp> Array_<_Tp> GetInterpolatedGradient(
            const Array_<_Tp> &gradient, const Array_<int> &pos, const Array_<int> &dirs, double epsilon,
            int interp
    ) {
        if (pos.size != 2)
            throw "GetInterpolatedGradient: pos should contain only 2 elements";
        _Tp value = _Tp(MAX_VAL);
        switch (interp) {
            case I_LINEAR:
                Array_<_Tp> f0 = GetGradientFromGrid(gradient, pos);
                Array_<int> p1 = pos + dirs;
                Array_<_Tp> f1 = (gradient.contains(p1)) ? GetGradientFromGrid(gradient, p1) : f0;
                value = get_linear_score(f0, f1, epsilon);
                break;
            case I_QUADATRIC:
                Array_<_Tp> fs[3];
                Array_<int> p = pos - dirs;
                int init = (gradient.contains(p)) ? -1 : 0, end = 3 + init;
                p = pos + dirs;
                if (!gradient.contains(p))
                    value = GetGradientFromGrid(gradient, pos);
                else {
                    for (int i=init; i < end; i++) {
                        p = pos + i * dirs;
                        fs[i - init] = (gradient.contains(p)) ? GetGradientFromGrid(gradient, p1) : fs[i - init - 1];
                    }
                    _Tp m[3];
                    get_quad_coeffs(fs, m, -init, 1-init);
                    value = get_quad_score(m[0], m[1], m[2], epsilon);
                }
                break;
            default: //I_HERMITE and I_PCHIP
                Array_<_Tp> fs[4];
                Array_<int> p = pos - dirs;
                int init = (gradient.contains(p)) ? -1 : 0, end = 4 + init;
                for (int i=init; i < end; i++) {
                    p = pos + i * dirs;
                    if (!gradient.contains(p)) {
                        end = i;
                        break;
                    }
                    fs[i - init] = GetGradientFromGrid(gradient, p);
                }
                _Tp m[2];
                if (interp == I_SPLINE) {
                    get_spline_coeffs(fs, m, 4, -init, 1-init);
                    value = get_spline_score(fs[-init], fs[1-init], m[0], m[1], epsilon);
                }
                else if (interp == I_HERMITE) {
                    get_hermite_coeffs(fs, m, 4, -init, 2-init);
                    value = get_hermite_score(fs[-init], fs[1-init], m[0], m[1], epsilon);
                }
                else {
                    get_pchip_coeffs(fs, m, 4, -init, 2-init);
                    value = get_pchip_score(fs[-init], fs[1-init], m[0], m[1], epsilon);
                }
                break;
        }
        return value;
    }

    template<typename _Tp> _Tp GetInterpValue(Segment_<_Tp> &segment, double epsilon, int interp) {

        double value = _Tp(MAX_VAL);
        switch (interp) {
            case I_LINEAR:
                value = (1.0 - epsilon)*segment.vs[0] + epsilon*segment.vs[1];
                break;
            case wmm::I_QUADATRIC:
                value = get_quad_score(segment.ms[0], segment.ms[1], segment.ms[2], epsilon);
                break;
            case wmm::I_SPLINE:
                value = get_spline_score(segment.vs[0], segment.vs[1], segment.ms[0], segment.vs[1], epsilon);
                break;
            case wmm::I_HERMITE:
                value = get_hermite_score(segment.vs[0], segment.vs[1], segment.ms[0], segment.vs[1], epsilon);
                break;
            default: //I_PCHIP
                value = get_pchip_score(segment.vs[0], segment.vs[1], segment.ms[0], segment.vs[1], epsilon);
        }
        if ((value < segment.vs[0] && value < segment.vs[1]) || (value >= segment.vs[0] && value >= segment.vs[1]))
            return GetInterpValue(segment, epsilon, I_LINEAR);
        return value;
    }


    template<typename _Tp> Array_<_Tp> GetInterpolatedIntegral(
            const Array_<_Tp> &gradient, const Array_<int> &pos, const Array_<int> &dirs, const Array_<double> &neighD,
            double epsilon, int interp
    ) {

        if (pos.size != 2)
            throw "GetInterpolatedGradient: pos should contain only 2 elements";
        auto posD = pos + epsilon * dirs;
        auto diff = neighD - posD;
        _Tp value = _Tp(MAX_VAL);
        switch (interp) {
            case I_LINEAR:
                Array_<_Tp> f0 = GetGradientFromGrid(gradient, posD);
                Array_<_Tp> f1 = GetGradientFromGrid(gradient, neighD);
                value = get_linear_integral_score(f0, f1);
                break;
            case I_QUADATRIC:
                Array_<_Tp> fs[3];
                auto p = posD - diff;
                int init = (gradient.contains(p)) ? -1 : 0, end = 3 + init;
                p = posD + (end - 1)*diff;
                if (!gradient.contains(p))
                    value = GetInterpolatedIntegral(gradient, pos, dirs, neighD, epsilon, I_LINEAR);
                else {
                    for (int i=init; i < end; i++) {
                        p = pos + i * dirs;
                        fs[i - init] = (gradient.contains(p)) ? GetGradientFromGrid(gradient, p1) : fs[i - init - 1];
                    }
                    _Tp m[3];
                    get_quad_coeffs(fs, m, -init, 1-init);
                    value = get_quad_integral_score(m[0], m[1], m[2], epsilon);
                }
                break;
            default: //I_HERMITE and I_PCHIP
                Array_<_Tp> fs[4];
                auto p = posD - diff;
                int init = (gradient.contains(p)) ? -1 : 0, end = 4 + init;
                for (int i=init; i < end; i++) {
                    p = posD + i * diff;
                    if (!gradient.contains(p)) {
                        end = i;
                        break;
                    }
                    fs[i - init] = GetGradientFromGrid(gradient, p);
                }
                _Tp m[2];
                if (interp == I_SPLINE) {
                    get_spline_coeffs(fs, m, end-init, -init, 1-init);
                    value = get_spline_integral_score(fs[-init], fs[1-init], m[0], m[1], epsilon);
                }
                else if (interp == I_HERMITE) {
                    get_hermite_coeffs(fs, m, end-init, -init, 2-init);
                    value = get_hermite_integral_score(fs[-init], fs[1-init], m[0], m[1], epsilon);
                }
                else {
                    get_pchip_coeffs(fs, m, end-init, -init, 2-init);
                    value = get_pchip_integral_score(fs[-init], fs[1-init], m[0], m[1], epsilon);
                }
                break;
        }
        return value;
    }

    template<typename _Tp> _Tp GetSurfaceValue(
            const Array_<_Tp> &gradient, const Wavefront_<_Tp> &wavefront, const Array_<double> &neighD, int interp,
            int mode
    ) {

        if (pos.size != 2)
            throw "GetInterpolatedGradient: pos should contain only 2 elements";
        auto posD = pos + epsilon * dirs;
        auto diff = neighD - posD;
        _Tp value = _Tp(MAX_VAL);
        switch (interp) {
            case I_LINEAR:
                Array_<_Tp> f0 = GetGradientFromGrid(gradient, posD);
                Array_<_Tp> f1 = GetGradientFromGrid(gradient, neighD);
                value = get_linear_integral_score(f0, f1);
                break;
            case I_QUADATRIC:
                Array_<_Tp> fs[3];
                auto p = posD - diff;
                int init = (gradient.contains(p)) ? -1 : 0, end = 3 + init;
                p = posD + (end - 1)*diff;
                if (!gradient.contains(p))
                    value = GetInterpolatedIntegral(gradient, pos, dirs, neighD, epsilon, I_LINEAR);
                else {
                    for (int i=init; i < end; i++) {
                        p = pos + i * dirs;
                        fs[i - init] = (gradient.contains(p)) ? GetGradientFromGrid(gradient, p1) : fs[i - init - 1];
                    }
                    _Tp m[3];
                    get_quad_coeffs(fs, m, -init, 1-init);
                    value = get_quad_integral_score(m[0], m[1], m[2], epsilon);
                }
                break;
            default: //I_HERMITE and I_PCHIP
                Array_<_Tp> fs[4];
                auto p = posD - diff;
                int init = (gradient.contains(p)) ? -1 : 0, end = 4 + init;
                for (int i=init; i < end; i++) {
                    p = posD + i * diff;
                    if (!gradient.contains(p)) {
                        end = i;
                        break;
                    }
                    fs[i - init] = GetGradientFromGrid(gradient, p);
                }
                _Tp m[2];
                if (interp == I_SPLINE) {
                    get_spline_coeffs(fs, m, end-init, -init, 1-init);
                    value = get_spline_integral_score(fs[-init], fs[1-init], m[0], m[1], epsilon);
                }
                else if (interp == I_HERMITE) {
                    get_hermite_coeffs(fs, m, end-init, -init, 2-init);
                    value = get_hermite_integral_score(fs[-init], fs[1-init], m[0], m[1], epsilon);
                }
                else {
                    get_pchip_coeffs(fs, m, end-init, -init, 2-init);
                    value = get_pchip_integral_score(fs[-init], fs[1-init], m[0], m[1], epsilon);
                }
                break;
        }
        return value;
    }

}
#endif	/* UTILS_H */

