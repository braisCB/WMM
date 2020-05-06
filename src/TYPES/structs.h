#include <stdlib.h>
#include <math.h>

#ifndef STRUCTS_H
#define	STRUCTS_H

namespace structs {

    const int yarray[8] = {-1, -1, -1, 0, 1, 1,  1,  0};
    const int xarray[8] = {-1,  0,  1, 1, 1, 0, -1, -1};


    const int M_GRADIENT = 0;
    const int M_HOPFLAX = 1;
    const int M_GOLDENSEARCH = 2;
    const int S_RIGHT = 0;
    const int S_LEFT = 1;
    const unsigned char P_FAR = 0;
    const unsigned char P_ALIVE = 1;
    const unsigned char P_TRIAL = 2;
    const unsigned char P_INTERMEDIATE = 3;
    const unsigned char P_DEFINITE_TRIAL = 2;
    const double TAU = 1e-05;
    const double EPSILON = 1e-12;
    const double MAX_VAL = 100000000;
    const double RESPHI = (sqrt(5.0)-1.0)/2.0;

    template<typename _Tp> inline _Tp* initArray(int anum, _Tp v) {
        _Tp *out = (_Tp*) malloc(anum*sizeof(_Tp));
        for (int i=0; i < anum; i++) {
            out[i] = _Tp(v);
        }
        return out;
    }

    template<typename _Tp> struct Array_ {
        char mode;
        int* shape;
        int ndim;
        int size;
        _Tp* data;
        
        inline __init_info(int *_shape, int _ndim, char _mode='C') {
            this->shape = new int[_ndim];
            this->ndim = _ndim;
            this->size = 1;
            for (int i=0; i < _ndim; i++) {
                this->ndim[i] = _shape[i];
                this->size *= _shape[i];
            }
        }

        inline Array_() : ndim(_Tp(0)), size(_Tp(0)), mode('C'), data(NULL) {}
        inline Array_(int _size, char _mode='C') {
            int _shape[1] = { _size };
            this->__init_info(_shape, 1, _mode);
            this->data = new _Tp[this->size];
        }
        inline Array_(_Tp _v, int _size, char _mode='C') {
            int _shape[1] = { _size };
            this->__init_info(_shape, 1, _mode);
            this->data = initArray(this->size, _v);
        }
        inline Array_(_Tp* _data, int _size, char _mode='C') {
            int _shape[1] = { _size };
            this->__init_info(_shape, 1, _mode);
            this->data = data;
        }
        inline Array_(int *_shape, int _ndim, char _mode='C') {
            this->__init_info(_shape, _ndim, _mode);
            this->data = new _Tp[this->size];
        }
        inline Array_(_Tp* _data, int *_shape, int _ndim, char _mode='C') {
            this->__init_info(_shape, _ndim, _mode);
            this->data = data;
        }
        inline Array_(_Tp _v, int *_shape, int _ndim, char _mode='C') {
            this->__init_info(_shape, _ndim, _mode);
            this->data = initArray(this->size, _v);
        }

        template<class _cl> _Tp& operator[] (_cl index) {
            return this->at[index];
        }

        inline _Tp& at(int pos) {
            return data[pos];
        }
        
        inline _Tp& at(const Array_<int> &pos) {
            if (pos.size > this->ndim) {
                throw "at: position dimension is higher than array dimension"
            }
            if ((mode == 'C') || (mode == 'C')) {
                int position = pos[ndim - 1];
                for (int i=ndim - 2; i >= 0; i--) {
                    position = position * this->shape[i] + pos[i];
                }
                return position;
            }
            else if ((mode == 'C') || (mode == 'c')) {
                int position = pos[0];
                for (int i=1; i < pos.ndim; i++) {
                    position = position * this->shape[i] + pos[i];
                }
                return position;
            }
            throw "mode should be 'C' or 'C'";
        }

        template<typename _Tp2> inline bool contains(const Array_<_Tp2> &pos) {
            if (pos.size > this->ndim) {
                throw "contains: position dimension is higher than array dimension"
            }
            for (int i=0; i < pos.ndim; i++) {
                if (pos[i] < 0 || pos[i] >= this->shape[i])
                    return false;
            }
            return true;
        }

        template<typename _Tp2> inline bool has_same_shape_as(const Array_<_Tp2> &pos) {
            if (pos.size != this->size || pos.ndim != this->ndim) {
                return false;
            }
            for (int i=0; i < pos.ndim; i++) {
                if (pos.shape[i] != this->shape[i])
                    return false;
            }
            return true;
        }

        inline _Tp sum() {
            _Tp result = _Tp(0);
            for (int i = 0; i < this->size; i++) {
                result += this->data[i];
            }
            return result;
        }

        template<typename _Tp2> inline void copyTo(const Array_<_Tp2> &array) {
            array.__init_info(this->shape, this->ndim, this->mode);
            for (int i = 0; i < this->size; i++) {
                array.data[i] = _Tp2(this->data[i]);
            }
        }

    };

    typedef Array_<int> Array;
    typedef Array_<double> ArrayD;

    template<typename _Tp, typename _Tp2> inline auto operator+(const Array_<_Tp> &array1, const Array_<_Tp2> &array2) {
        if (!array1.has_same_shape_as(array2))
            throw "+: dimension mismatch";

        Array_<decltype(array1.data[0] + array2.data[0])> result(array1.shape, array1.ndim, array1.mode);
        for (int i=0; i < array1.size; i++) {
            result.data[i] = array1.data[i] + array2.data[i];
        }
        return result;
    }
    
    template<typename _Tp, typename _Tp2> inline auto operator+(const Array_<_Tp> &array1, const _Tp2 &v) {
        Array_<decltype(array1.data[0] + v)> result(array1.shape, array1.ndim, array1.mode);
        for (int i=0; i < array1.size; i++) {
            result.data[i] = array1.data[i] + v;
        }
        return result;
    }
    
    template<typename _Tp, typename _Tp2> inline auto operator+(const _Tp2 &v, const Array_<_Tp> &array1) {
        return array1 + v;
    }
    
    template<typename _Tp, typename _Tp2> inline auto operator*(const Array_<_Tp> &array1, const Array_<_Tp2> &array2) {
        if (!array1.has_same_shape_as(array2))
            throw "*: dimension mismatch";

        Array_<decltype(array1.data[0] * array2.data[0])> result(array1.shape, array1.ndim, array1.mode);
        for (int i=0; i < array1.size; i++) {
            result.data[i] = array1.data[i] * array2.data[i];
        }
        return result;
    }
    
    template<typename _Tp, typename _Tp2> inline auto operator*(const Array_<_Tp> &array1, const _Tp2 &v) {
        Array_<decltype(array1.data[0] * v)> result(array1.shape, array1.ndim, array1.mode);
        for (int i=0; i < array1.size; i++) {
            result.data[i] = array1.data[i] * v;
        }
        return result;
    }
    
    template<typename _Tp, typename _Tp2> inline auto operator*(const _Tp2 &v, const Array_<_Tp> &array1) {
        return array1 * v;
    }
    
    template<typename _Tp, typename _Tp2> inline auto operator-(const Array_<_Tp> &array1, const Array_<_Tp2> &array2) {
        if (!array1.has_same_shape_as(array2))
            throw "*: dimension mismatch";

        Array_<decltype(array1.data[0] - array2.data[0])> result(array1.shape, array1.ndim, array1.mode);
        for (int i=0; i < array1.size; i++) {
            result.data[i] = array1.data[i] - array2.data[i];
        }
        return result;
    }
    
    template<typename _Tp, typename _Tp2> inline auto operator-(const Array_<_Tp> &array1, const _Tp2 &v) {
        Array_<decltype(array1.data[0] - v)> result(array1.shape, array1.ndim, array1.mode);
        for (int i=0; i < array1.size; i++) {
            result.data[i] = array1.data[i] - v;
        }
        return result;
    }
    
    template<typename _Tp, typename _Tp2> inline auto operator-(const _Tp2 &v, const Array_<_Tp> &array1) {
        Array_<decltype(array1.data[0] - v)> result(array1.shape, array1.ndim, array1.mode);
        for (int i=0; i < array1.size; i++) {
            result.data[i] = v - array1.data[i];
        }
        return result;
    }
    
    template<typename _Tp, typename _Tp2> inline auto operator/(const Array_<_Tp> &array1, const Array_<_Tp2> &array2) {
        if (!array1.has_same_shape_as(array2))
            throw "*: dimension mismatch";

        Array_<decltype(array1.data[0] / array2.data[0])> result(array1.shape, array1.ndim, array1.mode);
        for (int i=0; i < array1.size; i++) {
            result.data[i] = array1.data[i] / array2.data[i];
        }
        return result;
    }
    
    template<typename _Tp, typename _Tp2> inline auto operator/(const Array_<_Tp> &array1, const _Tp2 &v) {
        Array_<decltype(array1.data[0] / v)> result(array1.shape, array1.ndim, array1.mode);
        for (int i=0; i < array1.size; i++) {
            result.data[i] = array1.data[i] / v;
        }
        return result;
    }
    
    template<typename _Tp, typename _Tp2> inline auto operator/(const _Tp2 &v, const Array_<_Tp> &array1) {
        Array_<decltype(array1.data[0] / v)> result(array1.shape, array1.ndim, array1.mode);
        for (int i=0; i < array1.size; i++) {
            result.data[i] = v / array1.data[i];
        }
        return result;
    }
    
    template<typename _Tp> inline double norm(const Array_<_Tp> &array1) {
        double result = 0.0;
        for (int i=0; i < array1.size; i++) {
            result += array1.data[i] * array1.data[i];
        }
        return sqrt(result);
    }

    template<typename _Tp> inline auto sqrt(const Array_<_Tp> &array1) {
        Array_<decltype(array1.data[0] / v)> result(array1.shape, array1.ndim, array1.mode);
        for (int i=0; i < array1.size; i++) {
            result.data[i] = sqrt(array1.data[i]);
        }
        return result;
    }

    template<typename _Tp> inline Array_<double> to_polar(const Array_<_Tp> &array1) {
        if (array1.size != 2)
            throw "to_polar: array should contain only 2 elements"
        double n = norm(array1);
        Array_<double> polar;
        if (n == 0.0)
            polar = Array_<double>(0, 3, array1.mode);
        else {
            double polar_values[3] = {array1.data[1] / n, array1.data[0] / n, n};
            polar = Array_<double>(polar_values, 3, array1.mode);
        }
        return polar;
    }

    template<typename _Tp> inline Array_<_Tp> to_cartesian(const Array_<double> &array1) {
        if (array1.size != 3)
            throw "to_cartesian: array should contain only 3 elements"
        double n = norm(polar);
        Array_<_Tp> cartesian;
        if (n == 0)
            cartesian = Array_<_Tp>(0, 2, array1.mode);
        else {
            _Tp cartesian_values[2] = {
                _Tp(array1.data[1] * array1.data[2]/ n),
                _Tp(array1.data[0] * array1.data[2] / n)};
            cartesian = Array_<_Tp>(cartesian_values, 3, array1.mode);
        }
        return cartesian;
    }

    template<typename _Tp> struct Segment_ {
        Array dirs;
        Array_<_Tp> vs;
        Array_<_Tp> ms;

        inline Segment_() {}
    }

    template<typename _Tp> struct Wavefront_ {
        int ndir;
        int nsegments;
        Segment_<_Tp> *segments;

        inline Wavefront_() : nsegments(0), ndir(-1), segments(NULL) {}
    }

    template<typename _Tp> struct AugWmm_ {
        Node p;
        _Tp* v;
        int dir;
        int child_dir;
        Node3D* fm;
        double* m;
        int N;
        int M;
        int gamma;

        inline AugWmm_() {}

        inline AugWmm_(int _N, int _M, int _gamma) : N(_N), M(_M), gamma(_gamma), v((_Tp*) calloc(2*_gamma + 1,sizeof(_Tp))),
                fm((Node3D*) calloc(_N,sizeof(Node3D))), m((_Tp*) calloc(_M,sizeof(_Tp))) {}

        inline AugWmm_(AugWmm_<_Tp> *wmm) : p(wmm->p), dir(wmm->dir), child_dir(wmm->child_dir), N(wmm->N), M(wmm->M),
                gamma(wmm->gamma), v((_Tp*) calloc(2*wmm->gamma + 1,sizeof(_Tp))),
                fm((Node3D*) calloc(wmm->N,sizeof(Node3D))), m((_Tp*) calloc(wmm->M,sizeof(_Tp))) {

            for (int i=0; i < wmm->M; i++) {
                this->m[i] = wmm->m[i];
            }
            for (int i=0; i < wmm->N; i++) {
                this->fm[i] = wmm->fm[i];
            }
            for (int i=0; i <= 2*wmm->gamma; i++) {
                this->v[i] = wmm->v[i];
            }

        }

    };

    template<typename _Tp> inline void del(AugWmm_<_Tp> *wmm) {
        free(wmm->v);
        free(wmm->fm);
        free(wmm->m);
    }

    typedef Grid_<double> Grid;
    typedef AugWmm_<double> AugWmm;

    template<typename _Tp> struct Wmm0_ {
        Node p;
        _Tp v[7];
        wmm::NodeD f[7];
        int dir;
        int angle;
    };

    template<typename _Tp, int N> struct Wmm_ : public Wmm0_<_Tp> {
        _Tp m[N];
        _Tp fm[N];
    };

    template<typename _Tp> struct Wmm_<_Tp, 0> : public Wmm0_<_Tp> {
        _Tp m[1];
        _Tp fm[1];
    };

    template<typename _Tp, int N> struct AniWmm_ : public Wmm0_<_Tp> {
        _Tp m[N];
        Node3D fm[N];
    };

    template<typename _Tp> struct AniWmm_<_Tp, 0> : public Wmm0_<_Tp> {
        _Tp m[1];
        Node3D fm[1];
    };

    typedef Wmm_<double, 4> WmmNode;
    typedef AniWmm_<double, 4> AniWmmNode;

}
#endif	/* STRUCTS_H */

