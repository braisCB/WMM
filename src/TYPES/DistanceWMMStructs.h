#include <stdlib.h>
#include <math.h>

#ifndef MPPSTRUCTS_H
#define	MPPSTRUCTS_H

namespace wmm {

    const int yarray[8] = {-1, -1, -1, 0, 1, 1,  1,  0};
    const int xarray[8] = {-1,  0,  1, 1, 1, 0, -1, -1};

    const int I_LINEAR = 0;
    const int I_HERMITE = 1;
    const int I_PCHIP = 2;
    const int I_SPLINE = 3;
    const int I_QUADATRIC = 4;
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

    template<typename _Tp> struct Node_ {
        _Tp y;
        _Tp x;

        inline Node_(_Tp _y = 0, _Tp _x = 0) : y(_y), x(_x) {}
    };
    
    template<typename _Tp> struct Node3_ {
        _Tp x;
        _Tp y;
        _Tp z;

        inline Node3_(_Tp _x = 0, _Tp _y = 0, _Tp _z = 0) : x(_x), y(_y), z(_z) {}

    };
    
    template<typename _Tp, typename _Tp2> inline auto operator+(const Node_<_Tp> &p1, const Node_<_Tp2> &p2) {
        Node_<decltype(p1.y + p2.y)> p(p1.y + p2.y, p1.x + p2.x);
        return p;
    }

    template<typename _Tp, typename _Tp2> inline auto operator-(const Node_<_Tp> &p1, const Node_<_Tp2> &p2) {
        Node_<decltype(p1.y - p2.y)> p(p1.y - p2.y, p1.x - p2.x);
        return p;
    }

    template<typename _Tp, typename _Tp2> inline auto operator*(const _Tp &v, const Node_<_Tp2> &p) {
        Node_<decltype(v * p.y)> np(v*p.y, v*p.x);
        return np;
    }
    
    template<typename _Tp, typename _Tp2> inline auto operator*(const Node_<_Tp> &p1, const Node_<_Tp2> &p2) {
        Node_<decltype(p1.y * p2.y)> p(p1.y * p2.y, p1.x * p2.x);
        return p;
    }

    template<typename _Tp, typename _Tp2> inline auto operator/(const Node_<_Tp> &p1, const Node_<_Tp2> &p2) {
        Node_<decltype(p1.y / p2.y)> p(p1.y / p2.y, p1.x / p2.x);
        return p;
    }

    template<typename _Tp, typename _Tp2> inline auto operator/(const Node_<_Tp> &p1, const _Tp2 &p2) {
        Node_<decltype(p1.y / p2)> p(p1.y / p2, p1.x / p2);
        return p;
    }

    template<typename _Tp> inline double norm(const Node_<_Tp> &p) {
        return sqrt((double) (p.x*p.x + p.y*p.y));
    }

    template<typename _Tp> inline auto operator==(const Node_<_Tp> &p1, const Node_<_Tp> &p2) {
        return (p1.x == p2.x && p1.y == p2.y);
    }
    
    template<typename _Tp> inline Node3_<_Tp> operator+(const Node3_<_Tp> &p1, const Node3_<_Tp> &p2) {
        Node3_<_Tp> p(p1.x + p2.x, p1.y + p2.y, p1.z + p2.z);
        return p;
    }

    template<typename _Tp> inline Node3_<_Tp> operator+(const _Tp &v, const Node3_<_Tp> &p2) {
        Node3_<_Tp> p(v + p2.x, v + p2.y, v + p2.z);
        return p;
    }

    template<typename _Tp> inline Node3_<_Tp> operator+(const Node3_<_Tp> &p2, const _Tp &v) {
        Node3_<_Tp> p(v + p2.x, v + p2.y, v + p2.z);
        return p;
    }

    template<typename _Tp> inline Node3_<_Tp> operator-(const Node3_<_Tp> &p1, const Node3_<_Tp> &p2) {
        Node3_<_Tp> p(p1.x - p2.x, p1.y - p2.y, p1.z - p2.z);
        return p;
    }
    
    template<typename _Tp> inline Node3_<_Tp> operator*(const Node3_<_Tp> &p1, const Node3_<_Tp> &p2) {
        Node3_<_Tp> p(p1.x * p2.x, p1.y * p2.y, p1.z * p2.z);
        return p;
    }

    template<typename _Tp, typename _Tp2> inline auto operator*(const Node3_<_Tp> &p, const _Tp2 &v) {
        Node3_<decltype(v * p.y)> np(v*p.x, v*p.y, v*p.z);
        return np;
    }

    template<typename _Tp, typename _Tp2> inline auto operator*(const _Tp2 &v, const Node3_<_Tp> &p) {
        return p * v;
    }

    template<typename _Tp, typename _Tp2> inline auto operator/(const Node3_<_Tp> &p1, const Node3_<_Tp2> &p2) {
        Node3_<decltype(p1.y / p2.y)> p(p1.x / p2.x, p1.y / p2.y, p1.z / p2.z);
        return p;
    }

    template<typename _Tp> inline Node3_<_Tp> operator/(const _Tp &v, const Node3_<_Tp> &p) {
        Node3_<_Tp> np(v/p.x, v/p.y, v/p.z);
        return np;
    }

    template<typename _Tp> inline Node3_<_Tp> operator/(const Node3_<_Tp> &p, const _Tp &v) {
        Node3_<_Tp> np(p.x/v, p.y/v, p.z/v);
        return np;
    }

    template<typename _Tp> inline double norm(const Node3_<_Tp> &p) {
        return sqrt((double) (p.x*p.x + p.y*p.y + p.z*p.z));
    }

    template<typename _Tp> inline Node3_<_Tp> node_sqrt(const Node3_<_Tp> &p) {
        Node3_<_Tp> np(sqrt((double) p.x), sqrt((double) p.y), sqrt((double) p.z));
        return np;
    }
    
    template<typename _Tp> inline Node3_<_Tp> to_polar(const Node_<_Tp> &p) {
        _Tp n = wmm::norm(p);
        Node3_<_Tp> polar;
        if (n == 0.0)
            polar = Node3_<_Tp>(0.0, 0.0, 0.0);
        else
            polar = Node3_<_Tp>(p.x/n, p.y/n, n);
        return polar;
    }
    
    template<typename _Tp> inline Node_<_Tp> to_cartesian(const Node3_<_Tp> &p) {
        _Tp n = sqrt(p.x*p.x + p.y*p.y);
        Node_<_Tp> cartesian;
        if (n == 0)
            cartesian = Node_<_Tp>(0.0, 0.0);
        else
            cartesian = Node_<_Tp>(p.z*p.y/n, p.z*p.x/n);
        return cartesian;
    }

    typedef Node_<int> Node;
    typedef Node3_<int> Node3;
    typedef Node_<double> NodeD;
    typedef Node3_<double> Node3D;

    template<typename _Tp> inline _Tp* initArray(int anum, _Tp v) {
        _Tp *out = (_Tp*) malloc(anum*sizeof(_Tp));
        for (int i=0; i < anum; i++) {
            out[i] = v;
        }
        return out;
    }

    template<typename _Tp> struct Grid_ {
        int rows;
        int cols;
        int channels;
        char mode;
        _Tp* data;

        inline Grid_() : rows(0), cols(0), channels(0), mode('F'), data(NULL) {}
        inline Grid_(int _rows, int _cols, int _channels=1, char _mode='F') : rows(_rows), cols(_cols), channels(_channels), mode(_mode), data((_Tp*) calloc(_rows*_cols*_channels,sizeof(_Tp))) {}
        inline Grid_(_Tp* _data, int _rows, int _cols, int _channels=1, char _mode='F') : rows(_rows), cols(_cols), channels(_channels), mode(_mode), data(_data) {}
        inline Grid_(_Tp _v, int _rows, int _cols, int _channels, char _mode='F') : rows(_rows), cols(_cols), channels(_channels), mode(_mode), data(initArray(_rows*_cols*_channels, _v)) {}

        inline _Tp& at(Node &p, int dim=0) {
            if ((mode == 'F') || (mode == 'f')) {
                return data[p.y + rows*(p.x + cols*dim)];
            }
            else if ((mode == 'C') || (mode == 'c')) {
                return data[dim + channels*(p.x + cols*p.y)];
            }
            throw "mode should be 'C' or 'F'";
        }

        template<typename _Tp2> inline bool contains(Node_<_Tp2> &p) {
            return (p.x >= 0 && p.y >= 0 && p.x < cols && p.y < rows);
        }

    };
    
    /*template<typename _Tp> struct AugWmm_ {
        Node p;
        int dir;
        int child_dir;
        int N;
        int M;
        int gamma;
        _Tp* v;
        Node3D* fm;
        double* m;

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
    }*/

    template<typename _Tp> struct AugWmm_aux_ {
        _Tp* v;
        _Tp* d;
        Node dir;
        Node3D* fm;
    };

    template<typename _Tp> struct AugWmm_ {
        Node p;
        Node father;
        AugWmm_aux_<_Tp> planes[2];
        int N;
        int gamma;
        int ndirs;
        _Tp v;
        _Tp d;

        inline AugWmm_() {}

        inline AugWmm_(int _N, int _gamma, int _ndirs=2) : N(_N), gamma(_gamma), ndirs(_ndirs) {
            for (int i=0; i < 2; i++) {
                this->planes[i].v = (_Tp*) calloc((_gamma + 1),sizeof(_Tp));
                this->planes[i].d = (_Tp*) calloc((_gamma + 1),sizeof(_Tp));
                this->planes[i].fm = (Node3D*) calloc(_N,sizeof(Node3D));
            }
        }

    };

    typedef AugWmm_<double> AugWmm;

    template<typename _Tp> inline void del(AugWmm_<_Tp> *wmm) {
        for (int i=0; i < 2; i++) {
            free(wmm->planes[i].v);
            free(wmm->planes[i].d);
            free(wmm->planes[i].fm);
        }
    }

    typedef Grid_<double> Grid;
    typedef AugWmm_<double> AugWmm;

    template<typename _Tp> struct Wmm0_ {
        Node p;
        _Tp v[7];
        _Tp d[7];
        wmm::NodeD f[7];
        int dir;
        int angle;
    };

    template<typename _Tp, int N> struct Wmm_ : public Wmm0_<_Tp> {
        _Tp m[N];
        _Tp dm[N];
        _Tp fm[N];
    };

    template<typename _Tp> struct Wmm_<_Tp, 0> : public Wmm0_<_Tp> {
        _Tp m[1];
        _Tp dm[1];
        _Tp fm[1];
    };
    
    template<typename _Tp, int N> struct AniWmm_ : public Wmm0_<_Tp> {
        _Tp m[N];
        _Tp dm[N];
        Node3D fm[N];
    };

    template<typename _Tp> struct AniWmm_<_Tp, 0> : public Wmm0_<_Tp> {
        _Tp m[1];
        _Tp dm[1];
        Node3D fm[1];
    };

    typedef Wmm_<double, 4> WmmNode;
    typedef AniWmm_<double, 4> AniWmmNode;

}
#endif	/* MPPSTRUCTS_H */

