#include <stdlib.h>
#include <math.h>

#ifndef MPPSTRUCTS3D_H
#define	MPPSTRUCTS3D_H

namespace wmm3D {

    const int yarray[8] = {-1, -1, -1, 0, 1, 1, 1, 0};
    const int xarray[8] = {-1, 0, 1, 1, 1, 0, -1, -1};

    const int I_BILINEAR = 0;
    const int I_HERMITE = 1;
    const int I_PCHIP = 2;
    const int I_SPLINE = 3;
    const int I_QUADATRIC = 4;
    const int I_CUBIC = 9;
    const int I_BICUBIC = 11;
    const int M_GRADIENT = 0;
    const int M_HOPFLAX = 1;
    const int M_GOLDENSEARCH = 2;
    const unsigned char P_ALIVE = 1;
    const unsigned char P_TRIAL = 2;
    double TAU = 1e-03;
    double MAX_VAL = 100000000;
    double RESPHI = (sqrt(5.0)-1.0)/2.0;

    template<typename _Tp> struct Node3_ {
        _Tp x;
        _Tp y;
        _Tp z;

        inline Node3_(_Tp _x = 0, _Tp _y = 0, _Tp _z = 0) : x(_x), y(_y), z(_z) {}

    };
    
    template<typename _Tp> struct Node5_ {
        _Tp v;
        _Tp w;
        _Tp x;
        _Tp y;
        _Tp z;

        inline Node5_(_Tp _v = 0, _Tp _w = 0, _Tp _x = 0, _Tp _y = 0, _Tp _z = 0) : v(_v), w(_w), x(_x), y(_y), z(_z) {}

    };

    template<typename _Tp, typename _Tp2> inline auto operator+(const Node3_<_Tp> &p1, const Node3_<_Tp2> &p2) {
        Node3_<decltype(p1.x + p2.x)> p(p1.x + p2.x, p1.y + p2.y, p1.z + p2.z);
        return p;
    }

    template<typename _Tp, typename _Tp2> inline auto operator-(const Node3_<_Tp> &p1, const Node3_<_Tp2> &p2) {
        Node3_<decltype(p1.x - p2.x)> p(p1.x - p2.x, p1.y - p2.y, p1.z - p2.z);
        return p;
    }

    template<typename _Tp, typename _Tp2> inline auto operator*(const _Tp2 &v, const Node3_<_Tp> &p) {
        Node3_<decltype(v * p.x)> np(v*p.x, v*p.y, v*p.z);
        return np;
    }

    template<typename _Tp, typename _Tp2> inline auto operator*(const Node3_<_Tp> &p1, const Node3_<_Tp2> &p2) {
        Node3_<decltype(p1.x * p2.x)> p(p1.x * p2.x, p1.y * p2.y, p1.z * p2.z);
        return p;
    }

    template<typename _Tp> inline double norm(const Node3_<_Tp> &p) {
        return sqrt((double) (p.x*p.x + p.y*p.y + p.z*p.z));
    }
    
    template<typename _Tp> inline Node5_<_Tp> operator+(const Node5_<_Tp> &p1, const Node5_<_Tp> &p2) {
        Node5_<_Tp> p(p1.v + p2.v, p1.w + p2.w, p1.x + p2.x, p1.y + p2.y, p1.z + p2.z);
        return p;
    }

    template<typename _Tp> inline Node5_<_Tp> operator-(const Node5_<_Tp> &p1, const Node5_<_Tp> &p2) {
        Node5_<_Tp> p(p1.v - p2.v, p1.w - p2.w, p1.x - p2.x, p1.y - p2.y, p1.z - p2.z);
        return p;
    }

    template<typename _Tp, typename _Tp2> inline Node5_<_Tp> operator*(const _Tp2 &v, const Node5_<_Tp> &p) {
        Node5_<_Tp> np(v*p.v, v*p.w, v*p.x, v*p.y, v*p.z);
        return np;
    }

    template<typename _Tp> inline Node3_<_Tp> operator/(const Node3_<_Tp> &p1, const Node3_<_Tp> &p2) {
        Node3_<_Tp> p(p1.x / p2.x, p1.y / p2.y, p1.z / p2.z);
        return p;
    }

    template<typename _Tp, typename _Tp2> inline Node3_<_Tp> operator/(const _Tp2 &v, const Node3_<_Tp> &p) {
        Node3_<_Tp> np(v/p.x, v/p.y, v/p.z);
        return np;
    }

    template<typename _Tp, typename _Tp2> inline Node3_<_Tp> operator/(const Node3_<_Tp> &p, const _Tp2 &v) {
        Node3_<_Tp> np(p.x/v, p.y/v, p.z/v);
        return np;
    }

    template<typename _Tp> inline Node5_<_Tp> operator/(const Node5_<_Tp> &p1, const Node5_<_Tp> &p2) {
        Node5_<_Tp> p(p1.v / p2.v, p1.w / p2.w, p1.x / p2.x, p1.y / p2.y, p1.z / p2.z);
        return p;
    }

    template<typename _Tp, typename _Tp2> inline Node5_<_Tp> operator/(const _Tp2 &v, const Node5_<_Tp> &p) {
        Node5_<_Tp> np(v/p.v, v/p.w, v/p.x, v/p.y, v/p.z);
        return np;
    }

    template<typename _Tp, typename _Tp2> inline Node5_<_Tp> operator/(const Node5_<_Tp> &p, const _Tp2 &v) {
        Node5_<_Tp> np(p.v/v, p.w/v, p.x/v, p.y/v, p.z/v);
        return np;
    }

    typedef Node3_<int> Node3;
    typedef Node3_<double> Node3D;
    typedef Node5_<double> Node5D;
    
    template<typename _Tp> inline Node5_<_Tp> to_spherical(const Node3_<_Tp> &p) {
        _Tp n = norm(p);
        Node5_<_Tp> spherical;
        if (n == 0.0)
            spherical = Node5_<_Tp>(0.0, 0.0, 0.0, 0.0, 0.0);
        else {
            double phi = acos(p.z/n);
            double rho = atan2(p.y, p.x);
            spherical = Node5_<_Tp>(cos(rho), sin(rho), cos(phi), sin(phi), n);
        }
        return spherical;
    }
    
    template<typename _Tp> inline Node3_<_Tp> to_cartesian(const Node5_<_Tp> &p) {
        _Tp nxy = sqrt(p.x*p.x + p.y*p.y);
        _Tp nvw = sqrt(p.v*p.v + p.w*p.w);
        Node3_<_Tp> cartesian;
        if (nvw == 0 || nxy == 0)
            cartesian = Node3_<_Tp>(0.0, 0.0, 0.0);
        else
            cartesian = Node3_<_Tp>(p.z*(p.v/nvw)*(p.y/nxy), p.z*(p.w/nvw)*(p.y/nxy), p.z*p.x/nxy);
        return cartesian;
    }
    
    template<typename _Tp> struct AugWmm3_aux_ {
        _Tp** v;
        Node3 dirs[2];
        Node5D* fm;
    };
    
    template<typename _Tp> struct AugWmm3_ {
        Node3 p;
        Node3 father;
        AugWmm3_aux_<_Tp> *planes;
        int N;
        int gamma;
        int ndirs;
        _Tp v;
        
        inline AugWmm3_() {}

        inline AugWmm3_(int _N, int _gamma, int _ndirs) : N(_N), gamma(_gamma), ndirs(_ndirs) {
            if (_ndirs > 0) {
                this->planes = (AugWmm3_aux_<_Tp> *) malloc(_ndirs*sizeof(AugWmm3_aux_<_Tp>));
                for (int i=0; i < _ndirs; i++) {
                    this->planes[i].v = (_Tp**) calloc((_gamma + 1),sizeof(_Tp*));
                    for (int j=0; j <= _gamma; j++) {
                        this->planes[i].v[j] = (_Tp*) calloc((_gamma + 1),sizeof(_Tp));
                    }
                    this->planes[i].fm = (Node5D*) calloc(_N,sizeof(Node5D));
                }
            }
        }

    };
    
    typedef AugWmm3_<double> AugWmm3;
    
    template<typename _Tp> inline void del(AugWmm3_<_Tp> *wmm) {
        if (wmm->ndirs > 0) {
            for (int i=0; i < wmm->ndirs; i++) {
                for (int j=0; j <= wmm->gamma; j++) {
                    free(wmm->planes[i].v[j]);
                }
                free(wmm->planes[i].v);
                free(wmm->planes[i].fm);
            }
            free(wmm->planes);
        }
    }
    
    template<typename _Tp> inline _Tp* initArray(int anum, _Tp v) {
        _Tp *out = (_Tp*) malloc(anum*sizeof(_Tp));
        for (int i=0; i < anum; i++) {
            out[i] = v;
        }
        return out;
    }

    template<typename _Tp> struct Grid3_ {
        int rows;
        int cols;
        int channels;
        int dim;
        char mode;
        _Tp* data;

        inline Grid3_(int _rows, int _cols, int _channels, int _dim=1, char _mode='F') :
                rows(_rows), cols(_cols), channels(_channels), dim(_dim), mode(_mode), data((_Tp*) calloc(_rows*_cols*_channels*_dim,sizeof(_Tp))) {}
        inline Grid3_(_Tp* _data, int _rows, int _cols, int _channels, int _dim=1, char _mode='F') :
                rows(_rows), cols(_cols), channels(_channels), dim(_dim), mode(_mode), data(_data) {}
        inline Grid3_(_Tp _v, int _rows, int _cols, int _channels, int _dim=1, char _mode='F') :
                rows(_rows), cols(_cols), channels(_channels), dim(_dim), mode(_mode), data(initArray(_rows*_cols*_channels*_dim, _v)) {}

        inline _Tp& at(Node3 &p, int d=0) {
            if ((mode == 'F') || (mode == 'f')) {
                return data[p.x + rows*(p.y + cols*(p.z + channels*d))];
            }
            else if ((mode == 'C') || (mode == 'c')) {
                return data[d + dim*(p.z + channels*(p.y + cols*p.x))];
            }
            throw "mode should be 'C' or 'F'";
        }

        inline bool contains(Node3 &p) {
            return (p.x >= 0 && p.y >= 0 && p.z >= 0 && p.x < rows && p.y < cols && p.z < channels);
        }
        
        inline bool contains(Node3D &p) {
            return (p.x >= 0 && p.y >= 0 && p.z >= 0 && p.x < rows && p.y < cols && p.z < channels);
        }

        template<typename _Tp2> inline bool is_border(Node3_<_Tp2> &p) {
            return (p.x == 0 || p.y == 0 || p.z == 0 || p.x == (rows - 1) || p.y == (cols - 1) || p.z == (channels - 1));
        }

        template<typename _Tp2> inline Node3_<_Tp2> get_closest_node(Node3_<_Tp2> &p) {
            Node3_<_Tp2> &p2 = p;
            if (p2.x < 0) p2.x = 0;
            if (p2.y < 0) p2.y = 0;
            if (p2.z < 0) p2.z = 0;
            if (p2.x >= cols) p2.x = cols - 1;
            if (p2.y >= rows) p2.y = rows - 1;
            if (p2.z >= channels) p2.z = channels - 1;
            return p2;
        }

    };

    typedef Grid3_<double> Grid3;

    template<typename _Tp, int N> struct wmm30_aux_ {
        _Tp v[3];
        Node3 dirs[2];
        _Tp m[N];
        _Tp fm[N];
    };

    template<typename _Tp, int N> struct wmm3_ {
        Node3 p;
        Node3 father;
        _Tp v;
        int ndirs;
        wmm30_aux_<_Tp, N> planes[4];
    };

    template<typename _Tp> struct wmm3_<_Tp, 0> {
        Node3 p;
        Node3 father;
        _Tp v;
        int ndirs;
        wmm30_aux_<_Tp, 1> planes[4];
    };

}
#endif	/* MPPSTRUCTS3D_H */

