#pragma once

using Float = double;
using Long = long;

#if PY == 1

#include <pybind11/numpy.h>
namespace py = pybind11;

#endif


#include <cmath>
#include "pcg32.h"

#include <complex>

static constexpr auto pi = 3.141592653589793238462643383279502884197169399375105820974944592307816;
using namespace std::complex_literals;

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

struct Entry {
    Entry(Long idx, Float t) : t(t), idx(idx) {} 
    Float t;
    Long idx;
    bool operator<(const Entry& rhs) const {return t < rhs.t;}
    friend std::ostream& operator<<(std::ostream& out, const Entry& obj);
};

template<class CollisionTask> 
struct TaskParticle {
    TaskParticle() : x(0), v(0), t(0), nCollisions(0), index(0), sheet(0) {} 
    TaskParticle(Float x, Float v, int index, short sheet) : x(x), v(v), nCollisions(0), index(index), sheet(sheet){}
    Float x, v, t;
    Long nCollisions;
    int index;
    short sheet;
    CollisionTask task;
};

class PowerSpectrum {
public:
    virtual ~PowerSpectrum() {}
    virtual Float eval(Float k) const  = 0;
};

class PowerLaw : public PowerSpectrum {
public: 
    Float A = 0.1;
    virtual ~PowerLaw() {}
    virtual Float eval(Float k) const { 
        //return A/(0.1*k+1); 
        return A*k*k/(100+0*k*k*k*k);
    }
};

