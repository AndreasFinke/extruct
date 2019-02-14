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

template<class CollisionTask> 
struct TaskParticle {
    TaskParticle() : x(0), v(0), t(0) {} 
    TaskParticle(Float x, Float v) : x(x), v(v) {}
    Float x, v, t;
    Float collTime = -1;
    CollisionTask task;
};

class PowerSpectrum {
public:
    virtual ~PowerSpectrum() {}
    virtual Float eval(Float k) const  {}
};

class PowerLaw : public PowerSpectrum {
public: 
    virtual ~PowerLaw() {}
    virtual Float eval(Float k) const { 
        return 1/(0.1*k+1); 
    }
};

