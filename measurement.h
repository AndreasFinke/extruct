#pragma once

#include <vector>
#include <iostream>
#include <cstring>

#if PY == 1
#include <pybind11/numpy.h>
#endif

class Measurement {
public:
    Measurement(int bytes) : bytes(bytes) {
        data = new char[bytes];
    }

    virtual ~Measurement() { 
        delete[] data;
    }

    void reset() {
        std::memset(data, 0, bytes);
    }

    virtual void measure(const Universe& universe, int N) = 0;

    auto getResult() {
#if PY == 1 
        py::array_t<float> ret({rows,cols});
        std::memcpy(data, ret.mutable_data(), bytes) ;
        return ret;
#else
        return data;
#endif
    }

    void saveResult() {

    }

protected:

    int nMax = 0;
    int nCurr = 0;
    //void* data = nullptr;
    char* data = nullptr;
    int bytes;
    int dataTypeSize = 4;
    int rows = 0;
    int cols = 0;

};

#include <complex>
using namespace std::complex_literals;

class PowerSpectrumObs : public Measurement {


public:

    PowerSpectrumObs(int method, int res, Float kmin, Float kmax) : Measurement(4*res*2), method(method), res(res) {
        dataTypeSize = 4;
        rows = 2;
        cols = res;
        reset();
        datap = (float *) data;
        for (int i = 0; i < res; ++i) {
           datap[i] = kmin + (kmax-kmin)*i/Float(res); 
        }
    }
    
    virtual ~PowerSpectrumObs() {} 

    virtual void measure(const Universe& universe, int N) {
        if (nCurr == 0)
            nMax = N;
        if (++nCurr > nMax)  { 
            std::cout << "Resetted measurements after more than the requested " << nMax << " calls. Average begins afresh." << std::endl; 
            nCurr = 1; 
            nMax = N;
        }

        float A = 1.f / N; 

        if (method == 0) 
        {
            for (int i = 0; i < res; ++i) {
                std::complex<float> c = 0;
                for (int j = 0; j < universe.nParticles; ++j) {
                    c += A*std::exp(1if * float(datap[i] * universe.get_particle_pos(j)));
                    std::cout << c; 
                }
                datap[res+i] += std::abs(c);
            }
        }


    }
private:

    int method;
    int res;

    float* datap; 


};
