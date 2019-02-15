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
        reset();
    }

    virtual ~Measurement() { 
        delete[] data;
    }

    virtual void reset() {
        std::memset(data, 0, bytes);
    }

    virtual void measure(const Universe& universe, int N) {
        if (nCurr == 0)
            nMax = N;
        if (++nCurr > nMax)  { 
            std::cout << "Resetted measurements after more than the requested " << nMax << " calls. Average begins afresh." << std::endl; 
            nCurr = 1; 
            nMax = N;
            reset();
        }
    }

    auto getResult() {
#if PY == 1 
        py::array_t<float> ret({rows,cols});
        float * ret_ptr = ret.mutable_data();
        //std::cout << std::endl << std::endl << "Density = ";
        //for (int i = 0; i <100; ++i) {
            //std::cout << ((float*)data)[i] << " " << ((float*)data)[100+i] << std::endl;
        //}
        //std::cout <<std::endl << std::endl;
        std::memcpy(ret_ptr, (float*) data, bytes) ;
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

class DensityObs : public Measurement {


public:

    DensityObs(int res) : Measurement(4*2*res), res(res) {
        dataTypeSize = 4;
        rows = 2;
        cols = res;
        datap = (float *) data;
        reset();
    }
   
    virtual void reset() { 
        for (int i = 0; i < res; ++i) {
           datap[i] = 0 +  1*i/Float(res); 
           datap[i+res] = 0;
        }
    }

    virtual ~DensityObs() {} 

    virtual void measure(const Universe& universe, int N) {
        Measurement::measure(universe, N);

        float A = 1.f / N; 

        for (int i = 0; i < universe.nParticles; ++i) {
            Float x = universe.get_particle_pos(i);
            int idx = x * res;
            if (idx > res - 1) idx = res - 1;
            if (idx < 0) idx = 0;
            datap[res+idx] += A;
        }

    }
private:

    int res;
    float* datap; 


};
#include <complex>
using namespace std::complex_literals;

class PowerSpectrumObs : public Measurement {


public:

    PowerSpectrumObs(int method, int res, Float kmin, Float kmax) : Measurement(4*res*2), method(method), res(res), kmin(kmin), kmax(kmax) {
        dataTypeSize = 4;
        rows = 2;
        cols = res;
        datap = (float *) data;
        reset();
    }
    
    virtual void reset() { 
        for (int i = 0; i < res; ++i) {
           datap[i] = kmin + (kmax-kmin)*i/Float(res); 
           datap[i+res] = 0;
        }
    }

    virtual ~PowerSpectrumObs() {} 

    virtual void measure(const Universe& universe, int N) {
        Measurement::measure(universe, N);

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

    float kmin;
    float kmax;
    float* datap; 


};
