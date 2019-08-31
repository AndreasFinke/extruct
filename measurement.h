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
            //std::cout << "Resetting measurements after more than the requested " << nMax << " calls. Average will begin afresh." << std::endl; 
            nCurr = 1; 
            nMax = N;
            reset();
        }
    }

    auto getResult() {
#if PY == 1 
        py::array_t<float> ret({rows,cols});
        float * ret_ptr = ret.mutable_data();
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
            Float x = universe.get_particle_pos_standardized(i);
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

class PhaseSpaceDensityObs : public Measurement {


public:

    PhaseSpaceDensityObs(int res) : Measurement(2*4*res+4*res*res), res(res) {
        dataTypeSize = 4;
        rows = 2+res;
        cols = res;
        datap = (float *) data;
        reset();
    }
   
    virtual void reset() { 
        for (int i = 0; i < res; ++i) {
           datap[i] = 0 +  1*i/Float(res); 
           datap[res+i] = 0;
           for (int j = 0; j < res; ++j) 
               datap[2*res + i*res + j] = 0;
        }
    }

    virtual ~PhaseSpaceDensityObs() {} 

    virtual void measure(const Universe& universe, int N) {
        Measurement::measure(universe, N);

        float A = 1.f / N; 

        float minV = 0;
        float maxV = 0;

        for (int i = 0; i < universe.nParticles; ++i) {
            if (universe.get_particle_vel(i) > maxV)
                maxV = universe.get_particle_vel(i);
            if (universe.get_particle_vel(i) < minV)
                minV = universe.get_particle_vel(i);
        }
        for (int i = 0; i < res; ++i)
            datap[res+i] = minV + float(i)/(res-1)*(maxV-minV);
        for (int i = 0; i < universe.nParticles; ++i) {
            Float x = universe.get_particle_pos_standardized(i);
            Float v = universe.get_particle_vel(i);
            int idx = x * res;
            int idx_v = (v-minV)/(maxV-minV) * res;
            if (idx > res - 1) idx = res - 1;
            if (idx < 0) idx = 0;
            if (idx_v > res - 1) idx_v = res - 1;
            if (idx_v < 0) idx_v = 0;
            datap[2*res + idx*res + idx_v] += 1;
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

    PowerSpectrumObs(int method, int res, int skip) : Measurement(4*res*2), method(method), res(res), skip(skip){
        dataTypeSize = 4;
        rows = 2;
        cols = res;
        datap = (float *) data;
        reset();
    }
    
    virtual void reset() { 
        for (int i = 0; i < res; ++i) {
           datap[i] = 2*pi*i*skip; 
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
                //std::complex<float> c = 0;
                double rs = 0, rc = 0;
                for (int j = 0; j < universe.nParticles; ++j) {
                    //c += A*std::exp(1if * float(datap[i] * universe.get_particle_pos(j)));
                    rs += std::sin(datap[i] * universe.get_particle_pos(j));
                    rc += std::cos(datap[i] * universe.get_particle_pos(j));
                }
                //datap[res+i] += std::abs(c);
                datap[res+i] += A*std::sqrt(rs*rs+rc*rc);
            }
        }


    }
private:

    int method;
    int res;
    int skip;

    float* datap; 


};
