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
            idx = std::max(0, idx);
            idx = std::min(res-1, idx);

            datap[res+idx] += A;
        }

    }
private:

    int res;
    float* datap; 


};

class CollisionObs : public Measurement {

public:

    CollisionObs(int res) : Measurement(4*2*res), res(res) {
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

    virtual ~CollisionObs() {} 

    virtual void measure(const Universe& universe, int N) {
        Measurement::measure(universe, N);


        for (int i = 0; i < universe.nParticles; ++i) {
            Float x = universe.get_particle_pos_standardized(i);
            int idx = x * res;
            idx = std::max(0, idx);
            idx = std::min(res-1, idx);
            
            datap[res+idx] = std::max(float(universe.get_particle_collision_number(i)), datap[res+idx]);
        }

    }
private:

    int res;
    float* datap; 


};

class DensityObs2 : public Measurement {

    friend class PowerSpectrumObs;

public:

    DensityObs2(int res) : Measurement(4*2*res), res(res) {
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

    virtual ~DensityObs2() {} 

    virtual void measure(const Universe& universe, int N) {
        Measurement::measure(universe, N);

        float A = 1.f / N; 

        /* will hold indices of particles sorted by their id */
        std::vector<Long> sortedById(universe.nParticles);
        /* fill with 0, 1, ... nParticles-1 */  
        std::iota(std::begin(sortedById), std::end(sortedById), 0);

        //std::sort(ps.begin(), ps.end(), [](auto& lhs, auto& rhs) { return lhs.index < rhs.index;} );
        std::sort(sortedById.begin(), sortedById.end(), 
                [&](const Long& a, const Long& b) { 
                    return universe.get_particle_id(a) < universe.get_particle_id(b);
                }  
        );

        for (int i = 0; i < universe.nParticles; ++i) {

            Float xL = universe.get_real_particle_pos_standardized(sortedById[i]);
            /* wrap around at the end */
            Float xR = ( i+1 < universe.nParticles ) ? universe.get_real_particle_pos_standardized(sortedById[i+1])
                : (universe.get_real_particle_pos_standardized(sortedById[0]) + 1);
            xL *= zoom;
            xR *= zoom;
            //Float xL = ps[i].x/universe.L + 0.5 + ps[i].sheet;
            //Float xR = ps[(i+1)%nParticles].x/universe.L + 0.5 + ps[i].sheet;

            if (xL > xR) 
                std::swap(xL, xR);

            int idxL = std::floor(xL * res);
            int idxR = std::floor(xR * res);

            Float fracL = 1- (res*xL-idxL);
            Float fracR = 1- (res*xR-idxR);

            //if (i == universe.nParticles-1)
            //{
                //fracL = 1-fracL;
                //fracR = 1-fracR;
            //}

            //if (idxL > res - 1) idxL = res - 1;
            //if (idxR > res - 1) idxR = res - 1;
            //if (idxL < 0) idxR = 0;
            //if (idxL < 0) idxR = 0;


            /* this function always returns a positive number, unlike % */
            auto mod = [res=res](Long n) {return (n%res+res)%res;};

            if (idxL == idxR) { 
                //Float A2 = A*0.5;
                //
                datap[res+mod(idxL)] += A;
                //datap[res+idxL] += A2 * ( fracL + fracR );
                //datap[res+idxL+1] += A2 * ( 2  - fracL - fracR );
                //
                //std::cout << " added " << A << " to bin " << res + idxL << std::endl;

            }
            else if (idxR == idxL + 1) {

                datap[res+mod(idxL)] += A * fracL/(fracL + 1 - fracR);
                datap[res+mod(idxR)] += A * (1-fracR)/(fracL + 1 - fracR);
                //datap[res+idxL] += A2 * fracL;
                //datap[res+idxR+1] += A2 * ( 1 - fracR );
                //datap[res+idxR] += A - A2*fracL - A2*(1-fracR)
                //std::cout << " added " <<  A * fracL/(fracL + 1 - fracR) << " and " << (1-fracR)/(fracL + 1 - fracR) << " to bin " << res + idxL << " and " << res + idxR << std::endl;
            }
            else {

                int d = idxR - idxL;

                datap[res+mod(idxL)] += A * fracL/(fracL + d - fracR);
                datap[res+mod(idxR)] += A * (1-fracR)/(fracL + d - fracR);

                //std::cout << " added " <<  A * fracL/(fracL + d - fracR) << " and " << (1-fracR)/(fracL + d - fracR) << " to bin " << res + idxL << " and " << res + idxR;

                for (int k = idxL+1; k < idxR; ++k) { 
                        datap[res+mod(k)] += A/(fracL + d - fracR);

                        //std::cout << " and also " << A/(fracL + d - fracR) << " to " << res + k;
                }
                //std::cout << std::endl;

            }
            if  (fracL/(fracL + 1 - fracR) < 0 || (1-fracR)/(fracL + 1 - fracR) < 0 ) 
                std::cout << "fracs are " << fracL << " " << fracR << std::endl;

        }

    }

    void set_zoom(float z) { zoom = z; }
private:

    float zoom = 1;
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

    friend class CorrelationFunctionObs; 

public:

    PowerSpectrumObs(int method, int res, int skip) : Measurement(4*res*2), method(method), res(res), skip(skip){
        dataTypeSize = 4;
        rows = 2;
        cols = res;
        datap = (float *) data;
        reset();

        densobs = new DensityObs2(2*res);
    }
    
    virtual void reset() { 
        for (int i = 0; i < res; ++i) {
           datap[i] = 2*pi*i*skip; 
           datap[i+res] = 0;
        }
    }


    virtual ~PowerSpectrumObs() {
        delete densobs;
    } 

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
                    rs += std::sin(datap[i] * universe.get_particle_pos_standardized(j));
                    rc += std::cos(datap[i] * universe.get_particle_pos_standardized(j));
                }
                //datap[res+i] += std::abs(c);
                datap[res+i] += A*(rs*rs+rc*rc);
            }
        }

        if (method == 1) {
            /* will hold indices of particles sorted by their id */
            std::vector<Long> sortedById(universe.nParticles);
            /* fill with 0, 1, ... nParticles-1 */  
            std::iota(std::begin(sortedById), std::end(sortedById), 0);

            //std::sort(ps.begin(), ps.end(), [](auto& lhs, auto& rhs) { return lhs.index < rhs.index;} );
            std::sort(sortedById.begin(), sortedById.end(), 
                    [&](const Long& a, const Long& b) { 
                        return universe.get_particle_id(a) < universe.get_particle_id(b);
                    }  
            );

            for (int i = 1; i < res; ++i) {
                //std::complex<float> c = 0;
                double rs = 0, rc = 0;
                for (int j = 0; j < universe.nParticles; ++j) {
                    // FT of line segment
                    // int_a^b dx exp(i x k_j)  = 1/ik (exp(iak) - exp(ibk)) = exp(i (a+b)/2 k)/ik (exp(i (a-b)k/2) - exp(-i(a-b)k/2)) 
                    // but for density, there is another 1/(b-a) 
                    // so  -exp(i(a+b)k/2) sinc((b-a)k/2)) = (-cos((a+b)k/2) sinc((b-a)k/2), -sin((a+b)k/2) sinc((b-a)k/2))
                    Float xL = universe.get_real_particle_pos_standardized(sortedById[j]);
                    //xL = universe.get_particle_pos_standardized(j);
                    /* wrap around at the end */
                    Float xR = ( j+1 < universe.nParticles ) ? universe.get_real_particle_pos_standardized(sortedById[j+1])
                        : (universe.get_real_particle_pos_standardized(sortedById[0]) + 1);
                    //xR = universe.get_particle_pos_standardized()
                    if (xL > xR)
                        std::swap(xL,xR);
                    Float sincarg = (xR-xL)*datap[i]*0.5;
                    Float exparg = (xR+xL)*datap[i]*0.5;
                    rs += -std::sin(exparg)*std::sin(sincarg)/sincarg; 
                    rc += -std::cos(exparg)*std::sin(sincarg)/sincarg; 
                }
                //datap[res+i] += std::abs(c);
                datap[res+i] += A*(rs*rs+rc*rc);
            }
            /* we excluded k=0 above - this is the exact result */ 
            datap[res] += A*universe.nParticles*universe.nParticles; 
        }

        if (method == 2) {
            /* don't pass N: reset measurment after each and average later */
            densobs->set_zoom(skip);
            densobs->measure(universe, 1); 
            fftw_complex *out;
            double *in;
            fftw_plan p;
            in = (double*) fftw_malloc(sizeof(double) * 2 * res);
            
            for (int i = 0; i < 2*res; ++i) {
                //in[i] = (densobs->getData())[i];
                in[i] = densobs->datap[2*res+i];
            }


            out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * 2 * res);

        //p = fftw_plan_dft_1d(nGrid, reinterpret_cast<fftw_complex*>(&modes[0]), reinterpret_cast<fftw_complex*>(&field[0]), FFTW_BACKWARD, FFTW_ESTIMATE);
            p = fftw_plan_dft_r2c_1d(2 * res, in, out, FFTW_ESTIMATE);
            fftw_execute(p); 

            for (int i = 0; i < res; ++i) {
                datap[i+res] += A*(out[i][0]*out[i][0] + out[i][1]*out[i][1]);
            }

            fftw_destroy_plan(p);
            fftw_free(in); 
            fftw_free(out);
        }


    }
private:

    int method;
    int res;
    int skip;

    float* datap; 

    DensityObs2 * densobs;


};

class CorrelationFunctionObs : public Measurement {


public:

    CorrelationFunctionObs(int method, int res) : Measurement(4*res*2), method(method), res(res) {
        psobs = new PowerSpectrumObs(1, res, 1); 
        dataTypeSize = 4;
        rows = 2;
        cols = res;
        datap = (float *) data;
        reset();
    }
    
    virtual void reset() { 
        for (int i = 0; i < res; ++i) {
           datap[i] = 2*pi*i; 
           datap[i+res] = 0;
        }
    }

    virtual ~CorrelationFunctionObs() {
        delete psobs; 
    } 

    virtual void measure(const Universe& universe, int N) {
        Measurement::measure(universe, N);

        float A = 1.f / N; 

        psobs->measure(universe, 1);

        double *out;
        double *in;
        fftw_plan p;
        in = (double*) fftw_malloc(sizeof(double) * 1 * res);
        
        for (int i = 0; i < 1*res; ++i) {
            in[i] = psobs->datap[1*res+i];
        }

        out = (double*) fftw_malloc(sizeof(double) * 1 * res);

    //p = fftw_plan_dft_1d(nGrid, reinterpret_cast<fftw_complex*>(&modes[0]), reinterpret_cast<fftw_complex*>(&field[0]), FFTW_BACKWARD, FFTW_ESTIMATE);
        p = fftw_plan_r2r_1d(1*res, in, out, FFTW_REDFT00, FFTW_ESTIMATE);
        //p = fftw_plan_r2r_1d(1*res, in, out, FFTW_REDFT00);
        fftw_execute(p); 

        for (int i = 0; i < res; ++i) {
            datap[i+res] += A*out[i];
        }

        fftw_destroy_plan(p);
        fftw_free(in); 
        fftw_free(out);
    }
private:

    int method;
    int res;

    PowerSpectrumObs * psobs; 

    float* datap; 


};
