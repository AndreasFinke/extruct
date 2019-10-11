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

class DisplacementField : public Measurement {

public:

    DisplacementField(int resol) : Measurement(4*2*resol), res(resol) {
        dataTypeSize = 4;
        rows = 2;
        cols = res;
        datap = (float *) data;
        reset();
    }
   
    virtual void reset() { 
        for (int i = 0; i < 2*res; ++i) {
           datap[i] = 0;
        }
    }

    virtual ~DisplacementField() {} 

    virtual void measure(const Universe& universe, int N) {
        Measurement::measure(universe, N);

        float A = 1.f / N; 


        for (int i = 0; i < res; ++i) {
            datap[i] += A*universe.initDisplacement.get_field(Float(i)/res);
            datap[res+i] += A*universe.initDisplacement.get_displacement(Float(i)/res);
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

    //friend class CorrelationFunctionObs; 
    //friend class PowerSpectrum3DObs;

public:

    PowerSpectrumObs(int method, int res, Float L, int skip) : Measurement(4*res*2), method(method), res(res), skip(skip), L(L) {
        dataTypeSize = 4;
        rows = 2;
        cols = res;
        datap = (float *) data;
        reset();

    }
    
    const Float AVG_BIN_WIDTH = 0.01;

    Float L; 

    virtual void reset() { 
        for (int i = 0; i < res; ++i) {
           datap[i+res] = 0;
        }
        //Float kmax = skip*2*pi*res;
        //Float base = std::pow(kmax, Float(1)/res);
        datap[0] = 2*pi;
        Float base = 1.1;
        Float growth = 1;
        int last = 0;
        while (last < res-1) { 
            growth *= base;
            Float  prop = 2*pi*int(growth);

            if (std::fabs(prop-datap[last]) > 1) {
                ++last;
                datap[last] = prop/(L*hor);
            }
        }

            //datap[i] = datap[i-1]*base;
    }


    virtual ~PowerSpectrumObs() {
        //delete densobs;
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
                    Float x = universe.get_particle_pos_standardized(j);
                    Float w = std::fabs(Float(x-1/2));
                    w = w*w;
                    w = 1-w;
                    w = 1;
                    //c += A*std::exp(1if * float(datap[i] * universe.get_particle_pos(j)));
                    rs += std::sin(datap[i] * x)*w;
                    rc += std::cos(datap[i] * x)*w;
                    //rs += std::sin(datap[i] * x);
                    //rc += std::cos(datap[i] * x);
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
                //double rs = 0, rc = 0;
                //
                Float centralk = datap[i]*L*hor;
                Float deltak = AVG_BIN_WIDTH*centralk;
                deltak = std::min(deltak, Float(2*pi*50));
                //deltak = 0;
                /* how many other nontrivial k's fit in a symmetric bin from centralk - deltak ... centrak + deltak ? */
                int width = int(deltak/(2*pi));
                int nPoints = 2*width+1;
                Float kstep = 2*pi;
                //if (nPoints > 1) 
                    //kstep = 2*deltak/(nPoints-1);
                //else 
                    //deltak = 0;

                //std::cout << "Deltak " << deltak<< " npoints " << nPoints << std::endl;
                Float k = centralk - width*2*pi; //deltak; 
                for (int m = 0; m < nPoints; ++m) {

                    //std::cout << "k = " << k << " ";
                    Sum rs, rc;
                    for (int j = 0; j < universe.nParticles; ++j) {
                        // FS coeffs are FT sampled at right k / Volume
                        // We need something in between for the power spectrum: FT ones / sqrt(Volume) 
                        // FT of line segment
                        // int_a^b dx exp(i x k_j)  = 1/ik (exp(iak) - exp(ibk)) = exp(i (a+b)/2 k)/ik (exp(i (a-b)k/2) - exp(-i(a-b)k/2)) 
                        // but fixed mass is distributed into invervall a..b so there is an additional 1/(b-a)
                        // but we want  density contrast, which up to constant is density/mean density = mass / total mass
                        // and total mass is nParticles. 
                        //
                        // so  -exp(i(a+b)k/2) sinc((b-a)k/2)) = (-cos((a+b)k/2) sinc((b-a)k/2), -sin((a+b)k/2) sinc((b-a)k/2))
                        // and division by nParticles. 
                        //
                        // we can normalize positions going to x/L coordinates, picking up a factor of L in front. But we said we want only a sqrt(L) factor here 
                        // to get the right power spectrum after squaring 
                        Float xL = universe.get_real_particle_pos_standardized(sortedById[j]);
                        //xL = universe.get_particle_pos_standardized(j);
                        /* wrap around at the end */
                        Float xR = ( j+1 < universe.nParticles ) ? universe.get_real_particle_pos_standardized(sortedById[j+1])
                            : (universe.get_real_particle_pos_standardized(sortedById[0]) + 1);
                        //xR = universe.get_particle_pos_standardized()
                        if (xL > xR)
                            std::swap(xL,xR);
                        Float phasearg = (xR+xL)*k*Float(0.5);
                        Float sincarg = (xR-xL)*k*Float(0.5);
                        rs += -std::sin(phasearg)*std::sin(sincarg)/sincarg; 
                        rc += -std::cos(phasearg)*std::sin(sincarg)/sincarg; 

                        // NEW: do a gaussian smoothing here instead. FT of gaussian from a to b with sigma b-a has the same shifting phase factor, 
                        // but sinc -> exp(-(k*(b-a))^2) 
                        //Float garg = sincarg*sincarg;
                        //rs += -std::sin(phasearg)*std::exp(-garg*10); 
                        //rc += -std::cos(phasearg)*std::exp(-garg*10); 
                    }
                    //datap[res+i] += std::abs(c);
                    /* note that we use usual Mpc units */
                    double rsd = rs/universe.nParticles;
                    double rcd = rc/universe.nParticles;
                    datap[res+i] += A*(L*hor)*(rsd*rsd+rcd*rcd)/nPoints;

                    k += kstep;

                }
            }
            /* we excluded k=0 above - but for density contrast average is zero, so don't add into zero mode */ 
            //A*universe.nParticles*universe.nParticles; 
        }

        if (method == 2) {

            //for (int i = 0; i < res; ++i ) 
                //std::cout<< datap[i]<<" ";
            //std::cout << std::endl;
            //std::cout << "large int is " << datap[res-1]*L*hor/(2*pi) << std::endl;
            Float exponent = 1 + int(std::log(datap[res-1]*L*hor/(2*pi))/std::log(2.0));
            //std::cout << "exponent is " << exponent << std::endl;
            int res_fourier = std::round(std::pow(2, exponent));
            const int res_max = 4096*128*4;
            //std::cout << 2*res_fourier  << std::endl;
            skip = std::max(1, res_fourier/res_max);
            //std::cout << skip  << std::endl;
            res_fourier /= skip;

            //std::cout << 2*res_fourier  << std::endl;
            densobs = new DensityObs2(2*res_fourier);
            densobs->set_zoom(skip);
            /* don't pass N: reset measurment after each and average later */
            densobs->measure(universe, 1); 
            fftw_complex *out;
            double *in;
            fftw_plan p;
            in = (double*) fftw_malloc(sizeof(double) * 2 * res_fourier);

            //std::cout << "alive " << std::endl;
            Float avg = 0;
            Float l = densobs->datap[2*res_fourier];
            Float r = densobs->datap[2*res_fourier+2*res_fourier-1];
            for (int i = 0; i < 2*res_fourier; ++i) {
                //Float e = 1 - Float(i)/(2*res_fourier-1);
                //densobs->datap[2*res_fourier+i] -= e*l+ (1-e)*r;
                avg += densobs->datap[2*res_fourier+i];
            }
            //std::cout << "alive 2" << std::endl;
            avg /= 2*res_fourier;
            for (int i = 0; i < 2*res_fourier; ++i) {
                //in[i] = (densobs->getData())[i];
                Float w = std::fabs(Float(i-(res_fourier-1/2))/Float(res_fourier-1/2));
                w = w*w;
                w = 1 - w;
                w = 1;
                in[i] = w*(densobs->datap[2*res_fourier+i]-avg)/universe.nParticles;
            }

            //std::cout << "alive 3" << std::endl;
            delete densobs;

            out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * 2*res_fourier);

        //p = fftw_plan_dft_1d(nGrid, reinterpret_cast<fftw_complex*>(&modes[0]), reinterpret_cast<fftw_complex*>(&field[0]), FFTW_BACKWARD, FFTW_ESTIMATE);
            p = fftw_plan_dft_r2c_1d(2 * res_fourier, in, out, FFTW_ESTIMATE);
            fftw_execute(p); 

            //std::cout << "res_fourier " << res_fourier << std::endl;

            for (int i = 1; i < res; ++i) {
                Float centralk = datap[i]*L*hor;
                int centralid = std::round(centralk/(2*pi)/skip);
                Float deltak = AVG_BIN_WIDTH*centralk;
                //deltak = std::min(deltak, Float(2*pi*50));
                int width = int(deltak/(2*pi)/skip);
                int nPoints = 2*width+1;

                //std::cout << " i = " << i << "   ";
                for (int m = 0; m < nPoints; ++m) {
                    //std::cout << " m = " << m;
                    int idx = centralid - width + m;
                    idx = std::min(idx, res_fourier-1);
                    idx = std::max(idx, 0);
                    datap[i+res] += A*L*hor*(out[idx][0]*out[idx][0] + out[idx][1]*out[idx][1])/nPoints;

                }
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

//class CorrelationFunctionObs : public Measurement {

    //friend class PowerSpectrum3DObs;

//public:

    //CorrelationFunctionObs(int method, int res) : Measurement(4*res*2), res(res) {
        //psobs = new PowerSpectrumObs(method, res+1, 1); 
        //dataTypeSize = 4;
        //rows = 2;
        //cols = res;
        //datap = (float *) data;
        //reset();
    //}
    
    //virtual void reset() { 
        //for (int i = 0; i < res; ++i) {
           //datap[i] = 2*pi*i; 
           //datap[i+res] = 0;
        //}
    //}

    //virtual ~CorrelationFunctionObs() {
        //delete psobs; 
    //} 

    //virtual void measure(const Universe& universe, int N) {
        //Measurement::measure(universe, N);

        //float A = 1.f / N; 

        //psobs->measure(universe, 1);

        //double *out;
        //double *in;
        //fftw_plan p;
        //in = (double*) fftw_malloc(sizeof(double) * 1 * (res+1));
        
        //for (int i = 0; i < 1*res + 1; ++i) {
            //in[i] = psobs->datap[1*res+1+i];
        //}

        //[> remove delta at k=0, that is, constant offset of correlation fctn <]
        //in[0] = 0;

        //out = (double*) fftw_malloc(sizeof(double) * 1 * (res+1));

    ////p = fftw_plan_dft_1d(nGrid, reinterpret_cast<fftw_complex*>(&modes[0]), reinterpret_cast<fftw_complex*>(&field[0]), FFTW_BACKWARD, FFTW_ESTIMATE);
    ////
        //[> fast if N = 2(n-1) where n is the first argument has small factors! that is, res should be a power of 2!<]
        //p = fftw_plan_r2r_1d(1*res + 1, in, out, FFTW_REDFT00, FFTW_ESTIMATE); //can try 01 too (odd at right boundary... what does it mean? 
        ////p = fftw_plan_r2r_1d(1*res, in, out, FFTW_REDFT00);
        //fftw_execute(p); 

        //[> drop last <]

        //for (int i = 0; i < res; ++i) {
            //datap[i+res] += A*out[i];
        //}

        //fftw_destroy_plan(p);
        //fftw_free(in); 
        //fftw_free(out);
    //}
//private:

    //int res;

    //PowerSpectrumObs * psobs; 

    //float* datap; 


//};

//class PowerSpectrum3DObs : public Measurement {


//public:

    //PowerSpectrum3DObs(int method, int res) : Measurement(4*res*2), res(res), method(method) {
        //corrobs = new CorrelationFunctionObs(2, res);
        //psobs = new PowerSpectrumObs(2, res, 1);
        //dataTypeSize = 4;
        //rows = 2;
        //cols = res;
        //datap = (float *) data;
        //reset();
    //}
    
    //virtual void reset() { 
        //for (int i = 0; i < res; ++i) {
           //datap[i] = 2*pi*i; 
           //datap[i+res] = 0;
        //}
    //}

    //virtual ~PowerSpectrum3DObs() {
        //delete corrobs; 
        //delete psobs;
    //} 

    //virtual void measure(const Universe& universe, int N) {
        //Measurement::measure(universe, N);

        //float A = 1.f / N; 

        //if (method == 0) {

            //psobs->measure(universe, 1);

            //for (int i = 0; i < res-1; ++i) {
                //datap[i+res+1] += -A*(psobs->datap[res+i+1]-psobs->datap[res+i])/(i+1);
            //}
            //datap[res] = 0;
        //}
        //else if (method == 1) { 
            //corrobs->measure(universe, 1);

            //double *out;
            //double *in;
            //fftw_plan p;
            //in = (double*) fftw_malloc(sizeof(double) * (res-1));
           
            ////[> drop the for the odd function - note we add +1 in the indices and have one term less<] 
            //for (int i = 0; i < res-1; ++i) {
                //in[i] = corrobs->datap[res+i+1]*corrobs->datap[i+1]*(universe.L/res);
            //}

            //out = (double*) fftw_malloc(sizeof(double) * (res-1));

        ////p = fftw_plan_dft_1d(nGrid, reinterpret_cast<fftw_complex*>(&modes[0]), reinterpret_cast<fftw_complex*>(&field[0]), FFTW_BACKWARD, FFTW_ESTIMATE);
            ////p = fftw_plan_dft_r2c_1d(2*res, in, out, FFTW_ESTIMATE);
            ////[> symmetric about index -1 (which corresponds to the dropped first index) indeed! <]
            ////[> note this is fast if N=2(n+1) has small factors, so res should again be a power of 2! <]
            //p = fftw_plan_r2r_1d(1*res-1, in, out, FFTW_RODFT00, FFTW_ESTIMATE);
            //fftw_execute(p); 

            //for (int i = 0; i < res-1; ++i) {
                //datap[i+res+1] += A*8*pi*out[i]/datap[i];
            //}
            ////[> cannot obtain this one - zero (from oddness ) / zero... <]
            //datap[res] = 0;

            //fftw_destroy_plan(p);
            //fftw_free(in); 
            //fftw_free(out);
        //}
    //}
//private:

    //int method;
    //int res;

    //CorrelationFunctionObs * corrobs; 
    //PowerSpectrumObs * psobs; 

    //float* datap; 


//};
