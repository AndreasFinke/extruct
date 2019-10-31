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

    virtual void postprocess() { 
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

template<class T> class CorrelationFunction;
//class DensityObs2;


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
            datap[i] += A*universe.initDisplacement.get_field_at(Float(i)/res);
            datap[res+i] += A*universe.initDisplacement.get_displacement_at(Float(i)/res);
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
    friend class CorrelationFunction<DensityObs2>;

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
           datap[i] = i/Float(res); 
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


class DensityLin : public Measurement {

    friend class CorrelationFunction<DensityLin>;

public:

    DensityLin(int res) : Measurement(4*2*res), res(res) {
        dataTypeSize = 4;
        rows = 2;
        cols = res;
        datap = (float *) data;
        reset();
    }
   
    virtual void reset() { 
        for (int i = 0; i < res; ++i) {
           datap[i] = i/Float(res); 
           datap[i+res] = 0;
        }
    }

    virtual ~DensityLin() {} 

    virtual void measure(const Universe& universe, int N) {
        Measurement::measure(universe, N);

        float A = 1.f / N; 

        /* the field actually gives the density contrast modes. density = meandensity*(delta+1)
         * where meandensity is mostly arbitrary (to go back to delta / power etc) but we align it here with the convention for the other observables
         * which is total mass of nParticles */
        Float meandens = Float(universe.nParticles)/res;
        std::vector<Float> delta = universe.initDisplacement.get_field(res);
        Float g = universe.bg.getGrowth(universe.most_recent_particle_time());
        for (int i = 0; i < res; ++i) {
            datap[res+i] += A*(g*delta[i]+1)*meandens;
        }

    }

private:

    int res;
    float* datap; 


};

template <int nLoops>
class DensitySPT : public Measurement {

    friend class CorrelationFunction<DensitySPT<nLoops>>;

public:

    DensitySPT(int res) : Measurement(4*2*res), res(res) {
        dataTypeSize = 4;
        rows = 2;
        cols = res;
        datap = (float *) data;
        reset();
    }
   
    virtual void reset() { 
        for (int i = 0; i < res; ++i) {
           datap[i] = i/Float(res); 
           datap[i+res] = 0;
        }
    }

    virtual ~DensitySPT() {} 

    virtual void measure(const Universe& universe, int N) {
        Measurement::measure(universe, N);

        float A = 1.f / N; 

        /* work on grids of double resolution, do multiplications one at a time */
        /* can change to res ... more than 2*res should not change things */
        //int calcres = 2*res;
        int calcres = 4*res;

        /* this is fixed */
        int fourierres = calcres/2+1;
        int deltakres = res/2+1;

        std::vector<Float> psi = universe.initDisplacement.get_displacement(calcres, res);

        fftw_complex *out, *deltak;
        //char *convergedk;
        double *in, *delta;
        fftw_plan p1, p2, pdelta;

        in = (double*) fftw_malloc(sizeof(double) * calcres);
        out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fourierres);
        delta = (double*) fftw_malloc(sizeof(double) * res);
        deltak = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * deltakres);
        //convergedk = (char*) fftw_malloc(sizeof(char) * deltakres);

        for (int i = 0; i < deltakres; ++i) {
            deltak[i][0] = 0;
            deltak[i][1] = 0;
        }
        for (int i = 0; i < calcres; ++i) {
            in[i] = 1;
        }

             //FFT normalization 
            //for (int i = 0; i < calcres; ++i) {
                //in[i] /= calcres;
            //}

        p1 = fftw_plan_dft_r2c_1d(calcres, in, out, FFTW_ESTIMATE);
        p2 = fftw_plan_dft_c2r_1d(calcres, out, in, FFTW_ESTIMATE);
        pdelta = fftw_plan_dft_c2r_1d(res, deltak, delta, FFTW_ESTIMATE);

        Float g = universe.bg.getGrowth(universe.most_recent_particle_time());

        /* determine order of magnitude of psi factor */
        double avg = 0;
        for (int i = 0; i < calcres; ++i) {
            avg += std::fabs(psi[i]);
        }
        double typicalValue = avg*g/calcres/10; 
        int order = std::log(typicalValue)/std::log(10);
        double orderVal = std::pow(10, order);

        for (int loop = 1; loop <= nLoops; ++loop) {

            /* initialized "in" with 1. 
             * In a loop:
             * multiply "in" by psi 
             * do FFT, project (remove higher part of the modes; initially redundant since there are already only ~deltakres modes populated )
             * add result to deltak (with k factors)
             * transform back
             * repeat  ... 
             * one final fft is needed to get delta to real space
             */


            for (int i = 0; i < calcres; ++i) {
                in[i] *= g*psi[i]/loop/calcres / orderVal;
            }


            fftw_execute(p1);

            //std::cout << avg << " or with growth " << g*avg << " out is " << std::sqrt(out[1][0]*out[1][0] + out[1][1]*out[1][1]) << " " << std::sqrt
                //(out[deltakres/2][1]*out[deltakres/2][1]+out[deltakres/2][0]*out[deltakres/2][0]) << " orderval " << orderVal<<  "\n";
            /* project */
            double p = 0;
            for (int i = deltakres; i < fourierres; ++i) {
                p += out[i][0]*out[i][0] + out[i][1]*out[i][1];
                out[i][0] = 0;
                out[i][1] = 0;
            }
            std::cout << "killed aliasing power of " << p << " at loop " << loop << "\n";

            int nNotConverged = 0; 

            for (int i = 0; i < deltakres; ++i) {
                double k = 2*pi*i; 
                double kn = k*orderVal;

                // kn = k^loop 
                for (int j = 1; j < loop; ++j) 
                    kn*=k*orderVal;

                if ( (loop % 2) == 0) { 
                    int sign = 1;
                    if ( (loop % 4) != 0) sign = -1;
                    deltak[i][0] += sign*kn*out[i][0];
                    deltak[i][1] += sign*kn*out[i][1];
                }
                else { 
                    /* multiply by i (loop = 3,7,...) , or -i (loop = 1,5,... ) */
                    int sign = 1;
                    if ( ((loop+1) % 4) != 0) sign = -1;
                    deltak[i][0] -= sign*kn*out[i][1];
                    deltak[i][1] += sign*kn*out[i][0];
                }

                /* estimate convergence by size of last term wrt to sum */
                if (loop == nLoops && nLoops > 3) {
                    if (std::fabs(kn*out[i][1]) + std::fabs(kn*out[i][0]) > 1e-1*(std::fabs(deltak[i][0]) + std::fabs(deltak[i][1])) ) {
                        deltak[i][0] = 0;
                        deltak[i][1] = 0;
                        nNotConverged++;
                    }
                }
            }

            if (nNotConverged > 0) 
                std::cout << "SPT not converged to 10 percent for " << nNotConverged << " requested Fourier modes at loop " << nLoops << ", which were zeroed. (This check is only applied at orders higher than 3.)\n";

            /* note: c2r destroys input, i.e. out */
            fftw_execute(p2);



        }

        /* fourier to real space (normalize below) */
        fftw_execute(pdelta);

        /* from delta to density; note FFT normalization */
        for (int i = 0; i < res; ++i) {
            Float meandens = Float(universe.nParticles)/res;
            datap[res+i] += A*(delta[i]+1)*meandens;
        }


        fftw_free(in);
        fftw_free(out);
        fftw_free(deltak);
        fftw_free(delta);
        fftw_destroy_plan(p1);
        fftw_destroy_plan(p2);
        fftw_destroy_plan(pdelta);

    }

private:

    int res;
    float* datap; 


};
using CorrelationFunctionObs = CorrelationFunction<DensityObs2>;
using SPTCorrelationFunctionObs = CorrelationFunction<DensitySPT<100>>;
using LinCorrelationFunctionObs = CorrelationFunction<DensityLin>;

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

        if (method == 2) { 
            Float exponent = 1 + int(std::log(datap[res-1]*L*hor/(2*pi))/std::log(2.0));
            /* FFTW r2c algorithm is fastest for small factors in the size of the real input array - use power of 2 here */
            res_fourier = std::round(std::pow(2, exponent));
            //std::cout << 2*res_fourier  << std::endl;

            /* dividing two int powers of 2 gives int power of two */
            skip = std::max(1, res_fourier/res_max);
            /* so this is int power of two */
            res_fourier /= skip;

            densobs = new DensityObs2(2*res_fourier);
            //densobs = new DensityObs2(res_fourier);
            densobs->set_zoom(skip);

            std::cout << "res fourier : " << res_fourier << "\n";
            if (skip > 1)
                std::cout << "Warning: requested maximal Fourier mode is so large that modes are skipped after zooming " << skip << " times into density." << std::endl;
        }
    }
    
    const Float AVG_BIN_WIDTH = 0.01;
    const int res_max = 4096*128*16*4;

    int res_fourier; 

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

            if (std::fabs(prop-datap[last]*L*hor) > 1) {
                ++last;
                datap[last] = prop/(L*hor);
            }
        }

            //datap[i] = datap[i-1]*base;
    }


    virtual ~PowerSpectrumObs() {
        if (method == 2)
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

            /* don't pass N: reset measurment after each and average later */
            densobs->measure(universe, 1); 
            fftw_complex *out;
            //double *out;
            double *in;
            fftw_plan p;
            in = (double*) fftw_malloc(sizeof(double) * 2 * res_fourier);
            //in = (double*) fftw_malloc(sizeof(double) * res_fourier);

            Float avg = 0;
            ////Float l = densobs->datap[2*res_fourier];
            ////Float r = densobs->datap[2*res_fourier+2*res_fourier-1];
            for (int i = 0; i < 2*res_fourier; ++i) {
            //for (int i = 0; i < res_fourier; ++i) {
                ////Float e = 1 - Float(i)/(2*res_fourier-1);
                ////densobs->datap[2*res_fourier+i] -= e*l+ (1-e)*r;
                avg += densobs->datap[2*res_fourier+i];
                //avg += densobs->datap[res_fourier+i]/universe.nParticles;
            }
            avg /= 2*res_fourier;
            //avg /= res_fourier;
            for (int i = 0; i < 2*res_fourier; ++i) {
            //for (int i = 0; i < res_fourier; ++i) {
                Float w = std::fabs(Float(i-(res_fourier-1/2))/Float(res_fourier-1/2));
                w = w*w;
                w = 1 - w;
                w = 1;
                /* let's already normalize the FFT by 1/N = 1/(2*res_fourier) before*/
                in[i] = w*(densobs->datap[2*res_fourier+i]/avg - 1)/(2*res_fourier);
                //in[i] = w*(densobs->datap[res_fourier+i]/universe.nParticles - avg);
            }

            out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * 2*res_fourier);
            //out = (double*) fftw_malloc(sizeof(double) *res_fourier);

            p = fftw_plan_dft_r2c_1d(2 * res_fourier, in, out, FFTW_ESTIMATE);
            //p = fftw_plan_r2r_1d(res_fourier, in, out, FFTW_REDFT10, FFTW_ESTIMATE); //can try 01 too (odd at right boundary... what does it mean? 
            fftw_execute(p); 


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
                    //datap[i+res] += A*L*hor*(out[idx]*out[idx])/nPoints;

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

template<class Dens>
class CorrelationFunction : public Measurement {

public:

    CorrelationFunction(int res,/* int zoom, */Float L, bool pseudo3d) : Measurement(4*res*2), res(res), pseudo3d(pseudo3d) {
        dataTypeSize = 4;
        rows = 2;
        cols = res;
        datap = (float *) data;
        maxSize = L*hor/*/zoom*/;

        reset();
        densobs = new Dens(2*res);
        //densobs->set_zoom(zoom);
    }
   
    bool pseudo3d; 

    virtual void reset() { 
        for (int i = 0; i < res; ++i) {
            /* the factor 1/2 deserves an explantion: 
             * in the final FT from power spectrum to correlation we use variants of the FFTW r2r trafos. Since their modes are living on a mirrored space, they are "double as long".
             * this means that transforming from modes for k= i * 2 pi / L to the signal, the osciallations are actually half a sin/cos, one, 1.5, 2, ... 
             * such that the *last* mode corresponds to Nyquist (which is why with these trafos there is never the need to throw away some symmetric part 
             * Since our modes were obtained by DFT on an interval of length L, the reconstruction will correspond to only half that and have length L/2 
             * (even if the modes are double as long as usual, this is the way it goes - L was fixed and we only see half of the reconstruction) 
             */
            datap[i] = 0.5*(maxSize*i)/(res-1);
            datap[i+res] = 0;
        }
    }

    virtual ~CorrelationFunction() {
        delete densobs;
    } 

    virtual void measure(const Universe& universe, int N) {

        Measurement::measure(universe, N);

        float A = 1.f / N; 

        /* don't pass N: reset measurment after each and average later */
        densobs->measure(universe, 1); 
        fftw_complex *out;
        double *in;
        fftw_plan p;
        in = (double*) fftw_malloc(sizeof(double) * 2 * res);


        Float avg = 0;
        for (int i = 0; i < 2*res; ++i)
            avg += densobs->datap[2*res+i];
        avg /= (2*res);

        /* assume that average is nonzero, that is, density is not a density contrast! overall normalization of density drops when converting to the contrast. 
         * It's normalization in turn is not arbitrary. The conversion can clearly only be carried out once. 
           Again, the last factor 1/(2*res) is just to pre-deal with the missing FFTW normalization */

        for (int i = 0; i < 2*res; ++i) {
            in[i] = (densobs->datap[2*res+i]/avg-1)/(2*res);
        }

        /* remove linear trend. it is okay to introduce shift to nonzero mean here */ 

        //Float l = in[0];
        //Float r = in[2*res-1];

        //for (int i = 0; i < 2*res; ++i) 
            //in[i] -= Float(2*res-1-i)/Float(2*res-1)*l + Float(i)/Float(2*res-1)*r;

        out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * 2*res);

        p = fftw_plan_dft_r2c_1d(2 * res, in, out, FFTW_ESTIMATE);
        fftw_execute(p); 


        datap[res] = 0;

        /* no factors of V = L*hor needed, after postprocess they will always drop out */
        /* in this case, the result is P_1D(k) / (L*hor) */
        for (int i = 1; i < res; ++i) {
            Float fac = 1;
            if (pseudo3d) {
                /* divide by k^2/pi to get 3d power spectrum */
                /* in this case, it is P_3d / (L*hor)^3 */ 
                Float k = i*2*pi;
                fac = 1/(k*k/pi); 
            }
            datap[i+res] += A*fac*(out[i][0]*out[i][0] + out[i][1]*out[i][1]);
            //datap[i+res] += in[2*i]; //A*fac*(out[i][0]*out[i][0] + out[i][1]*out[i][1]);
        }
        
        fftw_destroy_plan(p);
        fftw_free(in); 
        fftw_free(out);

    }

    void postprocess() { 

        double *out;
        double *in;
        fftw_plan p;
        in = (double*) fftw_malloc(sizeof(double) * res);
        out = (double*) fftw_malloc(sizeof(double) * res);
       
        if (!pseudo3d) { 

            /* copy k=0 mode in the zero index */
            for (int i = 0; i < res; ++i) {
                in[i] = datap[res+i];//  - ((res-1.-i)/(res-1)*datap[res+1] + Float(i)/(res-1)*datap[2*res-1]);
            }

            /* remove delta at k=0, that is, constant offset of correlation fctn if there is any*/
            in[0] = in[1];

            /* cos trafo: even on first index (k=0) and odd on last (should work well if there function has fallen off and small derivative too) 
             * note that it does not fall of to k->0 for the equalCorrelation P_1d(k) which goes to const!) */
            p = fftw_plan_r2r_1d(res, in, out, FFTW_REDFT00, FFTW_ESTIMATE); //can try 01 too (odd at right boundary... what does it mean? 
            fftw_execute(p); 

            /* FFTW includes factor of 2 from the mirrored part */ 
            for (int i = 0; i < res; ++i) {
                datap[i+res] = out[i];
            }
        }
        else { 

            /* skip k=0 mode in the zero index */
            for (int i = 0; i < res-1; ++i) {
                in[i] = datap[res+i+1]*2*pi*i/8/pi/pi;
                //in[i] = datap[res+i]  - ((res-1.-i)/(res-1)*datap[res+1] + Float(i)/(res-1)*datap[2*res-1]);
            }
            /* 0 padding */
            in[res-1] = 0;

            /* sin trafo: odd around -1 (k=0) and even on last (again we hope it has fallen off) */
            p = fftw_plan_r2r_1d(res, in, out, FFTW_RODFT01, FFTW_ESTIMATE); //can try 01 too (odd at right boundary... what does it mean? 
            fftw_execute(p); 

            /* drop last */
            for (int i = 1; i < res; ++i) {
                datap[i+res] = out[i] / (Float(i)/(res-1));
            }
        }

        fftw_destroy_plan(p);
        fftw_free(in); 
        fftw_free(out);
    }

private:

    int res;

    //PowerSpectrumObs * psobs; 
    Dens * densobs;
    Float maxSize;

    float* datap; 


};

//using CorrelationFunctionObs = CorrelationFunction<DensityObs2>;


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
