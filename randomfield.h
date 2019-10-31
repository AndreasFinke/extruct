#pragma once 

#include "main.h"
#include <vector>
#include <complex>
#include <cmath>
#include <memory>
#include <iostream>
#include <limits>

#include "fftw3.h"
#include "ExactSum.h"

class RandomField {
public:

    RandomField(Long n) : nGrid(n) {

        /* allocate memory */ 
        field.reserve(nGrid);
        modes.reserve(nGrid);

        /* these will be deleted at the end of generate! */

        Long nModes = nGrid/2+1;
        out_ = (double*) fftw_malloc(sizeof(double) * nGrid);
        in_ = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nModes);
        /* planner is not thread save. We assume that RandomFields are always only constructed from a single thread
         * (the plan could also be global and passed as it should be possible to use the same plan in parallel, but this is ugly - I'd like to contain the FFTW stuff in here */

        p_ = fftw_plan_dft_c2r_1d(nGrid, in_, out_, FFTW_ESTIMATE);
        //p = fftw_plan_dft_1d(nGrid, reinterpret_cast<fftw_complex*>(&modes[0]), reinterpret_cast<fftw_complex*>(&field[0]), FFTW_BACKWARD, FFTW_ESTIMATE);
       
        /* this plan will not be deleted in generate, again, it is not clear this is thread-safe */

    }

    fftw_complex *in_;
    double *out_;
    fftw_plan p_;

    RandomField(RandomField&& other) : RandomField(other.nGrid) {
        field = std::move(other.field);
        modes = std::move(other.modes);
        in_ = other.in_;
        other.in_ = nullptr;
        out_ = other.out_;
        other.out_ = nullptr;
        p_ = other.p_;
        other.p_ = nullptr;
    }

    RandomField(const RandomField& other) = delete;

    ~RandomField() {
        if (p_ != nullptr)
            fftw_destroy_plan(p_);
    }

    /* this parallelizes but grid size has to be known */
    void generate(pcg32& pcg, PowerSpectrum* spec) {

        Long nModes = nGrid/2+1;
        modes.clear();
        modes.resize(nModes);
        
        for (int i = 1; i < nModes; ++i) {
            Float u1 = pcg.nextFloat();
            Float u2 = pcg.nextFloat();
            while (u1 <= std::numeric_limits<Float>::min())
                u1 = pcg.nextFloat();
            Float s, c;
            c = cospi(2*u2);
            Float normal = std::sqrt(-2 * std::log(u1)) * c;
            std::complex<double> mode;
            sincospi(pcg.nextFloat()*2, &s, &c);
            mode.real(s); 
            mode.imag(c); 
            Float k = 2*pi*i;
            modes[i] = mode*std::sqrt(spec->eval_dimless(i))*normal/**std::exp(-k*k*1e-8);*/;
        }
        modes[0] = 0;

        // to be absolutely sure it is safe we avoid reinterpret_cast<fftw_complex*>(&mode[0]) 

        for (int i = 1; i < nModes; ++i) {
            Float k = 2*pi*i;

            /* already divide by k here to get the displacement field */
            in_[i][0] = -modes[i].imag()/k;
            in_[i][1] = modes[i].real()/k;
        }

        in_[0][0] = 0;
        in_[0][1] = 0;
        /* let's not put stuff in nyquist */
        in_[nModes-1][1] = 0; 
        in_[nModes-1][0] = 0; 

        fftw_execute(p_); 

        displacementfield.clear(); 
        displacementfield.resize(nGrid);
        for (int i = 0; i < nGrid; ++i) {
            displacementfield[i] = out_[i];
        }

        /* c2r in FFTW destroys input array even for out-of-place trafo */

        for (int i = 1; i < nModes; ++i) {

            in_[i][0] = modes[i].real();
            in_[i][1] = modes[i].imag();
        }

        in_[0][0] = 0;
        in_[0][1] = 0;
        /* let's not put stuff in nyquist */
        in_[nModes-1][1] = 0; 
        in_[nModes-1][0] = 0; 

        fftw_execute(p_); 

        field.clear(); 
        field.resize(nGrid);
        for (int i = 0; i < nGrid; ++i) {
            field[i] = out_[i];
        }

        //fftw_destroy_plan(p_);
        fftw_free(in_); 
        fftw_free(out_);

        //std::cout << "field0 " << field[0] << " field1 " << field[1] << std::endl;

    }


    Float get_field_at(Float x) const {

        assert(x >= 0);
        Long idx = x*nGrid;
            
        idx = std::max(long(0), idx);
        idx = std::min(idx, nGrid-1);
        return field[idx];
        
        
            //Sum s;
        //for (int i = 1; i < modes.size(); ++i) {
            //double sn, cs;
            //sincospi(2.0*i*x, &sn, &cs);
            //s += 2*(sn*modes[i].imag() + cs*modes[i].real());
        //}
        //return s;
        
        //Sum r;
       
        /* assume nGrid is divisible by 2 */
        //for (int i = idx-nGrid/2; i < idx+nGrid/2; ++i) {
            //Float sincarg = x*nGrid-i;
            //Float sinc = 1;
            //if (std::fabs(sincarg) > 0.0000001)
                //sinc = sinpi(sincarg)/(pi*sincarg);
            //else
                //sinc = 1 - pi*pi*sincarg*sincarg/6;

            //int j = i;
            //if (j < 0) j = j + nGrid;
            //if (j >= nGrid) j = j - nGrid;
            //r += field[j]*sinc;
        //}
        //return r;
    }

    /* this displacement still needs to be multiplied by L to be correct */
    Float get_displacement_at(Float x) const {

        assert(x >= 0);
        Long idx = x*nGrid;
            
        //Sum s;
        //for (int i = 1; i < modes.size(); ++i) {
            //double sn, cs;
            //sincospi(2*i*x, &sn, &cs);
            //Float k = 2*pi*i;
            //s += 2*(-cs*modes[i].imag() + sn*modes[i].real())/k;
        //}
        //return s;

        
        //if (std::fabs(x*nGrid-std::round(x*nGrid)) < 1e-4) {
            idx = std::max(long(0), idx);
            idx = std::min(idx, nGrid-1);
            return displacementfield[idx];
        //}

        /* else, interpolate with Whittaker-Shannon: this does not include higher frequency components but gives the unique 
         * band-limited signal of highest band-limit that can exactly be represented by some discrete samples */

        //Sum r;
       
        //[> assume nGrid is divisible by 2 <]
        //for (int i = idx-nGrid/2; i < idx+nGrid/2; ++i) {
            //Float sincarg = x*nGrid-i;
            //Float sinc = 1;
            //if (std::fabs(sincarg) > 0.0000001)
                //sinc = sinpi(sincarg)/(pi*sincarg);
            //else
                //sinc = 1 - pi*pi*sincarg*sincarg/6;

            //int j = i;
            //if (j < 0) j = j + nGrid;
            //if (j >= nGrid) j = j - nGrid;
            //r += displacementfield[j]*sinc;
        //}
        //return r;

    }

    /* does not parallelize! */
    std::vector<Float> get_field(int res, int bandlimit = 0) const {

        if (bandlimit == 0) 
            bandlimit = res;

        std::vector<Float> ret; 
        ret.resize(res);

        /* return native field if limit exceeds it */
        if (bandlimit >= nGrid) {
            for (int i = 0; i < res; ++i) 
                ret[i] = get_field_at(Float(i)/Float(res));
            return ret;
        }

        /* else, we create an aliasing problem when simply subsampling the field. the best solution is to recompute the FFT with less modes */

        double * out = (double*) fftw_malloc(sizeof(double) * res);
        Long nModes = res/2+1;
        fftw_complex * in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nModes);

        fftw_plan pl;
        pl = fftw_plan_dft_c2r_1d(res, in, out, FFTW_ESTIMATE);

        for (int i = 1; i < nModes; ++i) {
            Float k = 2*pi*i;

            if (i < bandlimit/2+1) {
                in[i][0] = modes[i].real(); 
                in[i][1] = modes[i].imag(); 
            }
            else {
                in[i][0] = 0;
                in[i][1] = 0;
            }
        }

        in[0][0] = 0;
        in[0][1] = 0;
        //[> let's not put stuff in nyquist <]
        in[nModes-1][1] = 0; 
        in[nModes-1][0] = 0; 

        fftw_execute(pl); 

        for (int i = 0; i < res; ++i) {
            ret[i] = out[i];
        }

        fftw_destroy_plan(pl);
        fftw_free(in); 
        fftw_free(out);
        return ret;
    }


    /* does not parallelize! */
    std::vector<Float> get_displacement(int res, int bandlimit = 0) const {

        if (bandlimit == 0) 
            bandlimit = res;

        std::vector<Float> ret; 
        ret.resize(res);

        /* return native field if limit exceeds it */
        if (bandlimit >= nGrid) {
            for (int i = 0; i < res; ++i) 
                ret[i] = get_displacement_at(Float(i)/Float(res));
            return ret;
        }

        /* else, we have an aliasing problem simply subsampling the field. the best solution is to recompute the FFT with less modes */

        double * out = (double*) fftw_malloc(sizeof(double) * res);
        Long nModes = res/2+1;
        fftw_complex * in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nModes);

        fftw_plan pl = fftw_plan_dft_c2r_1d(res, in, out, FFTW_ESTIMATE);

        for (int i = 1; i < nModes; ++i) {
            Float k = 2*pi*i;

            if (i < bandlimit/2+1) {
                in[i][0] = -modes[i].imag()/k; 
                in[i][1] = modes[i].real()/k; 
            }
            else {
                in[i][0] = 0;
                in[i][1] = 0;
            }
        }

        in[0][0] = 0;
        in[0][1] = 0;
        /* let's not put stuff in nyquist */
        in[nModes-1][1] = 0; 
        in[nModes-1][0] = 0; 

        fftw_execute(pl); 

        for (int i = 0; i < res; ++i) {
            ret[i] = out[i];
        }

        fftw_destroy_plan(pl);

        fftw_free(in); 
        fftw_free(out);
        return ret;
    }
    /* Assuming to be in the linear regime, the density is related to the displacement of an Eulerian grid in a simple way. Similarly, this allows us to get the velocities from the time derivative. The result will be rather uniform, which is also good, as averaging different nonlinear realizations is not a good strategy (this Monte Carlo strategy averaging results for delta-function ICs that average to the right ones makes sense only for linear evolution). Zeldovich holds also at later times in 1D, but in the nonlinear regime it becomes more difficult to get the right initial displacements that give a given nonlinear density. One could find any other method to put particles according to the random field at any time and get the displacements from comparing to the nearest free Eulerian grid position... but all of it is simple in the linear regime. */


private:
    std::vector<double> field;
    std::vector<double> displacementfield;
    std::vector<std::complex<double>> modes;
    const Long nGrid;
};

