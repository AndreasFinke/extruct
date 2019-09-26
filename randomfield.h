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
    }

    RandomField(RandomField&& other) : RandomField(other.nGrid) {
        field = std::move(other.field);
        modes = std::move(other.modes);
    }
    
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

        fftw_complex *in;
        double *out;
        fftw_plan p;
        in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nModes);
        for (int i = 1; i < nModes; ++i) {
            Float k = 2*pi*i;

            /* already divide by k here to get the displacement field */
            in[i][0] = modes[i].real()/k;
            in[i][1] = modes[i].imag()/k;
        }

        in[0][0] = 0;
        in[0][1] = 0;
        /* let's not put stuff in nyquist */
        in[nModes-1][1] = 0; 
        in[nModes-1][0] = 0; 

        out = (double*) fftw_malloc(sizeof(double) * nGrid);

        //p = fftw_plan_dft_1d(nGrid, reinterpret_cast<fftw_complex*>(&modes[0]), reinterpret_cast<fftw_complex*>(&field[0]), FFTW_BACKWARD, FFTW_ESTIMATE);
        p = fftw_plan_dft_c2r_1d(nGrid, in, out, FFTW_ESTIMATE);
        fftw_execute(p); 

        field.clear(); 
        field.resize(nGrid);
        for (int i = 0; i < nGrid; ++i) {
            field[i] = out[i];
        }

        fftw_destroy_plan(p);
        fftw_free(in); 
        fftw_free(out);


    }

    Float get_field(Float x) const {

        Sum s;
        for (int i = 1; i < modes.size(); ++i) {
            double sn, cs;
            sincospi(2.0*i*x, &sn, &cs);
            s += 2*(sn*modes[i].imag() + cs*modes[i].real());
        }
        return s;
    }

    /* this displacement still needs to be multiplied by L to be correct */
    Float get_displacement(Float x) const {

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
            idx = std::max(long(1), idx);
            idx = std::min(idx, nGrid-1);
            return field[idx];
        //}

        /* else, interpolate with Whittaker-Shannon: this does not include higher frequency components but gives the unique 
         * band-limited signal of highest band-limit that can exactly be represented by some discrete samples */

        Sum r; 
        for (int i = 0; i < nGrid; ++i) {
            Float sincarg = nGrid*x-i;
            r += field[i]*sinpi(sincarg)/(pi*sincarg);
        }
        return r;

    }
    /* Assuming to be in the linear regime, the density is related to the displacement of an Eulerian grid in a simple way. Similarly, this allows us to get the velocities from the time derivative. The result will be rather uniform, which is also good, as averaging different nonlinear realizations is not a good strategy (this Monte Carlo strategy averaging results for delta-function ICs that average to the right ones makes sense only for linear evolution). Zeldovich holds also at later times in 1D, but in the nonlinear regime it becomes more difficult to get the right initial displacements that give a given nonlinear density. One could find any other method to put particles according to the random field at any time and get the displacements from comparing to the nearest free Eulerian grid position... but all of it is simple in the linear regime. */


private:
    std::vector<double> field;
    std::vector<std::complex<double>> modes;
    const Long nGrid;
};

