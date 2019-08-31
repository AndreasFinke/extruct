#pragma once 

#include "main.h"
#include <vector>
#include <complex>
#include <cmath>
#include <memory>
#include <iostream>

#include "fftw3.h"

class RandomField {
public:

    RandomField(const Long nGrid) : nGrid(nGrid) {

        /* allocate memory */ 
        field.reserve(nGrid);
        modes.reserve(nGrid);
    }

    
    void generate(pcg32& pcg, PowerSpectrum* spec) {


        Long nModes = nGrid/2+1;
        modes.clear();
        modes.resize(nModes);

        for (int i = 0; i < nModes; ++i) {
            Float phase = pcg.nextFloat()*2*pi;
            std::complex<double> mode;
            mode.real(std::sin(phase)); 
            mode.imag(std::cos(phase)); 
            Float k = 2*pi*i;
            modes[i] = mode*spec->eval(k) / (1+k*k*k*k*1e-4);
        }
        modes[0] *= 0;


        // to be absolutely sure it is safe we avoid reinterpret_cast<fftw_complex*>(&mode[0]) 

        fftw_complex *in;
        double *out;
        fftw_plan p;
        in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nModes);
         for (int i = 0; i < nModes; ++i) {
            in[i][0] = modes[i].real();
            in[i][1] = modes[i].imag();
        }
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

    Float get_field(Float x) {
        assert(x >= 0);
        Long idx = std::min(Long(x*nGrid), nGrid);
        return field[idx];
    }

    /* Assuming to be in the linear regime, the density is related to the displacement of an Eulerian grid in a simple way. Similarly, this allows us to get the velocities from the time derivative. The result will be rather uniform, which is also good, as averaging different nonlinear realizations is not a good strategy (this Monte Carlo strategy averaging results for delta-function ICs that average to the right ones makes sense only for linear evolution). Zeldovich holds also at later times in 1D, but in the nonlinear regime it becomes more difficult to get the right initial displacements that give a given nonlinear density. One could find any other method to put particles according to the random field at any time and get the displacements from comparing to the nearest free Eulerian grid position... but all of it is simple in the linear regime. */


private:
    std::vector<double> field;
    std::vector<std::complex<double>> modes;
    const Long nGrid = 256;
};

