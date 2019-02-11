#pragma once 

#include "main.h"
#include <vector>
#include <complex>
#include <cmath>

#include "fftw3.h"

class RandomField {
public:

    RandomField(const Long nGrid) : nGrid(nGrid) {

        // set up random number generator

        /* allocate memory */ 
        field.reserve(nGrid);
        modes.reserve(nGrid);
    }

    
    void generate(pcg32& pcg, PowerSpectrum const * spec) {

        modes.clear();
        modes.reserve(nGrid);

        for (int i = 0; i < nGrid; ++i) {
            Float phase = pcg.nextFloat()*2*pi;
            std::complex<double> mode;
            mode.real(std::sin(phase)); 
            mode.imag(std::cos(phase)); 
            Float k = 2*pi*i;
            double A = 0.01; 
            modes.push_back(mode*spec->eval(k)*A);
        }
        modes[0] *= 0;

        field.clear();
        field.resize(nGrid);

        //fftw_complex *in, *out;
        fftw_plan p;
        //in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nGrid);
        //out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nGrid);
        p = fftw_plan_dft_1d(nGrid, reinterpret_cast<fftw_complex*>(&modes[0]), reinterpret_cast<fftw_complex*>(&field[0]), FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_execute(p); /* repeat as needed */
        fftw_destroy_plan(p);
        //fftw_free(in); 
        //fftw_free(out);


    }

    Float get_field(Float x) {
        Long idx = x*nGrid;
        assert(idx < nGrid && idx >= 0);
        return field[idx];
    }

    /* Assuming to be in the linear regime, the density is related to the displacement of an Eulerian grid in a simple way. Similarly, this allows us to get the velocities from the time derivative. The result will be rather uniform, which is also good, as averaging different nonlinear realizations is not a good strategy (this Monte Carlo strategy averaging results for delta-function ICs that average to the right ones makes sense only for linear evolution). Zeldovich holds also at later times in 1D, but in the nonlinear regime it becomes more difficult to get the right initial displacements that give a given nonlinear density. One could find any other method to put particles according to the random field at any time and get the displacements from comparing to the nearest free Eulerian grid position... but all of it is simple in the linear regime. */

    //std::vector<Particle> sample(Float DH)  { 
        //std::vector<Particle> p(nParticles);
        //for (int i = 0; i < nParticles; ++i) {
            //Float psi = 1;
            //p.push_back(Particle(psi, DH*psi));
        //}
        //return p;
    //}

private:
    std::vector<double> field;
    std::vector<std::complex<double>> modes;
    const Long nGrid = 256;
};

