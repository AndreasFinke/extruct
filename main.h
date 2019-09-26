#pragma once

using Float = double;
using Long = long;

#if PY == 1

#include <pybind11/numpy.h>

namespace py = pybind11;

#endif


#include <cmath>
#include "pcg32.h"
#include "fftw3.h"
#include "background.h"
#include "sum.h"

#include <complex>
#include <array>
#include <iostream>


static constexpr Float pi = 3.141592653589793238462643383279502884197169399375105820974944592307816;

/* [100 km/s /Mpc (/c)]^-1 in Mpc */
static constexpr Float hor = 2997.92348; 

using namespace std::complex_literals;

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

struct Entry {
    Entry(Long idx, Float t) : t(t), idx(idx) {} 
    Float t;
    Long idx;
    bool operator<(const Entry& rhs) const {return t < rhs.t;}
    friend std::ostream& operator<<(std::ostream& out, const Entry& obj);
};

template<class CollisionTask> 
struct TaskParticle {
    TaskParticle() : x(0), v(0), t(0), nCollisions(0), index(0), sheet(0) {} 
    TaskParticle(Float x, Float v, int index, short sheet) : x(x), v(v), nCollisions(0), index(index), sheet(sheet){}
    Float x, v, t;
    Long nCollisions;
    int index;
    short sheet;
    CollisionTask task;
};

class PowerSpectrum {
public:
    virtual ~PowerSpectrum() {}
    virtual Float eval_dimless(int k) const  = 0;
};

class PowerLaw : public PowerSpectrum {
public: 
    Float A = 0.1;
    virtual ~PowerLaw() {}
    virtual Float eval_dimless(int k) const { 
        Float kp = 2*pi*k;
        return A*kp*kp/(100+0*kp*kp*kp*kp);
    }
};

class BBKS : public PowerSpectrum {

public:

    BBKS(const Background& bg, Float L, bool forceEqualCorrelation = true) : L(L) {

        pk.reserve(nModes);
        Omh2 = bg.h*bg.h*bg.Om;

        const Float n = 0.9656;
        const Float a1 = 2.205, a2 = 4.05,  a3 = 18.3,  a4 = 8.725,  a5 = 8; 

        double * in; 
        in = (double*) fftw_malloc(sizeof(double) * (nModes-1 + 2));

        /* for now we don't use the first and last element in the array */

        for (int i = 1; i < nModes; ++i) { 

            /* go from 2pi*i k's to 1/Mpc k's : divide by L in units of Mpc instead of inverse of 
             * H0 = h 100 km/s/Mpc /c = h 1/ (2997.92458 Mpc ) */
            Float k = 2*pi*i;
            Float kphys = k / (L*hor); 

            /*Float q = k./OmH2; %q formula actually also contains /Mpc^-1 so it is dimless (see BBKS App 7)  */
            Float q = kphys / Omh2;  
            Float b1 = a1*q;
            Float b2 = a2*q;
            Float b3 = a3*q;
            Float b4 = a4*q;
            Float b5 = a5*q;

            Float pre = std::log(1+b1)/b1;
            pre = pre*pre;

            b5 = b5*b5;

            in[i] = std::pow(kphys, n) *pre/std::sqrt(1 + b2 + b3*b3 + b4*b4*b4 + b5*b5);

            if (forceEqualCorrelation) {
                Float dkphys = 2*pi/(L*hor); 
                in[i] *= kphys / (8*pi) *dkphys;
            }
            else
                in[i] *= kphys*kphys/pi;

        }

/* the simplest way to go from 3d to 1d is to force sigma(R) = 1/(2pi)^n int^1/R d^n k P(k) equal for all R, that is, the integrands in 1D and 3D:
 * sigma(0) may be divergent, but formally is the correlation function at 0, so the inverse FT of P(k) without exponential. 
 *
 *  1/2pi P_1D(k) = 1/(2pi)^3 4 pi k^2 P(k)
 *
 * so that 
 *
 * P_1d(k) = k^2 P(k) / pi 
 *
 * Alternatively, we can put the same correlation function: 
 *
 * <delta_k delta_k'*> = int int d^n r1 d^n r2 exp(i k r1) exp(-i k' r2) C(|r1-r2|) 
 * Let r = r1-r2, R = (r1+r2)/2
 * _____
 * in 3d:
 * del(r,R)/del(r1,r2) = |det(  id3 -id3   ) |
 *                       |   ( id3/2 id3/2 ) | =  |-1| = 1
 *
 * in 1d: |det( [ [1, -1,], [0.5, 0.5] ] ) | = 1 
 * _____ 
 *
 *   = int int d^n r d^n R exp(i k (r/2+R) ) exp(-i k' (-r/2 + R) ) C(|r|) 
 *   =                     exp(i (k-k')R ) exp(i (k+k') r/2) C(|r|) 
 *   = kron(k-k') int d^n r exp(i k r) C(|r|) 
 *
 *  in 1D: simple FT between P(k) = <delta_k delta_k*> and C(|r|) 
 *         note that P(k) = P(-k) can be extended to negative k and is real symmetric like C(r) == C(|r|) 
 *  in 3D: z = cos theta, then k r -> k r z (before k and r were vetors, now they are standing for |k|, |r| ) 
 *
 *  P(k) = 4 pi int dr r^2 int dz exp(i k r z) C(|r|) 
 *
 *       = 4 pi int dr r^2 (-i) (k r)^(-1) (exp(ikr) - exp(-ikr)) C(r)
 *        
 *       = p(k) + p*(k) where
 *
 *       p(k) = -i 4 pi/k int dr r exp(ikr) C(r)
 *
 *  so P(k) = 2 Re p(k) or 
 *
 *     P(k) = Im h(k) where
 *
 *       h(k) = 8 pi/k int dr r exp(ikr) C(r)
 *
 *     since r C(r) is real-odd, h(k) is purely imaginary and we can also do a sine transform to get it 
 *
 *  but also, 
 *
 *  h(k) = 8pi/k -i \partial_k int dr exp(ikr) C(r) = 8pi/k -i \partial_k P_1D (k)  
 *
 *  so 
 *
 *  P(k) k = - 8pi P_1D'(k)  !!! 
 *
 *  P_1D(k) = const - 1/8pi int dk P(k) k 
 *
 *  const such that for large k P_1D is zero? 
 */
      
        pk.push_back(0);
        
        if (forceEqualCorrelation) {
            for (int i = 1; i < nModes; ++i) {
                pk.push_back(pk.back() - in[i]);
            }
            for (int i = 1; i < nModes; ++i) {
                pk[i] = pk[i] - pk[nModes-1];  
                //pk[i] = A*in[i];
            }
        }
        else { 
            for (int i = 1; i < nModes; ++i) { 
                pk.push_back(in[i]);
            }
        }
        

        //double * out; 
        //fftw_plan p;
        //out = (double*)fftw_malloc(sizeof(double)*(nModes-1 + 2));

        //p = fftw_plan_r2r_1d(nModes-1, in+1, out+1, FFTW_RODFT00, FFTW_ESTIMATE);
        //fftw_execute(p);

        //[> out now contains C(r) r - and we want to transform C(r) back <]

        //for (int i = 0; i < nModes-1; ++i) {
            //in[i+1] = out[i+1]/((i+1)/double(nModes));
        //}
        
        /* now we try to guess value of correlation function at 0 and at the last index. 
         * For zero, we use a quadratic fit with zero slope at 0 to the next two points.
         * far out we set it to zero. */

        //in[nModes+1] = 0;
        //in[0] = (4*in[1] - in[2])/3;

        //fftw_destroy_plan(p);
        //p = fftw_plan_r2r_1d(nModes+1, in, out, FFTW_REDFT00, FFTW_ESTIMATE);
        //fftw_execute(p);

        //for (int i = 0; i < nModes; ++i) 
            //pk[i] = out[i];

        fftw_free(in); 
        //fftw_free(out);


    }

    /* returns the dimensionless power of the dimensionless fourier series mode, P(k)/L ! */
    virtual Float eval_dimless(int k) const  {


        if (k >= 0 && k < nModes)
            return A*pk[k]/L;
        else {
            std::cout << "Warning: evaluating BBKS power spectrum for higher mode number " << k << " than the maximal computed one " << nModes << " - returning 0." << std::endl;
            return 0;
        }

    }

    virtual ~BBKS() {
    }

    Float A = 0.1;
    Float Omh2 = 0.15;
    Float L = 1;
    static constexpr int nModes = 1000000*2;

private:
    std::vector<Float> pk; 
};


class File : public PowerSpectrum {

public:

    File(Float L) : L(L) {

        pk.reserve(nModes);
        pk.push_back(0);
        
        for (int i = 1; i < nModes; ++i) { 
            //pk.push_back(in[i]);
        }
    }
    /* returns the dimensionless power of the dimensionless fourier series mode, P(k)/L ! */
    virtual Float eval_dimless(int k) const  {
        if (k >= 0 && k < nModes)
            return A*pk[k]/L;
        else {
            std::cout << "Warning: evaluating BBKS power spectrum for higher mode number " << k << " than the maximal computed one " << nModes << " - returning 0." << std::endl;
            return 0;
        }

    }

    virtual ~File() {
    }

    Float A = 0.1;
    Float L = 1;
    static constexpr int nModes = 1000000*2;

private:
    std::vector<Float> pk; 
};
#include <cmath>
#include <cstdint>

/* Credits to User njuffa. 
 * https://stackoverflow.com/questions/42792939/implementation-of-sinpi-and-cospi-using-standard-c-math-library
 *
   "Writes result sine result sin(πa) to the location pointed to by sp
   Writes result cosine result cos(πa) to the location pointed to by cp

   In extensive testing, no errors > 0.97 ulp were found in either the sine
   or cosine results, suggesting the results returned are faithfully rounded."

*/
inline void sincospi (double a, double *sp, double *cp)
{
    double c, r, s, t, az;
    int64_t i;

    az = a * 0.0; // must be evaluated with IEEE-754 semantics
    /* for |a| >= 2**53, cospi(a) = 1.0, but cospi(Inf) = NaN */
    a = (std::fabs (a) < 9.0071992547409920e+15) ? a : az;  // 0x1.0p53
    /* reduce argument to primary approximation interval (-0.25, 0.25) */
    r = std::nearbyint (a + a); // must use IEEE-754 "to nearest" rounding
    i = (int64_t)r;
    t = std::fma (-0.5, r, a);
    /* compute core approximations */
    s = t * t;
    /* Approximate cos(pi*x) for x in [-0.25,0.25] */
    r =            -1.0369917389758117e-4;
    r = std::fma (r, s,  1.9294935641298806e-3);
    r = std::fma (r, s, -2.5806887942825395e-2);
    r = std::fma (r, s,  2.3533063028328211e-1);
    r = std::fma (r, s, -1.3352627688538006e+0);
    r = std::fma (r, s,  4.0587121264167623e+0);
    r = std::fma (r, s, -4.9348022005446790e+0);
    c = std::fma (r, s,  1.0000000000000000e+0);
    /* Approximate sin(pi*x) for x in [-0.25,0.25] */
    r =             4.6151442520157035e-4;
    r = std::fma (r, s, -7.3700183130883555e-3);
    r = std::fma (r, s,  8.2145868949323936e-2);
    r = std::fma (r, s, -5.9926452893214921e-1);
    r = std::fma (r, s,  2.5501640398732688e+0);
    r = std::fma (r, s, -5.1677127800499516e+0);
    s = s * t;
    r = r * s;
    s = std::fma (t, 3.1415926535897931e+0, r);
    /* map results according to quadrant */
    if (i & 2) {
        s = 0.0 - s; // must be evaluated with IEEE-754 semantics
        c = 0.0 - c; // must be evaluated with IEEE-754 semantics
    }
    if (i & 1) { 
        t = 0.0 - s; // must be evaluated with IEEE-754 semantics
        s = c;
        c = t;
    }
    /* IEEE-754: sinPi(+n) is +0 and sinPi(-n) is -0 for positive integers n */
    if (a == std::floor (a)) s = az;
    *sp = s;
    *cp = c;
}


inline double sinpi (double a) {
    double s, c;
    sincospi(a, &s, &c);
    return s;
}

inline double cospi (double a) {
    double s, c;
    sincospi(a, &s, &c);
    return c;
}


