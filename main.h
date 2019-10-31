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
static constexpr Float e = 2.718281828459045235360287471352662497757247093699959574966967627724076630353;

/* [100 km/s /Mpc (/c)]^-1 in Mpc */
static constexpr Float hor = 2997.92458; 

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

    BBKS(const Background& bg, Float L, bool useEisensteinHu, bool forceEqualCorrelation = true, Float omegab = 0.0224) : L(L), EisensteinHu(useEisensteinHu) {

        pk.reserve(nModes);
        Omh2 = bg.h*bg.h*bg.Om;
        Om = bg.Om;

        /* Planck 2018 */ 
        Obh2 = omegab;
        Ob = Obh2 / bg.h / bg.h;

        h = bg.h;
        //Ob = 0.00001;
        //Obh2 = Ob*bg.h*bg.h;

        double * in; 
        in = (double*) fftw_malloc(sizeof(double) * (nModes-1 + 2));

        /* for now we don't use the first and last element in the array */

        for (int i = 1; i < nModes; ++i) { 

            /* go from 2pi*i k's to 1/Mpc k's : divide by L in units of Mpc instead of inverse of 
             * H0 = h 100 km/s/Mpc /c = h 1/ (2997.92458 Mpc ) */
            Float k = 2*pi*i;
            Float kphys = k / (L*hor); 

            in[i] = eval(kphys);

            if (forceEqualCorrelation) {
                Float dkphys = 2*pi/(L*hor); 
                in[i] *= kphys / (4*pi) *dkphys;
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
 *       h(k) = 8 pi/k int dr_0^infty r exp(ikr) C(r)
 *
 *  the imaginart part, all we care about, selects the sin from the exp(ikr) only. If we extend the integrand r C(r) by antisymmetry to negative r, 
 *  we just double that part and the real part will come out to be zero then. Let's adopt this and actually define h(k) (where still P(k) = Im h(k)) 
 *
 *      h(k) = 4 pi/k int dr_-infty^infty r C(r) exp(ikr) 
 *
 *   but then 
 *      
 *      h(k) = 4pi/k -i \partial_k int dr exp(ikr) C(r) = 4pi/k -i \partial_k P_1D (k)  
 *
 *  so P(k) = Im h(k) implies 
 *
 *  P(k) k = - 4pi P_1D'(k)  !!! 
 *
 *  or 
 *
 *      P_1D(k) = const - 1/4pi int dk P(k) k 
 *
 *  const such that for large k P_1D is zero? 
 *
 *********
 *  
 *  As an aside, how is C(r) obtained from P(k)?
 *
 *  Since h(k) is purely imaginary, it is fully determined by P(k) from P(k) = Im h(k) as h(k) = i P(k)
 *
 *  It is also clear from its formula that h(k) k is hermitian so P(-k) for positive k and negative arguments must be = P(k) in the following when we write P(k) 
 *  where there was h(k)
 *
 *  From the expression for h(k), 
 *      
 *      i P(k) k / 4pi = int dr r C(r) exp(ikr) 
 *
 *  and this means that 
 *
 *      C(r) r =  1/8pi^2 int dk i P(k) k exp(-ikr)  where P(-k) = P(k), which is real and antisymmetric as it should be 
 *
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

    /* returns the 3d fits in Mpc^3 expecting k in 1/Mpc without little h! */
    Float eval(Float kphys) {

        if (!EisensteinHu) { /*BBKS fit */

            const Float a1 = 2.205, a2 = 4.05,  a3 = 18.3,  a4 = 8.725,  a5 = 8; 
            /*Float q = k./OmH2; %q formula actually also contains /Mpc^-1 so it is dimless (see BBKS App 7)  */
            /* I think it expects units of h/Mpc though, so to convert there by dividing our k by h */
            Float q = (kphys/h) / Omh2;  
            Float b1 = a1*q;
            Float b2 = a2*q;
            Float b3 = a3*q;
            Float b4 = a4*q;
            Float b5 = a5*q;

            Float pre = std::log(1+b1)/b1;
            pre = pre*pre;

            b5 = b5*b5;

            return A*std::pow(kphys/h, n) *pre/std::sqrt(1 + b2 + b3*b3 + b4*b4*b4 + b5*b5);
        }
        else {

            Float Thetainv = Float(2.7)/Tcmb;
            Float Thetainvsq = Thetainv*Thetainv;
            Float zeq = Float(2.5e4)*Omh2*Thetainvsq*Thetainvsq;
            Float keq = Float(7.46e-2)*Omh2*Thetainvsq;

            Float b1 = Float(0.313)*std::pow(Omh2, Float(-0.419))*(1+Float(0.607)*std::pow(Omh2, Float(0.674)));
            Float b2 = Float(0.238)*std::pow(Omh2, Float(0.223));
            Float zd = Float(1291)*std::pow(Omh2, Float(0.251))/(1+Float(0.659)*std::pow(Omh2, Float(0.828))) * (1+b1*std::pow(Obh2, b2));

            Float Req = Float(31.5)*Obh2*Thetainvsq*Thetainvsq/(Float(1e-3)*zeq);
            Float Rd =  Float(31.5)*Obh2*Thetainvsq*Thetainvsq/(Float(1e-3)*zd);

            Float s = 2/(3*keq)*std::sqrt(6/Req)*std::log((std::sqrt(1+Rd)+std::sqrt(Rd+Req))/(1+std::sqrt(Req)));
            Float ksilk = Float(1.6)*std::pow(Obh2,Float(0.52))*std::pow(Omh2, Float(0.73))*(1+std::pow(Float(10.4)*Omh2, Float(-0.95)));

            Float q = kphys / (Float(13.41)*keq);
            Float a1 = std::pow(Float(46.9)*Omh2, Float(0.670))*(1+std::pow(Float(32.1*Omh2),Float(-0.532)));
            Float a2 = std::pow(Float(12.0)*Omh2, Float(0.424))*(1+std::pow(Float(45.0*Omh2),Float(-0.582)));
                  b1 = Float(0.944)/(1+std::pow(458*Omh2, Float(-0.708)));
                  b2 = std::pow(Float(0.395)*Omh2, Float(-0.0266));

            Float ObOm = Ob/Om;
            Float alphac = std::pow(a1, -ObOm)*std::pow(a2, -ObOm*ObOm*ObOm);
            Float betac = 1/( 1+b1*(std::pow((Om-Ob)/Om,b2)-1) );

            auto C = [&] (Float aa) {return Float(14.2)/aa + 386/(1+Float(69.9)*std::pow(q, Float(1.08)));};
            auto T0 = [&] (Float aa, Float bb) {return std::log(e + Float(1.8)*bb*q)/(std::log(e + Float(1.8)*bb*q) + C(aa)*q*q);};
            Float x = kphys*s/Float(5.4);
            x *= x;
            x *= x;
            Float f = 1/(1+x);
            //f = 1;
            
            Float Tc = f*T0(1,betac) + (1-f)*T0(alphac, betac);


            Float y = (1+zeq)/(1+zd);
            Float G = y*(-6*std::sqrt(1+y) + (2+3*y)*std::log((std::sqrt(1+y)+1)/(std::sqrt(1+y)-1)));
            Float alphab = Float(2.07)*keq*s*std::pow(1+Rd, Float(-0.75))*G;
            Float betab = Float(0.5) + Ob/Om + (3-2*Ob/Om)*std::sqrt(Float(17.2*17.2)*Omh2*Omh2+1);
            Float betanode = Float(8.41)*std::pow(Omh2, Float(0.435));
            Float t = betanode/(kphys*s);
            t = t*t*t;
            Float stilde = s/std::pow(1+t, 1/Float(3));

            t = betab/(kphys*s);
            t = t*t*t;
            x = kphys*s/Float(5.2);
            x *= x;
            Float sincarg = kphys*stilde;
            Float Tb = (T0(1,1)/(1+x) + alphab/(1+t)*std::exp( - std::pow(kphys/ksilk, Float(1.4)))) * std::sin(sincarg)/sincarg;
           
            Float T = Ob/Om * Tb + (Om - Ob)/Om * Tc; 

            Float deltaH = 0;
             //Bunn & White, 96, "The Four-Year COBE Normalization and Large-Scale Structure" 
            //if (bg.isIntegrated) {
                //deltaH = As*bg.getGrowth(0)/bg.Om;  
            //}
            //else {
                Float nbar = n-1;
                /* see also (A1) - (A3) in Eisenstein-Hu 97 "Baryonic Features ... "*/
                deltaH = Float(1.94e-5)*std::pow(Om, -Float(0.785)-Float(0.05)*std::log(Om))*std::exp(-Float(0.95)*nbar - Float(0.169)*nbar*nbar);
                //std::cout << "Background not integrated yet when computing power spectrum. Using old normalization fit based on COBE and ignoring As." << std::endl;
            //}

            return A*deltaH*deltaH*2*pi*pi*hor*hor*hor/(h*h*h)*std::pow(kphys*hor/h, n) * T * T;
        }
    }
    /* returns the dimensionless power of the dimensionless fourier series mode, P(k)/L ! */
    virtual Float eval_dimless(int k) const  {

        /* since pk array got filled in constructor when A=1, we need to multiply by A here, too */

        if (k >= 0 && k < nModes)
            return A*pk[k]/(L*hor);
        else {
            std::cout << "Warning: evaluating BBKS power spectrum for higher mode number " << k << " than the maximal computed one " << nModes << " - returning 0." << std::endl;
            return 0;
        }

    }

    void set_dimless(int k, Float p)  {

        if (k >= 0 && k < nModes)
            pk[k] = p*L*hor/A;
        else {
            std::cout << "Warning: setting BBKS power spectrum for higher mode number " << k << " than the maximal computed one " << nModes << " - returning 0." << std::endl;
        }

    }

    int numModes() const {
        return nModes;
    }

    virtual ~BBKS() {
    }
  
    Float L;

    Float A = 1;
    Float n = 0.965;
    Float Tcmb = 2.728;

    Float Om, Ob, Omh2, Obh2, h;

    bool EisensteinHu = false;

    static constexpr int nModes = 1000000*8;

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


