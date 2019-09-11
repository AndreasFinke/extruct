
#pragma once

#include "spline.h"
#include "main.h"
#include <utility>
#include <array>

struct Background {

    friend class Universe; 
    // default is standard LCDM without radiation but curvature 

    Background() {} //use standard parameters (see in class definition)

    Background(Float zin, Float zfin, Float Om, Float h, Float hf, Float Oc) : Om(Om), h(h), hf(hf), Oc(Oc), zin(zin), zfin(zfin) { 
    }

    Background(Float zin, Float zfin, Float Om, Float h, Float Oc) : Om(Om), h(h), hf(h), Oc(Oc), zin(zin), zfin(zfin) { 
    }
    // for other cosmologies, instead load the bg functions from file 
    Background(std::string& file) {
        //TODO read from file
        isIntegrated = true;
    }

    bool isIntegrated = false;

    //Background(Background&& orig) : Background(orig) {}

    //Background& operator=()
    //Background(const Background& orig) 
    //: Om(orig.Om), h(orig.h), hf(orig.hf), Oc(orig.Oc), zin(orig.zin), zfin(orig.zfin) {

        //isIntegrated = orig.isIntegrated;

        //if (isIntegrated) {

            //taufin = orig.taufin;
            //dtau = orig.dtau;

            //a = orig.a;
            //D1 = orig.D1;
            //D2 = orig.D2;
            //D1d = orig.D1d;
            //D2d = orig.D2d;
            //Pec = orig.Pec;
            //Pecd = orig.Pecd;
            //t = orig.t;
            ////memcpy(a, orig.a, NTABLE*sizeof(Float));  // copy data
            ////memcpy(D1, orig.D1, NTABLE*sizeof(Float));  // copy data
            ////memcpy(D2, orig.D2, NTABLE*sizeof(Float));  // copy data
            ////memcpy(D1d, orig.D1d, NTABLE*sizeof(Float));  // copy data
            ////memcpy(D2d, orig.D2d, NTABLE*sizeof(Float));  // copy data
            ////memcpy(Pec, orig.Pec, NTABLE*sizeof(Float));  // copy data
            ////memcpy(Pecd, orig.Pecd, NTABLE*sizeof(Float));  // copy data
            ////memcpy(t, orig.t, NTABLE*sizeof(Float));  // copy data
        //}
        ////else
            ////integrate();
    //}


private:
    Float zin = 100, zfin = 0, Om = 0.3, h = 0.7, hf = 0.7, Oc = 0;
    Float taufin = 0;
    Float dtau = 1;

    static constexpr int NTABLE = 2048*8;
    
    // the actual integration proceeds in finer (regular) steps - for this we require some temporary variables 
    static constexpr Float dtau_fine = 0.00001;

    std::array<Float, NTABLE> a; 
    std::array<Float, NTABLE> D1; 
    std::array<Float, NTABLE> D2; 
    std::array<Float, NTABLE> D1d; 
    std::array<Float, NTABLE> D2d; 
    std::array<Float, NTABLE> Pec; 
    std::array<Float, NTABLE> Pecd; 
    std::array<Float, NTABLE> t; 
    //Float a[NTABLE];
    //Float D1[NTABLE];
    //Float D2[NTABLE];
    //Float D1d[NTABLE];
    //Float D2d[NTABLE];
    //Float Pec[NTABLE];
    //Float Pecd[NTABLE];
    //Float t[NTABLE];

    // this is the Hubble rate (up to a missing factor of H0) given the scale factor, Om, and Oc 
    Float E(Float a) { 
        Float ainv = 1/a;
        Float asqinv = ainv*ainv;
        return sqrt(Om*asqinv*ainv + Oc*asqinv + (1-Oc-Om));
    };

    // returns pair of index of time bin to the left of argument as well as that bins time 
    auto getTimeBin(Float tau){
        int idx = std::min(NTABLE-1, int(tau/dtau));
        return std::make_pair(idx, idx*dtau);
    }

public:

    Float getOm() {return Om;} 
    Float getFinalTau() { return taufin; }

    Float getScaleFactor(Float tau) {
        int idx;
        Float taul;
        std::tie(idx, taul) = getTimeBin(tau);

        if (idx == NTABLE-1)
            return a[NTABLE-1];

        return Spline::y(a[idx], a[idx]*a[idx]*a[idx]*E(a[idx])*h, a[idx+1], a[idx+1]*a[idx+1]*a[idx+1]*E(a[idx+1])*h, taul, dtau, tau);
    }

    Float getTauOfZ(Float z) {
        int idxR = 1;
        for (; a[idxR] < 1/(1+z); ++idxR) {
            if (idxR == NTABLE)
                return taufin;
        }
        int idxL = idxR - 1;

        Float yL = a[idxL];
        Float yR = a[idxR];
        Float DyL = a[idxL]*a[idxL]*a[idxL]*E(a[idxL])*h;
        Float DyR = a[idxR]*a[idxR]*a[idxR]*E(a[idxR])*h;

        return Spline::find_zero(1/(1+z)-yL, -DyL, 1/(1+z)-yR, -DyR, idxL*dtau, dtau);
    }

    Float getAOfZ(Float z) { return getScaleFactor(getTauOfZ(z)); }

    Float getPhysTime(Float tau) {
        int idx;
        Float taul;
        std::tie(idx, taul) = getTimeBin(tau);

        if (idx == NTABLE-1)
            return t[NTABLE-1];

        return Spline::y(t[idx], a[idx]*a[idx], t[idx+1], a[idx+1]*a[idx+1], taul, dtau, tau);
    }

    Float getD1(Float tau) {
        int idx;
        Float taul;
        std::tie(idx, taul) = getTimeBin(tau);

        if (idx == NTABLE-1)
            return D1[NTABLE-1];

        return Spline::y(D1[idx], D1d[idx], D1[idx+1], D1d[idx+1], taul, dtau, tau);
    }
    
    Float getD2(Float tau) {
        int idx;
        Float taul;
        std::tie(idx, taul) = getTimeBin(tau);

        if (idx == NTABLE-1)
            return D2[NTABLE-1];

        return Spline::y(D2[idx], D2d[idx], D2[idx+1], D2d[idx+1], taul, dtau, tau);
    }

    Float getD1d(Float tau) {
        int idx;
        Float taul;
        std::tie(idx, taul) = getTimeBin(tau);

        if (idx == NTABLE-1)
            return D1d[NTABLE-1];

        Float m = 1.5*Om*h*h;
        //return Spline::DDy(D1[idx], D1d[idx], D1[idx+1], D1d[idx+1], taul, dtau, tau);
        return Spline::y(D1d[idx], m*a[idx]*D1[idx], D1d[idx+1], m*a[idx+1]*D1[idx+1], taul, dtau, tau);
    }

    Float getD2d(Float tau) {
        int idx;
        Float taul;
        std::tie(idx, taul) = getTimeBin(tau);

        if (idx == NTABLE-1)
            return D2d[NTABLE-1];

        Float m = 1.5*Om*h*h;
        //return Spline::Dy(D1d[idx], m*a[idx]*D1[idx], D1d[idx+1], m*a[idx+1]*D1[idx+1], taul, dtau, tau);
        return Spline::y(D2d[idx], m*a[idx]*D2[idx], D2d[idx+1], m*a[idx+1]*D2[idx+1], taul, dtau, tau);
    }
    
    Float getPec(Float tau) {
        int idx;
        Float taul;
        std::tie(idx, taul) = getTimeBin(tau);

        if (idx == NTABLE-1)
            return Pec[NTABLE-1];

        return Spline::y(Pec[idx], Pecd[idx], Pec[idx+1], Pecd[idx+1], taul, dtau, tau);
    }

    Float getPecd(Float tau) {
        int idx;
        Float taul;
        std::tie(idx, taul) = getTimeBin(tau);

        if (idx == NTABLE-1)
            return Pecd[NTABLE-1];

        Float m = 1.5*Om*h*h;
        return Spline::y(Pecd[idx], m*a[idx]*(Pec[idx]-1), Pecd[idx+1], m*a[idx+1]*(Pec[idx+1]-1), taul, dtau, tau);
    }

    void integrate() {

        std::cout << "Integrating bg... " << std::endl;


        isIntegrated = true;


        std::vector<Float> afine, D1fine, Pecfine, S2, tfine;
        afine.push_back(1/(zin+1));
        D1fine.push_back(0);
        Pecfine.push_back(0);
        S2.push_back(0);
        tfine.push_back(0);

        Float S1 = 0;


        // proceed until we are past the end point AND we have enough points to divide the set of points into NTABLE-1 equal intervals with NTABLE boundary points
        // e.g. 0 1 2 3 4 and NTABLE=3: can be divided into 2 intervals of length l=2 in the sense that (size-1)/(3-1) = 4/2 = 2, and then we get 3 points starting at index 0 and adding l=2: 0 2 4
        // (below, after the loop, this length l is called steplen)

        while ( (afine.back() < 1/(zfin + 1)) || ((afine.size()-1)%(NTABLE-1)) ) {

            Float acurr = afine.back();

            // note that dot refers to superconformal time, so adot = da/dtau = a^2 da/dt = a^3 H 
            Float adot = acurr*acurr*acurr*h*E(acurr);
            Float amid = acurr + 0.5*dtau_fine*adot;

            adot       = amid*amid*amid*h*E(amid);
            afine.push_back(acurr +  dtau_fine*adot);

            tfine.push_back(tfine.back() + dtau_fine*amid*amid);
            Float Einvsq = 1/E(amid);
            Einvsq = Einvsq*Einvsq;

            // S1 = int_tauin^tauout 1/E^2 dtau
            // S2 = int_tauin^tauout (a^-1 - ain^-1)/E^2 dtau
            // apply midpoint rule to differential equations obtained from differentiating these integrals 

            S1 += Einvsq*dtau_fine;
            S2.push_back(S2.back() + dtau_fine*Einvsq*(1/amid-1/afine[0]));

            //obtain D1 and peculiar solution at time step by multiplying S1 and S2 by constants and E(a) (no midpoint rule here!) 
            D1fine.push_back(S1*E(afine.back()));
            Pecfine.push_back(Float(1.5)*Om*hf*hf/h*S2.back()*E(afine.back()));
        }

        int steplen = (afine.size()-1)/(NTABLE-1);
        dtau = steplen * dtau_fine;

        int ind = 0;
        Float taucurr = 0;

        for (int i = 0; i < NTABLE; ++i) {
            //tau[i] = taucurr;

            a[i] = afine[ind];
            t[i] = tfine[ind];
            D2[i] = E(a[i]);
            // from differentiating Friedmann 
            D2d[i] = -h*(Float(1.5)*Om/a[i] + Oc); 
            D1[i] = D1fine[ind];
            // from conserved Wronskian 
            D1d[i] = (1+D1[i]*D2d[i])/D2[i];

            Pec[i] = Pecfine[ind];
            // from direct differentiation
            Pecd[i] = Float(1.5)*Om*hf*hf/h*(D2d[i]*S2[ind] + (1/a[i]-1/a[0])/D2[i]);

            taufin = taucurr; 
            taucurr += dtau;

            ind += steplen;
        }
    }
};

