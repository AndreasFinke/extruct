#pragma once

#include "main.h"

class Spline {
public:

    static Float y(Float yl, Float Dyl, Float yr, Float Dyr, Float xl, Float dx, Float x) {
        Float t = x - xl;
        Float tt = t/dx;
        Float ttt = t*tt;
        Float tttt = tt*tt;
        Float dy = yr - yl;

        return yl + t*Dyl - ttt*(2*Dyl + Dyr) + tttt*( (3-2*tt)*dy + t*(Dyl + Dyr) );
    }

    static Float Dy(Float yl, Float Dyl, Float yr, Float Dyr, Float xl, Float dx, Float x) {
        Float t = x - xl;
        Float tt = t/dx;
        Float dy = yr - yl;

        return Dyl + tt*(-2*(2*Dyl + Dyr) + 6*dy/dx*(1-tt) + 3*tt*(Dyl+Dyr) );
    }

    static Float DDy(Float yl, Float Dyl, Float yr, Float Dyr, Float xl, Float dx, Float x) {
        Float t = x - xl;
        Float tt = t/dx;
        Float dy = yr - yl;

        Float dxinv = 1/dx;
        Float dydx = dy*dxinv;
        return ( -2*(2*Dyl + Dyr) + 6*dydx + tt*(-12*dydx + 6*(Dyl+Dyr)))*dxinv;
    }

    static Float find_zero(Float yl, Float Dyl, Float DDyl, Float yr, Float Dyr, Float DDyr, Float xl, Float dx) {

        /* first find zero of a quadratic fit (from the left) as a proxy - 
         * it would be fine to estimate the second derivative DDyl as finite diff, 
         * but since in the present application we have it exactly, it is passed as a parameter here */ 

        Float Dylinv = 1/Dyl;
        Float D = 1 - 2*yl*DDyl*Dylinv*Dylinv;
        assert(D >= 0);
        
        Float x = xl - 2*yl*Dylinv/(1 + std::sqrt(D) );

        /* find zero of spline zero using another quadratic fit (2nd order Taylor) based on the spline evaluated at the proxy.
         * Alternatively, we could have just do one step of Newton (1st order Taylor) or two steps (which would be very good, but require a lot of calculation for the extra step of Spline evaluations. A single quadratic fit is slightly faster and good. TODO: is this even true? need to evaluate second spline derivative, which isn't great either - probably quite comparable to two steps of Newton after all. 
         * Of course, it is safe to Taylor expand above formula for zero of quadratic (inverse and the square root). The third order expansion in y (which is small close to the zero) is almost for free compared to the second order, so we use it.  */


        /* spline is more accurate than derivative of spline (first has continuous second derivatives, second has only continuous first derivatives). We need interpolation of first derivative, and can get it from interpolating that directly. Similarly, 
         * first spline derivative is more accurate (has continuous first derivatives still) than second spline derivative.
         * We need interpolation of second derivative of y, 
         * but we can get it from a first spline derivative interpolating the first derivative of y by spline */

        Float Dyinv = 1/Dy(yl, Dyl, yr, Dyr, xl, dx, x);
        //Float Dyinv = 1/y(Dyl, DDyl, Dyr, DDyr, xl, dx, x);
        Float foo = y(yl, Dyl, yr, Dyr, xl, dx, x)*Dyinv;
        Float foo2 = foo*Dyinv*DDy(yl, Dyl, yr, Dyr, xl, dx, x);
        //Float foo2 = foo*Dyinv*Dy(Dyl, DDyl, Dyr, DDyr, xl, dx, x);
        return x - foo*(1 + 0.5*foo2*(1+foo2));
    }
};
