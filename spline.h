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

    /* this function assumes l and r bracket a sign change from positive to negative. it can deal with the case that the left bracket 
     * sits on a second zero from neg. to pos. or even to the left of that, as long as the function is similar to a quadratic opening downwards in this case. It will then return the second zero (in the application, the left zero is often the current time after a collision, and one wants the next collision  */
    static Float find_zero(Float yl, Float Dyl, Float DDyl, Float yr, Float Dyr, Float DDyr, Float xl, Float dx) {

        /* first find zero of a quadratic fit (from the left) as a proxy - 
         * it would be fine to estimate the second derivative DDyl as finite diff (a function for this is included below)
         * but since in the present application we have it exactly, it is passed as a parameter here */ 

        /* the zeros of a quadratic fit at any point xl are 
         *  
         *   x0 = xl -y'/y'' +- [(y'/y'')^2 - 2y/y'']^(1/2)
         *
         *  when y is small catastrophic cancellation can happen for the + sign (right zero) if sign(y' y'') = + and for the - sign (left zero) if sign(y' y'') = - 
         *
         *  if we are close to the correct zero, like below when refining the first estimate of the zero, this is always the case. 
         *             
         *  then we should use the alternative solution formula (obtained from expanding with  (+- []^(1/2) + y'/y'')/(+- []^(1/2) + y'/y'') )
         *  
         *  x0 = xl - 2y/y'' 1/( y'/y'' +- [(y'/y'')^2 - 2 y/y'']^1/2)
         *
         *  which can be written with Q = 2 y y'' / y'^2  (small) as 
         *
         *  x0 = xl - 2y/y' / (1 +- sign(y' y'') sqrt(1 - Q) )  (*) 
         *
         *  and again, close to the correct zero, +- sign(y' y'') is actually always +1 so 
         *
         *  x0 = xl - 2y/y' / (1 + sqrt(1 - Q) )                (A)    
         *
         *  and we use a taylor expansion of this in the refinement step below.
         *
         *  First, however, we estimate the position of the zero more roughly.
         *  In the case of two zeros and negative y'' the correct zero for our application is the RIGHT one (+ sign). 
         *  This is always true for negative y'', because there is a guarantee we do bracket a positive region and we want the right zero crossing from positive to negative. 
         *  Conversely, if y'' is positive, because of the guarantee crossing from positive to negative at the zero of interest 
         *  we want the *left* zero. 
         *
         *  Let us consider the y'' < 0 case first. Often, y ~ 0 (can be positive or negative) but we are close to the left zero, where y' > 0, so sign(y' y'') = -1 and we want the + sign. In this case, the original formula is much better than the "improved" one! 
         *
         *   x0 = xl - y'/y'' + [(y'/y'')^2 - 2y/y'']^(1/2)     (B) 
         *
         *   We may also be closer to the right zero, but typically y is not small anymore, so this formula is still fine. 
         *   Finally, very close the the correct right zero, the formula starts becoming inaccurate but maybe this would not be dramatic. 
         *   An exception could maybe be if y'' -> 0 also, which may really break the square root completely. 
         *   To cover this case, when y'' < 0, we use (B) if y' > 0 and switch to (A) as soon as y' < 0
         *
         *   But what if y '' < 0, y' > 0, and y'' -> 0? In this case it is normal that something explodes; the right zero is far away.
         *   y'' cannot be zero, since then there would not be any zero crossing from + to - as was guaranteed. 
         *   So no worries about the y'' small case! But an assert(y'' != 0)... 
         *
         *  If y'' > 0, we are left of the left zero, and y' < 0. This is well dealt with by the formula (*) where we choose the negative sifgn and sign(y' y'') = -1, which again leads (A) and entertains the limit y'' -> 0.
         *
         */ 

        Float x = xl;

        if ((DDyl <= 0) && (Dyl > 0))  {
            std::cout << "other branch " << std::endl;
            assert(DDyl != 0);
            Float DDylinv = 1/DDyl;
            Float foo = Dyl * DDylinv; 
            Float D = foo*foo - 2*yl*DDylinv;
                std::cout << "First case" << std::endl;
                std::cout << "D = " << D << " is negative." << std::endl;
                std::cout << "yl " << yl << " Dyl " << Dyl << " DDyl " << DDyl << " yr " << yr << " Dyr " << Dyr << " DDyr " << DDyr << " xl " << xl << " dx " << dx << std::endl << std::endl; 
            if (D >= 0)
                x -= foo + sqrt(D);
            else {
                std::cout << "First case" << std::endl;
                std::cout << "D = " << D << " is negative." << std::endl;
                std::cout << "yl " << yl << " Dyl " << Dyl << " DDyl " << DDyl << " yr " << yr << " Dyr " << Dyr << " DDyr " << DDyr << " xl " << xl << " dx " << dx << std::endl << std::endl; 
            }
        }
        else {

            Float Dylinv = 1/Dyl;
            Float D = 1 - 2*yl*DDyl*Dylinv*Dylinv;


                //std::cout << "Second case" << std::endl;
                //std::cout << "D = " << D << " is negative." << std::endl;
                //std::cout << "yl " << yl << " Dyl " << Dyl << " DDyl " << DDyl << " yr " << yr << " Dyr " << Dyr << " DDyr " << DDyr << " xl " << xl << " dx " << dx << std::endl << std::endl; 
            if (D >= 0)
                x -= - 2*yl*Dylinv/(1 + std::sqrt(D) );
            else {
                std::cout << "Second case" << std::endl;
                std::cout << "D = " << D << " is negative." << std::endl;
                std::cout << "yl " << yl << " Dyl " << Dyl << " DDyl " << DDyl << " yr " << yr << " Dyr " << Dyr << " DDyr " << DDyr << " xl " << xl << " dx " << dx << std::endl << std::endl; 
            }

            assert(D>=0);
        }


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
        //return x - foo*(1 + 0.5*foo2*(1+foo2));
        x -= foo*(1 + 0.5*foo2*(1+foo2));

        for (int k = 0; k < 1; ++k) { 
            Float Dyinv = 1/Dy(yl, Dyl, yr, Dyr, xl, dx, x);
            //Float Dyinv = 1/y(Dyl, DDyl, Dyr, DDyr, xl, dx, x);
            Float foo = y(yl, Dyl, yr, Dyr, xl, dx, x)*Dyinv;
            x -= foo;
        }

        return x;


    }

    static Float find_zero(Float yl, Float Dyl, Float yr, Float Dyr, Float xl, Float dx) {

        Float DD = Dyr - Dyl;
        DD /= dx;

        Float Dylinv = 1/Dyl;
        Float D = 1 - 2*yl*DD*Dylinv*Dylinv;
        assert(D >= 0);
        
        Float x = xl - 2*yl*Dylinv/(1 + std::sqrt(D) );

        Float Dyinv = 1/Dy(yl, Dyl, yr, Dyr, xl, dx, x);
        Float foo = y(yl, Dyl, yr, Dyr, xl, dx, x)*Dyinv;
        Float foo2 = foo*Dyinv*DDy(yl, Dyl, yr, Dyr, xl, dx, x);
        return x - foo*(1 + 0.5*foo2*(1+foo2));
    }
};
