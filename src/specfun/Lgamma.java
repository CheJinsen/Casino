package specfun;

import specfun.detail.Cospi;
import specfun.detail.Lgammacor;
import specfun.detail.Sign;

public class Lgamma {
    private final static double M_LN_SQRT_2PI	= 0.918938533204672741780329736406;	// log(sqrt(2*pi))
    private final static double M_LN_SQRT_PId2	= 0.225791352644727432363097614947;	// log(sqrt(pi/2))

    public static double lammafn(double x)
    {
        Sign sgn = new Sign(0);
        return lgammafn_sign(x, sgn);
    }

    public static double lgammafn_sign(double x, Sign sgn)
    {
        double ans, y, sinpiy;


        /* For IEEE double precision DBL_EPSILON = 2^-52 = 2.220446049250313e-16 :
           xmax  = DBL_MAX / log(DBL_MAX) = 2^1024 / (1024 * log(2)) = 2^1014 / log(2)
           dxrel = sqrt(DBL_EPSILON) = 2^-26 = 5^26 * 1e-26 (is *exact* below !)
         */
        final double xmax = 2.5327372760800758e+305;
        final double dxrel = 1.490116119384765625e-8;

        if (sgn.val != 0) sgn.val = 1;

        if(Double.isNaN(x)) return x;

        if (sgn.val != 0 && x < 0 && fmod(Math.floor(-x), 2.0) == 0) {
            sgn.val = -1;
        }

        if (x <= 0 && x == (int)x) { /* Negative integer argument */
            // No warning: this is the best answer; was  ML_WARNING(ME_RANGE, "lgamma");
            return Double.POSITIVE_INFINITY;/* +Inf, since lgamma(x) = log|gamma(x)| */
        }

        y = Math.abs(x);

        if (y < 1e-306) return - Math.log(y); // denormalized range, R change
        if (y <= 10) return Math.log(Math.abs(Gamma.gammafn(x)));
    /*
      ELSE  y = |x| > 10 ---------------------- */

        if (y > xmax) {
            // No warning: +Inf is the best answer
            return Double.POSITIVE_INFINITY;
        }

        if (x > 0) { /* i.e. y = x > 10 */
            double ret = M_LN_SQRT_2PI + (x - 0.5) * Math.log(x) - x;
            if(x > 1e17)
                return(x*(Math.log(x) - 1.));
            else if(x > 4934720.)
                return ret;
            else
                return ret + Lgammacor.lgammacor(x);
        }
        /* else: x < -10; y = -x */
        sinpiy = Math.abs(Cospi.sinpi(y));

        if (sinpiy == 0) { /* Negative integer argument ===
			  Now UNNECESSARY: caught above */
            System.out.printf(" ** should NEVER happen! *** [lgamma.c: Neg.int, y=%g]\n",y);
            return Double.NaN;
        }

        ans = M_LN_SQRT_PId2 + (x - 0.5) * Math.log(y) - x - Math.log(sinpiy) - Lgammacor.lgammacor(y);

        if(Math.abs((x - (int)(x - 0.5)) * ans / x) < dxrel) {

            /* The answer is less than half precision because
             * the argument is too near a negative integer. */

            System.out.println("Full precision may not have been achieved in lgamma");
        }

        return ans;
    }

    private static double fmod(double x1, double x2)
    {
        double q = x1 / x2;
        return x1 - Math.floor(q) * x2;
    }
}
