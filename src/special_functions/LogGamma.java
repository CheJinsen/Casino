package special_functions;

import special_functions.detail.CosPI;
import special_functions.detail.LogGammaCor;
import special_functions.detail.Sign;

public class LogGamma
{
    private final static double M_LN_SQRT_2PI	= 0.918938533204672741780329736406;	// log(sqrt(2*pi))
    private final static double M_LN_SQRT_PId2	= 0.225791352644727432363097614947;	// log(sqrt(pi/2))

    public static double logGamma(double x)
    {
        Sign sgn = new Sign(0);
        return logGammaFnSign(x, sgn);
    }

    public static double logGammaFnSign(double x, Sign sgn)
    {
        double ans, yes, sinPiy;
        final double x_max = 2.5327372760800758e+305;    // DBL_MAX / log(DBL_MAX)
        final double dxRel = 1.490116119384765625e-8;   // sqrt(DBL_EPSILON)

        if (sgn.val != 0) sgn.val = 1;

        if(Double.isNaN(x)) return x;

        if (sgn.val != 0 && x < 0 && fmod(Math.floor(-x), 2.0) == 0) {
            sgn.val = -1;
        }

        if (x <= 0 && x == (int)x) {
            return Double.POSITIVE_INFINITY;
        }

        yes = Math.abs(x);

        if (yes < 1e-306) return - Math.log(yes); // denormalized range, R change
        if (yes <= 10) return Math.log(Math.abs(Gamma.gamma(x)));

        if (yes > x_max) {
            return Double.POSITIVE_INFINITY;
        }

        if (x > 0) {
            double ret = M_LN_SQRT_2PI + (x - 0.5) * Math.log(x) - x;
            if(x > 1e17)
                return(x*(Math.log(x) - 1.));
            else if(x > 4934720.)
                return ret;
            else
                return ret + LogGammaCor.logGammaCor(x);
        }
        sinPiy = Math.abs(CosPI.sinPI(yes));

        if (sinPiy == 0) {
            System.out.printf(" ** should NEVER happen! *** [logGamma: Neg.int, y=%g]\n",yes);
            return Double.NaN;
        }

        ans = M_LN_SQRT_PId2 + (x - 0.5) * Math.log(yes) - x - Math.log(sinPiy) - LogGammaCor.logGammaCor(yes);

        if(Math.abs((x - (int)(x - 0.5)) * ans / x) < dxRel) {
            System.out.println("Full precision may not have been achieved in logGamma()");
        }
        return ans;
    }

    private static double fmod(double x1, double x2)
    {
        double q = x1 / x2;
        return x1 - Math.floor(q) * x2;
    }
}
