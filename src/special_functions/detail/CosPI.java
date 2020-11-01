package special_functions.detail;

import distributions.detail.Dpq;

public class CosPI
{
    // sin(pi * x)  -- exact when x = k/2  for all integer k
    public static double sinPI(double x) {
        if (Double.isNaN(x)) return x;
        if(Double.isInfinite(x)) {
            System.out.println("Argument out of domain in chebyshev_eval");
            return Double.NaN;
        }

        x = fmod(x, 2.0);
        if(x <= -1)
            x += 2.0;
        else if (x > 1.0)
            x -= 2.0;

        if(x == 0.0 || x == 1.0) return 0.;
        if(x ==  0.5)	return  1.0;
        if(x == -0.5)	return -1.0;

        return Math.sin(Math.PI * x);
    }

    public static double tanPI(double x)
    {
        if (Double.isNaN(x)) {
            return x;
        }
        if (Double.isInfinite(x)) {
            return Dpq.nanWarn();
        }

        x = fmod(x, 1.0);
        if (x <= -0.5) {
            x++;
        } else if (x > 0.5) {
            x--;
        }
        return x == 0.0 ? 0.0 : (x == 0.5 ? Double.NaN : Math.tan(Math.PI * x));
    }

    private static double fmod(double x1, double x2)
    {
        double q = x1 / x2;
        return x1 - Math.floor(q) * x2;
    }
}
