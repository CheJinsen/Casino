package specfun.detail;

public class Cospi {
    // sin(pi * x)  -- exact when x = k/2  for all integer k
    public static double sinpi(double x) {
        if (Double.isNaN(x)) return x;
        if(Double.isInfinite(x)) {
            System.out.println("Argument out of domain in chebyshev_eval");
            return Double.NaN;
        }

        x = fmod(x, 2.); // sin(pi(x + 2k)) == sin(pi x)  for all integer k
        // map (-2,2) --> (-1,1] :
        if(x <= -1) x += 2.; else if (x > 1.) x -= 2.;
        if(x == 0. || x == 1.) return 0.;
        if(x ==  0.5)	return  1.;
        if(x == -0.5)	return -1.;
        // otherwise
        return Math.sin(Math.PI * x);
    }

    private static double fmod(double x1, double x2)
    {
        double q = x1 / x2;
        return x1 - Math.floor(q) * x2;
    }
}
