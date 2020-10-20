package specfun.detail;

public class Chebyshev {
    public static int chebyshev_init(double[] dos, int nos, double eta)
    {
        int i, ii;
        double err;

        if (nos < 1)
            return 0;

        err = 0.0;
        i = 0;			/* just to avoid compiler warnings */
        for (ii=1; ii<=nos; ii++) {
            i = nos - ii;
            err += Math.abs(dos[i]);
            if (err > eta) {
                return i;
            }
        }
        return i;
    }

    public static double chebyshev_eval(double x, double[] a, int n)
    {
        double b0, b1, b2, twox;

        if (n < 1 || n > 1000) {
            System.out.println("Argument out of domain in chebyshev_eval");
            return Double.NaN;
        }

        if (x < -1.1 || x > 1.1) {
            System.out.println("Argument out of domain in chebyshev_eval");
            return Double.NaN;
        }

        twox = x * 2;
        b2 = b1 = 0;
        b0 = 0;
        for (int i = 1; i <= n; i++) {
            b2 = b1;
            b1 = b0;
            b0 = twox * b1 - b2 + a[n - i];
        }
        return (b0 - b2) * 0.5;
    }
}
