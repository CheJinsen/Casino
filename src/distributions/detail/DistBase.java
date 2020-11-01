package distributions.detail;

import random.UniformRand;
import special_functions.LogGamma;
import special_functions.detail.Catherine;

public class DistBase
{
    protected static double poissonPdfRaw(double x, double lambda, boolean give_log)
    {
        if (lambda == 0.0) {
            return (x == 0.0) ? Dpq.D1(give_log) : Dpq.D0(give_log);
        }
        if (Double.isInfinite(lambda)) {
            return Dpq.D0(give_log);
        }
        if (x < 0.0) {
            return Dpq.D0(give_log);
        }
        if (x <= lambda * Double.MIN_VALUE) {
            return Dpq.DExp(-lambda, give_log);
        }
        if (lambda < x * Double.MIN_VALUE) {
            if (Double.isInfinite(x))
                return Dpq.D0(give_log);
            return Dpq.DExp(-lambda + x * Math.log(lambda) - LogGamma.logGamma(x + 1.0), give_log);
        }
        return Dpq.DFExp(2.0 * Math.PI * x, -Catherine.stirlingError(x) - Catherine.bd0(x, lambda), give_log);
    }

    protected static double binomialPdfRaw(double x, double n, double p, double q, boolean give_log)
    {
        double lf, lc;
        if (p == 0.0) {
            return x == 0.0 ? Dpq.D1(give_log) : Dpq.D0(give_log);
        }
        if (q == 0.0) {
            return x == n ? Dpq.D1(give_log) : Dpq.D0(give_log);
        }
        if (x == 0.0) {
            if (n == 0.0)
                return Dpq.D1(give_log);
            lc = p < 0.1 ? -Catherine.bd0(n, n * q) - n * p : n * Math.log(q);
            return Dpq.DExp(lc, give_log);
        }
        if (x == n) {
            lc = q < 0.1 ? -Catherine.bd0(n, n * p) - n * q : n * Math.log(p);
            return Dpq.DExp(lc, give_log);
        }
        if (x < 0.0 || x > n) {
            return Dpq.D0(give_log);
        }

        lc = Catherine.stirlingError(n) - Catherine.stirlingError(x) - Catherine.stirlingError(n - x) -
                Catherine.bd0(x, n * p) - Catherine.bd0(n - x, n * q);
        lf = Dpq.M_LN_2PI + Math.log(x) + Math.log1p(-x / n);
        return Dpq.DExp(lc - 0.5 * lf, give_log);
    }

    protected static double logSpaceAdd(double log_x, double log_y)
    {
        return Math.max(log_x, log_y) + Math.log1p(Math.exp(-Math.abs(log_x - log_y)));
    }

    protected static double expRand()
    {
        double[] q = {
                0.6931471805599453,
                0.9333736875190459,
                0.9888777961838675,
                0.9984959252914960,
                0.9998292811061389,
                0.9999833164100727,
                0.9999985691438767,
                0.9999998906925558,
                0.9999999924734159,
                0.9999999995283275,
                0.9999999999728814,
                0.9999999999985598,
                0.9999999999999289,
                0.9999999999999968,
                0.9999999999999999,
                1.0000000000000000
        };

        double a = 0.0;
        double u = UniformRand.rand();

        while (u <= 0.0 || u >= 1.0) {
            u = UniformRand.rand();
        }
        for (;;) {
            u += u;
            if (u > 1.0)
                break;
            a += q[0];
        }

        u -= 1.0;
        if (u <= q[0]) {
            return a + u;
        }

        int i = 0;
        double uStar = UniformRand.rand();
        double uMin = uStar;
        do {
            uStar = UniformRand.rand();
            if (uMin > uStar)
                uMin = uStar;
            i++;
        } while (u > q[i]);
        return a + uMin * q[0];
    }
}
