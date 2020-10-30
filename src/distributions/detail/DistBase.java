package distributions.detail;

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
}
