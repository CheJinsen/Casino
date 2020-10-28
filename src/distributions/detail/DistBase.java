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
}
