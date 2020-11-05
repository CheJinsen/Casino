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

    protected static double nonCentralBetaCdf2(double x, double o_x, double a, double b,
                               double ncp, boolean lower_tail, boolean log_p)
    {
        double ans = nonCentralBetaCdfRaw(x, o_x, a,b, ncp);
        if (lower_tail) {
            return log_p ? Math.log(ans) : ans;
        } else {
            if (ans > 1.0 - 1e-10) {
                System.out.println("Full precision may not have been achieved in NBeta.cdf()");
            }
            if (ans > 1.0) {
                ans = 1.0;
            }
            return log_p ? Math.log1p(-ans) : (1.0 - ans);
        }
    }

    private static double nonCentralBetaCdfRaw(double x, double oXes, double a, double b, double ncp)
    {
        final double errMax = 1.0e-9;
        final int itrMax = 10000;

        if (ncp < 0.0 || a <= 0.0 || b <= 0.0) {
            return Dpq.nanWarn();
        }
        if (x < 0.0 || oXes > 1.0 || (x == 0.0 && oXes == 1.0)) {
            return 0.0;
        }
        if (x > 1.0 || oXes < 0.0 || (x == 1.0 && oXes == 0.0)) {
            return 1.0;
        }

        double c = ncp / 2.0;
        double x0 = Math.floor(Math.max(c - 7.0 * Math.sqrt(c), 0.0));
        double a0 = a + x0;
        double logBeta = LogGamma.logGamma(a0) + LogGamma.logGamma(b) - LogGamma.logGamma(a0 + b);

        RefDouble temp = new RefDouble(0.0);
        RefDouble tmp_c = new RefDouble(0.0);
        RefInt i_err = new RefInt(0);
        Toms.betaRatio(a0, b, x, oXes, temp, tmp_c, i_err, false);

        double gx = Math.exp(a0 * Math.log(x) + b * (x < 0.5 ? Math.log1p(-x) : Math.log(oXes)) -
                logBeta - Math.log(a0));
        double q;
        if (a0 > a) {
            q = Math.exp(-c + x0 * Math.log(c) - LogGamma.logGamma(x0 + 1.));
        } else {
            q = Math.exp(-c);
        }

        double ax, err_bd;
        double sumQ = 1.0 - q;
        double ans = q * temp.val;

        // recurse over subsequent terms until convergence is achieved
        double j = Math.floor(x0); // x0 could be billions, and is in package EnvStats
        do {
            j++;
            temp.val -= gx;
            gx *= x * (a + b + j - 1.) / (a + j);
            q *= c / j;
            sumQ -= q;
            ax = temp.val * q;
            ans += ax;
            err_bd = (temp.val - gx) * sumQ;
        } while (err_bd > errMax && j < itrMax + x0);

        if (err_bd > errMax) {
            System.out.println("Full precision may not have been achieved in NBeta.cdf()");
        }
        if (j >= itrMax + x0) {
            System.out.println("Convergence failed in NBeta.cdf()");
        }

        return ans;
    }
}
