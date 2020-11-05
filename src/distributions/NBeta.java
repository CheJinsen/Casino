package distributions;

import distributions.detail.*;

public class NBeta extends DistBase  // Non-central Beta Distribution
{
    public static double pdf(double x, double a, double b, double ncp)
    {
        return pdf(x, a, b, ncp, false);
    }

    public static double pdf(double x, double a, double b, double ncp, boolean give_log)
    {
        if (Double.isNaN(x) || Double.isNaN(a) || Double.isNaN(b) || Double.isNaN(ncp)) {
            return x + a + b + ncp;
        }
        if (ncp < 0 || a <= 0 || b <= 0) {
            return Dpq.nanWarn();
        }
        if (Double.isInfinite(a) || Double.isInfinite(b) || Double.isInfinite(ncp)) {
            return Dpq.nanWarn();
        }
        if (x < 0 || x > 1) {
            return Dpq.D0(give_log);
        }
        if(ncp == 0) {
            return Beta.pdf(x, a, b, give_log);
        }

        int kMax;
        double ncp2 = 0.5 * ncp;
        double dx2 = ncp2 * x;
        double d = (dx2 - a - 1)/2;
        double D = d * d + dx2 * (a + b) - a;
        final double eps = 1.e-15;

        if (D <= 0) {
            kMax = 0;
        } else {
            D = Math.ceil(d + Math.sqrt(D));
            kMax = (D > 0) ? (int)D : 0;
        }

        double term = Beta.pdf(x, a + kMax, b, true);
        double p_k = poissonPdfRaw(kMax, ncp2, true);
        if (x == 0.0 || Double.isInfinite(term) || Double.isInfinite(p_k)) {
            return Dpq.DExp(p_k + term, give_log);
        }

        p_k += term;
        double sum = term = 1.0;
        double k = kMax;
        double q;
        while (k > 0 && term > sum * eps) {
            k--;
            q = (k + 1) * (k + a) / (k + a + b) / dx2;
            term *= q;
            sum += term;
        }

        term = 1.0;
        k = kMax;
        do {
            q = dx2 * (k + a + b) / (k + a) / (k + 1);
            k++;
            term *= q;
            sum += term;
        } while (term > sum * eps);

        return Dpq.DExp(p_k + Math.log(sum), give_log);
    }

    public static double cdf(double x, double a, double b, double ncp)
    {
        return cdf(x, a, b, ncp, true, false);
    }

    public static double cdf(double x, double a, double b, double ncp, boolean lower_tail)
    {
        return cdf(x, a, b, ncp, lower_tail, false);
    }

    public static double cdf(double x, double a, double b, double ncp, boolean lower_tail, boolean log_p)
    {
        if (Double.isNaN(x) || Double.isNaN(a) || Double.isNaN(b) || Double.isNaN(ncp)) {
            return x + a + b + ncp;
        }
        if (x <= 0.0) {
            return Dpq.DT0(lower_tail, log_p);
        }
        if (x >= 1.0) {
            return Dpq.DT1(lower_tail, log_p);
        }
        return nonCentralBetaCdf2(x, 1.0 - x, a, b, ncp, lower_tail, log_p);
    }

    public static double quantile(double p, double a, double b, double ncp)
    {
        return quantile(p, a, b, ncp, true, false);
    }

    public static double quantile(double p, double a, double b, double ncp, boolean lower_tail)
    {
        return quantile(p, a, b, ncp, lower_tail, false);
    }

    public static double quantile(double p, double a, double b, double ncp, boolean lower_tail, boolean log_p)
    {
        final double accu = 1e-15;
        final double Eps = 1e-14;
        double ux, lx, nx, pp;

        if (Double.isNaN(p) || Double.isNaN(a) || Double.isNaN(b) || Double.isNaN(ncp)) {
            return p + a + b + ncp;
        }
        if (Double.isInfinite(a)) {
            return Dpq.nanWarn();
        }
        if (ncp < 0.0 || a <= 0.0 || b <= 0.0) {
            return Dpq.nanWarn();
        }

        double l = 0.0;
        double r = 1.0;
        if (log_p) {
            if (p > 0.0) {
                return Dpq.nanWarn();
            }
            if (p == 0.0) {
                return lower_tail ? r : l;
            }
            if (p == Double.NEGATIVE_INFINITY) {
                return lower_tail ? l : r;
            }
        } else {
            if (p < 0.0 || p > 1.0) {
                return Dpq.nanWarn();
            }
            if (p == 0.0) {
                return lower_tail ? l : r;
            }
            if (p == 1.0) {
                return lower_tail ? r : l;
            }
        }

        p = Dpq.DTqIv(p, lower_tail, log_p);
        if (p > 1 - Dpq.DBL_EPSILON) {
            return 1.0;
        }

        pp = Math.min(1.0 - Dpq.DBL_EPSILON, p * (1.0 + Eps));
        ux = 0.5;
        while (ux < 1 - Dpq.DBL_EPSILON && cdf(ux, a, b, ncp, true, false) < pp) {
            ux = 0.5 * (1 + ux);
        }

        pp = p * (1 - Eps);
        lx = 0.5;
        while (lx > Double.MIN_VALUE && cdf(lx, a, b, ncp, true, false) > pp) {
            lx *= 0.5;
        }

        do {
            nx = 0.5 * (lx + ux);
            if (cdf(nx, a, b, ncp, true, false) > p) {
                ux = nx;
            } else {
                lx = nx;
            }
        } while ((ux - lx) / nx > accu);

        return 0.5 * (ux + lx);
    }

    public static double rand(double shape1, double shape2, double ncp)
    {
        if (Double.isNaN(shape1) || Double.isInfinite(shape2) || Double.isInfinite(ncp)) {
            return Dpq.nanWarn();
        }
        if (shape1 < 0.0 || shape2 < 0.0 || ncp < 0.0) {
            return Dpq.nanWarn();
        }
        double x = NChiSquared.rand(2.0 * shape1, ncp);
        return x / (x + ChiSquared.rand(2.0 * shape2));
    }

    public static void main(String[] args)
    {
        System.out.println(pdf(0.56, 6, 8, 7));
        System.out.println(cdf(0.756, 6, 8, 70));
        System.out.println(quantile(0.756, 6, 8, 70));
        System.out.println();

        for (int i = 0; i < 10; i++) {
            System.out.println(rand(6, 8, 70));
        }
    }
}
