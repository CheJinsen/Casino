package distributions;

import distributions.detail.DistBase;
import distributions.detail.Dpq;

public class Geometric extends DistBase
{
    public static double pdf(double x, double p)
    {
        return pdf(x, p, false);
    }

    public static double pdf(double x, double p, boolean give_log)
    {
        if (Double.isNaN(x) || Double.isNaN(p)) {
            return Dpq.nanWarn();
        }
        if (p <= 0.0 || p > 1.0) {
            return Dpq.nanWarn();
        }
        if (Dpq.nonInt(x)) {
            System.out.println("non-integer x = " + x);
            return Dpq.D0(give_log);
        }
        if (x < 0.0 || Double.isInfinite(x) || p == 0.0) {
            return Dpq.D0(give_log);
        }

        x = Math.round(x);
        double prob = binomialPdfRaw(0.0, x, p, 1.0 - p, give_log);
        return give_log ? Math.log(p) + prob : p * prob;
    }

    public static double cdf(double x, double p)
    {
        return cdf(x, p, true, false);
    }

    public static double cdf(double x, double p, boolean lower_tail)
    {
        return cdf(x, p, lower_tail, false);
    }

    public static double cdf(double x, double p, boolean lower_tail, boolean log_p)
    {
        if (Double.isNaN(x) || Double.isNaN(p)) {
            return x + p;
        }
        if (p <= 0.0 || p > 1.0) {
            return Dpq.nanWarn();
        }
        if (x < 0.0) {
            return Dpq.DT0(lower_tail, log_p);
        }
        if (Double.isInfinite(x)) {
            return Dpq.DT1(lower_tail, log_p);
        }

        x = Math.floor(x + 1e-7);
        if (p == 1.0) {
            x = lower_tail ? 1: 0;
            return log_p ? Math.log(x) : x;
        }

        x = Math.log1p(-p) * (x + 1);
        if (log_p) {
            return Dpq.DTCLog(x, lower_tail, true);
        } else {
            return lower_tail ? -Math.expm1(x) : Math.exp(x);
        }
    }

    public static double quantile(double p, double prob)
    {
        return quantile(p, prob, true, false);
    }

    public static double quantile(double p, double prob, boolean lower_tail)
    {
        return quantile(p, prob, lower_tail, false);
    }

    public static double quantile(double p, double prob, boolean lower_tail, boolean log_p)
    {
        if (Double.isNaN(p) || Double.isNaN(prob)) {
            return p + prob;
        }
        if (prob <= 0.0 || prob > 1.0) {
            return Dpq.nanWarn();
        }
        if ((log_p && p > 0) || (!log_p && (p < 0 || p > 1))) {
            return Dpq.nanWarn();
        }
        if (prob == 1.0) {
            return 0.0;
        }

        double l = 0.0;
        double r = Double.POSITIVE_INFINITY;
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

        return Math.max(0.0, Math.ceil(Dpq.DTCLog(p, lower_tail, log_p) / Math.log1p(-prob) - 1 - 1e-12));
    }

    public static double rand(double p)
    {
        if (Double.isInfinite(p) || p <= 0.0 || p > 1.0) {
            return Dpq.nanWarn();
        }
        return Poisson.rand(expRand() * ((1.0 - p) / p));
    }

    public static void main(String[] args)
    {
        System.out.println(pdf(79, 0.5));
        System.out.println(cdf(9, 0.651));
        System.out.println(quantile(0.956, 0.025));

        for (int i = 0; i < 10; i++) {
            System.out.println(rand(0.025));
        }
    }
}
