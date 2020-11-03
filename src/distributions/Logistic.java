package distributions;

import distributions.detail.Dpq;
import random.UniformRand;

public class Logistic
{
    public static double pdf(double x)
    {
        return pdf(x, 0.0, 1.0, false);
    }

    public static double pdf(double x, double location, double scale)
    {
        return pdf(x, location, scale, false);
    }

    public static double pdf(double x, double location, double scale, boolean give_log)
    {
        if (Double.isNaN(x) || Double.isNaN(location) || Double.isNaN(scale)) {
            return x + location + scale;
        }
        if (scale <= 0.0) {
            return Dpq.nanWarn();
        }

        x = Math.abs((x - location) / scale);
        double e = Math.exp(-x);
        double f = 1.0 + e;
        return give_log ? -(x + Math.log(scale * f * f)) : e / (scale * f * f);
    }

    public static double cdf(double x)
    {
        return cdf(x, 0.0, 1.0, true, false);
    }

    public static double cdf(double x, double location, double scale)
    {
        return cdf(x, location, scale, true, false);
    }

    public static double cdf(double x, double location, double scale, boolean lower_tail)
    {
        return cdf(x, location, scale, lower_tail, false);
    }

    public static double cdf(double x, double location, double scale, boolean lower_tail, boolean log_p)
    {
        if (Double.isNaN(x) || Double.isNaN(location) || Double.isNaN(scale)) {
            return x + location + scale;
        }
        if (scale <= 0.0) {
            return Dpq.nanWarn();
        }

        x = (x - location) / scale;
        if (Double.isNaN(x)) {
            return Dpq.nanWarn();
        }
        if (Double.isInfinite(x)) {
            if (x > 0)  return Dpq.DT1(lower_tail, log_p);
            else        return Dpq.DT0(lower_tail, log_p);
        }
        if (log_p) {
            return -log1pExp(lower_tail ? -x : x);
        } else {
            return 1 / (1 + Math.exp(lower_tail ? -x : x));
        }
    }

    public static double quantile(double p)
    {
        return quantile(p, 0.0, 1.0, true, false);
    }

    public static double quantile(double p, double location, double scale)
    {
        return quantile(p, location, scale, true, false);
    }

    public static double quantile(double p, double location, double scale, boolean lower_tail)
    {
        return quantile(p, location, scale, lower_tail, false);
    }

    public static double quantile(double p, double location, double scale, boolean lower_tail, boolean log_p)
    {
        if (Double.isNaN(p) || Double.isNaN(location) || Double.isNaN(scale)) {
            return p + location + scale;
        }
        double l = Double.NEGATIVE_INFINITY;
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

        if (scale <	 0.0) {
            return Dpq.nanWarn();
        }
        if (scale == 0.0) {
            return location;
        }

        if (log_p) {
            if (lower_tail)
                p = p - Dpq.log1Exp(p);
            else
                p = Dpq.log1Exp(p) - p;
        } else {
            p = Math.log(lower_tail ? (p / (1. - p)) : ((1. - p) / p));
        }

        return location + scale * p;
    }

    public static double rand()
    {
        return rand(0.0, 1.0);
    }

    public static double rand(double location, double scale)
    {
        if (Double.isNaN(location) || Double.isInfinite(scale)) {
            return Dpq.nanWarn();
        }
        if (scale == 0.0 || Double.isInfinite(location)) {
            return location;
        } else {
            double u = UniformRand.rand();
            return location + scale * Math.log(u / (1.0 - u));
        }
    }

    private static double log1pExp(double x)
    {
        if (x <= 18.0)
            return Math.log1p(Math.exp(x));
        if (x > 33.3)
            return x;
        return x + Math.exp(-x);
    }

    public static void main(String[] args)
    {
        System.out.println(pdf(6.976));
        System.out.println(cdf(0.975));
        System.out.println(quantile(0.975));

        for (int i = 0; i < 10; i++) {
            System.out.println(rand());
        }
    }
}
