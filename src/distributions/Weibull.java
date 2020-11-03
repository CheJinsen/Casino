package distributions;

import distributions.detail.Dpq;
import random.UniformRand;

public class Weibull
{
    public static double pdf(double x, double shape)
    {
        return pdf(x, shape, 1.0, false);
    }

    public static double pdf(double x, double shape, double scale)
    {
        return pdf(x, shape, scale, false);
    }

    public static double pdf(double x, double shape, double scale, boolean give_log)
    {
        if (Double.isNaN(x) || Double.isNaN(shape) || Double.isNaN(scale)) {
            return x + shape + scale;
        }
        if (shape <= 0 || scale <= 0) {
            return Dpq.nanWarn();
        }
        if (x < 0) {
            return Dpq.D0(give_log);
        }
        if (Double.isInfinite(x)) {
            return Dpq.D0(give_log);
        }
        if (x == 0 && shape < 1) {
            return Double.POSITIVE_INFINITY;
        }

        double tmp1 = Math.pow(x / scale, shape - 1);
        double tmp2 = tmp1 * (x / scale);
        return  give_log ?
                -tmp2 + Math.log(shape * tmp1 / scale) :
                shape * tmp1 * Math.exp(-tmp2) / scale;
    }

    public static double cdf(double x, double shape)
    {
        return cdf(x, shape, 1.0, true, false);
    }

    public static double cdf(double x, double shape, double scale)
    {
        return cdf(x, shape, scale, true, false);
    }

    public static double cdf(double x, double shape, double scale, boolean lower_tail)
    {
        return cdf(x, shape, scale, lower_tail, false);
    }

    public static double cdf(double x, double shape, double scale, boolean lower_tail, boolean log_p)
    {
        if (Double.isNaN(x) || Double.isNaN(shape) || Double.isNaN(scale)) {
            return x + shape + scale;
        }
        if (shape <= 0.0 || scale <= 0.0) {
            return Dpq.nanWarn();
        }
        if (x <= 0.0) {
            return Dpq.DT0(lower_tail, log_p);
        }

        x = -Math.pow(x / scale, shape);
        return lower_tail
                ? log_p ? Dpq.log1Exp(x) : -Math.expm1(x)
                : Dpq.DExp(x, log_p);
    }

    public static double quantile(double p, double shape)
    {
        return quantile(p, shape, 1.0, true, false);
    }

    public static double quantile(double p, double shape, double scale)
    {
        return quantile(p, shape, scale, true, false);
    }

    public static double quantile(double p, double shape, double scale, boolean lower_tail)
    {
        return quantile(p, shape, scale, lower_tail, false);
    }

    public static double quantile(double p, double shape, double scale, boolean lower_tail, boolean log_p)
    {
        if (Double.isNaN(p) || Double.isNaN(shape) || Double.isNaN(scale)) {
            return p + shape + scale;
        }
        if (shape <= 0.0 || scale <= 0.0) {
            return Dpq.nanWarn();
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

        return scale * Math.pow(-Dpq.DTCLog(p, lower_tail, log_p), 1.0 / shape);
    }

    public static double rand(double shape)
    {
        return rand(shape, 1.0);
    }

    public static double rand(double shape, double scale)
    {
        if (Double.isInfinite(shape) || Double.isInfinite(scale) || shape <= 0.0 || scale <= 0.0) {
            if (scale == 0.0)
                return 0.0;
            return Dpq.nanWarn();
        }
        return scale * Math.pow(-Math.log(UniformRand.rand()), 1.0 / shape);
    }

    public static void main(String[] args)
    {
        System.out.println(pdf(0.8765, 1.0));
        System.out.println(pdf(0.8765, 1.0, 9));
        System.out.println(cdf(0.8765, 10));
        System.out.println(cdf(0.8765, 10, 20));
        System.out.println(cdf(0.8765, 10, 20, false));
        System.out.println(quantile(0.8765, 10));

        for (int i = 0; i < 10; i++) {
            System.out.println(rand(9.88));
        }
    }
}
