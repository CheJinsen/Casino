package distributions;

import distributions.detail.DistBase;
import distributions.detail.Dpq;

public class Exponential extends DistBase
{
    public static double pdf(double x)
    {
        return pdf(x, 1.0, false);
    }

    public static double pdf(double x, double scale)
    {
        return pdf(x, scale, false);
    }

    public static double pdf(double x, double scale, boolean give_log)
    {
        if (Double.isNaN(x) || Double.isNaN(scale)) {
            return x + scale;
        }
        if (scale <= 0.0) {
            return Dpq.nanWarn();
        }
        if (x < 0.0) {
            return Dpq.D0(give_log);
        }
        return give_log ? (-x / scale) - Math.log(scale) : Math.exp(-x / scale) / scale;
    }

    public static double cdf(double x)
    {
        return cdf(x, 1.0, true, false);
    }

    public static double cdf(double x, double scale)
    {
        return cdf(x, scale, true, false);
    }

    public static double cdf(double x, double scale, boolean lower_tail)
    {
        return cdf(x, scale, lower_tail, false);
    }

    public static double cdf(double x, double scale, boolean lower_tail, boolean log_p)
    {
        if (Double.isNaN(x) || Double.isNaN(scale)) {
            return x + scale;
        }
        if (scale < 0.0) {
            return Dpq.nanWarn();
        }
        if (scale <= 0.0) {
            return Dpq.nanWarn();
        }
        if (x <= 0.0) {
            return Dpq.DT0(lower_tail, log_p);
        }

        x = -(x / scale);
        return lower_tail
                ? (log_p ? Dpq.log1Exp(x) : -Math.expm1(x))
                : Dpq.DExp(x, log_p);
    }

    public static double quantile(double p, double scale, boolean lower_tail, boolean log_p)
    {
        if (Double.isNaN(p) || Double.isNaN(scale)) {
            return p + scale;
        }
        if (scale < 0.0) {
            return Dpq.nanWarn();
        }

        if ((log_p && p > 0) || (!log_p && (p < 0 || p > 1))) {
            return Dpq.nanWarn();
        }
        if (p == Dpq.DT0(lower_tail, log_p)) {
            return 0.0;
        }

        return -scale * Dpq.DTCLog(p, lower_tail, log_p);
    }

    public static double rand()
    {
        return rand(1.0);
    }

    public static double rand(double scale)
    {
        if (Double.isInfinite(scale) || scale <= 0.0) {
            if (scale == 0.0)
                return 0.0;
            return Dpq.nanWarn();
        }
        return scale * expRand();
    }

    public static void main(String[] args)
    {
        System.out.println(pdf(0.987));
        System.out.println(pdf(0.987, 10, false));
        System.out.println(cdf(0.025));
        System.out.println(cdf(0.025, 5));
        System.out.println(cdf(0.025, 5, false, true));

        for (int i = 0; i < 10; i++) {
            System.out.println(rand());
        }
    }
}
