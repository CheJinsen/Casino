package distributions;

import distributions.detail.Dpq;
import random.UniformRand;
import special_functions.detail.CosPI;

public class Cauchy
{
    public static double pdf(double x)
    {
        return pdf(x, 0.0, 1.0, false);
    }

    public static double pdf(double x, double location, double scale, boolean give_log)
    {
        if (Double.isNaN(x) || Double.isNaN(location) || Double.isNaN(scale)) {
            return x + location + scale;
        }
        if (scale <= 0.0) {
            return Dpq.nanWarn();
        }

        double y = (x - location) / scale;
        return give_log ? -Math.log(Math.PI * scale * (1.0 + y * y)) : 1.0 / (Math.PI * scale * (1.0 + y * y));
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
            if (x < 0.0)
                return Dpq.DT0(lower_tail, log_p);
            else
                return Dpq.DT1(lower_tail, log_p);
        }
        if (!lower_tail) {
            x = -x;
        }

        if (Math.abs(x) > 1.0) {
            double y = Math.atan(1.0 / x) / Math.PI;
            return x > 0.0 ? Dpq.DCLog(y, log_p) : Dpq.DVal(-y, log_p);
        } else {
            return Dpq.DVal(0.5 + Math.atan(x) / Math.PI, log_p);
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
        if ((log_p && p > 0) || (!log_p && (p < 0 || p > 1))) {
            return Dpq.nanWarn();
        }
        if (scale <= 0 || Double.isInfinite(scale)) {
            if (scale == 0)
                return location;
            return Dpq.nanWarn();
        }

        final double my_INF = location + (lower_tail ? scale : -scale) * Double.POSITIVE_INFINITY;
        if (log_p) {
            if (p > -1) {
                if (p == 0.0)
                    return my_INF;
                lower_tail = !lower_tail;
                p = -Math.expm1(p);
            } else {
                p = Math.exp(p);
            }
        } else {
            if (p > 0.5) {
                if (p == 1.0)
                    return my_INF;
                p = 1 - p;
                lower_tail = !lower_tail;
            }
        }

        if (p == 0.5) return location;
        if (p == 0.) return location + (lower_tail ? scale : -scale) * Double.NEGATIVE_INFINITY;
        return location + (lower_tail ? -scale : scale) / CosPI.tanPI(p);
    }

    public static double rand()
    {
        return rand(0.0, 1.0);
    }

    public static double rand(double location, double scale)
    {
        if (Double.isNaN(location) || Double.isInfinite(scale) || scale < 0.0) {
            return Dpq.nanWarn();
        }
        if (scale == 0.0 || Double.isInfinite(location)) {
            return location;
        } else {
            return location + scale * Math.tan(Math.PI * UniformRand.rand());
        }
    }

    public static void main(String[] args)
    {
        System.out.println(pdf(0.9987));
        System.out.println(cdf(9.8765));
        System.out.println(quantile(0.98765));
        System.out.println(quantile(0.98765, 50, 60));
        System.out.println(quantile(0.98765, 50, 60, false));

        System.out.println(rand());
    }
}
