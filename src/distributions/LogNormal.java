package distributions;

import distributions.detail.Dpq;

public class LogNormal
{
    public static double pdf(double x)
    {
        return pdf(x, 0.0, 1.0, false);
    }

    public static double pdf(double x, double meanLog, double sdLog, boolean give_log)
    {
        if (Double.isNaN(x) || Double.isNaN(meanLog) || Double.isNaN(sdLog)) {
            return x + meanLog + sdLog;
        }
        if (sdLog < 0.0) {
            return Dpq.nanWarn();
        }
        if (Double.isInfinite(x) && Math.log(x) == meanLog) {
            return Double.NaN;
        }
        if (sdLog == 0.0) {
            return Math.log(x) == meanLog ? Double.POSITIVE_INFINITY : Dpq.D0(give_log);
        }
        if (x <= 0.0) {
            return Dpq.D0(give_log);
        }

        double y = (Math.log(x) - meanLog) / sdLog;
        return give_log
                ? -(Dpq.M_LN_SQRT_2PI + 0.5 * y * y + Math.log(x * sdLog))
                : Dpq.M_1_SQRT_2PI * Math.exp(-0.5 * y * y) / (x * sdLog);
    }

    public static double cdf(double x)
    {
        return cdf(x, 0.0, 1.0, true, false);
    }

    public static double cdf(double x, double meanLog, double sdLog)
    {
        return cdf(x, meanLog, sdLog, true, false);
    }

    public static double cdf(double x, double meanLog, double sdLog, boolean lower_tail)
    {
        return cdf(x, meanLog, sdLog, lower_tail, false);
    }

    public static double cdf(double x, double meanLog, double sdLog, boolean lower_tail, boolean log_p)
    {
        if (Double.isNaN(x) || Double.isNaN(meanLog) || Double.isNaN(sdLog)) {
            return x + meanLog + sdLog;
        }
        if (sdLog < 0.0) {
            return Dpq.nanWarn();
        }
        if (x > 0.0) {
            return Normal.cdf(Math.log(x), meanLog, sdLog, lower_tail, log_p);
        }
        return Dpq.DT0(lower_tail, log_p);
    }

    public static double quantile(double p)
    {
        return quantile(p, 0.0, 1.0, true, false);
    }

    public static double quantile(double p, double meanLog, double sdLog)
    {
        return quantile(p, meanLog, sdLog, true, false);
    }

    public static double quantile(double p, double meanLog, double sdLog, boolean lower_tail)
    {
        return quantile(p, meanLog, sdLog, lower_tail, false);
    }

    public static double quantile(double p, double meanLog, double sdLog, boolean lower_tail, boolean log_p)
    {
        if (Double.isNaN(p) || Double.isNaN(meanLog) || Double.isNaN(sdLog)) {
            return p + meanLog + sdLog;
        }
        double r = Double.POSITIVE_INFINITY;
        double l = 0.0;
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

        return Math.exp(Normal.quantile(p, meanLog, sdLog, lower_tail, log_p));
    }

    public static double rand()
    {
        return rand(0.0, 1.0);
    }

    public static double rand(double meanLog, double sdLog)
    {
        if (Double.isNaN(meanLog) || Double.isNaN(sdLog) || sdLog < 0.0) {
            return Dpq.nanWarn();
        }
        return Math.exp(Normal.rand(meanLog, sdLog));
    }

    public static void main(String[] args)
    {
        System.out.println(pdf(9.8765));
        System.out.println(cdf(9.8765, 2.0, 5, false));
        System.out.println(quantile(0.12345));

        for (int i = 0; i < 10; i++) {
            System.out.println(rand());
        }
    }
}
