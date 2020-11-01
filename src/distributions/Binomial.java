package distributions;

import distributions.detail.DistBase;
import distributions.detail.Dpq;
import distributions.detail.RefDouble;

public class Binomial extends DistBase
{
    public static double pdf(double x, double n, double p)
    {
        return pdf(x, n, p, false);
    }

    public static double pdf(double x, double n, double p, boolean give_log)
    {
        if (Double.isNaN(x) || Double.isNaN(n) || Double.isNaN(p)) {
            return x + n + p;
        }
        if (p < 0.0 || p > 1.0 || Dpq.negINonInt(n)) {
            return Dpq.nanWarn();
        }
        if (x < 0.0 || Double.isInfinite(x)) {
            return Dpq.D0(give_log);
        }

        n = Math.round(n);
        x = Math.round(x);
        return binomialPdfRaw(x, n, p, 1.0 - p, give_log);
    }

    public static double cdf(double x, double n, double p)
    {
        return cdf(x, n, p, true, false);
    }

    public static double cdf(double x, double n, double p, boolean lower_tail)
    {
        return cdf(x, n, p, lower_tail, false);
    }

    public static double cdf(double x, double n, double p, boolean lower_tail, boolean log_p)
    {
        if (Double.isNaN(x) || Double.isNaN(n) || Double.isNaN(p)) {
            return x + n + p;
        }
        if (Double.isInfinite(n) || Double.isInfinite(p)) {
            return Dpq.nanWarn();
        }
        if (Dpq.nonInt(n)) {
            System.out.println("Non-integer n = " + n);
            return Dpq.nanWarn();
        }

        n = Math.round(n);
        if (n < 0.0 || p < 0.0 || p > 1.0) {
            return Dpq.nanWarn();
        }
        if (x < 0.0) {
            return Dpq.DT0(lower_tail, log_p);
        }

        x = Math.floor(x + 1e-7);
        if (n <= x) {
            return Dpq.DT1(lower_tail, log_p);
        }
        return Beta.cdf(p, x + 1.0, n - x, !lower_tail, log_p);
    }

    public static double quantile(double p, double n, double pr)
    {
        return quantile(p, n, pr, true, false);
    }

    public static double quantile(double p, double n, double pr, boolean lower_tail)
    {
        return quantile(p, n, pr, lower_tail, false);
    }

    public static double quantile(double p, double n, double pr, boolean lower_tail, boolean log_p)
    {
        if (Double.isNaN(p) || Double.isNaN(n) || Double.isNaN(pr)) {
            return p + n + pr;
        }
        if (Double.isInfinite(n) || Double.isInfinite(pr)) {
            return Dpq.nanWarn();
        }
        if (Double.isInfinite(p) && !log_p) {
            return Dpq.nanWarn();
        }
        if (n != Math.floor(n + 0.5)) {
            return Dpq.nanWarn();
        }
        if (pr < 0.0 || pr > 1.0 || n < 0.0) {
            return Dpq.nanWarn();
        }

        double l = 0.0;
        if (log_p) {
            if (p > 0.0) {
                return Dpq.nanWarn();
            }
            if (p == 0.0) {
                return lower_tail ? n : l;
            }
            if (p == Double.NEGATIVE_INFINITY) {
                return lower_tail ? l : n;
            }
        } else {
            if (p < 0.0 || p > 1.0) {
                return Dpq.nanWarn();
            }
            if (p == 0.0) {
                return lower_tail ? l : n;
            }
            if (p == 1.0) {
                return lower_tail ? n : l;
            }
        }

        if (pr == 0.0 || n == 0.0) {
            return 0.0;
        }

        double q = 1.0 - pr;
        if (q == 0.0) {
            return n;
        }

        double mu = n * pr;
        double sigma = Math.sqrt(n * pr * q);
        double gamma = (q - pr) / sigma;

        if (!lower_tail || log_p) {
            p = Dpq.DTqIv(p, lower_tail, log_p);
            if (p == 0.0) return 0.0;
            if (p == 1.0) return n;
        }
        if (p + 1.01 * Dpq.DBL_EPSILON >= 1.0) {
            return n;
        }

        RefDouble z = new RefDouble(0.0);
        z.val = Normal.quantile(p, 0.0, 1.0, true, false);
        double yes = Math.floor(mu + sigma * (z.val + gamma * (z.val * z.val - 1.0) / 6) + 0.5);

        if (yes > n) {
            yes = n;
        }

        z.val = cdf(yes, n, pr, true, false);
        p *= 1.0 - 64.0 * Dpq.DBL_EPSILON;

        if (n < 1e5) {
            return doSearch(yes, z, p, n, pr, 1.0);
        }

        {
            double incR = Math.floor(n * 0.001), oldInc;
            do {
                oldInc = incR;
                yes = doSearch(yes, z, p, n, pr, incR);
                incR = Math.max(1.0, Math.floor(incR / 100.0));
            } while (oldInc > 1.0 && incR > n * 1e-15);
            return yes;
        }
    }

    private static double doSearch(double yes, RefDouble z, double p, double n, double pr, double incR)
    {
        if (z.val >= p) {
            for (;;) {
                double newZ;
                if (yes == 0.0 || (newZ = cdf(yes - incR, n, pr, true, false)) < p)
                    return yes;
                yes = Math.max(0.0, yes - incR);
                z.val = newZ;
            }
        } else {
            for (;;) {
                yes = Math.min(yes + incR, n);
                if (yes == n || (z.val = cdf(yes, n, pr, true, false)) >= p)
                    return yes;
            }
        }
    }

    public static void main(String[] args)
    {
        System.out.println(pdf(75, 1500, 0.05));
        System.out.println(cdf(74, 1500, 0.05));
        System.out.println(cdf(74, 1500, 0.05, false));

        System.out.println(quantile(0.9568, 1500, 0.87));
    }
}
