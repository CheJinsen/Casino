package distributions;

import distributions.detail.Dpq;
import random.UniformRand;

import java.util.Vector;

public class SignRank
{
    private static final Vector<Double> w = new Vector<>();

    public static double pdf(double x, double n)
    {
        return pdf(x, n, false);
    }

    public static double pdf(double x, double n, boolean give_log)
    {
        if (Double.isNaN(x) || Double.isNaN(n)) {
            return (x + n);
        }

        n = Math.round(n);
        if (n <= 0) {
            return Dpq.nanWarn();
        }
        if (Math.abs(x - Math.round(x)) > 1e-7) {
            return Dpq.D0(give_log);
        }

        x = Math.round(x);
        if ((x < 0) || (x > (n * (n + 1) / 2))) {
            return Dpq.D0(give_log);
        }

        int nn = (int)n;
        wInit(nn);
        return Dpq.DExp(Math.log(cSignRank((int) x, nn)) - n * Dpq.M_LN2, give_log);
    }

    public static double cdf(double x, double n)
    {
        return cdf(x, n, true, false);
    }

    public static double cdf(double x, double n, boolean lower_tail)
    {
        return cdf(x, n, lower_tail, false);
    }

    public static double cdf(double x, double n, boolean lower_tail, boolean log_p)
    {
        if (Double.isNaN(x) || Double.isNaN(n)) {
            return x + n;
        }
        if (Double.isInfinite(n)) {
            return Dpq.nanWarn();
        }

        n = Math.round(n);
        if (n <= 0) {
            return Dpq.nanWarn();
        }

        x = Math.round(x + 1e-7);
        if (x < 0.0) {
            return Dpq.DT0(lower_tail, log_p);
        }
        if (x >= n * (n + 1) / 2) {
            return Dpq.DT1(lower_tail, log_p);
        }

        int nn = (int)n;
        wInit(nn);
        double f = Math.exp(-n * Dpq.M_LN2);
        double p = 0.0;
        if (x <= (n * (n + 1) / 4)) {
            for (int i = 0; i <= x; i++)
                p += cSignRank(i, nn) * f;
        } else {
            x = n * (n + 1) / 2 - x;
            for (int i = 0; i < x; i++)
                p += cSignRank(i, nn) * f;
            lower_tail = !lower_tail; /* p = 1 - p; */
        }

        return Dpq.DTVal(p, lower_tail, log_p);
    }

    public static double quantile(double x, double n)
    {
        return quantile(x, n, true, false);
    }

    public static double quantile(double x, double n, boolean lower_tail)
    {
        return quantile(x, n, lower_tail, false);
    }

    public static double quantile(double x, double n, boolean lower_tail, boolean log_p)
    {
        if (Double.isNaN(x) || Double.isNaN(n)) {
            return x + n;
        }
        if (Double.isInfinite(x) || Double.isInfinite(n)) {
            return Dpq.nanWarn();
        }
        if ((log_p && x > 0) || (!log_p && (x < 0 || x > 1))) {
            return Dpq.nanWarn();
        }

        n = Math.round(n);
        if (n <= 0) {
            return Dpq.nanWarn();
        }
        if (x == Dpq.DT0(lower_tail, log_p)) {
            return 0;
        }
        if (x == Dpq.DT1(lower_tail, log_p)) {
            return n * (n + 1) / 2;
        }
        if (log_p || !lower_tail) {
            x = Dpq.DTqIv(x, lower_tail, log_p);
        }

        int nn = (int)n;
        wInit(nn);
        double f = Math.exp(-n * Dpq.M_LN2);
        double p = 0.0;
        int q = 0;
        if (x <= 0.5) {
            x = x - 10 * Dpq.DBL_EPSILON;
            for (;;) {
                p += cSignRank(q, nn) * f;
                if (p >= x)
                    break;
                q++;
            }
        } else {
            x = 1 - x + 10 * Dpq.DBL_EPSILON;
            for (;;) {
                p += cSignRank(q, nn) * f;
                if (p > x) {
                    q = (int)(n * (n + 1) / 2 - q);
                    break;
                }
                q++;
            }
        }
        return q;
    }

    public static double rand(double n)
    {
        if (Double.isNaN(n)) {
            return n;
        }

        n = Math.round(n);
        if (n < 0.0) {
            return Dpq.nanWarn();
        }
        if (n == 0.0) {
            return 0.0;
        }

        double r = 0.0;
        int k = (int)n;
        for (int i = 0; i < k;) {
            r += (++i) * Math.floor(UniformRand.rand() + 0.5);
        }
        return r;
    }

    private static void wInit(int n)
    {
        int u = n * (n + 1) / 2;
        int c = (u / 2);

        w.clear();
        w.setSize(c + 1);
        for (int i = 0; i < c + 1; i++) {
            w.set(i, 0.0);
        }
    }

    private static double cSignRank(int k, int n)
    {
        int u = n * (n + 1) / 2;
        int c = u / 2;

        if (k < 0 || k > u) {
            return 0;
        }
        if (k > c) {
            k = u - k;
        }

        if (n == 1) {
            return 1.0;
        }
        if (w.get(0) == 1.0) {
            return w.get(k);
        }

        w.set(0, 1.0);
        w.set(1, 1.0);
        for(int j = 2; j < n+1; ++j) {
            int end = Math.min(j * (j + 1) / 2, c);
            for(int i = end; i >= j; --i)
                w.set(i, w.get(i) + w.get(i - j));
        }

        return w.get(k);
    }

    public static void main(String[] args)
    {
        System.out.println(pdf(76, 89));
        System.out.println(cdf(76, 89));
        System.out.println(quantile(0.76, 89));

        for (int i = 0; i < 10; i++) {
            System.out.println(rand(76));
        }
    }
}
