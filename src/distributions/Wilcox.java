package distributions;

import distributions.detail.Dpq;
import random.UniformRand;
import special_functions.Choose;

import java.util.Vector;

public class Wilcox
{
    private static final Vector<Vector<Vector<Double>>> w = new Vector<>();

    public static double pdf(double x, double m, double n)
    {
        return pdf(x, m, n, false);
    }

    public static double pdf(double x, double m, double n, boolean give_log)
    {
        if (Double.isNaN(x) || Double.isNaN(m) || Double.isNaN(n)) {
            return x + m + n;
        }

        m = Math.round(m);
        n = Math.round(n);
        if (m <= 0.0 || n <= 0.0) {
            return Dpq.nanWarn();
        }

        if (Math.abs(x - Math.round(x)) > 1e-7) {
            return Dpq.D0(give_log);
        }

        x = Math.round(x);
        if (x < 0.0 || x > m * n) {
            return Dpq.D0(give_log);
        }

        int mm = (int)m, nn = (int)n, xx = (int)x;
        wInit(mm, nn);
        return give_log
                ? Math.log(cWilcox(xx, mm, nn)) - Choose.logChoose(m + n, n)
                : cWilcox(xx, mm, nn) / Choose.choose(m + n, n);
    }

    public static double cdf(double q, double m, double n)
    {
        return cdf(q, m, n, true, false);
    }

    public static double cdf(double q, double m, double n, boolean lower_tail)
    {
        return cdf(q, m, n, lower_tail, false);
    }

    public static double cdf(double q, double m, double n, boolean lower_tail, boolean log_p)
    {
        if (Double.isNaN(q) || Double.isNaN(m) || Double.isNaN(n)) {
            return q + m + n;
        }
        if (Double.isInfinite(m) || Double.isInfinite(n)) {
            return Dpq.nanWarn();
        }
        m = Math.round(m);
        n = Math.round(n);
        if (m <= 0 || n <= 0) {
            return Dpq.nanWarn();
        }

        q = Math.floor(q + 1e-7);
        if (q < 0.0) {
            return Dpq.DT0(lower_tail, log_p);
        }
        if (q >= m * n) {
            return Dpq.DT1(lower_tail, log_p);
        }

        int mm = (int)m, nn = (int) n;
        wInit(mm, nn);
        double c = Choose.choose(m + n, n);
        double p = 0.0;
        if (q <= (m * n / 2)) {
            for (int i = 0; i <= q; i++)
                p += cWilcox(i, mm, nn) / c;
        } else {
            q = m * n - q;
            for (int i = 0; i < q; i++) {
                p += cWilcox(i, mm, nn) / c;
            }
            lower_tail = !lower_tail;
        }

        return Dpq.DTVal(p, lower_tail, log_p);
    }

    public static double quantile(double x, double m, double n)
    {
        return quantile(x, m, n, true, false);
    }

    public static double quantile(double x, double m, double n, boolean lower_tail)
    {
        return quantile(x, m, n, lower_tail, false);
    }

    public static double quantile(double x, double m, double n, boolean lower_tail, boolean log_p)
    {
        if (Double.isNaN(x) || Double.isNaN(m) || Double.isNaN(n)) {
            return x + m + n;
        }
        if (Double.isInfinite(x) || Double.isInfinite(m) || Double.isInfinite(n)) {
            return Dpq.nanWarn();
        }
        if ((log_p && x > 0) || (!log_p && (x < 0 || x > 1))) {
            return Dpq.nanWarn();
        }

        m = Math.round(m);
        n = Math.round(n);
        if (m <= 0 || n <= 0) {
            return Dpq.nanWarn();
        }
        if (x == Dpq.DT0(lower_tail, log_p)) {
            return 0;
        }
        if (x == Dpq.DT1(lower_tail, log_p)) {
            return m * n;
        }
        if (log_p || !lower_tail) {
            x = Dpq.DTqIv(x, lower_tail, log_p);
        }

        int mm = (int) m, nn = (int) n;
        wInit(mm, nn);
        double c = Choose.choose(m + n, n);
        double p = 0;
        int q = 0;
        if (x <= 0.5) {
            x = x - 10 * Dpq.DBL_EPSILON;
            for (;;) {
                p += cWilcox(q, mm, nn) / c;
                if (p >= x)
                    break;
                q++;
            }
        } else {
            x = 1 - x + 10 * Dpq.DBL_EPSILON;
            for (;;) {
                p += cWilcox(q, mm, nn) / c;
                if (p > x) {
                    q = (int)(m * n - q);
                    break;
                }
                q++;
            }
        }
        return q;
    }

    public static double rand(double m, double n)
    {
        if (Double.isNaN(m) || Double.isNaN(n)) {
            return Dpq.nanWarn();
        }

        m = Math.round(m);
        n = Math.round(n);
        if (m < 0.0 || n < 0.0) {
            return Dpq.nanWarn();
        }
        if (m == 0.0 || n == 0.0) {
            return 0.0;
        }

        double r = 0.0;
        int k = (int)(m + n);
        Vector<Integer> x = new Vector<>();
        x.setSize(k);

        for (int i = 0; i < k; i++) {
            x.set(i, i);
        }
        for (int i = 0; i < n; i++) {
            int j = (int)uniformIndex(k);
            r += x.get(j);
            x.set(j, x.get(--k));
        }

        return r - n * (n - 1) / 2;
    }

    // This counts the number of choices with statistic = k
    private static double cWilcox(int k, int m, int n)
    {
        int u = m * n;
        if (k < 0 || k > u) {
            return 0;
        }

        int c = u / 2;
        if (k > c) {
            k = u - k;
        }

        int i, j;
        if (m < n) {
            i = m; j = n;
        } else {
            i = n; j = m;
        }

        if (j == 0) {
            return k == 0 ? 1.0 : 0.0;
        }
        if (j > 0 && k < j) {
            return cWilcox(k, i, k);
        }
        if (w.get(i).get(j) == null) {
            Vector<Double> temp = new Vector<>();
            temp.setSize(c + 1);
            w.get(i).set(j, temp);
            for (int l = 0; l <= c; l++) {
                w.get(i).get(j).set(l, -1.0);
            }
        }
        if (w.get(i).get(j).get(k) < 0.0) {
            w.get(i).get(j).set(k, cWilcox(k - j, i - 1, j) + cWilcox(k, i, j - 1));
        }

        return w.get(i).get(j).get(k);
    }

    private static void wInit(int m, int n)
    {
        if (m > n) {
            int temp = n; n = m; m = temp;
        }

        w.clear();
        m = Math.max(m, 50);
        n = Math.max(n, 50);
        w.setSize(m + 1);

        for (int i = 0; i <= m; i++) {
            Vector<Vector<Double>> temp = new Vector<>();
            temp.setSize(n + 1);
            w.set(i, temp);
        }
    }

    private static double rBits(int bits)
    {
        long v = 0L;
        for (int n = 0; n <= bits; n += 16) {
            int v1 = (int)Math.floor(UniformRand.rand() * 65536);
            v = 65536 * v + v1;
        }

        final long one64 = 1L;
        return (double)(v & ((one64 << bits) - 1));
    }

    private static double uniformIndex(double dn)
    {
        if (dn <= 0.0) {
            return 0.0;
        }
        int bits = (int)Math.ceil(Math.log(dn) / Dpq.M_LN2);
        double dv;
        do {
            dv = rBits(bits);
        } while (dn <= dv);
        return dv;
    }

    public static void main(String[] args)
    {
        //System.out.println(pdf(5806, 1570, 175)); // OutOfMemoryError
        System.out.println(pdf(806, 157, 170));
        System.out.println(cdf(806, 157, 170));
        System.out.println(quantile(0.806, 15, 17));
        System.out.println();

        for (int i = 0; i < 10; i++) {
            System.out.println(rand(157, 170));
        }
    }
}
