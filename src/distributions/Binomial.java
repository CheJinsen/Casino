package distributions;

import distributions.detail.DistBase;
import distributions.detail.Dpq;
import distributions.detail.RefDouble;
import random.UniformRand;

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

    public static double rand(double nin, double pp)
    {
        double c = 0.0, fm = 0.0, npq = 0.0, p1 = 0.0, p2 = 0.0, p3 = 0.0, p4 = 0.0, qn = 0.0;
        double xl = 0.0, xll = 0.0, xlr = 0.0, xm = 0.0, xr = 0.0;

        double p_save = -1.0;
        int n_save = -1;
        int m = 0;

        double f, f1, f2, u, v, w, w2, x, x1, x2, z, z2;
        double p, q, np, g, r, al, alv, a_max_p, ffm, y_norm;
        int i, ix, k, n;

        if (Double.isInfinite(nin)) {
            return Dpq.nanWarn();
        }
        r = Math.round(nin);
        if (r != nin) {
            return Dpq.nanWarn();
        }
        if (Double.isInfinite(pp) || r < 0 || pp < 0.0 || pp > 1.0)	{
            return Dpq.nanWarn();
        }
        if (r == 0 || pp == 0.) {
            return 0;
        }
        if (pp == 1.) {
            return r;
        }
        if (r >= Integer.MAX_VALUE) {
            return quantile(UniformRand.rand(), r, pp, false, false);
        }

        n = (int)r;

        p = Math.min(pp, 1.0 - pp);
        q = 1.0 - p;
        np = n * p;
        r = p / q;
        g = r * (n + 1);

        if (pp != p_save || n != n_save) {
            p_save = pp;
            n_save = n;
            if (np < 30.0) {
                qn = Math.pow(q, n);
                //L_np_small:
                for (;;) {
                    ix = 0;
                    f = qn;
                    u = UniformRand.rand();
                    for (;;) {
                        if (u < f) {
                            if (p_save > 0.5)
                                ix = n - ix;
                            return ix;
                        }
                        if (ix > 110)
                            break;
                        u -= f;
                        ix++;
                        f *= (g / ix - r);
                    }
                }
            } else {
                ffm = np + p;
                m = (int)ffm;
                fm = m;
                npq = np * q;
                p1 = (int)(2.195 * Math.sqrt(npq) - 4.6 * q) + 0.5;
                xm = fm + 0.5;
                xl = xm - p1;
                xr = xm + p1;
                c = 0.134 + 20.5 / (15.3 + fm);
                al = (ffm - xl) / (ffm - xl * p);
                xll = al * (1.0 + 0.5 * al);
                al = (xr - ffm) / (xr * q);
                xlr = al * (1.0 + 0.5 * al);
                p2 = p1 * (1.0 + c + c);
                p3 = p2 + c / xll;
                p4 = p3 + c / xlr;
            }
        } else {
            if (np < 30.0) {
                //L_np_small:
                for (;;) {
                    ix = 0;
                    f = qn;
                    u = UniformRand.rand();
                    for (;;) {
                        if (u < f) {
                            return ix;
                        }
                        if (ix > 110)
                            break;
                        u -= f;
                        ix++;
                        f *= (g / ix - r);
                    }
                }
            }
        }

        for (;;) {
            u = UniformRand.rand() * p4;
            v = UniformRand.rand();

            if (u <= p1) {
                ix = (int)(xm - p1 * v + u);
                if (p_save > 0.5)
                    ix = n - ix;
                return ix;
            }

            if (u <= p2) {
                x = xl + (u - p1) / c;
                v = v * c + 1.0 - Math.abs(xm - x) / p1;
                if (v > 1.0 || v <= 0.)
                    continue;
                ix = (int) x;
            } else {
                if (u > p3) {
                    ix = (int)(xr - Math.log(v) / xlr);
                    if (ix > n)
                        continue;
                    v = v * (u - p3) * xlr;
                } else {
                    ix = (int)(xl + Math.log(v) / xll);
                    if (ix < 0)
                        continue;
                    v = v * (u - p2) * xll;
                }
            }

            k = Math.abs(ix - m);
            if (k <= 20 || k >= npq / 2 - 1) {
                f = 1.0;
                if (m < ix) {
                    for (i = m + 1; i <= ix; i++)
                        f *= (g / i - r);
                } else if (m != ix) {
                    for (i = ix + 1; i <= m; i++)
                        f /= (g / i - r);
                }
                if (v <= f) {
                    if (p_save > 0.5)
                        ix = n - ix;
                    return ix;
                }
            } else {
                a_max_p = (k / npq) * ((k * (k / 3. + 0.625) + 0.1666666666666) / npq + 0.5);
                y_norm = -k * k / (2.0 * npq);
                alv = Math.log(v);
                if (alv < y_norm - a_max_p) {
                    if (p_save > 0.5)
                        ix = n - ix;
                    return ix;
                }
                if (alv <= y_norm + a_max_p) {
                    x1 = ix + 1;
                    f1 = fm + 1.0;
                    z = n + 1 - fm;
                    w = n - ix + 1.0;
                    z2 = z * z;
                    x2 = x1 * x1;
                    f2 = f1 * f1;
                    w2 = w * w;
                    double temp = xm * Math.log(f1 / x1) + (n - m + 0.5) * Math.log(z / w) + (ix - m) *
                            Math.log(w * p / (x1 * q)) + (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / f2) /
                            f2) / f2) / f2) / f1 / 166320.0 + (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / z2) /
                            z2) / z2) / z2) / z / 166320.0 + (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / x2) /
                            x2) / x2) / x2) / x1 / 166320.0 + (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / w2) /
                            w2) / w2) / w2) / w / 166320.0;
                    if (alv <= temp) {
                        if (p_save > 0.5)
                            ix = n - ix;
                        return ix;
                    }
                }
            }
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

        for (int i = 0; i < 10; i++) {
            System.out.println(rand(1500, 0.87));
        }
    }
}
