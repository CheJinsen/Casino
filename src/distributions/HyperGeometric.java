package distributions;

import distributions.detail.DistBase;
import distributions.detail.Dpq;
import random.UniformRand;
import special_functions.Choose;

public class HyperGeometric extends DistBase
{
    public static double pdf(double x, double r, double b, double n)
    {
        return pdf(x, r, b, n, false);
    }

    public static double pdf(double x, double r, double b, double n, boolean give_log)
    {
        if (Double.isNaN(x) || Double.isNaN(r) || Double.isNaN(b) || Double.isNaN(n)) {
            return x + r + b + n;
        }
        if (Dpq.negINonInt(r) || Dpq.negINonInt(b) || Dpq.negINonInt(n) || n > r+b) {
            return Dpq.nanWarn();
        }
        if(x < 0) {
            return Dpq.D0(give_log);
        }
        if (Dpq.nonInt(x)) {
            System.out.printf("non-integer x = %f", x);
            return Dpq.D0(give_log);
        }

        x = Math.round(x);
        r = Math.round(r);
        b = Math.round(b);
        n = Math.round(n);

        if (n < x || r < x || n - x > b) {
            return Dpq.D0(give_log);
        }
        if (n == 0) {
            return (x == 0) ? Dpq.D1(give_log) : Dpq.D0(give_log);
        }

        double p = (n) / (r + b);
        double q = (r + b - n) / (r + b);

        double p1 = binomialPdfRaw(x, r, p, q, give_log);
        double p2 = binomialPdfRaw(n - x,b, p, q, give_log);
        double p3 = binomialPdfRaw(n,r + b, p, q, give_log);

        return (give_log) ? p1 + p2 - p3 : p1 * p2 / p3;
    }

    public static double cdf(double x, double NR, double NB, double n)
    {
        return cdf(x, NR, NB, n, true, false);
    }

    public static double cdf(double x, double NR, double NB, double n, boolean lower_tail)
    {
        return cdf(x, NR, NB, n, lower_tail, false);
    }

    public static double cdf(double x, double NR, double NB, double n, boolean lower_tail, boolean log_p)
    {
        if (Double.isNaN(x) || Double.isNaN(NR) || Double.isNaN(NB) || Double.isNaN(n)) {
            return x + NR + NB + n;
        }

        x = Math.floor(x + 1e-7);
        NR = Math.round(NR);
        NB = Math.round(NB);
        n  = Math.round(n);

        if (NR < 0 || NB < 0 || Double.isInfinite(NR + NB) || n < 0 || n > NR + NB) {
            return Dpq.nanWarn();
        }

        if (x * (NR + NB) > n * NR) {
            double oldNB = NB;
            NB = NR;
            NR = oldNB;
            x = n - x - 1;
            lower_tail = !lower_tail;
        }

        if (x < 0 || x < n - NB) {
            return Dpq.DT0(lower_tail, log_p);
        }
        if (x >= NR || x >= n) {
            return Dpq.DT1(lower_tail, log_p);
        }

        double d  = pdf(x, NR, NB, n, log_p);
        if ((!log_p && d == 0.) || (log_p && d == Double.NEGATIVE_INFINITY)) {
            return Dpq.DT0(lower_tail, log_p);
        }

        double pd = pdHyper(x, NR, NB, n, log_p);
        return log_p ? Dpq.DT_Log(d + pd, lower_tail) : Dpq.DLVal(d * pd, lower_tail);
    }

    public static double quantile(double p, double NR, double NB, double n)
    {
        return quantile(p, NR, NB, n, true, false);
    }

    public static double quantile(double p, double NR, double NB, double n, boolean lower_tail)
    {
        return quantile(p, NR, NB, n, lower_tail, false);
    }

    public static double quantile(double p, double NR, double NB, double n, boolean lower_tail, boolean log_p)
    {
        if (Double.isNaN(p) || Double.isNaN(NR) || Double.isNaN(NB) || Double.isNaN(n)) {
            return p + NR + NB + n;
        }
        if (Double.isInfinite(p) || Double.isInfinite(NR) || Double.isInfinite(NB) || Double.isInfinite(n)) {
            return Dpq.nanWarn();
        }

        NR = Math.round(NR);
        NB = Math.round(NB);
        double N = NR + NB;
        n = Math.round(n);
        if (NR < 0 || NB < 0 || n < 0 || n > N) {
            return Dpq.nanWarn();
        }

        double x_start = Math.max(0, n - NB);
        double x_end = Math.min(n, NR);

        if (log_p) {
            if (p > 0.0) {
                return Dpq.nanWarn();
            }
            if (p == 0.0) {
                return lower_tail ? x_end : x_start;
            }
            if (p == Double.NEGATIVE_INFINITY) {
                return lower_tail ? x_start : x_end;
            }
        } else {
            if (p < 0.0 || p > 1.0) {
                return Dpq.nanWarn();
            }
            if (p == 0.0) {
                return lower_tail ? x_start : x_end;
            }
            if (p == 1.0) {
                return lower_tail ? x_end : x_start;
            }
        }

        double xr = x_start;
        double xb = n - xr;
        boolean small_N = (N < 1000);
        double term = Choose.logFastChoose(NR, xr) + Choose.logFastChoose(NB, xb) - Choose.logFastChoose(N, n);
        if (small_N) {
            term = Math.exp(term);
        }

        NR -= xr;
        NB -= xb;
        if(!lower_tail || log_p) {
            p = Dpq.DTqIv(p, lower_tail, log_p);
        }
        p *= 1 - 1000 * Dpq.DBL_EPSILON;
        double sum = small_N ? term : Math.exp(term);

        while (sum < p && xr < x_end) {
            xr++;
            NB++;
            if (small_N)
                term *= (NR / xr) * (xb / NB);
            else
                term += Math.log((NR / xr) * (xb / NB));
            sum += small_N ? term : Math.exp(term);
            xb--;
            NR--;
        }
        return xr;
    }

    public static double rand(double nn1in, double nn2in, double kk_in)
    {
        int nn1, nn2, kk;
        int ix; // return value (coerced to double at the very end)
        boolean setup1, setup2;

        int ks = -1, n1s = -1, n2s = -1;
        int m = 0, min_jx = 0, max_jx = 0;
        int k = 0, n1 = 0, n2 = 0; // <- not allowing larger integer par
        double N = 0.0, w;
        double a, d, s, xl, xr, kl, kr;
        double lam_dl, lam_dr, p1, p2, p3;

        if(Double.isInfinite(nn1in) || Double.isInfinite(nn2in) || Double.isInfinite(kk_in)) {
            return Dpq.nanWarn();
        }

        nn1in = Math.round(nn1in);
        nn2in = Math.round(nn2in);
        kk_in  = Math.round(kk_in);

        if (nn1in < 0 || nn2in < 0 || kk_in < 0 || kk_in > nn1in + nn2in) {
            return Dpq.nanWarn();
        }
        if (nn1in >= Integer.MAX_VALUE || nn2in >= Integer.MAX_VALUE || kk_in >= Integer.MAX_VALUE) {
            if(kk_in == 1.0) {
                return Binomial.rand(kk_in, nn1in / (nn1in + nn2in));
            }
            return quantile(UniformRand.rand(), nn1in, nn2in, kk_in, false, false);
        }
        nn1 = (int)nn1in;
        nn2 = (int)nn2in;
        kk  = (int)kk_in;

        if (nn1 != n1s || nn2 != n2s) { // n1 | n2 is changed: setup all
            setup1 = true;	setup2 = true;
        } else if (kk != ks) { // n1 & n2 are unchanged: setup 'k' only
            setup1 = false;	setup2 = true;
        } else { // all three unchanged ==> no setup
            setup1 = false;	setup2 = false;
        }
        if (setup1) {
            N = nn1 + (double)nn2; // avoid int overflow
            if (nn1 <= nn2) {
                n1 = nn1; n2 = nn2;
            } else { // nn2 < nn1
                n1 = nn2; n2 = nn1;
            }
            // now have n1 <= n2
        }
        if (setup2) { // k
            if ((double)kk + kk >= N) { // this could overflow
                k = (int)(N - kk);
            } else {
                k = kk;
            }
        }
        if (setup1 || setup2) {
            m = (int) ((k + 1.) * (n1 + 1.) / (N + 2.)); // m := floor(adjusted mean E[.])
            min_jx = Math.max(0, k - n2);
            max_jx = Math.min(n1, k);

        }

        if (min_jx == max_jx) {
            ix = max_jx;
	        //goto L_finis; // return appropriate variate
            if ((double)kk + kk >= N) {
                if (nn1 > nn2) {
                    ix = kk - nn2 + ix;
                } else {
                    ix = nn1 - ix;
                }
            } else if (nn1 > nn2) {
                ix = kk - ix;
            }
            return ix;
        } else if (m - min_jx < 10) {
	        final double scale = 1e25;
	        final double con = 57.5646273248511421;
            double lw;
            if (k < n2) {
                lw = afc(n2) + afc(n1 + n2 - k) - afc(n2 - k) - afc(n1 + n2);
            } else {
                lw = afc(n1) + afc(     k     ) - afc(k - n2) - afc(n1 + n2);
            }
            w = Math.exp(lw + con);
            double p, u;

            L10:
            do {
                p = w;
                ix = min_jx;
                u = UniformRand.rand() * scale;

                while (u > p) {
                    u -= p;
                    p *= ((double) n1 - ix) * (k - ix);
                    ix++;
                    p = p / ix / (n2 - k + ix);

                    if (ix > max_jx)
		                continue L10;
                }
            } while (false);

        } else {

            double u,v;

            s = Math.sqrt((N - k) * k * n1 * n2 / (N - 1) / N / N);

            d = (int) (1.5 * s) + .5;
            xl = m - d + .5;
            xr = m + d + .5;
            a = afc(m) + afc(n1 - m) + afc(k - m) + afc(n2 - k + m);
            kl = Math.exp(a - afc((int) (xl)) - afc((int) (n1 - xl)) - afc((int) (k - xl))
                    - afc((int) (n2 - k + xl)));
            kr = Math.exp(a - afc((int) (xr - 1)) - afc((int) (n1 - xr + 1)) - afc((int) (k - xr + 1))
                    - afc((int) (n2 - k + xr - 1)));
            lam_dl = -Math.log(xl * (n2 - k + xl) / (n1 - xl + 1) / (k - xl + 1));
            lam_dr = -Math.log((n1 - xr + 1) * (k - xr + 1) / xr / (n2 - k + xr));
            p1 = d + d;
            p2 = p1 + kl / lam_dl;
            p3 = p2 + kr / lam_dr;

            int n_uv = 0;
            L30:
            do {
                u = UniformRand.rand() * p3;
                v = UniformRand.rand();
                n_uv++;
                if (n_uv >= 10000) {
                    System.out.printf("HyperGeometric.rand(*, n1=%d, n2=%d, k=%d): branch III: " +
                            "giving up after %d rejections\n", nn1, nn2, kk, n_uv);
                    return Dpq.nanWarn();
                }


                if (u < p1) {
                    ix = (int) (xl + u);
                } else if (u <= p2) {
                    ix = (int) (xl + Math.log(v) / lam_dl);
                    if (ix < min_jx)
                        continue;
                    v = v * (u - p1) * lam_dl;
                } else {
                    ix = (int) (xr - Math.log(v) / lam_dr);
                    if (ix > max_jx)
		                continue;
                    v = v * (u - p2) * lam_dr;
                }

                boolean reject = true;

                if (m < 100 || ix <= 50) {

                    int i;
                    double f = 1.0;
                    if (m < ix) {
                        for (i = m + 1; i <= ix; i++)
                            f = f * (n1 - i + 1) * (k - i + 1) / (n2 - k + i) / i;
                    } else if (m > ix) {
                        for (i = ix + 1; i <= m; i++)
                            f = f * i * (n2 - k + i) / (n1 - i + 1) / (k - i + 1);
                    }
                    if (v <= f) {
                        reject = false;
                    }
                } else {

                    final double del_tal = 0.0078;
                    final double del_tau = 0.0034;

                    double e, g, r, t, y;
                    double de, dg, dr, ds, dt, gl, gu, nk, nm, ub;
                    double xk, xm, xn, y1, ym, yn, yk, alv;

                    y = ix;
                    y1 = y + 1.0;
                    ym = y - m;
                    yn = n1 - y + 1.0;
                    yk = k - y + 1.0;
                    nk = n2 - k + y1;
                    r = -ym / y1;
                    s = ym / yn;
                    t = ym / yk;
                    e = -ym / nk;
                    g = yn * yk / (y1 * nk) - 1.0;
                    dg = 1.0;
                    if (g < 0.0)
                        dg = 1.0 + g;
                    gu = g * (1.0 + g * (-0.5 + g / 3.0));
                    gl = gu - .25 * (g * g * g * g) / dg;
                    xm = m + 0.5;
                    xn = n1 - m + 0.5;
                    xk = k - m + 0.5;
                    nm = n2 - k + xm;
                    ub = y * gu - m * gl + del_tau
                            + xm * r * (1. + r * (-0.5 + r / 3.0))
                            + xn * s * (1. + s * (-0.5 + s / 3.0))
                            + xk * t * (1. + t * (-0.5 + t / 3.0))
                            + nm * e * (1. + e * (-0.5 + e / 3.0));
                    // test against upper bound
                    alv = Math.log(v);
                    if (alv > ub) {
                        reject = true;
                    } else {
                        // test against lower bound
                        dr = xm * (r * r * r * r);
                        if (r < 0.0)
                            dr /= (1.0 + r);
                        ds = xn * (s * s * s * s);
                        if (s < 0.0)
                            ds /= (1.0 + s);
                        dt = xk * (t * t * t * t);
                        if (t < 0.0)
                            dt /= (1.0 + t);
                        de = nm * (e * e * e * e);
                        if (e < 0.0)
                            de /= (1.0 + e);
                        if (alv < ub - 0.25 * (dr + ds + dt + de)
                                + (y + m) * (gl - gu) - del_tal) {
                            reject = false;
                        } else {
                            reject = !(alv <= (a - afc(ix) - afc(n1 - ix) - afc(k - ix) - afc(n2 - k + ix)));
                        }
                    }
                }
                if (reject)
                    continue L30;
            } while (false);
        }


        //L_finis:

        if ((double)kk + kk >= N) {
            if (nn1 > nn2) {
                ix = kk - nn2 + ix;
            } else {
                ix = nn1 - ix;
            }
        } else if (nn1 > nn2) {
            ix = kk - ix;
        }

        return ix;
    }

    private static double pdHyper(double x, double NR, double NB, double n, boolean log_p)
    {
        /*
         * Calculate
         *
         *	    cdf(x, NR, NB, n, TRUE, FALSE)
         *   [log]  ----------------------------------
         *	       pdf(x, NR, NB, n, FALSE)
         *
         * without actually calling cdf.  This assumes that
         *
         *     x * (NR + NB) <= n * NR
         *
         */
        double sum = 0.0;
        double term = 1.0;

        while (x > 0 && term >= Dpq.DBL_EPSILON * sum) {
            term *= x * (NB - n + x) / (n + 1 - x) / (NR + 1 - x);
            sum += term;
            x--;
        }

        double ss = sum;
        return log_p ? Math.log1p(ss) : 1.0 + ss;
    }

    private static double afc(int i)
    {
        double[] al = {
            0.0, // ln(0!) = ln(1)
            0.0, // ln(1!) = ln(1)
            0.69314718055994530941723212145817, // ln(2)
            1.79175946922805500081247735838070, // ln(6)
            3.17805383034794561964694160129705, // ln(24)
            4.78749174278204599424770093452324,
            6.57925121201010099506017829290394,
            8.52516136106541430016553103634712
        };

        if (i < 0) {
            System.out.printf("HyperGeometric.rand(): afc(i), i=%d < 0 -- SHOULD NOT HAPPEN!\n", i);
            return -1;
        }
        if (i <= 7) {
            return al[i];
        }
        double i2 = (double)i * (double)i;
        return ((double)i + 0.5) * Math.log(i) - (double)i + Dpq.M_LN_SQRT_2PI +
                (0.0833333333333333 - 0.00277777777777778 / i2) / (double)i;
    }

    public static void main(String[] args)
    {
        System.out.println(pdf(8, 10, 7, 8));
        System.out.println(cdf(2, 10, 7, 8));
        System.out.println(quantile(0.975, 10, 7, 8));

        for (int i = 0; i < 10; i++) {
            System.out.println(rand(10, 7, 8));
        }
    }
}
