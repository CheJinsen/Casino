package distributions;

import distributions.detail.DistBase;
import distributions.detail.Dpq;
import special_functions.LogGamma;

// Non central chis squared distribution
public class NChiSquared extends DistBase
{
    public static double pdf(double x, double df, double ncp)
    {
        return pdf(x, df, ncp, false);
    }

    public static double pdf(double x, double df, double ncp, boolean give_log)
    {
        final double eps = 5e-15;

        double i, ncp2, q, mid, df_mid = 0.0, imax;
        double sum, term;

        if (Double.isNaN(x) || Double.isNaN(df) || Double.isNaN(ncp)) {
            return x + df + ncp;
        }
        if (Double.isInfinite(df) || Double.isInfinite(ncp) || ncp < 0.0 || df < 0.0) {
            return Dpq.nanWarn();
        }
        if(x < 0) {
            return Dpq.D0(give_log);
        }
        if(x == 0 && df < 2.0) {
            return Double.POSITIVE_INFINITY;
        }
        if(ncp == 0) {
            return (df > 0) ? ChiSquared.pdf(x, df, give_log) : Dpq.D0(give_log);
        }
        if(x == Double.POSITIVE_INFINITY) {
            return Dpq.D0(give_log);
        }

        ncp2 = 0.5 * ncp;

        /* find max element of sum */
        imax = Math.ceil((-(2.0 + df) + Math.sqrt((2.0 - df) * (2.0 - df) + 4.0 * ncp * x)) / 4.0);
        if (imax < 0) {
            imax = 0;
        }
        if (Double.isFinite(imax)) {
            df_mid = df + 2 * imax;
            mid = poissonPdfRaw(imax, ncp2, false) * ChiSquared.pdf(x, df_mid, false);
        } else {
            mid = 0;
        }

        if(mid == 0) {
            if (give_log || ncp > 1000.) {
                double nl = df + ncp, ic = nl/(nl + ncp);
                return ChiSquared.pdf(x * ic, nl * ic, give_log);
            } else {
                return Dpq.D0(give_log);
            }
        }

        sum = mid;

        // upper tail
        term = mid; df = df_mid; i = imax;
        double x2 = x * ncp2;
        do {
            i++;
            q = x2 / i / df;
            df += 2;
            term *= q;
            sum += term;
        } while (q >= 1 || term * q > (1-q)*eps || term > 1e-10*sum);
        // lower tail
        term = mid; df = df_mid; i = imax;
        while (i != 0) {
            df -= 2;
            q = i * df / x2;
            i--;
            term *= q;
            sum += term;
            if (q < 1 && term * q <= (1-q)*eps)
                break;
        }
        return Dpq.DVal(sum, give_log);
    }

    public static double cdf(double x, double df, double ncp)
    {
        return cdf(x, df, ncp, true, false);
    }

    public static double cdf(double x, double df, double ncp, boolean lower_tail)
    {
        return cdf(x, df, ncp, lower_tail, false);
    }

    public static double cdf(double x, double df, double ncp, boolean lower_tail, boolean log_p)
    {
        if (Double.isNaN(x) || Double.isNaN(df) || Double.isNaN(ncp)) {
            return x + df + ncp;
        }
        if (Double.isInfinite(df) || Double.isInfinite(ncp)) {
            return Dpq.nanWarn();
        }
        if (df < 0.0 || ncp < 0.0) {
            return Dpq.nanWarn();
        }

        double ans = cdfRaw(x, df, ncp, 1e-12, 8.0 * Dpq.DBL_EPSILON,
                1000000, lower_tail, log_p);

        if (x <= 0.0 || x == Double.POSITIVE_INFINITY) {
            return ans;
        }
        if (ncp >= 80.0) {
            if (lower_tail) {
                ans = Math.min(ans, Dpq.D1(log_p));
            } else {
                if (ans < (log_p ? (-10.0 * Math.log(10.0)) : 1e-10))
                    System.out.println("Full precision may not have been achieved in NChiSquared.cdf()");
                if (!log_p && ans < 0.0)
                    ans = 0.0;
            }
        }

        if (!log_p || ans < -1e-8) {
            return ans;
        } else {
            ans = cdfRaw(x, df, ncp, 1e-12, 8.0 * Dpq.DBL_EPSILON,
                    1000000, !lower_tail, false);
            return Math.log1p(-ans);
        }
    }

    public static double quantile(double p, double df, double ncp)
    {
        return quantile(p, df, ncp, true, false);
    }

    public static double quantile(double p, double df, double ncp, boolean lower_tail)
    {
        return quantile(p, df, ncp, lower_tail, false);
    }

    public static double quantile(double p, double df, double ncp, boolean lower_tail, boolean log_p)
    {
        final double accu = 1e-13;
        final double r_acc = 4.0 * Dpq.DBL_EPSILON;
        final double eps = 1e-11;
        final double r_eps = 1e-10;

        double ux, lx, ux0, nx, pp;

        if (Double.isNaN(p) || Double.isNaN(df) || Double.isNaN(ncp)) {
            return p + df + ncp;
        }
        if (Double.isInfinite(df)) {
            return Dpq.nanWarn();
        }
        if (df < 0 || ncp < 0) {
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

        pp = Dpq.DqIv(p, log_p);
        if (pp > 1 - Dpq.DBL_EPSILON) {
            return lower_tail ? Double.POSITIVE_INFINITY : 0.0;
        }

        /* Invert cdf(.) :
         * 1. finding an upper and lower bound */
        {
       /* This is Pearson's (1959) approximation,
          which is usually good to 4 figs or so.  */
            double b, c, ff;
            b = (ncp * ncp) / (df + 3 * ncp);
            c = (df + 3 * ncp) / (df + 2 * ncp);
            ff = (df + 2 * ncp) / (c * c);
            ux = b + c * ChiSquared.quantile(p, ff, lower_tail, log_p);
            if(ux < 0) ux = 1;
            ux0 = ux;
        }

        if (!lower_tail && ncp >= 80) {
            if(pp < 1e-10) {
                System.out.println("full precision may not have been achieved in NChiSquared.cdf()");
            }
            p = log_p ? -Math.expm1(p) : (0.5 - (p) + 0.5);
            lower_tail = true;
        } else {
            p = pp;
        }

        pp = Math.min(1 - Dpq.DBL_EPSILON, p * (1 + eps));
        if (lower_tail) {
            while (ux < Double.MAX_VALUE &&
                    cdfRaw(ux, df, ncp, eps, r_eps, 10000, true, false) < pp) {
                ux *= 2.0;
            }
            pp = p * (1 - eps);
            lx = Math.min(ux0, Double.MAX_VALUE);
            while (lx > Double.MIN_VALUE &&
                        cdfRaw(lx, df, ncp, eps, r_eps, 10000, true, false) > pp) {
                lx *= 0.5;
            }
        } else {
            while (ux < Double.MAX_VALUE &&
                    cdfRaw(ux, df, ncp, eps, r_eps, 10000, false, false) > pp) {
                ux *= 2.0;
            }
            pp = p * (1 - eps);
            lx = Math.min(ux0, Double.MAX_VALUE);
            while (lx > Double.MIN_VALUE &&
                        cdfRaw(lx, df, ncp, eps, r_eps, 10000, false, false) < pp) {
                lx *= 0.5;
            }
        }

        if (lower_tail) {
            do {
                nx = 0.5 * (lx + ux);
                if (cdfRaw(nx, df, ncp, accu, r_acc, 100000, true, false) > p)
                    ux = nx;
                else
                    lx = nx;
            } while ((ux - lx) / nx > accu);
        } else {
            do {
                nx = 0.5 * (lx + ux);
                if (cdfRaw(nx, df, ncp, accu, r_acc, 100000, false, false) < p)
                    ux = nx;
                else
                    lx = nx;
            } while ((ux - lx) / nx > accu);
        }
        return 0.5 * (ux + lx);
    }

    public static double rand(double df, double lambda)
    {
        if (Double.isNaN(df) || Double.isInfinite(lambda) || df < 0.0 || lambda < 0.0) {
            return Dpq.nanWarn();
        }
        if (lambda == 0.0) {
            return df == 0.0 ? 0.0 : Gamma.rand(df / 2.0, 2.0);
        } else {
            double r = Poisson.rand(lambda / 2.0);
            if (r > 0.0)
                r = ChiSquared.rand(2.0 * r);
            if (df > 0.0)
                r += Gamma.rand(df / 2.0, 2.0);
            return r;
        }
    }

    private static double cdfRaw(double x, double f, double theta, double err_max, double rel_tol,
                                 int itr_max, boolean lower_tail, boolean log_p)
    {
        double lam, x2, f2, term, bound, f_x_2n, f_2n;
        double l_lam = -1.0, l_x = -1.0;
        boolean lamSml, tSml;
        double ans, u, v, t, lt, lu =-1;
        final double _L = (-0.5 * theta); // = -lambda
        final double _dbl_min_exp = Dpq.M_LN2 * Double.MIN_EXPONENT;

        if (x <= 0.0) {
            if(x == 0.0 && f == 0.0) { // chi^2_0(.) has point mass at zero
                return lower_tail ? Dpq.DExp(_L, log_p) : (log_p ? Dpq.log1Exp(_L) : -Math.expm1(_L));
            }
            return Dpq.DT0(lower_tail, log_p);
        }
        if (Double.isInfinite(x)) {
            return Dpq.DT1(lower_tail, log_p);
        }

        if (theta < 80) {
            if (lower_tail && f > 0.0 &&
                    Math.log(x) < Dpq.M_LN2 + 2 / f * (LogGamma.logGamma(f / 2.0 + 1.0) + _dbl_min_exp)) {

                double lambda = 0.5 * theta;
                double sum, sum2, pr = -lambda;
                sum = sum2 = Double.NEGATIVE_INFINITY;
                for(int i = 0; i < 110;  pr += Math.log(lambda) - Math.log(++i)) {
                    sum2 = logSpaceAdd(sum2, pr);
                    sum = logSpaceAdd(sum, pr + ChiSquared.cdf(x, f + 2 * i, true, true));
                    if (sum2 >= -1e-15)
                        break;
                }
                ans = sum - sum2;
                return log_p ? ans : Math.exp(ans);
            } else {
                double lambda = 0.5 * theta; // < 40
                double sum = 0, sum2 = 0, pr = Math.exp(-lambda);
                for(int i = 0; i < 110;  pr *= lambda/++i) {
                    sum2 += pr;
                    sum += pr * ChiSquared.cdf(x, f + 2 * i, lower_tail, false);
                    if (sum2 >= 1 - 1e-15)
                        break;
                }
                ans = sum / sum2;
                return log_p ? Math.log(ans) : ans;
            }
        }

        lam = 0.5 * theta; // = lambda = ncp/2
        lamSml = (-lam < _dbl_min_exp);
        if (lamSml) {
            u = 0;
            lu = -lam;
            l_lam = Math.log(lam);
        } else {
            u = Math.exp(-lam);
        }

        v = u;
        x2 = 0.5 * x;
        f2 = 0.5 * f;
        f_x_2n = f - x;

        if (f2 * Dpq.DBL_EPSILON > 0.125 && Math.abs(t = x2 - f2) < Math.sqrt(Dpq.DBL_EPSILON) * f2) {
            lt = (1 - t) * (2 - t / (f2 + 1)) - Dpq.M_LN_SQRT_2PI - 0.5 * Math.log(f2 + 1);
        } else {
            lt = f2 * Math.log(x2) -x2 - LogGamma.logGamma(f2 + 1.0);
        }

        tSml = (lt < _dbl_min_exp);
        if (tSml) {
            if (x > f + theta +  5* Math.sqrt(2 * (f + 2 * theta))) {
                return Dpq.DT1(lower_tail, log_p);
            }
            l_x = Math.log(x);
            ans = term = 0.; t = 0;
        } else {
            t = Math.exp(lt);
            ans = term = (v * t);
        }

        int n;
        for (n = 1, f_2n = f + 2.0, f_x_2n += 2.0; n <= itr_max ; n++, f_2n += 2, f_x_2n += 2) {
            if (f_x_2n > 0) {
                bound = t * x / f_x_2n;
                if (((bound <= err_max) && (term <= rel_tol * ans))) {
                    break;
                }
            }
            if (lamSml) {
                lu += l_lam - Math.log(n);
                if(lu >= _dbl_min_exp) {
                    v = u = Math.exp(lu);
                    lamSml = false;
                }
            } else {
                u *= lam / n;
                v += u;
            }
            if (tSml) {
                lt += l_x - Math.log(f_2n);
                if (lt >= _dbl_min_exp) {
                    t = Math.exp(lt);
                    tSml = false;
                }
            } else {
                t *= x / f_2n;
            }
            if (!lamSml && !tSml) {
                term = v * t;
                ans += term;
            }
        }

        if (n > itr_max) {
            System.out.printf("NChiSquared.cdf(x = %g, f = %g, theta = %g, ..): not converged in %d iterator.",
                    x, f, theta, itr_max);
        }

        return Dpq.DTVal(ans, lower_tail, log_p);
    }

    public static void main(String[] args)
    {
        System.out.println(pdf(12.87, 10, 20));
        System.out.println(cdf(12.87, 5.6, 0));
        System.out.println(quantile(0.887, 0.300, 200));

        System.out.println(rand(10, 3.4));
    }
}
