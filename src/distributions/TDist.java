package distributions;

import distributions.detail.Dpq;
import random.NormalRand;
import special_functions.detail.Catherine;
import special_functions.detail.CosPI;

public class TDist
{
    public static double pdf(double x, double n)
    {
        return pdf(x, n, false);
    }

    public static double pdf(double x, double n, boolean give_log)
    {
        if (Double.isNaN(x) || Double.isNaN(n)) {
            return x + n;
        }
        if (n <= 0) {
            return Dpq.nanWarn();
        }
        if (Double.isInfinite(x)) {
            return Dpq.D0(give_log);
        }
        if (Double.isInfinite(n)) {
            return Normal.pdf(x, 0.0, 1.0, give_log);
        }

        double u;
        double t = -Catherine.bd0(n / 2.0,(n + 1.0) / 2.0) +
                Catherine.stirlingError((n + 1.0) / 2.0) - Catherine.stirlingError(n / 2.0);
        double x2n = x * x / n;
        double ax = 0.0;
        double l_x2n;
        boolean lrg_x2n = (x2n > 1.0 / Dpq.DBL_EPSILON);

        if (lrg_x2n) {
            ax = Math.abs(x);
            l_x2n = Math.log(ax) - Math.log(n) / 2.0;
            u = n * l_x2n;
        } else if (x2n > 0.2) {
            l_x2n = Math.log(1.0 + x2n) / 2.0;
            u = n * l_x2n;
        } else {
            l_x2n = Math.log1p(x2n) / 2.0;
            u = -Catherine.bd0(n / 2.0,(n + x * x) / 2.0) + x * x / 2.0;
        }

        if (give_log) {
            return t - u - (Dpq.M_LN_SQRT_2PI + l_x2n);
        }

        double I_sqrt_ = (lrg_x2n ? Math.sqrt(n)/ax : Math.exp(-l_x2n));
        return Math.exp(t-u) * Dpq.M_1_SQRT_2PI * I_sqrt_;
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
        double val, nx;

        if (Double.isNaN(x) || Double.isNaN(n)) {
            return x + n;
        }
        if (n <= 0.0) {
            return Dpq.nanWarn();
        }
        if (Double.isInfinite(x)) {
            return (x < 0) ? Dpq.DT0(lower_tail, log_p) : Dpq.DT1(lower_tail, log_p);
        }
        if (Double.isInfinite(n)) {
            return Normal.cdf(x, 0.0, 1.0, lower_tail, log_p);
        }

        nx = 1 + (x/n)*x;
        if (nx > 1e100) {
            double l_val = -0.5 * n * (2 * Math.log(Math.abs(x)) - Math.log(n))
                    -special_functions.Beta.logBeta(0.5 * n, 0.5) - Math.log(0.5 * n);
            val = log_p ? l_val : Math.exp(l_val);
        } else {
            val = (n > x * x)
                    ? Beta.cdf(x * x / (n + x * x), 0.5, n / 2.0, false, log_p)
                    : Beta.cdf (1.0 / nx,n / 2.0, 0.5, true, log_p);
        }
        if (x <= 0.) {
            lower_tail = !lower_tail;
        }
        if (log_p) {
            if(lower_tail) {
                return Math.log1p(-0.5 * Math.exp(val));
            } else {
                return val - Dpq.M_LN2;
            }
        } else {
            val /= 2.0;
            return Dpq.DCVal(val, lower_tail);
        }
    }

    public static double quantile(double p, double ndf)
    {
        return quantile(p, ndf, true, false);
    }

    public static double quantile(double p, double ndf, boolean lower_tail)
    {
        return quantile(p, ndf, lower_tail, false);
    }

    public static double quantile(double p, double ndf, boolean lower_tail, boolean log_p)
    {
        final double eps = 1.e-12;
        double P, q;

        if (Double.isNaN(p) || Double.isNaN(ndf)) {
            return p + ndf;
        }

        double r = Double.POSITIVE_INFINITY;
        double l = Double.NEGATIVE_INFINITY;
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

        if (ndf <= 0) {
            return Dpq.nanWarn();
        }

        if (ndf < 1) {
            final double accu = 1e-13;
	        final double Eps = 1e-11;

            double ux, lx, nx, pp;

            int it = 0;

            p = Dpq.DTqIv(p, lower_tail, log_p);

            if (p > 1 - Dpq.DBL_EPSILON) {
                return Double.POSITIVE_INFINITY;
            }

            pp = Math.min(1.0 - Dpq.DBL_EPSILON, p * (1.0 + Eps));
            ux = 1.0;
            while (ux < Double.MAX_VALUE && cdf(ux, ndf, true, false) < pp) {
                ux *= 2.0;
            }

            pp = p * (1 - Eps);
            lx = -1.0;
            while (lx > -Double.MAX_VALUE && cdf(lx, ndf, true, false) > pp) {
                lx *= 2.0;
            }

            do {
                nx = 0.5 * (lx + ux);
                if (cdf(nx, ndf, true, false) > p)
                    ux = nx;
                else
                    lx = nx;
            } while ((ux - lx) / Math.abs(nx) > accu && ++it < 1000);

            if (it >= 1000) {
                System.out.println("Full precision may not have been achieved in TDist.cdf()");
            }

            return 0.5 * (lx + ux);
        }

        if (ndf > 1e20) {
            return Normal.quantile(p, 0.0, 1.0, lower_tail, log_p);
        }

        P = Dpq.DqIv(p, log_p);

        boolean neg = (!lower_tail || P < 0.5) && (lower_tail || P > 0.5);
        boolean is_neg_lower = (lower_tail == neg);

        if (neg) {
            P = 2 * (log_p ? (lower_tail ? P : -Math.expm1(p)) : Dpq.DLVal(p, lower_tail));
        } else {
            P = 2 * (log_p ? (lower_tail ? -Math.expm1(p) : P) : Dpq.DCVal(p, lower_tail));
        }

        if (Math.abs(ndf - 2) < eps) {
            if (P > Double.MIN_VALUE) {
                if (3 * P < Dpq.DBL_EPSILON)
                    q = 1 / Math.sqrt(P);
                else if (P > 0.9)
                    q = (1 - P) * Math.sqrt(2 /(P * (2 - P)));
                else
                    q = Math.sqrt(2 / (P * (2 - P)) - 2);
            } else {
                if(log_p)
                    q = is_neg_lower ? Math.exp(- p/2) / Math.sqrt(2.0) : 1 / Math.sqrt(-Math.expm1(p));
                else
                    q = Double.POSITIVE_INFINITY;
            }
        } else if (ndf < 1 + eps) {
            if(P == 1.0) {
                q = 0;
            } else if(P > 0) {
                q = 1 / CosPI.tanPI(P / 2.0);
            } else {
                if (log_p)
                    q = is_neg_lower ? 1.0 / Math.PI * Math.exp(-p) : -1.0 / (Math.PI * Math.expm1(p));
                else
                    q = Double.POSITIVE_INFINITY;
            }
        } else {
            double x = 0.0, y = 0.0, log_P2 = 0.0;
            double a = 1 / (ndf - 0.5);
            double b = 48 / (a * a);
            double c = ((20700 * a / b - 98) * a - 16) * a + 96.36;
            double d = ((94.5 / (b + c) - 3) / b + 1) * Math.sqrt(a * Math.PI / 2.0) * ndf;

            boolean P_ok1 = P > Double.MIN_VALUE || !log_p,  P_ok = P_ok1;
            if (P_ok1) {
                y = Math.pow(d * P, 2.0 / ndf);
                P_ok = (y >= Dpq.DBL_EPSILON);
            }
            if (!P_ok) {
                log_P2 = is_neg_lower ? Dpq.DLog(p, log_p) : Dpq.DLExp(p, log_p);
                x = (Math.log(d) + Dpq.M_LN2 + log_P2) / ndf;
                y = Math.exp(2 * x);
            }

            if ((ndf < 2.1 && P > 0.5) || y > 0.05 + a) {
                if(P_ok) {
                    x = Normal.quantile(0.5 * P, 0.0, 1.0, true, false);
                } else {
                    x = Normal.quantile(log_P2, 0.0, 1.0, lower_tail, true);
                }

                y = x * x;
                if (ndf < 5) {
                    c += 0.3 * (ndf - 4.5) * (x + 0.6);
                }

                c = (((0.05 * d * x - 5) * x - 7) * x - 2) * x + b + c;
                y = (((((0.4 * y + 6.3) * y + 36) * y + 94.5) / c - y - 3) / b + 1) * x;
                y = Math.expm1(a * y * y);
                q = Math.sqrt(ndf * y);

            } else if (!P_ok && x < - Dpq.M_LN2 * 53) {
                q = Math.sqrt(ndf) * Math.exp(-x);
            } else {
                y = ((1 / (((ndf + 6) / (ndf * y) - 0.089 * d - 0.822) * (ndf + 2) * 3) + 0.5 / (ndf + 4))
                        * y - 1) * (ndf + 1) / (ndf + 2) + 1 / y;
                q = Math.sqrt(ndf * y);
            }

            if (P_ok1) {
                int it = 0;
                while (it++ < 10 && (y = pdf(q, ndf, false)) > 0 &&
                        Double.isFinite(x = (cdf(q, ndf, false, false) - P / 2) / y) &&
                        Math.abs(x) > 1e-14 * Math.abs(q)) {
                    q += x * (1.0 + x * q * (ndf + 1) / (2 * (q * q + ndf)));
                }
            }
        }
        if (neg) q = -q;
        return q;
    }

    public static double rand(double df)
    {
        if (Double.isNaN(df) || df <= 0.0) {
            return Dpq.nanWarn();
        }
        if (Double.isInfinite(df)) {
            return NormalRand.rand();
        } else {
            double num = NormalRand.rand();
            return num / Math.sqrt(ChiSquared.rand(df) / df);
        }
    }

    public static void main(String[] args)
    {
        System.out.println(pdf(-0.98, 56));
        System.out.println(pdf(9, 7.8, true));

        System.out.println(cdf(-3, 199) * 2);
        System.out.println(cdf(1.58, 9, false));

        System.out.println(quantile(0.05, 99));
        System.out.println(quantile(0.95, 23));

        System.out.println("---------random number--------");
        for (int i = 0; i < 10; i++) {
            System.out.println(rand(0.975));
        }
    }
}
