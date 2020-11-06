package distributions;

import distributions.detail.Dpq;
import special_functions.LogGamma;

public class NTDist
{
    public static double pdf(double x, double df, double ncp)
    {
        return pdf(x, df, ncp, false);
    }

    public static double pdf(double x, double df, double ncp, boolean give_log)
    {
        if (Double.isNaN(x) || Double.isNaN(df)) {
            return x + df;
        }
        if (df <= 0.0) {
            return Dpq.nanWarn();
        }
        if (ncp == 0.0) {
            return TDist.pdf(x, df, give_log);
        }
        if (Double.isInfinite(x)) {
            return Dpq.D0(give_log);
        }
        if (Double.isInfinite(df) || df > 1e8) {
            return Normal.pdf(x, ncp, 1., give_log);
        }

        double u;
        if (Math.abs(x) > Math.sqrt(df * Dpq.DBL_EPSILON)) {
            u = Math.log(df) - Math.log(Math.abs(x)) +
                    Math.log(Math.abs(cdf(x * Math.sqrt((df+2)/df), df+2, ncp, true, false) -
                            cdf(x, df, ncp, true, false)));
        } else {
            u = LogGamma.logGamma((df + 1.0) / 2.0) - LogGamma.logGamma(df / 2.0) -
                    (Dpq.M_LN_SQRT_PI + 0.5 * (Math.log(df) + ncp * ncp));
        }

        return give_log ? u : Math.exp(u);
    }

    public static double cdf(double t, double df, double ncp)
    {
        return cdf(t, df, ncp, true, false);
    }

    public static double cdf(double t, double df, double ncp, boolean lower_tail)
    {
        return cdf(t, df, ncp, lower_tail, false);
    }

    public static double cdf(double t, double df, double ncp, boolean lower_tail, boolean log_p)
    {
        final int itrMax = 1000;
        final double errMax = 1.e-12;

        if (df <= 0.0) {
            return Dpq.nanWarn();
        }
        if(ncp == 0.0) {
            return TDist.cdf(t, df, lower_tail, log_p);
        }
        if (Double.isInfinite(t)) {
            return t < 0 ? Dpq.DT0(lower_tail, log_p) : Dpq.DT1(lower_tail, log_p);
        }

        boolean negDel;
        double tt, del;
        if (t >= 0.) {
            negDel = false; tt = t;	 del = ncp;
        } else {
            if (ncp > 40 && (!log_p || !lower_tail))
                return Dpq.DT0(lower_tail, log_p);
            negDel = true;	tt = -t; del = -ncp;
        }

        double s;
        if (df > 4e5 || del * del > 2 * Dpq.M_LN2 * (-Double.MIN_EXPONENT)) {
            s = 1.0 / (4.0 * df);
            return Normal.cdf(tt * (1.0 - s), del, Math.sqrt(1.0 + tt * tt * 2.0 * s),
                    lower_tail != negDel, log_p);
        }

        double tnc;
        double x = t * t;
        double rxb = df / (x + df);
        x = x / (x + df);

        if (x > 0.0) {
            double lambda = del * del;
            double p = 0.5 * Math.exp(-0.5 * lambda);

            if (p == 0.0) {
                System.out.println("Underflow occurred in NTDist.cdf()");
                System.out.println("Value out of range in NTDist.cdf()");
                return Dpq.DT0(lower_tail, log_p);
            }

            double q = Dpq.M_SQRT_2dPI * p * del;
            s = 0.5 - p;
            if (s < 1e-7) {
                s = -0.5 * Math.expm1(-0.5 * lambda);
            }
            double a = 0.5;
            double b = 0.5 * df;

            rxb = Math.pow(rxb, b);
            double alBeta = Dpq.M_LN_SQRT_PI + LogGamma.logGamma(b) - LogGamma.logGamma(0.5 + b);
            double xOdd = Beta.cdf(x, a, b, true, false);
            double gOdd = 2.0 * rxb * Math.exp(a * Math.log(x) - alBeta);
            tnc = b * x;
            double xEven = (tnc < Dpq.DBL_EPSILON) ? tnc : 1.0 - rxb;
            double gEven = tnc * rxb;
            tnc = p * xOdd + q * xEven;

            for(int it = 1; it <= itrMax; it++) {
                a += 1.0;
                xOdd  -= gOdd;
                xEven -= gEven;
                gOdd  *= x * (a + b - 1.0) / a;
                gEven *= x * (a + b - 0.5) / (a + 0.5);
                p *= lambda / (2 * it);
                q *= lambda / (2 * it + 1);
                tnc += p * xOdd + q * xEven;
                s -= p;

                if(s < -1.e-10) {
                    System.out.println("Full precision may not have been achieved in NTDist.cdf()");
                    //finis:
                    tnc += Normal.cdf(-del, 0.0, 1.0, true, false);
                    lower_tail = lower_tail != negDel;
                    if(tnc > 1 - 1e-10 && lower_tail) {
                        System.out.println("Full precision may not have been achieved in NTDist.cdf()");
                    }
                    return Dpq.DTVal(Math.min(tnc, 1.0), lower_tail, log_p);
                }
                if (s <= 0 && it > 1) {
                    //finis:
                    tnc += Normal.cdf(-del, 0.0, 1.0, true, false);
                    lower_tail = lower_tail != negDel;
                    if(tnc > 1 - 1e-10 && lower_tail) {
                        System.out.println("Full precision may not have been achieved in NTDist.cdf()");
                    }
                    return Dpq.DTVal(Math.min(tnc, 1.0), lower_tail, log_p);
                }

                double err_bd = 2.0 * s * (xOdd - gOdd);
                if (Math.abs(err_bd) < errMax) {
                    //finis:
                    tnc += Normal.cdf(-del, 0.0, 1.0, true, false);
                    lower_tail = lower_tail != negDel;
                    if(tnc > 1 - 1e-10 && lower_tail) {
                        System.out.println("Full precision may not have been achieved in NTDist.cdf()");
                    }
                    return Dpq.DTVal(Math.min(tnc, 1.0), lower_tail, log_p);
                }
            }
            System.out.println("Convergence failed in NTDist.cdf()");
        } else {
            tnc = 0.0;
        }

        //finis:
        tnc += Normal.cdf(-del, 0.0, 1.0, true, false);
        lower_tail = lower_tail != negDel;
        if(tnc > 1 - 1e-10 && lower_tail) {
            System.out.println("Full precision may not have been achieved in NTDist.cdf()");
        }
        return Dpq.DTVal(Math.min(tnc, 1.0), lower_tail, log_p);
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
        final double Eps = 1e-11;

        double ux, lx, nx, pp;

        if (Double.isNaN(p) || Double.isNaN(df) || Double.isNaN(ncp)) {
            return p + df + ncp;
        }
        if (df <= 0.0) {
            return Dpq.nanWarn();
        }
        if (ncp == 0.0 && df >= 1.0) {
            return TDist.quantile(p, df, lower_tail, log_p);
        }

        double r = Double.NEGATIVE_INFINITY;
        double l = Double.POSITIVE_INFINITY;
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

        if (Double.isInfinite(df)) {
            return Normal.quantile(p, ncp, 1.0, lower_tail, log_p);
        }

        p = Dpq.DTqIv(p, lower_tail, log_p);
        if (p > 1 - Dpq.DBL_EPSILON) {
            return Double.POSITIVE_INFINITY;
        }
        pp = Math.min(1 - Dpq.DBL_EPSILON, p * (1 + Eps));
        ux = Math.max(1.0, ncp);
        while (ux < Double.MAX_VALUE && cdf(ux, df, ncp, true, false) < pp) {
            ux *= 2.0;
        }
        pp = p * (1 - Eps);
        lx = Math.min(-1.0, -ncp);
        while (lx > -Double.MAX_VALUE && cdf(lx, df, ncp, true, false) > pp) {
            lx *= 2;
        }

        do {
            nx = 0.5 * (lx + ux); // could be zero
            if (cdf(nx, df, ncp, true, false) > p)
                ux = nx;
            else
                lx = nx;
        } while ((ux - lx) > accu * Math.max(Math.abs(lx), Math.abs(ux)));

        return 0.5 * (lx + ux);
    }

    public static double rand(double df, double ncp)
    {
        if (Double.isInfinite(df) || Double.isInfinite(ncp) || df < 0.0 || ncp < 0.0) {
            return Dpq.nanWarn();
        }
        return Normal.rand(ncp, 1.0) / Math.sqrt(ChiSquared.rand(df) / df);
    }
    public static void main(String[] args)
    {
        System.out.println(pdf(89.07, 56, 78));
        System.out.println(cdf(89.07, 56, 78));
        System.out.println(quantile(0.8907, 56, 78));
        System.out.println();

        for (int i = 0; i < 10; i++) {
            System.out.println(rand(56, 78));
        }
    }
}
