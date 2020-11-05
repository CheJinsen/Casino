package distributions;

import distributions.detail.DistBase;
import distributions.detail.Dpq;

public class NFDist extends DistBase // Non-central F Distribution
{
    public static double pdf(double x, double df1, double df2, double ncp)
    {
        return pdf(x, df1, df2, ncp, false);
    }

    public static double pdf(double x, double df1, double df2, double ncp, boolean give_log)
    {
        if (Double.isNaN(x) || Double.isNaN(df1) || Double.isNaN(df2) || Double.isNaN(ncp)) {
            return x + df2 + df1 + ncp;
        }
        if (df1 <= 0.0 || df2 <= 0.0 || ncp < 0.0) {
            return Dpq.nanWarn();
        }
        if (x < 0.0) {
            return Dpq.D0(give_log);
        }
        if (Double.isInfinite(ncp)) {
            return Dpq.nanWarn();
        }
        if (Double.isInfinite(df1) && Double.isInfinite(df2)) {
            if(x == 1.0)
                return Double.POSITIVE_INFINITY;
            else
                return Dpq.D0(give_log);
        }
        if (Double.isInfinite(df2)) {
            return df1 * NChiSquared.pdf(x * df1, df1, ncp, give_log);
        }

        double z;
        if (df1 > 1e14 && ncp < 1e7) {
            double f = 1 + ncp/df1;
            z = Gamma.pdf(1.0 / x / f, df2 / 2.0, 2.0 / df2, give_log);
            return give_log ? z - 2 * Math.log(x) - Math.log(f) : z / (x * x) / f;
        }

        double y = (df1 / df2) * x;
        z = NBeta.pdf(y / (1.0 + y), df1 / 2.0, df2 / 2.0, ncp, give_log);
        return  give_log ?
                z + Math.log(df1) - Math.log(df2) - 2.0 * Math.log1p(y) :
                z * (df1 / df2) / (1.0 + y) / (1.0 + y);
    }

    public static double cdf(double x, double df1, double df2, double ncp)
    {
        return cdf(x, df1, df2, ncp, true, false);
    }

    public static double cdf(double x, double df1, double df2, double ncp, boolean lower_tail)
    {
        return cdf(x, df1, df2, ncp, lower_tail, false);
    }

    public static double cdf(double x, double df1, double df2, double ncp, boolean lower_tail, boolean log_p)
    {
        if (Double.isNaN(x) || Double.isNaN(df1) || Double.isNaN(df2) || Double.isNaN(ncp)) {
            return x + df2 + df1 + ncp;
        }
        if (df1 <= 0.0 || df2 <= 0.0 || ncp < 0.0) {
            return Dpq.nanWarn();
        }
        if (Double.isInfinite(ncp)) {
            return Dpq.nanWarn();
        }
        if (Double.isInfinite(df1) && Double.isInfinite(df2)) {
            return Dpq.nanWarn();
        }
        if (x <= 0.0) {
            return Dpq.DT0(lower_tail, log_p);
        }
        if (x >= Double.POSITIVE_INFINITY) {
            return Dpq.DT1(lower_tail, log_p);
        }

        if (df2 > 1e8) {
            return NChiSquared.cdf(x * df1, df1, ncp, lower_tail, log_p);
        }

        double y = (df1 / df2) * x;
        return nonCentralBetaCdf2(y/(1. + y), 1./(1. + y), df1 / 2., df2 / 2.,
                ncp, lower_tail, log_p);
    }

    public static double quantile(double p, double df1, double df2, double ncp)
    {
        return quantile(p, df1, df2, ncp, true, false);
    }

    public static double quantile(double p, double df1, double df2, double ncp, boolean lower_tail)
    {
        return quantile(p, df1, df2, ncp, lower_tail, false);
    }

    public static double quantile(double p, double df1, double df2, double ncp, boolean lower_tail, boolean log_p)
    {
        if (Double.isNaN(p) || Double.isNaN(df1) || Double.isNaN(df2) || Double.isNaN(ncp)) {
            return p + df1 + df2 + ncp;
        }
        if (df1 <= 0.0 || df2 <= 0.0 || ncp < 0.0) {
            return Dpq.nanWarn();
        }
        if (Double.isInfinite(ncp)) {
            return Dpq.nanWarn();
        }
        if (Double.isInfinite(df1) && Double.isInfinite(df2)) {
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

        if (df2 > 1e8) {
            return NChiSquared.quantile(p, df1, ncp, lower_tail, log_p) / df1;
        }

        double y = NBeta.quantile(p, df1 / 2.0, df2 / 2.0, ncp, lower_tail, log_p);
        return y/(1-y) * (df2/df1);
    }

    public static double rand(double df1, double df2, double ncp)
    {
        if (Double.isNaN(df1) || Double.isInfinite(df2) || Double.isInfinite(ncp)) {
            return Dpq.nanWarn();
        }
        if (df1 < 0.0 || df2 < 0.0 || ncp < 0.0) {
            return Dpq.nanWarn();
        }

        return (NChiSquared.rand(df1, ncp) / df1) / (ChiSquared.rand(df2) / df2);
    }

    public static void main(String[] args)
    {
        System.out.println(pdf(8.67, 5, 67, 88));
        System.out.println(cdf(8.67, 5, 67, 8));
        System.out.println(quantile(0.867, 5, 67, 8));
        System.out.println();

        for (int i = 0; i < 10; i++) {
            System.out.println(rand(5, 67, 8));
        }
    }
}
