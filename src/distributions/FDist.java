package distributions;

import distributions.detail.DistBase;
import distributions.detail.Dpq;

public class FDist extends DistBase
{
    public static double pdf(double x, double m, double n)
    {
        return pdf(x, m, n, false);
    }

    public static double pdf(double x, double m, double n, boolean give_log)
    {
        double dens;
        if (Double.isNaN(x) || Double.isNaN(m) || Double.isNaN(n)) {
            return x + m + n;
        }
        if (m <= 0.0 || n <= 0.0) {
            return Dpq.nanWarn();
        }
        if (x < 0.0) {
            return Dpq.D0(give_log);
        }
        if (x == 0.0) {
            return m > 2 ? Dpq.D0(give_log) : (m == 2.0 ? Dpq.D1(give_log) : Double.POSITIVE_INFINITY);
        }
        if (Double.isInfinite(m) && Double.isInfinite(n)) {
            if (x == 1.0)
                return Double.POSITIVE_INFINITY;
            else
                return Dpq.D0(give_log);
        }
        if (Double.isInfinite(n)) {
            return Gamma.pdf(x, m / 2.0, 2.0 / m, give_log);
        }
        if (m > 1e14) {
            dens = Gamma.pdf(1.0 / x, n/ 2.0, 2.0 / n, give_log);
            return give_log ? dens - 2.0 * Math.log(x) : dens / (x * x);
        }

        double f = 1.0 / (n + x * m);
        double q = n * f;
        double p = x * m * f;

        if (m >= 2.0) {
            f = m * q / 2.0;
            dens = binomialPdfRaw((m - 2.0) / 2.0, (m + n - 2.0) / 2.0, p, q, give_log);
        } else {
            f = m * m * q / (2.0 * p * (m + n));
            dens = binomialPdfRaw(m / 2.0, (m + n) / 2.0, p, q, give_log);
        }
        return give_log ? Math.log(f) + dens : f * dens;
    }

    public static double cdf(double x, double df1, double df2)
    {
        return cdf(x, df1, df2, true, false);
    }

    public static double cdf(double x, double df1, double df2, boolean lower_tail)
    {
        return cdf(x, df1, df2, lower_tail, false);
    }

    public static double cdf(double x, double df1, double df2, boolean lower_tail, boolean log_p)
    {
        if (Double.isNaN(x) || Double.isNaN(df1) || Double.isNaN(df2)) {
            return x + df1 + df2;
        }
        if (df1 <= 0.0 || df2 <= 0.0) {
            return Dpq.nanWarn();
        }

        double xMin = 0.0;
        double xMax = Double.POSITIVE_INFINITY;
        if (x <= xMin)  return Dpq.DT0(lower_tail, log_p);
        if (x >= xMax)  return Dpq.DT1(lower_tail, log_p);

        if (df2 == Double.POSITIVE_INFINITY) {
            if (df1 == Double.POSITIVE_INFINITY) {
                if (x < 1.0) return Dpq.DT0(lower_tail, log_p);
                if (x == 1.0) return log_p ? -Dpq.M_LN2 : 0.5;
                if (x > 1.0) return Dpq.DT1(lower_tail, log_p);
            }
        }

        if (df1 == Double.POSITIVE_INFINITY) {
            return ChiSquared.cdf(df2 / x, df2, !lower_tail, log_p);
        }

        if (df1 * x > df2) {
            x = Beta.cdf(df2 / (df2 + df1 * x), df2 / 2.0, df1 / 2.0, !lower_tail, log_p);
        } else {
            x = Beta.cdf(df1 * x / (df2 + df1 * x) , df1 / 2.0, df2 / 2.0, lower_tail, log_p);
        }
        return !Double.isNaN(x) ? x : Double.NaN;
    }

    public static double quantile(double p, double df1, double df2)
    {
        return quantile(p, df1, df2, true, false);
    }

    public static double quantile(double p, double df1, double df2, boolean lower_tail)
    {
        return quantile(p, df1, df2, lower_tail, false);
    }

    public static double quantile(double p, double df1, double df2, boolean lower_tail, boolean log_p)
    {
        if (Double.isNaN(p) || Double.isNaN(df1) || Double.isNaN(df2)) {
            return p + df1 + df2;
        }
        if (df1 <= 0.0 || df2 <= 0.0) {
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

        if (df1 <= df2 && df2 > 4e5) {
            if (Double.isInfinite(df1))
                return 1.0;
            return ChiSquared.quantile(p, df1, lower_tail, log_p) / df1;
        }
        if (df1 > 4e5) {
            return df2 / ChiSquared.quantile(p, df2, !lower_tail, log_p);
        }

        p = (1.0 / Beta.quantile(p, df2 / 2.0, df1 / 2.0, !lower_tail, log_p) - 1.0) * (df2 / df1);
        return !Double.isNaN(p) ? p : Double.NaN;
    }

    public static double rand(double df1, double df2)
    {
        if (Double.isNaN(df1) || Double.isNaN(df2) || df1 <= 0.0 || df2 <= 0.0) {
            return Dpq.nanWarn();
        }
        double v1 = Double.isFinite(df1) ? ChiSquared.rand(df1) / df1 : 1.0;
        double v2 = Double.isFinite(df2) ? ChiSquared.rand(df2) / df2 : 1.0;
        return v1 / v2;
    }

    public static void main(String[] args)
    {
        System.out.println(pdf(0.975, 80, 90));
        System.out.println(pdf(8, 10, 40, true));
        System.out.println(cdf(0.975, 80, 90));
        System.out.println(cdf(10.975, 80, 90));

        System.out.println(quantile(0.975, 97, 78));
        System.out.println(quantile(0.975, 97, 78, false, false));

        for (int i = 0; i < 10; i++) {
            System.out.println(rand(8, 9));
        }
    }
}
