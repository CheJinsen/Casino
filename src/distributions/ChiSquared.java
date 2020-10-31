package distributions;

import distributions.detail.Dpq;

public class ChiSquared
{
    public static double pdf(double x, double df)
    {
        return pdf(x, df, false);
    }

    public static double pdf(double x, double df, boolean give_log)
    {
        return Gamma.pdf(x, df / 2.0, 2.0, give_log);
    }

    public static double cdf(double x, double df)
    {
        return cdf(x, df, true, false);
    }

    public static double cdf(double x, double df, boolean lower_tail)
    {
        return cdf(x, df, lower_tail, false);
    }

    public static double cdf(double x, double df, boolean lower_tail, boolean log_p)
    {
        return Gamma.cdf(x, df / 2.0, 2.0, lower_tail, log_p);
    }

    public static double quantile(double p, double df)
    {
        return quantile(p, df, true, false);
    }

    public static double quantile(double p, double df, boolean lower_tail)
    {
        return quantile(p, df, lower_tail, false);
    }

    public static double quantile(double p, double df, boolean lower_tail, boolean log_p)
    {
        return Gamma.quantile(p, 0.5 * df, 2.0, lower_tail, log_p);
    }

    public static double rand(double df)
    {
        if (Double.isInfinite(df) || df < 0.0) {
            return Dpq.nanWarn();
        }
        return Gamma.rand(df / 2.0, 2.0);
    }

    public static void main(String[] args)
    {
        System.out.println(pdf(10.12345, 10));
        System.out.println(cdf(10.12345, 10));
        System.out.println(quantile(0.92345, 0.5));

        for (int i = 0; i < 10; i++) {
            System.out.println(rand(10));
        }
    }
}
