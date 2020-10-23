package distributions;

import distributions.detail.Dpq;
import random.UniformRand;

public class Uniform
{
    public static double pdf(double x, double a, double b)
    {
        return pdf(x, a, b, false);
    }

    public static double pdf(double x, double a, double b, boolean give_log)
    {
        if (Double.isNaN(x) || Double.isNaN(a) || Double.isNaN(b)) {
            return x + a + b;
        }
        if (b <= a) {
            return Dpq.nanWarn();
        }
        if (a <= x && x <= b) {
            return give_log ? -Math.log(b - a) : 1.0 / (b - a);
        }
        return Dpq.D0(give_log);
    }

    public static double cdf(double x, double a, double b)
    {
        return cdf(x, a, b, true, false);
    }

    public static double cdf(double x, double a, double b, boolean lower_tail)
    {
        return cdf(x, a, b, lower_tail, false);
    }

    public static double cdf(double x, double a, double b, boolean lower_tail, boolean log_p)
    {
        if (Double.isNaN(x) || Double.isNaN(a) || Double.isNaN(b)) {
            return x + a + b;
        }
        if (b < a) {
            return Dpq.nanWarn();
        }
        if (Double.isInfinite(a) || Double.isInfinite(b)) {
            return Dpq.nanWarn();
        }

        if (x >= b) {
            return Dpq.DT1(lower_tail, log_p);
        }
        if (x <= a) {
            return Dpq.DT0(lower_tail, log_p);
        }
        if (lower_tail) {
            return Dpq.DVal((x - a) / (b - a), log_p);
        } else {
            return Dpq.DVal((b - x) / (b - a), log_p);
        }
    }

    public static double quantile(double p, double a, double b)
    {
        return quantile(p, a, b, true, false);
    }

    public static double quantile(double p, double a, double b, boolean lower_tail)
    {
        return quantile(p, a, b, lower_tail, false);
    }

    public static double quantile(double p, double a, double b, boolean lowe_tail, boolean log_p)
    {
        if (Double.isNaN(p) || Double.isNaN(a) || Double.isNaN(b)) {
            return p + a + b;
        }
        if ((log_p && p > 0) || (!log_p && (p < 0 || p > 1))) {
            return Dpq.nanWarn();
        }
        if (Double.isInfinite(a) || Double.isInfinite(b)) {
            return Dpq.nanWarn();
        }
        if (b < a) {
            return Dpq.nanWarn();
        }
        if (b == a) {
            return a;
        }
        return a + Dpq.DTqIv(p, lowe_tail, log_p) * (b - a);
    }

    public static double rand()
    {
        return rand(0.0, 1.0);  // default: min = 0.0, max = 1.0;
    }

    public static double rand(double a, double b)
    {
        if (Double.isInfinite(a) || Double.isInfinite(b) || b < a) {
            return Dpq.nanWarn();
        }

        if (a == b) {
            return a;
        } else {
            double u;
            do {
                u = UniformRand.rand();
            } while (u <= 0.0 || u >= 1.0); // always true, but protect against user-supplied ones.
            return a + (b - a) * u;
        }
    }

    public static void main(String[] args)
    {
        int length = 1000;
        double sum = 0.0;
        for (int i = 0; i < length; i++) {
            double temp = rand();
//            System.out.println(temp);
            sum += temp;
        }
        System.out.println("Avg: " + sum / length);

        System.out.println(pdf(20, 10, 50));
        System.out.println(cdf(20, 10, 50));
        System.out.println(quantile(0.975, 3, 2));
        System.out.println(quantile(0.987, 2, 30.321));
    }
}
