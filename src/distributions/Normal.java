package distributions;

import distributions.detail.Dpq;
import distributions.detail.RefDouble;
import random.NormalRand;

public class Normal
{
    public static double pdf(double x)
    {
        return pdf(x, 0.0, 1.0, false);
    }

    public static double pdf(double x, double mu, double sigma)
    {
        return pdf(x, mu, sigma, false);
    }

    public static double pdf(double x, double mu, double sigma, boolean give_log)
    {
        if (Double.isNaN(x) || Double.isNaN(mu) || Double.isNaN(sigma)) {
            return x + mu + sigma;
        }
        if (sigma < 0.0) {
            return Dpq.nanWarn();
        }
        if (Double.isInfinite(sigma)) {
            return Dpq.D0(give_log);
        }
        if (Double.isInfinite(x) && mu == x) {
            return Double.NaN;
        }
        if (sigma == 0.0) {
            return (x == mu) ? Double.POSITIVE_INFINITY : Dpq.D0(give_log);
        }

        x = (x - mu) / sigma;
        if (Double.isInfinite(x)) {
            return Dpq.D0(give_log);
        }

        x = Math.abs(x);
        if (x >= 2.0 * Math.sqrt(Double.MAX_VALUE)) {
            return Dpq.D0(give_log);
        }
        if (give_log) {
            return -(Math.log(Math.sqrt(2.0 * Math.PI)) + 0.5 * x * x + Math.log(sigma));
        }
        return 1.0 / Math.sqrt(2.0 * Math.PI) * Math.exp(-0.5 * x * x) / sigma;
    }

    public static double cdf(double x)
    {
        return cdf(x, 0.0, 1.0, true, false);
    }

    public static double cdf(double x, double mu, double sigma)
    {
        return cdf(x, mu, sigma, true, false);
    }

    public static double cdf(double x, double mu, double sigma, boolean lower_tail)
    {
        return cdf(x, mu, sigma, lower_tail, false);
    }

    public static double cdf(double x, double mu, double sigma, boolean lower_tail, boolean log_p)
    {
        if (Double.isNaN(x) || Double.isNaN(mu) || Double.isNaN(sigma)) {
            return x + mu + sigma;
        }
        if (Double.isInfinite(x) && mu == x) {
            return Double.NaN;
        }
        if (sigma <= 0.0) {
            if (sigma < 0.0)
                return Dpq.nanWarn();
            return (x < mu) ? Dpq.DT0(lower_tail, log_p) : Dpq.DT1(lower_tail, log_p);
        }

        RefDouble p = new RefDouble(0.0);
        RefDouble cp = new RefDouble(0.0);
        p.val = (x - mu) / sigma;
        if (Double.isInfinite(p.val)) {
            return (x < mu) ? Dpq.DT0(lower_tail, log_p) : Dpq.DT1(lower_tail, log_p);
        }

        x = p.val;
        cdfBoth(x, p, cp, (lower_tail ? 0 : 1), log_p);
        return lower_tail ? p.val : cp.val;
    }

    public static void cdfBoth(double x, RefDouble cum, RefDouble c_cum, int i_tail, boolean log_p)
    {
        double[] a = {
            2.2352520354606839287,
            161.02823106855587881,
            1067.6894854603709582,
            18154.981253343561249,
            0.065682337918207449113
        };
        double[] b = {
            47.20258190468824187,
            976.09855173777669322,
            10260.932208618978205,
            45507.789335026729956
        };
        double[] c = {
            0.39894151208813466764,
            8.8831497943883759412,
            93.506656132177855979,
            597.27027639480026226,
            2494.5375852903726711,
            6848.1904505362823326,
            11602.651437647350124,
            9842.7148383839780218,
            1.0765576773720192317e-8
        };
        double[] d = {
            22.266688044328115691,
            235.38790178262499861,
            1519.377599407554805,
            6485.558298266760755,
            18615.571640885098091,
            34900.952721145977266,
            38912.003286093271411,
            19685.429676859990727
        };
        double[] p = {
            0.21589853405795699,
            0.1274011611602473639,
            0.022235277870649807,
            0.001421619193227893466,
            2.9112874951168792e-5,
            0.02307344176494017303
        };
        double[] q = {
            1.28426009614491121,
            0.468238212480865118,
            0.0659881378689285515,
            0.00378239633202758244,
            7.29751555083966205e-5
        };

        if(Double.isNaN(x)) {
            cum.val = x;
            c_cum.val = x;
            return;
        }

        double eps = Dpq.DBL_EPSILON * 0.5;
        boolean lower = i_tail != 1;
        boolean upper = i_tail != 0;

        double xsq, x_num, x_den, del;
        double yes = Math.abs(x);
        if (yes <= 0.67448975) {
            if (yes > eps) {
                xsq = x * x;
                x_num = a[4] * xsq;
                x_den = xsq;
                for (int i = 0; i < 3; ++i) {
                    x_num = (x_num + a[i]) * xsq;
                    x_den = (x_den + b[i]) * xsq;
                }
            } else {
                x_num = x_den = 0.0;
            }

            double temp = x * (x_num + a[3]) / (x_den + b[3]);
            if (lower)  cum.val = 0.5 + temp;
            if (upper) c_cum.val = 0.5 - temp;
            if (log_p) {
                if (lower)  cum.val = Math.log(cum.val);
                if (upper) c_cum.val = Math.log(c_cum.val);
            }
        } else if (yes <= Math.sqrt(32.0)) {
            x_num = c[8] * yes;
            x_den = yes;
            for (int i = 0; i < 7; ++i) {
                x_num = (x_num + c[i]) * yes;
                x_den = (x_den + d[i]) * yes;
            }

            double temp = (x_num + c[7]) / (x_den + d[7]);
            xsq = (int)(yes * 16) / 16.0;
            del = (yes - xsq) * (yes + xsq);
            if (log_p) {
                cum.val = (-xsq * xsq * 0.5) + (-del * 0.5) + Math.log(temp);
                if ((lower && x > 0.) || (upper && x <= 0.0))
                    c_cum.val = Math.log1p(-Math.exp(-xsq * xsq * 0.5) * Math.exp(-del * 0.5) * temp);
            } else {
                cum.val = Math.exp(-xsq * xsq * 0.5) * Math.exp(-del * 0.5) * temp;
                c_cum.val = 1.0 - cum.val;
            }

            if (x > 0.0) {
                temp = cum.val;
                if(lower)
                    cum.val = c_cum.val;
                c_cum.val = temp;
            }
        } else if ((log_p && yes < 1e170) || (lower && -37.5193 < x  &&  x < 8.2924)
                || (upper && -8.2924  < x  &&  x < 37.5193)) {

            xsq = 1.0 / (x * x);
            x_num = p[5] * xsq;
            x_den = xsq;
            for (int i = 0; i < 4; ++i) {
                x_num = (x_num + p[i]) * xsq;
                x_den = (x_den + q[i]) * xsq;
            }
            double temp = xsq * (x_num + p[4]) / (x_den + q[4]);
            temp = (1.0 / Math.sqrt(2.0 * Math.PI) - temp) / yes;

            xsq = (int)(x * 16) / 16.0;
            del = (x - xsq) * (x + xsq);
            if (log_p) {
                cum.val = (-xsq * xsq * 0.5) + (-del * 0.5) + Math.log(temp);
                if ((lower && x > 0.) || (upper && x <= 0.0))
                    c_cum.val = Math.log1p(-Math.exp(-xsq * xsq * 0.5) * Math.exp(-del * 0.5) * temp);
            } else {
                cum.val = Math.exp(-xsq * xsq * 0.5) * Math.exp(-del * 0.5) * temp;
                c_cum.val = 1.0 - cum.val;
            }

            if (x > 0.0) {
                temp = cum.val;
                if(lower)
                    cum.val = c_cum.val;
                c_cum.val = temp;
            }
        } else {
            if (x > 0) {
                cum.val = Dpq.D1(log_p);
                c_cum.val = Dpq.D0(log_p);
            } else {
                cum.val = Dpq.D0(log_p);
                c_cum.val = Dpq.D1(log_p);
            }
        }
    }

    public static double quantile(double p)
    {
        return quantile(p, 0.0, 1.0, true, false);
    }

    public static double quantile(double p, double mu, double sigma)
    {
        return quantile(p, mu, sigma, true, false);
    }

    public static double quantile(double p, double mu, double sigma, boolean lower_tail)
    {
        return quantile(p, mu, sigma, lower_tail, false);
    }

    public static double quantile(double p, double mu, double sigma, boolean lower_tail, boolean log_p)
    {
        if (Double.isNaN(p) || Double.isNaN(mu) || Double.isNaN(sigma)) {
            return p + mu + sigma;
        }

        // Q_P01_boundaries
        if (log_p) {
            if (p > 0.0) {
                return Dpq.nanWarn();
            }
            if (p == 0.0) {
                return lower_tail ? Double.POSITIVE_INFINITY : Double.NEGATIVE_INFINITY;
            }
            if (p == Double.NEGATIVE_INFINITY) {
                return lower_tail ? Double.NEGATIVE_INFINITY : Double.POSITIVE_INFINITY;
            }
        } else {
            if (p < 0.0 || p > 1.0) {
                return Dpq.nanWarn();
            }
            if (p == 0.0) {
                return lower_tail ? Double.NEGATIVE_INFINITY : Double.POSITIVE_INFINITY;
            }
            if (p == 1.0) {
                return lower_tail ? Double.POSITIVE_INFINITY : Double.NEGATIVE_INFINITY;
            }
        }

        if (sigma < 0.0) {
            return Dpq.nanWarn();
        }
        if (sigma == 0.0) {
            return mu;
        }

        double p_ = Dpq.DTqIv(p, lower_tail, log_p);
        double q = p_ - 0.5;

        double r, val;
        if (Math.abs(q) <= .425) {
            r = 0.180625 - q * q;
            val = q * (((((((r * 2509.0809287301226727 +
                33430.575583588128105) * r + 67265.770927008700853) * r +
                45921.953931549871457) * r + 13731.693765509461125) * r +
                1971.5909503065514427) * r + 133.14166789178437745) * r +
                3.387132872796366608)
                / (((((((r * 5226.495278852854561 +
                28729.085735721942674) * r + 39307.89580009271061) * r +
                21213.794301586595867) * r + 5394.1960214247511077) * r +
                687.1870074920579083) * r + 42.313330701600911252) * r + 1.);
        } else {
            if (q > 0)
                r = Dpq.DTCIv(p, lower_tail, log_p);
            else
                r = p_;

            r = Math.sqrt(- ((log_p &&
                    ((lower_tail && q <= 0) || (!lower_tail && q > 0))) ?
                    p : Math.log(r)));

            if (r <= 5.) {
                r += -1.6;
                val = (((((((r * 7.7454501427834140764e-4 +
                    .0227238449892691845833) * r + .24178072517745061177) *
                    r + 1.27045825245236838258) * r +
                    3.64784832476320460504) * r + 5.7694972214606914055) *
                    r + 4.6303378461565452959) * r +
                    1.42343711074968357734)
                    / (((((((r *
                    1.05075007164441684324e-9 + 5.475938084995344946e-4) *
                    r + .0151986665636164571966) * r +
                    .14810397642748007459) * r + .68976733498510000455) *
                    r + 1.6763848301838038494) * r +
                    2.05319162663775882187) * r + 1.);
            } else {
                r += -5.;
                val = (((((((r * 2.01033439929228813265e-7 +
                    2.71155556874348757815e-5) * r +
                    .0012426609473880784386) * r + .026532189526576123093) *
                    r + .29656057182850489123) * r +
                    1.7848265399172913358) * r + 5.4637849111641143699) *
                    r + 6.6579046435011037772)
                    / (((((((r *
                    2.04426310338993978564e-15 + 1.4215117583164458887e-7)*
                    r + 1.8463183175100546818e-5) * r +
                    7.868691311456132591e-4) * r + .0148753612908506148525)
                    * r + .13692988092273580531) * r +
                    .59983220655588793769) * r + 1.);
            }

            if(q < 0.0)
                val = -val;
        }
        return mu + sigma * val;
    }

    public static double rand()
    {
        return rand(0.0, 1.0);
    }

    public static double rand(double mu, double sigma)
    {
        if (Double.isNaN(mu) || Double.isInfinite(sigma) || sigma < 0.0) {
            return Dpq.nanWarn();
        }
        if (sigma == 0.0 || Double.isInfinite(mu)) {
            return mu;
        } else {
            return mu + sigma * NormalRand.rand();
        }
    }

    public static void main(String[] args)
    {
        System.out.println(pdf(10, 20, 30));
        System.out.println(pdf(0.025));

        System.out.println(cdf(0.975));
        System.out.println(cdf(100.2, 200.3, 300.4));
        System.out.println(cdf(10.975, 20, 10, false));

        System.out.println(quantile(0.98756));
        System.out.println(quantile(10.98756, 5, 6.7));
        System.out.println(quantile(0.98756, 15, 6.7));
        System.out.println(quantile(0.98756, 15, 6.7));
        System.out.println(quantile(0.98756, 15, 6.7, false));

        System.out.println(rand());
    }
}
