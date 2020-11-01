package distributions;

import distributions.detail.DistBase;
import distributions.detail.Dpq;
import random.NormalRand;
import random.UniformRand;
import special_functions.LogGamma;

public class Gamma extends DistBase
{
    private static final double M_cutoff = Dpq.M_LN2 * Double.MAX_EXPONENT / Dpq.DBL_EPSILON;
    private static final double scaleFactor = sqr(sqr(sqr(4294967296.0)));
    private static double sqr(double x) { return x * x; }

    public static double pdf(double x, double shape, double scale)
    {
        return pdf(x, shape, scale, false);
    }

    public static double pdf(double x, double shape, double scale, boolean give_log)
    {
        if (Double.isNaN(x) || Double.isNaN(shape) || Double.isNaN(scale)) {
            return x + shape + scale;
        }
        if (shape < 0.0 || scale <= 0.0) {
            return Dpq.nanWarn();
        }
        if (x < 0.0) {
            return Dpq.D0(give_log);
        }
        if (shape == 0.0) {
            return (x == 0.0) ? Double.POSITIVE_INFINITY : Dpq.D0(give_log);
        }
        if (x == 0.0) {
            if (shape < 1.0)    return Double.POSITIVE_INFINITY;
            if (shape > 1.0)    return Dpq.D0(give_log);
            return give_log ? -Math.log(scale) : 1.0 / scale;
        }

        double pr;
        if (shape < 1.0) {
            pr = poissonPdfRaw(shape, x / scale, give_log);
            return give_log
                    ? pr + ((Double.isFinite(shape / x)) ? Math.log(shape / x) : Math.log(shape) - Math.log(x))
                    : pr * shape / x;
        }
        pr = poissonPdfRaw(shape - 1.0, x / scale, give_log);
        return give_log ? pr - Math.log(scale) : pr / scale;
    }

    public static double cdf(double x, double alpha, double scale)
    {
        return cdf(x, alpha, scale, true, false);
    }

    public static double cdf(double x, double alpha, double scale, boolean lower_tail, boolean log_p)
    {
        if (Double.isNaN(x) || Double.isNaN(alpha) || Double.isNaN(scale)) {
            return x + alpha + scale;
        }
        if (alpha < 0.0 || scale <= 0.0) {
            return Dpq.nanWarn();
        }

        x /= scale;
        if (Double.isNaN(x)) {
            return x;
        }

        if (alpha == 0.0) {
            return (x <= 0.0) ? Dpq.DT0(lower_tail, log_p) : Dpq.DT1(lower_tail, log_p);
        }
        return cdfRaw(x, alpha, lower_tail, log_p);
    }

    public static double quantile(double p, double alpha, double scale)
    {
        return quantile(p, alpha, scale, true, false);
    }

    public static double quantile(double p, double alpha, double scale, boolean lower_tail, boolean log_p)
    {
        final double EPS1 = 1E-2;
        final double EPS2 = 5E-7;
        final double EPS_N = 1E-15;
        //final double LN_EPS = -36.043653389117156;
        final double MAX_IT = 1000;
        final double pMIN = 1e-100;
        final double pMAX = 1.0 - 1e-14;

        final double i420 = 1.0 / 420.0;
        final double i2520 = 1.0 / 2520.0;
        final double i5040 = 1.0 / 5040.0;

        double p_, a, b, c, g, ch, ch0, p1;
        double p2, q, s1, s2, s3, s4, s5, s6, t, x;
        int max_it_Newton = 1;

        if (Double.isNaN(p) || Double.isNaN(alpha) || Double.isNaN(scale)) {
            return p + scale + alpha;
        }
        if (log_p) {
            if (p > 0.0)
                return Dpq.nanWarn();
            if (p == 0.0)
                return lower_tail ? Double.POSITIVE_INFINITY : 0.0;
            if (p == Double.NEGATIVE_INFINITY)
                return lower_tail ? 0.0 : Double.POSITIVE_INFINITY;
        } else {
            if (p < 0.0 || p > 1)
                return Dpq.nanWarn();
            if (p == 0.0)
                return lower_tail ? 0.0 : Double.POSITIVE_INFINITY;
            if (p == 1.0)
                return lower_tail ? Double.POSITIVE_INFINITY : 0.0;
        }

        if (alpha < 0.0 || scale <= 0.0) {
            return Dpq.nanWarn();
        }
        if (alpha == 0.0) {
            return 0.0;
        }
        if (alpha < 1e-10) {
            max_it_Newton = 7;
        }

        p_ = Dpq.DTqIv(p, lower_tail, log_p);
        g = LogGamma.logGamma(alpha);
        ch = quantileChisApp(p, 2.0 * alpha, g, lower_tail, log_p, EPS1);

        if (Double.isInfinite(ch)) {
            return 0.5 * scale * ch;
        }

        if (ch < EPS2) {
            max_it_Newton = 20;
            x = 0.5 * scale * ch;
            if (!log_p) {
                p = Math.log(p);
                log_p = true;
            }
            if (x == 0.0) {
                final double _1_p = 1.0 + 1e-7;
                final double _1_m = 1.0 - 1e-7;
                x = Double.MIN_VALUE;
                p_ = cdf(x, alpha, scale, lower_tail, log_p);
                if ((lower_tail && p_ > p * _1_p) || (!lower_tail && p_ < p * _1_m)) {
                    return 0.0;
                }
            } else {
                p_ = cdf(x, alpha, scale, lower_tail, log_p);
            }

            if (p_ == Double.NEGATIVE_INFINITY) {
                return 0.0;
            }

            for (int i = 1; i <= max_it_Newton; i++) {
                p1 = p_ - p;

                if(Math.abs(p1) < Math.abs(EPS_N * p)) {
                    break;
                }

                if((g = pdf(x, alpha, scale, log_p)) == Dpq.D0(log_p)) {
                    break;
                }

                t = log_p ? p1 * Math.exp(p_ - g) : p1 / g;
                t = lower_tail ? x - t : x + t;
                p_ = cdf(t, alpha, scale, lower_tail, log_p);
                if (Math.abs(p_ - p) > Math.abs(p1) || (i > 1 && Math.abs(p_ - p) == Math.abs(p1))) {
                    break;
                }
                x = t;
            }
            return x;
        }

        if (p_ > pMAX || p_ < pMIN) {
            max_it_Newton = 20;
            x = 0.5 * scale * ch;
            if (!log_p) {
                p = Math.log(p);
                log_p = true;
            }
            if (x == 0.0) {
                final double _1_p = 1.0 + 1e-7;
                final double _1_m = 1.0 - 1e-7;
                x = Double.MIN_VALUE;
                p_ = cdf(x, alpha, scale, lower_tail, log_p);
                if ((lower_tail && p_ > p * _1_p) || (!lower_tail && p_ < p * _1_m)) {
                    return 0.0;
                }
            } else {
                p_ = cdf(x, alpha, scale, lower_tail, log_p);
            }

            if (p_ == Double.NEGATIVE_INFINITY) {
                return 0.0;
            }

            for (int i = 1; i <= max_it_Newton; i++) {
                p1 = p_ - p;

                if(Math.abs(p1) < Math.abs(EPS_N * p)) {
                    break;
                }

                if((g = pdf(x, alpha, scale, log_p)) == Dpq.D0(log_p)) {
                    break;
                }

                t = log_p ? p1 * Math.exp(p_ - g) : p1 / g;
                t = lower_tail ? x - t : x + t;
                p_ = cdf(t, alpha, scale, lower_tail, log_p);
                if (Math.abs(p_ - p) > Math.abs(p1) || (i > 1 && Math.abs(p_ - p) == Math.abs(p1))) {
                    break;
                }
                x = t;
            }
            return x;
        }

        c = alpha - 1.0;
        s6 = (120.0 + c * (346 + 127 * c)) * i5040;
        ch0 = ch;

        for (int i = 1; i <= MAX_IT; i++) {
            q = ch;
            p1 = 0.5 * ch;
            p2 = p_ - cdfRaw(p1, alpha, true, false);

            if(Double.isInfinite(p2) || ch <= 0) {
                ch = ch0; max_it_Newton = 27; break;
            }

            t = p2 * Math.exp(alpha * Dpq.M_LN2 + g + p1 - c * Math.log(ch));
            b = t / ch;
            a = 0.5 * t - b * c;
            s1 = (210 + a * (140 + a * (105 + a * (84 + a * (70 + 60 * a))))) * i420;
            s2 = (420 + a * (735 + a * (966 + a * (1141 + 1278 * a)))) * i2520;
            s3 = (210 + a * (462 + a * (707 + 932 * a))) * i2520;
            s4 = (252 + a * (672 + 1182 * a) + c * (294 + a * (889 + 1740 * a))) * i5040;
            s5 = (84 + 2264 * a + c * (1175 + 606 * a)) * i2520;

            ch += t * (1 + 0.5 * t * s1 - b * c * (s1 - b * (s2 - b * (s3 - b * (s4 - b * (s5 - b * s6))))));

            if (Math.abs(q - ch) < EPS2 * ch) {
                break;
            }
            if (Math.abs(q - ch) > 0.1 * ch) {
                if(ch < q) ch = 0.9 * q; else ch = 1.1 * q;
            }
        }

        // END:
        x = 0.5 * scale * ch;
        if (!log_p) {
            p = Math.log(p);
            log_p = true;
        }
        if (x == 0.0) {
            final double _1_p = 1.0 + 1e-7;
            final double _1_m = 1.0 - 1e-7;
            x = Double.MIN_VALUE;
            p_ = cdf(x, alpha, scale, lower_tail, log_p);
            if ((lower_tail && p_ > p * _1_p) || (!lower_tail && p_ < p * _1_m)) {
                return 0.0;
            }
        } else {
            p_ = cdf(x, alpha, scale, lower_tail, log_p);
        }

        if (p_ == Double.NEGATIVE_INFINITY) {
            return 0.0;
        }

        for (int i = 1; i <= max_it_Newton; i++) {
            p1 = p_ - p;

            if(Math.abs(p1) < Math.abs(EPS_N * p)) {
                break;
            }

            if((g = pdf(x, alpha, scale, log_p)) == Dpq.D0(log_p)) {
                break;
            }

            t = log_p ? p1 * Math.exp(p_ - g) : p1 / g;
            t = lower_tail ? x - t : x + t;
            p_ = cdf(t, alpha, scale, lower_tail, log_p);
            if (Math.abs(p_ - p) > Math.abs(p1) || (i > 1 && Math.abs(p_ - p) == Math.abs(p1))) {
                break;
            }
            x = t;
        }
        return x;
    }

    public static double rand(double a, double scale)
    {
        final double sqrt32 = 5.656854;
        final double exp_m1 = 0.36787944117144232159;/* exp(-1) = 1/e */

        /*
         * Coefficients q[k] - for q0 = sum(q[k]*a^(-k))
         * Coefficients a[k] - for q = q0+(t*t/2)*sum(a[k]*v^k)
         * Coefficients e[k] - for exp(q)-1 = sum(e[k]*q^k)
         */
        final double q1 = 0.04166669;
        final double q2 = 0.02083148;
        final double q3 = 0.00801191;
        final double q4 = 0.00144121;
        final double q5 = -7.388e-5;
        final double q6 = 2.4511e-4;
        final double q7 = 2.424e-4;

        final double a1 = 0.3333333;
        final double a2 = -0.250003;
        final double a3 = 0.2000062;
        final double a4 = -0.1662921;
        final double a5 = 0.1423657;
        final double a6 = -0.1367177;
        final double a7 = 0.1233795;

        double aa = 0.0;
        double aaa = 0.0;
        double s = 0.0, s2 = 0.0, d = 0.0;    /* no. 1 (step 1) */
        double q0 = 0.0, b = 0.0, si = 0.0, c = 0.0;/* no. 2 (step 4) */

        double e, p, q, r, t, u, v, w, x, ret_val;

        if (Double.isNaN(a) || Double.isNaN(scale)) {
            return Dpq.nanWarn();
        }
        if (a <= 0.0 || scale <= 0.0) {
            if(scale == 0.0 || a == 0.0) {
                return 0.0;
            }
            return Dpq.nanWarn();
        }
        if(Double.isInfinite(a) || Double.isInfinite(scale)) {
            return Double.POSITIVE_INFINITY;
        }

        if (a < 1.0) {
            e = 1.0 + exp_m1 * a;
            for(;;) {
                p = e * UniformRand.rand();
                if (p >= 1.0) {
                    x = -Math.log((e - p) / a);
                    if (expRand() >= (1.0 - a) * Math.log(x))
                        break;
                } else {
                    x = Math.exp(Math.log(p) / a);
                    if (expRand() >= x)
                        break;
                }
            }
            return scale * x;
        }

        /* --- a >= 1 : GD algorithm --- */

        /* Step 1: Recalculations of s2, s, d if a has changed */
        if (a != aa) {
            aa = a;
            s2 = a - 0.5;
            s = Math.sqrt(s2);
            d = sqrt32 - s * 12.0;
        }
    /* Step 2: t = standard normal deviate,
               x = (s,1/2) -normal deviate. */

        /* immediate acceptance (i) */
        t = NormalRand.rand();
        x = s + 0.5 * t;
        ret_val = x * x;
        if (t >= 0.0)
            return scale * ret_val;

        /* Step 3: u = 0,1 - uniform sample. squeeze acceptance (s) */
        u = UniformRand.rand();
        if (d * u <= t * t * t)
            return scale * ret_val;

        /* Step 4: recalculations of q0, b, si, c if necessary */

        if (a != aaa) {
            aaa = a;
            r = 1.0 / a;
            q0 = ((((((q7 * r + q6) * r + q5) * r + q4) * r + q3) * r
                    + q2) * r + q1) * r;

            /* Approximation depending on size of parameter a */
            /* The constants in the expressions for b, si and c */
            /* were established by numerical experiments */

            if (a <= 3.686) {
                b = 0.463 + s + 0.178 * s2;
                si = 1.235;
                c = 0.195 / s - 0.079 + 0.16 * s;
            } else if (a <= 13.022) {
                b = 1.654 + 0.0076 * s2;
                si = 1.68 / s + 0.275;
                c = 0.062 / s + 0.024;
            } else {
                b = 1.77;
                si = 0.75;
                c = 0.1515 / s;
            }
        }
        /* Step 5: no quotient test if x not positive */

        if (x > 0.0) {
            /* Step 6: calculation of v and quotient q */
            v = t / (s + s);
            if (Math.abs(v) <= 0.25)
                q = q0 + 0.5 * t * t * ((((((a7 * v + a6) * v + a5) * v + a4) * v
                        + a3) * v + a2) * v + a1) * v;
            else
                q = q0 - s * t + 0.25 * t * t + (s2 + s2) * Math.log(1.0 + v);


            /* Step 7: quotient acceptance (q) */
            if (Math.log(1.0 - u) <= q)
                return scale * ret_val;
        }

        for(;;) {
            /* Step 8: e = standard exponential deviate
             *	u =  0,1 -uniform deviate
             *	t = (b,si)-double exponential (laplace) sample */
            e = expRand();
            u = UniformRand.rand();
            u = u + u - 1.0;
            if (u < 0.0)
                t = b - si * e;
            else
                t = b + si * e;
            /* Step	 9:  rejection if t < tau(1) = -0.71874483771719 */
            if (t >= -0.71874483771719) {
                /* Step 10:	 calculation of v and quotient q */
                v = t / (s + s);
                if (Math.abs(v) <= 0.25)
                    q = q0 + 0.5 * t * t *
                            ((((((a7 * v + a6) * v + a5) * v + a4) * v + a3) * v
                                    + a2) * v + a1) * v;
                else
                    q = q0 - s * t + 0.25 * t * t + (s2 + s2) * Math.log(1.0 + v);
                /* Step 11:	 hat acceptance (h) */
                /* (if q not positive go to step 8) */
                if (q > 0.0) {
                    w = Math.expm1(q);
                    /*  ^^^^^ original code had approximation with rel.err < 2e-7 */
                    /* if t is rejected sample again at step 8 */
                    if (c * Math.abs(u) <= w * Math.exp(e - 0.5 * t * t))
                        break;
                }
            }
        } /* repeat .. until  `t' is accepted */
        x = s + 0.5 * t;
        return scale * x * x;
    }

    private static double logCF(double x, double i, double d, double eps)
    {
        double c1 = 2 * d;
        double c2 = i + d;
        double c4 = c2 + d;
        double a1 = c2;
        double b1 = i * (c2 - i * x);
        double b2 = d * d * x;
        double a2 = c4 * c2 - b2;


        b2 = c4 * b1 - i * b2;

        while (Math.abs(a2 * b1 - a1 * b2) > Math.abs(eps * b1 * b2)) {
            double c3 = c2 * c2 * x;
            c2 += d;
            c4 += d;
            a1 = c4 * a2 - c3 * a1;
            b1 = c4 * b2 - c3 * b1;

            c3 = c1 * c1 * x;
            c1 += d;
            c4 += d;
            a2 = c4 * a1 - c3 * a2;
            b2 = c4 * b1 - c3 * b2;

            if (Math.abs (b2) > scaleFactor) {
                a1 /= scaleFactor;
                b1 /= scaleFactor;
                a2 /= scaleFactor;
                b2 /= scaleFactor;
            } else if (Math.abs (b2) < 1 / scaleFactor) {
                a1 *= scaleFactor;
                b1 *= scaleFactor;
                a2 *= scaleFactor;
                b2 *= scaleFactor;
            }
        }
        return a2 / b2;
    }

    private static double log1pmx(double x)
    {
        final double minLog1Value = -0.79149064;

        if (x > 1 || x < minLog1Value) {
            return Math.log1p(x) - x;
        } else {
            double r = x / (2 + x), yes = r * r;
            if (Math.abs(x) < 1e-2) {
                final double two = 2.0;
                return r * ((((two / 9 * yes + two / 7) * yes + two / 5) * yes + two / 3) * yes - x);
            } else {
                final double tol_logCF = 1e-14;
                return r * (2 * yes * logCF(yes, 3, 2, tol_logCF) - x);
            }
        }
    }

    private static double logGamma1p(double a)
    {
        final double euler_const = 0.5772156649015328606065120900824024;

        final int N = 40;
        final double[] coeFfs = {
                0.3224670334241132182362075833230126e-0,
                0.6735230105319809513324605383715000e-1,
                0.2058080842778454787900092413529198e-1,
                0.7385551028673985266273097291406834e-2,
                0.2890510330741523285752988298486755e-2,
                0.1192753911703260977113935692828109e-2,
                0.5096695247430424223356548135815582e-3,
                0.2231547584535793797614188036013401e-3,
                0.9945751278180853371459589003190170e-4,
                0.4492623673813314170020750240635786e-4,
                0.2050721277567069155316650397830591e-4,
                0.9439488275268395903987425104415055e-5,
                0.4374866789907487804181793223952411e-5,
                0.2039215753801366236781900709670839e-5,
                0.9551412130407419832857179772951265e-6,
                0.4492469198764566043294290331193655e-6,
                0.2120718480555466586923135901077628e-6,
                0.1004322482396809960872083050053344e-6,
                0.4769810169363980565760193417246730e-7,
                0.2271109460894316491031998116062124e-7,
                0.1083865921489695409107491757968159e-7,
                0.5183475041970046655121248647057669e-8,
                0.2483674543802478317185008663991718e-8,
                0.1192140140586091207442548202774640e-8,
                0.5731367241678862013330194857961011e-9,
                0.2759522885124233145178149692816341e-9,
                0.1330476437424448948149715720858008e-9,
                0.6422964563838100022082448087644648e-10,
                0.3104424774732227276239215783404066e-10,
                0.1502138408075414217093301048780668e-10,
                0.7275974480239079662504549924814047e-11,
                0.3527742476575915083615072228655483e-11,
                0.1711991790559617908601084114443031e-11,
                0.8315385841420284819798357793954418e-12,
                0.4042200525289440065536008957032895e-12,
                0.1966475631096616490411045679010286e-12,
                0.9573630387838555763782200936508615e-13,
                0.4664076026428374224576492565974577e-13,
                0.2273736960065972320633279596737272e-13,
                0.1109139947083452201658320007192334e-13
        };

        final double c = 0.2273736845824652515226821577978691e-12;

        if (Math.abs (a) >= 0.5)
            return LogGamma.logGamma(a + 1);

        double logGam = c * logCF(-a / 2, N + 2, 1, 1e-14);
        for (int i = N - 1; i >= 0; i--) {
            logGam = coeFfs[i] - a * logGam;
        }

        return (a * logGam - euler_const) * a - log1pmx(a);
    }

    private static double logSpaceSub(double log_x, double log_y)
    {
        return log_x + Dpq.log1Exp(log_y - log_x);
    }

    private static double logSpaceSum(double[] log_x, int n)
    {
        if (n == 0) return Double.NEGATIVE_INFINITY;
        if (n == 1) return log_x[0];
        if (n == 2) return logSpaceAdd(log_x[0], log_x[1]);

        double Mx = log_x[0];
        for (int i = 1; i < n; i++) {
            if (Mx < log_x[i])
                Mx = log_x[i];
        }

        double s = 0.0;
        for (int i = 0; i < n; i++) {
            s += Math.exp(log_x[i] - Mx);
        }
        return Mx + Math.log(s);
    }

    private static double poissonPdfWrap(double x_plus_1, double lambda, boolean give_log)
    {
        if (Double.isInfinite(lambda)) {
            return Dpq.D0(give_log);
        }
        if (x_plus_1 > 1.0) {
            return poissonPdfRaw(x_plus_1 - 1, lambda, give_log);
        }
        if (lambda > Math.abs(x_plus_1 - 1.0) * M_cutoff) {
            return Dpq.DExp(-lambda - LogGamma.logGamma(x_plus_1), give_log);
        } else {
            double d = poissonPdfRaw(x_plus_1, lambda, give_log);
            return give_log
                    ? d + Math.log(x_plus_1 / lambda)
                    : d * (x_plus_1 / lambda);
        }
    }

    private static double gammaCdfSmallX(double x, double alpha, boolean lower_tail, boolean log_p)
    {
        double sum = 0.0, c = alpha, n = 0.0, term;

        do {
            n++;
            c *= -x / n;
            term = c / (alpha + n);
            sum += term;
        } while (Math.abs(term) > Dpq.DBL_EPSILON * Math.abs(sum));

        if (lower_tail) {
            double f1 = log_p ? Math.log1p(sum) : 1 + sum;
            double f2;
            if (alpha > 1) {
                f2 = poissonPdfRaw(alpha, x, log_p);
                f2 = log_p ? f2 + x : f2 * Math.exp(x);
            } else if (log_p) {
                f2 = alpha * Math.log(x) - logGamma1p(alpha);
            } else {
                f2 = Math.pow(x, alpha) / Math.exp(logGamma1p(alpha));
            }
            return log_p ? f1 + f2 : f1 * f2;
        } else {
            double lf2 = alpha * Math.log(x) - logGamma1p(alpha);
            if (log_p) {
                return Dpq.log1Exp(Math.log1p(sum) + lf2);
            } else {
                double f2m1 = Math.expm1(lf2);
                return -(sum + f2m1 + sum * f2m1);
            }
        }
    }

    private static double pdUpperSeries(double x, double y, boolean log_p)
    {
        double term = x / y;
        double sum = term;

        do {
            y++;
            term *= x / y;
            sum += term;
        } while (term > sum * Dpq.DBL_EPSILON);
        return log_p ? Math.log(sum) : sum;
    }

    private static double pdLowerCF(double y, double d)
    {
        double f= 0.0, of, f0;
        final int max_it = 200000;
        if (y == 0) {
            return 0.0;
        }

        f0 = y / d;
        if (Math.abs(y - 1) < Math.abs(d) * Dpq.DBL_EPSILON) {
            return f0;
        }

        if (f0 > 1.0) {
            f0 = 1.0;
        }

        double c2 = y;
        double c4 = d;
        double a1 = 0;
        double b1 = 1;
        double a2 = y;
        double b2 = d;

        while (b2 > scaleFactor) {
            a1 /= scaleFactor;
            b1 /= scaleFactor;
            a2 /= scaleFactor;
            b2 /= scaleFactor;
        }

        double i = 0.0; of = -1.0;
        while (i < max_it) {

            i++;	c2--;
            double c3 = i * c2;	c4 += 2;
            a1 = c4 * a2 + c3 * a1;
            b1 = c4 * b2 + c3 * b1;

            i++;	c2--;	c3 = i * c2;	c4 += 2;
            a2 = c4 * a1 + c3 * a2;
            b2 = c4 * b1 + c3 * b2;

            if (b2 > scaleFactor) {
                a1 /= scaleFactor;
                b1 /= scaleFactor;
                a2 /= scaleFactor;
                b2 /= scaleFactor;
            }

            if (b2 != 0) {
                f = a2 / b2;
                /* convergence check: relative; "absolute" for very small f : */
                if (Math.abs (f - of) <= Dpq.DBL_EPSILON * Math.max(f0, Math.abs(f))) {
                    return f;
                }
                of = f;
            }
        }
        // should not happen...
        System.out.printf(" ** NON-convergence in Gamma.cdf()'s pdLowerCF() f = %g.\n", f);
        return f;
    }

    private static double pdLowerSeries(double lambda, double y)
    {
        double term = 1.0, sum = 0.0;

        while (y >= 1 && term > sum * Dpq.DBL_EPSILON) {
            term *= y / lambda;
            sum += term;
            y--;
        }

        if (y != Math.floor(y)) {
            double f = pdLowerCF(y, lambda + 1.0 - y);
            sum += term * f;
        }
        return sum;
    }

    private static double dp_norm(double x, boolean lower_tail, double lp)
    {
        if (x < 0) {
            x = -x;
            lower_tail = !lower_tail;
        }

        if (x > 10 && !lower_tail) {
            double term = 1 / x;
            double sum = term;
            double x2 = x * x;
            double i = 1;

            do {
                term *= -i / x2;
                sum += term;
                i += 2;
            } while (Math.abs (term) > Dpq.DBL_EPSILON * sum);

            return 1 / sum;
        } else {
            double d = pdf(x, 0.0, 1.0, false);
            return d / Math.exp(lp);
        }
    }

    private static double poissonAsymptotic(double x, double lambda, boolean lower_tail, boolean log_p)
    {
        final double[] coeFs_a = {
                -1e99, /* placeholder used for 1-indexing */
                2/3.,
                -4/135.,
                8/2835.,
                16/8505.,
                -8992/12629925.,
                -334144/492567075.,
                698752/1477701225.
        };

        final double[] coeFs_b = {
                -1e99, /* placeholder */
                1/12.0,
                1/288.0,
                -139/51840.0,
                -571/2488320.0,
                163879/209018880.0,
                5246819/75246796800.0,
                -534703531/902961561600.0
        };

        double dfm = lambda - x;
        // If lambda is large, the distribution is highly concentrated
        // about lambda.  So representation error in x or lambda can lead
        // to arbitrarily large values of pt_ and hence divergence of the
        // coefficients of this approximation.

        double pt_ = -log1pmx (dfm / x);
        double s2pt = Math.sqrt(2 * x * pt_);
        if (dfm < 0) {
            s2pt = -s2pt;
        }

        double res1_term;
        double res2_term;
        double res12 = 0;
        double res1_ig = res1_term = Math.sqrt(x);
        double res2_ig = res2_term = s2pt;
        for (int i = 1; i < 8; i++) {
            res12 += res1_ig * coeFs_a[i];
            res12 += res2_ig * coeFs_b[i];
            res1_term *= pt_ / i ;
            res2_term *= 2 * pt_ / (2 * i + 1);
            res1_ig = res1_ig / x + res1_term;
            res2_ig = res2_ig / x + res2_term;
        }

        double elFb = x;
        double elFb_term = 1;
        for (int i = 1; i < 8; i++) {
            elFb += elFb_term * coeFs_b[i];
            elFb_term /= x;
        }
        if (!lower_tail) {
            elFb = -elFb;
        }

        double f = res12 / elFb;
        double np = cdf(s2pt, 0.0, 1.0, !lower_tail, log_p);

        if (log_p) {
            double n_d_over_p = dp_norm(s2pt, !lower_tail, np);

            return np + Math.log1p(f * n_d_over_p);
        } else {
            double nd = pdf(s2pt, 0.0, 1.0, log_p);
            return np + f * nd;
        }
    }

    private static double cdfRaw(double x, double alpha, boolean lower_tail, boolean log_p)
    {
        double res;
        if (x <= 0.0) {
            return Dpq.DT0(lower_tail, log_p);
        }
        if (x >= Double.POSITIVE_INFINITY) {
            return Dpq.DT0(lower_tail, log_p);
        }

        if (x < 1) {
            res = gammaCdfSmallX(x, alpha, lower_tail, log_p);
        } else if (x <= alpha - 1 && x < 0.8 * (alpha + 50)) {
            double sum = pdUpperSeries(x, alpha, log_p);
            double d = poissonPdfWrap(alpha, x, log_p);

            if (!lower_tail) {
                res = log_p ? Dpq.log1Exp(d + sum) : 1.0 - d * sum;
            } else {
                res = log_p ? sum + d : sum * d;
            }
        } else if (alpha - 1 < x && alpha < 0.8 * (x + 50)) {
            double sum;
            double d = poissonPdfWrap(alpha, x, log_p);

            if (alpha < 1) {
                if (x * Dpq.DBL_EPSILON > 1 - alpha) {
                    sum = Dpq.D1(log_p);
                } else {
                    double f = pdLowerCF(alpha, x - (alpha - 1)) * x / alpha;
                    sum = log_p ? Math.log(f) : f;
                }
            } else {
                sum = pdLowerSeries(x, alpha - 1);
                sum = log_p ? Math.log1p(sum) : 1 + sum;
            }

            if (!lower_tail) {
                res = log_p ? sum + d : sum * d;
            } else {
                res = log_p ? Dpq.log1Exp(d + sum) : 1 - d * sum;
            }
        } else {
            res = poissonAsymptotic(alpha - 1, x, !lower_tail, log_p);
        }

        /*
         * We lose a fair amount of accuracy to underflow in the cases
         * where the final result is very close to DBL_MIN.	 In those
         * cases, simply redo via log space.
         */
        if (!log_p && res < Double.MIN_VALUE / Dpq.DBL_EPSILON) {
            return Math.exp (cdfRaw(x, alpha, lower_tail, true));
        } else {
            return res;
        }
    }

    private static double quantileChisApp(double p, double nu, double g,
                                          boolean lower_tail, boolean log_p, double tol)
    {
        final double C7	= 4.67;
        final double C8	= 6.66;
        final double C9	= 6.73;
        final double C10 = 13.32;

        if (Double.isNaN(p) || Double.isNaN(nu)) {
            return p + nu;
        }
        if ((log_p && p > 0) || (!log_p && (p < 0 || p > 1))) {
            return Dpq.nanWarn();
        }
        if (nu <= 0) {
            return Dpq.nanWarn();
        }

        double alpha, a, c, ch, p1;
        double p2, q, t, x;

        alpha = 0.5 * nu;
        c = alpha - 1;

        p1 = Dpq.DT_log(p, lower_tail, log_p);
        if (nu < -1.24 * p1) {
            double l_gam1pa = (alpha < 0.5) ? logGamma1p(alpha) : (Math.log(alpha) + g);
            ch = Math.exp((l_gam1pa + p1) / alpha + Dpq.M_LN2);

        } else if (nu > 0.32) {

            x = Normal.quantile(p, 0.0, 1.0, lower_tail, log_p);
            p1 = 2.0 / (9 * nu);
            ch = nu * Math.pow(x * Math.sqrt(p1) + 1 - p1, 3);

            if (ch > 2.2 * nu + 6) {
                ch = -2 * (Dpq.DTCLog(p, lower_tail, log_p) - c * Math.log(0.5 * ch) + g);
            }
        } else {
            ch = 0.4;
            a = Dpq.DTCLog(p, lower_tail, log_p) + g + c * Dpq.M_LN2;

            do {
                q = ch;
                p1 = 1.0 / (1 + ch * (C7 + ch));
                p2 = ch * (C9 + ch * (C8 + ch));
                t = -0.5 + (C7 + 2 * ch) * p1 - (C9 + ch * (C10 + 3 * ch)) / p2;
                ch -= (1- Math.exp(a + 0.5 * ch) * p2 * p1) / t;
            } while (Math.abs(q - ch) > tol * Math.abs(ch));
        }

        return ch;
    }

    public static void main(String[] args)
    {
        System.out.println(pdf(4, 1, 1));

        System.out.println(cdf(0.987, 1, 2));
        System.out.println(quantile(0.975, 10, 20, true, true));
        System.out.println(rand(1.0, 2.0));
    }
}