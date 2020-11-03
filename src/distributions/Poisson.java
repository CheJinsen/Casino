package distributions;

import distributions.detail.DistBase;
import distributions.detail.Dpq;
import distributions.detail.RefDouble;
import random.NormalRand;
import random.UniformRand;

public class Poisson extends DistBase
{
    public static double pdf(double x, double lambda)
    {
        return pdf(x, lambda, false);
    }

    public static double pdf(double x, double lambda, boolean give_log)
    {
        if (Double.isNaN(x) || Double.isNaN(lambda)) {
            return x + lambda;
        }
        if (lambda < 0.0) {
            return Dpq.nanWarn();
        }
        if (Dpq.nonInt(x)) {
            System.out.printf("non-integer x = %f", x);
            return Dpq.D0(give_log);
        }
        if (x < 0.0 || Double.isInfinite(x)) {
            return Dpq.D0(give_log);
        }

        x = Math.round(x);
        return poissonPdfRaw(x, lambda, give_log);
    }

    public static double cdf(double x, double lambda)
    {
        return cdf(x, lambda, true, false);
    }

    public static double cdf(double x, double lambda, boolean lower_tail)
    {
        return cdf(x, lambda, lower_tail, false);
    }

    public static double cdf(double x, double lambda, boolean lower_tail, boolean log_p)
    {
        if (Double.isNaN(x) || Double.isNaN(lambda)) {
            return x + lambda;
        }
        if (lambda < 0.0) {
            return Dpq.nanWarn();
        }
        if (x < 0.0) {
            return Dpq.DT0(lower_tail, log_p);
        }
        if (lambda == 0.0) {
            return Dpq.DT1(lower_tail, log_p);
        }
        if (Double.isInfinite(x)) {
            return Dpq.DT1(lower_tail, log_p);
        }

        x = Math.floor(x + 1e-7);
        return Gamma.cdf(lambda, x + 1.0, 1.0, !lower_tail, log_p);
    }

    public static double quantile(double p, double lambda)
    {
        return quantile(p, lambda, true, false);
    }

    public static double quantile(double p, double lambda, boolean lower_tail)
    {
        return quantile(p, lambda, lower_tail, false);
    }

    public static double quantile(double p, double lambda, boolean lower_tail, boolean log_p)
    {
        if (Double.isNaN(p) || Double.isNaN(lambda)) {
            return p + lambda;
        }
        if(Double.isInfinite(lambda)) {
            return Dpq.nanWarn();
        }
        if(lambda < 0) {
            return Dpq.nanWarn();
        }
        if ((log_p && p > 0) || (!log_p && (p < 0 || p > 1))) {
            return Dpq.nanWarn();
        }
        if (lambda == 0) {
            return 0;
        }
        if (p == Dpq.DT0(lower_tail, log_p)) {
            return 0;
        }
        if (p == Dpq.DT1(lower_tail, log_p)) {
            return Double.POSITIVE_INFINITY;
        }

        double sigma = Math.sqrt(lambda);
        double gamma = 1.0 / sigma;

        if (!lower_tail || log_p) {
            p = Dpq.DTqIv(p, lower_tail, log_p);
            if (p == 0.0) return 0;
            if (p == 1.0) return Double.POSITIVE_INFINITY;
        }
        if (p + 1.01 * Dpq.DBL_EPSILON >= 1.0) {
            return Double.POSITIVE_INFINITY;
        }

        RefDouble z = new RefDouble(0.0);
        z.val = Normal.quantile(p, 0.0, 1.0, true, false);
        double yes = Math.round(lambda + sigma * (z.val + gamma * (z.val * z.val - 1.0) / 6.0));

        z.val = cdf(yes, lambda, true, false);

        p *= 1 - 64 * Dpq.DBL_EPSILON;

        if(lambda < 1e5) {
            return doSearch(yes, z, p, lambda, 1.0);
        }

        {
            double incR = Math.floor(yes * 0.001), old_incR;
            do {
                old_incR = incR;
                yes = doSearch(yes, z, p, lambda, incR);
                incR = Math.max(1.0, Math.floor(incR / 100.0));
            } while (old_incR > 1.0 && incR > lambda * 1e-15);
            return yes;
        }
    }

    private static double doSearch(double yes, RefDouble z, double p, double lambda, double incR)
    {
        if (z.val >= p) {
            for (;;) {
                if(yes == 0 || (z.val = cdf(yes - incR, lambda, true, false)) < p)
                    return yes;
                yes = Math.max(0, yes - incR);
            }
        } else {
            for (;;) {
                yes = yes + incR;
                if ((z.val = cdf(yes, lambda, true, false)) >= p)
                    return yes;
            }
        }
    }

    public static double rand(double mu)
    {
        final double a0	= -0.5;
        final double a1	= 0.3333333;
        final double a2	= -0.2500068;
        final double a3	= 0.2000118;
        final double a4	= -0.1661269;
        final double a5	= 0.1421878;
        final double a6	= -0.1384794;
        final double a7	= 0.1250060;

        final double one_7 = 0.1428571428571428571;
        final double one_12 = 0.0833333333333333333;
        final double one_24 = 0.0416666666666666667;

        // Factorial Table (0:9)!
        final double[] fact = {1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0, 5040.0, 40320.0, 362880.0};

        int l = 0, m = 0;

        double b1, b2, c = 0.0, c0 = 0.0, c1 = 0.0, c2 = 0.0, c3 = 0.0;
        double[] pp = new double[36];
        double p0 = 0.0, p = 0.0, q = 0.0, s = 0.0, d = 0.0, omega = 0.0;
        double big_l = 0.0;
        double mu_prev = 0.0, mu_prev2 = 0.0;

        double del, dif_muk = 0.0, E, fk= 0.0, fx, fy, g, px, py, t, u= 0.0, v, x;
        double poi = -1.0;
        int k;
        boolean big_mu, new_big_mu = false;

        if (Double.isInfinite(mu) || mu < 0.0) {
            return Dpq.nanWarn();
        }

        if (mu <= 0.0) {
            return 0.0;
        }

        big_mu = mu >= 10.0;

        if (!(big_mu && mu == mu_prev)) {
            if (big_mu) {
                new_big_mu = true;
                s = Math.sqrt(mu);
                d = 6.0 * mu * mu;
                big_l = Math.floor(mu - 1.1484);
            } else {
                if (mu != mu_prev) {
                    m = Math.max(1, (int)mu);
                    l = 0;
                    q = p0 = p = Math.exp(-mu);
                }

                for (;;) {
                    u = UniformRand.rand();
                    if (u <= p0) {
                        return 0.0;
                    }
                    if (l != 0) {
                        for (k = (u <= 0.458) ? 1 : Math.min(l, m);  k <= l; k++) {
                            if (u <= pp[k])
                                return k;
                        }
                        continue;
                    }

                    l++;
                    for (k = l; k <= 35; k++) {
                        p *= mu / k;
                        q += p;
                        pp[k] = q;
                        if (u <= q) {
                            return k;
                        }
                    }
                    l = 35;
                }
            }
        }

        g = mu + s * NormalRand.rand();

        if (g >= 0.0) {
            poi = Math.floor(g);
            if (poi >= big_l)
                return poi;
            fk = poi;
            dif_muk = mu - fk;
            u = UniformRand.rand();
            if (d * u >= dif_muk * dif_muk * dif_muk)
                return poi;
        }

        if (new_big_mu || mu != mu_prev2) {
            omega = Dpq.M_1_SQRT_2PI / s;

            b1 = one_24 / mu;
            b2 = 0.3 * b1 * b1;
            c3 = one_7 * b1 * b2;
            c2 = b2 - 15. * c3;
            c1 = b1 - 6. * b2 + 45. * c3;
            c0 = 1. - b1 + 3. * b2 - 15. * c3;
            c = 0.1069 / mu;
        }

        if (g >= 0.0) {

            //goto Step_F;
            if (poi < 10) {
                px = -mu;
                py = Math.pow(mu, poi) / fact[(int)poi];
            } else {
                del = one_12 / fk;
                del = del * (1. - 4.8 * del * del);
                v = dif_muk / fk;
                if (Math.abs(v) <= 0.25) {
                    px = fk * v * v * (((((((a7 * v + a6) * v + a5) * v + a4) * v + a3) *
                            v + a2) * v + a1) * v + a0) - del;
                } else {
                    px = fk * Math.log(1.0 + v) - dif_muk - del;
                }
                py = Dpq.M_1_SQRT_2PI / Math.sqrt(fk);
            }
            x = (0.5 - dif_muk) / s;
            x *= x;
            fx = -0.5 * x;
            fy = omega * (((c3 * x + c2) * x + c1) * x + c0);
            if (fy - u * fy <= py * Math.exp(px - fx))
                return poi;

            for (;;) {
                E = expRand();
                u = 2 * UniformRand.rand() - 1.0;
                t = 1.8 + fSign(E, u);
                if (t > -0.6744) {
                    poi = Math.floor(mu + s * t);
                    fk = poi;
                    dif_muk = mu - fk;

                    //Step_F: // 'subroutine' F : calculation of px,py,fx,fy.

                    if (poi < 10) {
                        px = -mu;
                        py = Math.pow(mu, poi) / fact[(int)poi];
                    } else {
                        del = one_12 / fk;
                        del = del * (1. - 4.8 * del * del);
                        v = dif_muk / fk;
                        if (Math.abs(v) <= 0.25) {
                            px = fk * v * v * (((((((a7 * v + a6) * v + a5) * v + a4) * v + a3) *
                                    v + a2) * v + a1) * v + a0) - del;
                        } else {
                            px = fk * Math.log(1.0 + v) - dif_muk - del;
                        }
                        py = Dpq.M_1_SQRT_2PI / Math.sqrt(fk);
                    }
                    x = (0.5 - dif_muk) / s;
                    x *= x;
                    fx = -0.5 * x;
                    fy = omega * (((c3 * x + c2) * x + c1) * x + c0);
                    if (c * Math.abs(u) <= py * Math.exp(px + E) - fy * Math.exp(fx + E))
                        break;
                }
            }
            return poi;
        }

        for (;;) {
            E = expRand();
            u = 2 * UniformRand.rand() - 1.0;
            t = 1.8 + fSign(E, u);
            if (t > -0.6744) {
                poi = Math.floor(mu + s * t);
                fk = poi;
                dif_muk = mu - fk;

                //Step_F: // 'subroutine' F : calculation of px,py,fx,fy.

                if (poi < 10) {
                    px = -mu;
                    py = Math.pow(mu, poi) / fact[(int)poi];
                } else {
                    del = one_12 / fk;
                    del = del * (1. - 4.8 * del * del);
                    v = dif_muk / fk;
                    if (Math.abs(v) <= 0.25) {
                        px = fk * v * v * (((((((a7 * v + a6) * v + a5) * v + a4) * v + a3) *
                                v + a2) * v + a1) * v + a0) - del;
                    } else {
                        px = fk * Math.log(1.0 + v) - dif_muk - del;
                    }
                    py = Dpq.M_1_SQRT_2PI / Math.sqrt(fk);
                }
                x = (0.5 - dif_muk) / s;
                x *= x;
                fx = -0.5 * x;
                fy = omega * (((c3 * x + c2) * x + c1) * x + c0);
                if (c * Math.abs(u) <= py * Math.exp(px + E) - fy * Math.exp(fx + E))
                    break;
            }
        }
        return poi;
    }

    private static double fSign(double x, double y)
    {
        if (Double.isNaN(x) || Double.isNaN(y)) {
            return x + y;
        }
        return y >= 0.0 ? Math.abs(x) : -Math.abs(x);
    }

    public static void main(String[] args)
    {
        System.out.println(pdf(1, 2.3));
        System.out.println(cdf(5, 2.3));
        System.out.println(cdf(5, 2.3, false));
        System.out.println(quantile(0.987, 8.675));

        for (int i = 0; i < 10; i++) {
            System.out.print(rand(11) + " ");
        }
        System.out.println();
    }
}
