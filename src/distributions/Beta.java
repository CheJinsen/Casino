package distributions;

import distributions.detail.*;

public class Beta extends DistBase
{
    private static final double DBL_VERY_MIN = Double.MIN_VALUE / 4.0;
    private static final double DBL_LOG_V_MIN = Dpq.M_LN2 * (Double.MIN_EXPONENT - 2);
    private static final double DBL_1_EPS = 0x1.fffffffffffffp-1;

    private static final double USE_LOG_X_CUTOFF = -5.0;
    private static final int N_NEWTON_FREE = 4;

    private static final double fpu = 3e-308;
    private static final double acu_min = 1e-300;
    private static final double p_lo = fpu;
    private static final double p_hi = 1 - 2.22e-16;

    private static final double const1 = 2.30753;
    private static final double const2 = 0.27061;
    private static final double const3 = 0.99229;
    private static final double const4 = 0.04481;

    public static double pdf(double x, double a, double b)
    {
        return pdf(x, a, b, false);
    }

    public static double pdf(double x, double a, double b, boolean give_log)
    {
        if (Double.isNaN(x) || Double.isNaN(a) || Double.isNaN(b)) {
            return x + a + b;
        }
        if (a < 0.0 || b < 0.0) {
            return Dpq.nanWarn();
        }
        if (x < 0.0 || x > 1.0) {
            return Dpq.D0(give_log);
        }
        if (a == 0.0 || b == 0.0 || Double.isInfinite(a) || Double.isInfinite(b)) {
            if (a == 0.0 && b == 0.0) {
                if (x == 0.0 || x == 1.0)
                    return Double.POSITIVE_INFINITY;
                else
                    return Dpq.D0(give_log);
            }
            if (a == 0.0 || a / b == 0.0) {
                if (x == 0.0)
                    return Double.POSITIVE_INFINITY;
                else
                    return Dpq.D0(give_log);
            }
            if (b == 0.0 || b / a == 0.0) {
                if (x == 1.0)
                    return Double.POSITIVE_INFINITY;
                else
                    return Dpq.D0(give_log);
            }
            if (x == 0.0)
                return Double.POSITIVE_INFINITY;
            else
                return Dpq.D0(give_log);
        }

        if (x == 0.0) {
            if (a > 1.0)
                return Dpq.D0(give_log);
            if (a < 1.0)
                return Double.POSITIVE_INFINITY;
            return Dpq.DVal(b, give_log);
        }
        if (x == 1.0) {
            if (b > 1.0)
                return Dpq.D0(give_log);
            if (b < 1.0)
                return Double.POSITIVE_INFINITY;
            return Dpq.DVal(a, give_log);
        }

        double lVal;
        if (a <= 2.0 || b <= 2.0) {
            lVal = (a - 1.0) * Math.log(x) + (b - 1.0) * Math.log1p(-x) - special_functions.Beta.logBeta(a, b);
        } else {
            lVal = Math.log(a + b - 1) + binomialPdfRaw(a - 1, a + b - 2, x, 1 - x, true);
        }
        return Dpq.DExp(lVal, give_log);
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
        if (a < 0.0 || b < 0.0) {
            return Dpq.nanWarn();
        }
        if (x <= 0.0) {
            return Dpq.DT0(lower_tail, log_p);
        }
        if (x >= 1.0) {
            return Dpq.DT1(lower_tail, log_p);
        }

        return cdfRaw(x, a, b, lower_tail, log_p);
    }

    public static double quantile(double alpha, double p, double q)
    {
        return quantile(alpha, p, q, true, false);
    }

    public static double quantile(double alpha, double p, double q, boolean lower_tail, boolean log_p)
    {
        if (Double.isNaN(p) || Double.isNaN(q) || Double.isNaN(alpha)) {
            return p + q + alpha;
        }
        if (p < 0.0 || q < 0.0) {
            return Dpq.nanWarn();
        }

        double[] qBeta = new double[2];
        quantileRaw(alpha, p, q, lower_tail, log_p, qBeta);
        return qBeta[0];
    }

    private static void quantileRaw(double alpha, double p, double q,
                                    boolean lower_tail, boolean log_p, double[] qb)
    {
        boolean swap_tail;
        boolean log_;
        boolean use_log_x;
        boolean warned = false, add_N_step = true;

        int i_pb, i_inn;
        double a, la, log_beta, g, h, pp, p_, qq, r, s, t, w, y = -1.;
        double u, u_n = 0.0, xin_beta;

        if (alpha == Dpq.DT0(lower_tail, log_p)) {
            // return_q_0
            qb[0] = 0.0;
            qb[1] = 1.0;
            return;
        }

        if (alpha == Dpq.DT1(lower_tail, log_p)) {
            // return_q_1
            qb[0] = 1.0;
            qb[1] = 0.0;
            return;
        }

        if ((log_p && alpha > 0.0) || (!log_p && (alpha < 0.0 || alpha > 1.0))) {
            System.out.println("Argument out of domain");
            qb[0] = qb[1] = Double.NaN;
            return;
        }

        if (p == 0.0 || q == 0.0 || Double.isInfinite(p) || Double.isInfinite(q)) {
            if (p == 0.0 && q == 0.0) {
                if (alpha < Dpq.DHalf(log_p)) {
                    qb[0] = 0.0;
                    qb[1] = 1.0;
                    return;
                }
                if (alpha > Dpq.DHalf(log_p)) {
                    qb[0] = 1.0;
                    qb[1] = 0.0;
                    return;
                }
                // return_q_half
                qb[0] = qb[1] = 0.5;
                return;
            } else if (p == 0.0 || p / q == 0.0) {
                qb[0] = 0.0;
                qb[1] = 1.0;
                return;
            } else if (q == 0.0 || q / p == 0.0) {
                qb[0] = 1.0;
                qb[1] = 0.0;
                return;
            }
            qb[0] = qb[1] = 0.5;
            return;
        }

        p_ = Dpq.DTqIv(alpha, lower_tail, log_p);
        log_beta = special_functions.Beta.logBeta(p, q);
        swap_tail = p_ > 0.5;

        if (swap_tail) {
            a = Dpq.DTCIv(alpha, lower_tail, log_p);
            la = Dpq.DTCLog(alpha, lower_tail, log_p);
            pp = q; qq = p;
        } else {
            a = p_;
            la = Dpq.DT_log(alpha, lower_tail, log_p);
            pp = p; qq = q;
        }

        double acu = Math.max(acu_min, Math.pow(10.0, -13.0 - 2.5 / (pp * pp) - 0.5 / (a * a)));
        double tx, u0 = (la + Math.log(pp) + log_beta) / pp;
        double log_eps_c = Math.log(Dpq.DBL_EPSILON);
        r = pp * (1.0 - qq) / (pp + 1.0);
        t = 0.2;

        if (Dpq.M_LN2 * Double.MIN_EXPONENT < u0 && u0 < -0.01 &&
            u0 < (t * log_eps_c - Math.log(Math.abs(pp * (1.0 - qq) * (2.0 - qq) / (2.0 * (pp + 2.0))))) / 2.0) {
            r = r * Math.exp(u0);
            if (r > -1.0) {
                u = u0 - Math.log1p(r) / pp;
            } else {
                u = u0;
            }
            xin_beta = Math.exp(u);

            //goto L_Newton;
            boolean flag = true;

            // L_Newton
            L_converged:
            do {
                r = 1.0 - pp;
                t = 1.0 - qq;
                double w_prev = 0.0, prev = 1.0, adj = 1.0;

                for (i_pb = 0; i_pb < 1000; i_pb++) {
                    y = cdfRaw(xin_beta, pp, qq, true, true);
                    w = (y == Double.NEGATIVE_INFINITY)
                            ? 0. : (y - la) * Math.exp(y - u + log_beta + r * u + t * Dpq.log1Exp(u));
                    if (Double.isInfinite(w))
                        break;
                    if (i_pb >= Beta.N_NEWTON_FREE && w * w_prev <= 0.)
                        prev = Math.max(Math.abs(adj), fpu);

                    g = 1;
                    for (i_inn = 0; i_inn < 1000; i_inn++) {
                        adj = g * w;
                        if (Math.abs(adj) < prev) {
                            u_n = u - adj; // u_{n+1} = u_n - g*w
                            if (u_n <= 0.) {
                                if (prev <= acu || Math.abs(w) <= acu) {
                                    flag = false;
                                    break L_converged;
                                }
                                break;
                            }
                        }
                        g /= 3;
                    }
                    double D = Math.min(Math.abs(adj), Math.abs(u_n - u));

                    if (D <= 4e-16 * Math.abs(u_n + u)) {
                        flag = false;
                        break L_converged;
                    }
                    u = u_n;
                    xin_beta = Math.exp(u);
                    w_prev = w;
                } // for(i )

            } while (false);

            if (flag) {
                warned = true;
                System.out.println("Full precision may not have been achieved in Beta.quantile()");
            }


            // L_converged
            if (y == Double.NEGATIVE_INFINITY) {
                u_n = DBL_LOG_V_MIN;// = log(DBL_very_MIN)
                add_N_step = false; // not trying to do better anymore
            }
            else if(!warned && Math.abs(y - la) > 3) {
                if(!(y == Double.NEGATIVE_INFINITY && cdfRaw(DBL_1_EPS, pp, qq, true, true) > la + 2))
                    System.out.printf("Beta.quantile(a, *) =: x0 with | Beta.cdf(x0,*%s) - alpha| = %.5g " +
                                    "is not accurate", ", log_", Math.abs(y - la));
            }



            // L_return
            if(add_N_step) {
                xin_beta = Math.exp(u_n);
                y = cdfRaw(xin_beta, pp, qq, true, log_p);
                w = log_p
                        ? (y - la) * Math.exp(y + log_beta + r * Math.log(xin_beta) + t * Math.log1p(-xin_beta))
                        : (y - a)  * Math.exp(log_beta + r * Math.log(xin_beta) + t * Math.log1p(-xin_beta));
                tx = xin_beta - w;

            } else {
                if(swap_tail) {
                    qb[0] = -Math.expm1(u_n); qb[1] = Math.exp(u_n);
                } else {
                    qb[0] = Math.exp(u_n); qb[1] = -Math.expm1(u_n);
                }
                return;
            }
            if(swap_tail) {
                qb[0] = 1 - tx; qb[1] = tx;
            } else {
                qb[0] = tx;	qb[1] = 1 - tx;
            }
            return;

        }

        r = Math.sqrt(-2.8 * la);
        y = r - (const1 + const2 * r) / (1.0 + (const3 + const4 * r) * r);

        if (pp > 1.0 && qq > 1.0) {
            r = (y * y - 3.) / 6.;
            s = 1. / (pp + pp - 1.);
            t = 1. / (qq + qq - 1.);
            h = 2. / (s + t);
            w = y * Math.sqrt(h + r) / h - (t - s) * (r + 5. / 6. - 2. / (3. * h));

            if(w > 300) { // exp(w+w) is huge or overflows
                t = w + w + Math.log(qq) - Math.log(pp);
                u = (t <= 18) ? - Math.log1p(Math.exp(t)) : -t - Math.exp(-t);
                xin_beta = Math.exp(u);
            } else {
                xin_beta = pp / (pp + qq * Math.exp(w + w));
                u = -Math.log1p(qq / pp * Math.exp(w + w));
            }
        } else {
            r = qq + qq;
            t = 1. / (3. * Math.sqrt(qq));
            t = r * Math.pow(1.0 + t * (-t + y), 3);// = \chi^2_{alpha} of AS 64
            s = 4. * pp + r - 2.;// 4p + 2q - 2 = numerator of new t = (...) / chi^2

            if (t == 0 || (t < 0. && s >= t)) {
                double l1ma;
                if(swap_tail)
                    l1ma = Dpq.DT_log(alpha, lower_tail, log_p);
                else
                    l1ma = Dpq.DTCLog(alpha, lower_tail, log_p);

                double xx = (l1ma + Math.log(qq) + log_beta) / qq;
                if(xx <= 0.) {
                    xin_beta = -Math.expm1(xx);
                    u = Dpq.log1Exp(xx);
                } else { // xx > 0 ==> 1 - e^xx < 0 .. is nonsense

                    xin_beta = 0; u = Double.NEGATIVE_INFINITY;
                }
            } else {
                t = s / t;

                if (t <= 1.) {
                    u = u0;
                    xin_beta = Math.exp(u);
                } else { // (1+x0)/(1-x0) = t,  solved for x0 :
                    xin_beta = 1. - 2. / (t + 1.);
                    u = Math.log1p(-2. / (t + 1.));
                }
            }
        }

        if (swap_tail && u >= -Math.exp(Beta.USE_LOG_X_CUTOFF) ||
                !swap_tail && u >= -Math.exp(4.0 * Beta.USE_LOG_X_CUTOFF) && pp / qq < 1000.0) {
            swap_tail = !swap_tail;

            if(swap_tail) { // "swap now" (much less easily)
                a = Dpq.DTCIv(alpha, lower_tail, log_p); // needed ?
                la = Dpq.DTCLog(alpha, lower_tail, log_p);
                pp = q; qq = p;
            }
            else { // swap back :
                a = p_;
                la = Dpq.DT_log(alpha, lower_tail, log_p);
                pp = p; qq = q;
            }

            // we could redo computations above, but this should be stable
            u = Dpq.log1Exp(u);
            xin_beta = Math.exp(u);
        }

        use_log_x = u < Beta.USE_LOG_X_CUTOFF;

        boolean bad_u = Double.isInfinite(u);
        boolean bad_init = bad_u || xin_beta > p_hi;

        u_n = 1.0;
        tx = xin_beta;

        if (bad_u || u < Beta.USE_LOG_X_CUTOFF) {
            w = cdfRaw(DBL_VERY_MIN, pp, qq, true, log_p);
            if (w > (log_p ? la : a)) {

                if (log_p || Math.abs(w - a) < Math.abs(0 - a)) { // DBL_very_MIN is better than 0
                    tx   = DBL_VERY_MIN;
                    u_n  = DBL_LOG_V_MIN;// = log(DBL_very_MIN)
                } else {
                    tx   = 0.0;
                    u_n  = Double.NEGATIVE_INFINITY;
                }
                use_log_x = log_p;

                // L_return
                if(use_log_x) {
                    if(swap_tail) {
                        qb[0] = -Math.expm1(u_n); qb[1] = Math.exp(u_n);
                    } else {
                        qb[0] = Math.exp(u_n); qb[1] = -Math.expm1(u_n);
                    }
                    return;
                }
                if(swap_tail) {
                    qb[0] = 1 - tx; qb[1] = tx;
                } else {
                    qb[0] = tx;	qb[1] = 1 - tx;
                }
                return;
                // End L_return

            } else {
                if(u  < DBL_LOG_V_MIN) {
                    u = DBL_LOG_V_MIN;// = log(DBL_very_MIN)
                    xin_beta = DBL_VERY_MIN;
                }
            }
        }

        if (bad_init && !(use_log_x && tx > 0.0)) {
            if(u == Double.NEGATIVE_INFINITY) {
                u = Dpq.M_LN2 * Double.MIN_EXPONENT;
                xin_beta = Double.MIN_VALUE;
            } else {
                xin_beta = (xin_beta > 1.1) // i.e. "way off"
                        ? 0.5 // otherwise, keep the respective boundary:
                        : ((xin_beta < p_lo) ? Math.exp(u) : p_hi);
                if(bad_u)
                    u = Math.log(xin_beta);
                // otherwise: not changing "potentially better" u than the above
            }
        }

        boolean flag = true;

        // L_Newton
        L_converged:
        do {
            r = 1.0 - pp;
            t = 1.0 - qq;
            double w_prev = 0.0, prev = 1.0, adj = 1.0;

            if (use_log_x) {
                for (i_pb = 0; i_pb < 1000; i_pb++) {
                    y = cdfRaw(xin_beta, pp, qq, true, true);
                    w = (y == Double.NEGATIVE_INFINITY)
                            ? 0. : (y - la) * Math.exp(y - u + log_beta + r * u + t * Dpq.log1Exp(u));
                    if (Double.isInfinite(w))
                        break;
                    if (i_pb >= Beta.N_NEWTON_FREE && w * w_prev <= 0.)
                        prev = Math.max(Math.abs(adj), fpu);

                    g = 1;
                    for (i_inn = 0; i_inn < 1000; i_inn++) {
                        adj = g * w;
                        if (Math.abs(adj) < prev) {
                            u_n = u - adj; // u_{n+1} = u_n - g*w
                            if (u_n <= 0.) {
                                if (prev <= acu || Math.abs(w) <= acu) {
                                    flag = false;
                                    break L_converged;
                                }
                                break;
                            }
                        }
                        g /= 3;
                    }
                    double D = Math.min(Math.abs(adj), Math.abs(u_n - u));

                    if (D <= 4e-16 * Math.abs(u_n + u)) {
                        flag = false;
			            break L_converged;
                    }
                    u = u_n;
                    xin_beta = Math.exp(u);
                    w_prev = w;
                } // for(i )

            } else { // "normal scale" Newton
                for (i_pb = 0; i_pb < 1000; i_pb++) {
                    y = cdfRaw(xin_beta, pp, qq, true, log_p);

                    if (Double.isInfinite(y) && !(log_p && y == Double.NEGATIVE_INFINITY)) {
                        System.out.println("Argument out of domain");
                        qb[0] = qb[1] = Double.NaN;
                        return;
                    }

                    w = log_p
                            ? (y - la) * Math.exp(y + log_beta + r * Math.log(xin_beta) + t * Math.log1p(-xin_beta))
                            : (y - a) * Math.exp(log_beta + r * Math.log(xin_beta) + t * Math.log1p(-xin_beta));

                    if (i_pb >= Beta.N_NEWTON_FREE && w * w_prev <= 0.0)
                        prev = Math.max(Math.abs(adj), fpu);

                    g = 1;
                    for (i_inn = 0; i_inn < 1000; i_inn++) {
                        adj = g * w;
                        if (i_pb < Beta.N_NEWTON_FREE || Math.abs(adj) < prev) {
                            tx = xin_beta - adj; // x_{n+1} = x_n - g*w
                            if (0. <= tx && tx <= 1.) {
                                if (prev <= acu || Math.abs(w) <= acu) {
                                    flag = false;
                                    break L_converged;
                                }
                                if (tx != 0. && tx != 1)
                                    break;
                            }
                        }
                        g /= 3;
                    }

                    if (Math.abs(tx - xin_beta) <= 4e-16 * (tx + xin_beta)) {
				        flag = false;
                        break L_converged;
                    }

                    xin_beta = tx;
                    if (tx == 0) // "we have lost"
                        break;
                    w_prev = w;
                }
            }
        } while (false);

        if (flag) {
            warned = true;
            System.out.println("Full precision may not have been achieved in Beta.quantile()");
        }


        // L_converged
        log_ = log_p || use_log_x; // only for printing

        if ((log_ && y == Double.NEGATIVE_INFINITY) || (!log_ && y == 0)) {
            w = cdfRaw(DBL_VERY_MIN, pp, qq, true, log_);
            if (log_ || Math.abs(w - a) <= Math.abs(y - a)) {
                tx  = DBL_VERY_MIN;
                u_n = DBL_LOG_V_MIN;// = log(DBL_very_MIN)
            }
            add_N_step = false; // not trying to do better anymore
        }
        else if (!warned && (log_ ? Math.abs(y - la) > 3 : Math.abs(y - a) > 1e-4)) {
            if (!(log_ && y == Double.NEGATIVE_INFINITY &&
                    cdfRaw(DBL_1_EPS, pp, qq, true, true) > la + 2))
                System.out.printf("Beta.quantile(a, *) =: x0 with | Beta.cdf(x0,*%s) - alpha| = %.5g is " +
                        "not accurate", (log_ ? ", log_" : ""), Math.abs(y - (log_ ? la : a)));
        }

        // L_return
        if (use_log_x) {
            if (add_N_step) {
                xin_beta = Math.exp(u_n);
                y = cdfRaw(xin_beta, pp, qq, true, log_p);
                w = log_p
                        ? (y - la) * Math.exp(y + log_beta + r * Math.log(xin_beta) + t * Math.log1p(-xin_beta))
                        : (y - a)  * Math.exp(log_beta + r * Math.log(xin_beta) + t * Math.log1p(-xin_beta));
                tx = xin_beta - w;

            } else {
                if (swap_tail) {
                    qb[0] = -Math.expm1(u_n); qb[1] = Math.exp(u_n);
                } else {
                    qb[0] = Math.exp(u_n); qb[1] = -Math.expm1(u_n);
                }
                return;
            }
        }
        if (swap_tail) {
            qb[0] = 1 - tx; qb[1] = tx;
        } else {
            qb[0] = tx;	qb[1] = 1 - tx;
        }
    }

    private static double cdfRaw(double x, double a, double b, boolean lower_tail, boolean log_p)
    {
        if (a == 0.0 || b == 0.0 || Double.isInfinite(a) || Double.isInfinite(b)) {
            if (a == 0.0 && b == 0.0) {
                return log_p ? -Dpq.M_LN2 : 0.5;
            }
            if (a == 0.0 || a / b == 0.0) {
                return Dpq.DT1(lower_tail, log_p);
            }
            if (b == 0.0 || b / a == 0.0) {
                return Dpq.DT0(lower_tail, log_p);
            }
            return (x < 0.5) ? Dpq.DT0(lower_tail, log_p) : Dpq.DT1(lower_tail, log_p);
        }

        double xes1 = 0.5 - x + 0.5;
        RefInt i_err = new RefInt(0);
        RefDouble w = new RefDouble(0.0);
        RefDouble wc = new RefDouble(0.0);

        Toms.betaRatio(a, b, x, xes1, w, wc, i_err, log_p);
        if (i_err.val != 0 && i_err.val != 11 && i_err.val != 14) {
            System.out.printf(
                    "Beta.cdfRaw(%g, a = %g, b = %g, ..) -> betaRatio() gave error code %d",
                    x, a, b, i_err.val);
        }
        return lower_tail ? w.val : wc.val;
    }

    public static void main(String[] args)
    {
        System.out.println(pdf(0.3, 2, 10));
        System.out.println(Beta.cdf(0.975, 10, 2));
        System.out.println(Beta.cdf(0.005, 79, 50, false, true));
        System.out.println(Beta.quantile(0.992510, 0.910, 2.0));
    }
}
