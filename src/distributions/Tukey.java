package distributions;

import distributions.detail.Dpq;
import special_functions.LogGamma;

import java.math.BigDecimal;

public class Tukey  // Studentized Range Distribution
{
    public static double cdf(double q, double nRanges, double nMeans, double df)
    {
        return cdf(q, nRanges, nMeans, df, true, false);
    }

    public static double cdf(double q, double nRanges, double nMeans, double df, boolean lower_tail)
    {
        return cdf(q, nRanges, nMeans, df, lower_tail, false);
    }

    public static double cdf(double q, double nRanges, double nMeans, double df, boolean lower_tail, boolean log_p)
    {
        final int n_leg_q	= 16;
        final int i_half_q = 8;

        final double eps1 = -30.0;
        final double eps2 = 1.0e-14;
        final double d_haf  = 100.0;
        final double d_qua_r = 800.0;
        final double dei_gh = 5000.0;
        final double dl_arg = 25000.0;
        final double u_len1 = 1.0;
        final double u_len2 = 0.5;
        final double u_len3 = 0.25;
        final double u_len4 = 0.125;
        final double[] x_leg_q = {
            0.989400934991649932596154173450,
            0.944575023073232576077988415535,
            0.865631202387831743880467897712,
            0.755404408355003033895101194847,
            0.617876244402643748446671764049,
            0.458016777657227386342419442984,
            0.281603550779258913230460501460,
            0.950125098376374401853193354250e-1
        };
        final double[] a_leg_q = {
            0.271524594117540948517805724560e-1,
            0.622535239386478928628438369944e-1,
            0.951585116824927848099251076022e-1,
            0.124628971255533872052476282192,
            0.149595988816576732081501730547,
            0.169156519395002538189312079030,
            0.182603415044923588866763667969,
            0.189450610455068496285396723208
        };

        if (Double.isNaN(q) || Double.isNaN(nRanges) || Double.isNaN(nMeans) || Double.isNaN(df)) {
            return Dpq.nanWarn();
        }
        if (q <= 0) {
            return Dpq.DT0(lower_tail, log_p);
        }
        if (df < 2 || nRanges < 1 || nMeans < 2) {
            return Dpq.nanWarn();
        }
        if (Double.isInfinite(q)) {
            return Dpq.DT1(lower_tail, log_p);
        }
        if (df > dl_arg) {
            return Dpq.DTVal(wProb(q, nRanges, nMeans), lower_tail, log_p);
        }

        // calculate leading constant
        double f2 = df * 0.5;
        double f2lf = ((f2 * Math.log(df)) - (df * Dpq.M_LN2)) - LogGamma.logGamma(f2);
        double f21 = f2 - 1.0;
        double ff4 = df * 0.25;
        double uLen;
        if (df <= d_haf) {
            uLen = u_len1;
        } else if (df <= d_qua_r) {
            uLen = u_len2;
        } else if (df <= dei_gh) {
            uLen = u_len3;
        } else {
            uLen = u_len4;
        }

        f2lf += Math.log(uLen);

        // integrate over each sub-interval
        double ans = 0.0;
        double otSum = 0.0;
        for (int i = 1; i <= 50; i++) {
            otSum = 0.0;
            double twa1 = (2 * i - 1) * uLen;

            for (int jj = 1; jj <= n_leg_q; jj++) {
                int j;
                double t1;
                if (i_half_q < jj) {
                    j = jj - i_half_q - 1;
                    t1 = (f2lf + (f21 * Math.log(twa1 + (x_leg_q[j] * uLen)))) - (((x_leg_q[j] * uLen) + twa1) * ff4);
                } else {
                    j = jj - 1;
                    t1 = (f2lf + (f21 * Math.log(twa1 - (x_leg_q[j] * uLen)))) + (((x_leg_q[j] * uLen) - twa1) * ff4);
                }

                double q_sqz;
                if (t1 >= eps1) {
                    if (i_half_q < jj) {
                        q_sqz = q * Math.sqrt(((x_leg_q[j] * uLen) + twa1) * 0.5);
                    } else {
                        q_sqz = q * Math.sqrt(((-(x_leg_q[j] * uLen)) + twa1) * 0.5);
                    }
                    double wPrb = wProb(q_sqz, nRanges, nMeans);
                    double rotSum = (wPrb * a_leg_q[j]) * Math.exp(t1);
                    otSum += rotSum;
                }
            }

            if (i * uLen >= 1.0 && otSum <= eps2) {
                break;
            }
            ans += otSum;
        }

        if(otSum > eps2) {
            System.out.println("Full precision may not have been achieved in Tukey.cdf()");
        }
        if (ans > 1.0) {
            ans = 1.0;
        }
        return Dpq.DTVal(ans, lower_tail, log_p);
    }

    private static double wProb(double w, double rr, double cc)
    {
        final int nLeg = 12;
        final int iHalf = 6;

        final double C1 = -30.0;
        final double C2 = -50.0;
        final double C3 = 60.0;
        final double bb   = 8.0;
        final double wLar = 3.0;
        final double winCr1 = 2.0;
        final double winCr2 = 3.0;
        final double[] xLeg = {
            0.981560634246719250690549090149,
            0.904117256370474856678465866119,
            0.769902674194304687036893833213,
            0.587317954286617447296702418941,
            0.367831498998180193752691536644,
            0.125233408511468915472441369464
        };
        final double[] aLeg = {
            0.047175336386511827194615961485,
            0.106939325995318430960254718194,
            0.160078328543346226334652529543,
            0.203167426723065921749064455810,
            0.233492536538354808760849898925,
            0.249147045813402785000562436043
        };

        double qs_qz = w * 0.5;
        if (qs_qz >= bb) {
            return 1.0;
        }

        double pr_w = 2 * Normal.cdf(qs_qz, 0.0,1.0, true,false) - 1.0;
        if (pr_w >= Math.exp(C2 / cc)) {
            pr_w = Math.pow(pr_w, cc);
        } else {
            pr_w = 0.0;
        }

        double winCr;
        if (w > wLar) {
            winCr = winCr1;
        } else {
            winCr = winCr2;
        }

        double blb = qs_qz;
        double binC = (bb - qs_qz) / winCr;
        double bub = blb + binC;
        double einSum = 0.0;

        // integrate over each interval
        double cc1 = cc - 1.0;
        for (double wi = 1.0; wi <= winCr; wi++) {
            double elSum = 0.0;
            double a = 0.5 * (bub + blb);

            // legendre quadrature with order = nLeg
            double b = 0.5 * (bub - blb);
            for (int jj = 1; jj <= nLeg; jj++) {
                int j;
                double xx;
                if (iHalf < jj) {
                    j = (nLeg - jj) + 1;
                    xx = xLeg[j-1];
                } else {
                    j = jj;
                    xx = -xLeg[j-1];
                }

                double c = b * xx;
                double ac = a + c;
                double qExpo = ac * ac;
                if (qExpo > C3)
                    break;

                double pPlus = 2.0 * Normal.cdf(ac, 0.0, 1.0, true,false);
                double pMinus = 2.0 * Normal.cdf(ac, w,  1.0, true, false);
                double rinSum = (pPlus * 0.5) - (pMinus * 0.5);
                if (rinSum >= Math.exp(C1 / cc1)) {
                    rinSum = (aLeg[j-1] * Math.exp(-(0.5 * qExpo))) * Math.pow(rinSum, cc1);
                    elSum += rinSum;
                }
            }
            elSum *= (((2.0 * b) * cc) * Dpq.M_1_SQRT_2PI);
            einSum += elSum;
            blb = bub;
            bub += binC;
        }

        pr_w += einSum;
        if (pr_w <= Math.exp(C1 / rr)) {
            return 0.0;
        }
        pr_w = Math.pow(pr_w, rr);
        return Math.min(pr_w, 1.0);
    }

    public static void main(String[] args)
    {
        // Maybe have precision's problem after 15th digit point
        System.out.println(cdf(7, 10, 60, 5));
        System.out.println(cdf(70, 100, 60, 5));
    }
}
