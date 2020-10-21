package specfun;

import randist.detail.Dpq;
import specfun.detail.Lgammacor;

public class Beta {
    public static double logBeta(double a, double b)
    {
        if (Double.isNaN(a) || Double.isNaN(b)) {
            return a + b;
        }

        double p = a;
        double q = a;
        if (b < p) { p = b; }
        if (b > q) { q = b; }

        if (p < 0.0) {
            Dpq.nanWarn();
        } else if (p == 0.0) {
            return Double.POSITIVE_INFINITY;
        } else if (Double.isInfinite(q)) {
            return Double.NEGATIVE_INFINITY;
        }

        double corr;
        double ret = Math.log1p(-p / (p + q));
        if (p >= 10.0) {
            corr = Lgammacor.lgammacor(p) + Lgammacor.lgammacor(q) - Lgammacor.lgammacor(p + q);
            return Math.log(q) * -0.5 + Dpq.M_LN_SQRT_2PI + corr +
                    (p - 0.5) * Math.log(p / (p + q)) + q * ret;
        } else if (q >= 10.0) {
            corr = Lgammacor.lgammacor(q) - Lgammacor.lgammacor(p + q);
            return Lgamma.lammafn(p) + corr + p - p * Math.log(p + q) + (q - 0.5) * ret;
        } else {
            if (p < 1e-306) {
                return Lgamma.lammafn(p) + (Lgamma.lammafn(q) - Lgamma.lammafn(p + q));
            } else {
                return Math.log(Gamma.gammafn(p) * (Gamma.gammafn(q) / Gamma.gammafn(p + q)));
            }
        }
    }

    // test session
    public static void main(String[] args)
    {
        System.out.println(logBeta(0.31, 0.034));
    }
}
