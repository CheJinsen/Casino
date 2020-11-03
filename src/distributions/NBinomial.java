package distributions;

import distributions.detail.*;
import special_functions.LogGamma;

public class NBinomial extends DistBase  // Negative Binomial Distribution
{
    public static double pdf(double x, double size, double prob)
    {
        return pdf(x, size, prob, false);
    }

    public static double pdf(double x, double size, double prob, boolean give_log)
    {
        if (Double.isNaN(x) || Double.isNaN(size) || Double.isNaN(prob)) {
            return x + size + prob;
        }

        if (prob <= 0 || prob > 1 || size < 0) {
            return Dpq.nanWarn();
        }
        if (Dpq.nonInt(x)) {
            System.out.printf("non-integer x = %f", x);
            return Dpq.D0(give_log);
        }
        if (x < 0 || Double.isInfinite(x)) {
            return Dpq.D0(give_log);
        }
        if (x == 0 && size==0) {
            return Dpq.D1(give_log);
        }

        x = Math.round(x);
        if (Double.isInfinite(size)) {
            size = Double.MAX_VALUE;
        }

        double ans = binomialPdfRaw(size, x+size, prob, 1-prob, give_log);
        double p = size / (size + x);
        return (give_log) ? Math.log(p) + ans : p * ans;
    }

    public static double pdfMu(double x, double size, double mu)
    {
        return pdfMu(x, size, mu, false);
    }

    public static double pdfMu(double x, double size, double mu, boolean give_log)
    {
        if (Double.isNaN(x) || Double.isNaN(size) || Double.isNaN(mu)) {
            return x + size + mu;
        }
        if (mu < 0 || size < 0) {
            return Dpq.nanWarn();
        }
        if (Dpq.nonInt(x)) {
            System.out.printf("non-integer x = %f", x);
            return Dpq.D0(give_log);
        }
        if (x < 0 || Double.isInfinite(x)) {
            return Dpq.D0(give_log);
        }

        /* limiting case as size approaches zero is point mass at zero,
         * even if mu is kept constant. limit distribution does not
         * have mean mu, though.
         */
        if (x == 0 && size == 0) {
            return Dpq.D1(give_log);
        }
        x = Math.round(x);
        if (Double.isInfinite(size)) {
            return poissonPdfRaw(x, mu, give_log);
        }

        if (x == 0) {
            return Dpq.DExp(size * (size < mu
                    ? Math.log(size / (size + mu)) : Math.log1p(-mu / (size + mu))),
                    give_log);
        }
        if (x < 1e-10 * size) {
            double p = (size < mu ? Math.log(size / (1.0 + size/mu)) : Math.log(mu / (1.0 + mu/size)));
            return Dpq.DExp(x * p - mu - LogGamma.logGamma(x+1) +
                    Math.log1p(x * (x - 1) / (2 * size)), give_log);
        } else {
            double p = size / (size + x);
            double ans = binomialPdfRaw(size, x+size, size/(size+mu), mu/(size+mu), give_log);
            return (give_log) ? Math.log(p) + ans : p * ans;
        }
    }

    public static double cdf(double x, double size, double prob)
    {
        return cdf(x, size, prob, true, false);
    }

    public static double cdf(double x, double size, double prob, boolean lower_tail)
    {
        return cdf(x, size, prob, lower_tail, false);
    }

    public static double cdf(double x, double size, double prob, boolean lower_tail, boolean log_p)
    {
        if (Double.isNaN(x) || Double.isNaN(size) || Double.isNaN(prob)) {
            return x + size + prob;
        }
        if(Double.isInfinite(size) || Double.isInfinite(prob)) {
            return Dpq.nanWarn();
        }
        if (size < 0 || prob <= 0 || prob > 1) {
            return Dpq.nanWarn();
        }
        if (size == 0) {
            return (x >= 0) ? Dpq.DT1(lower_tail, log_p) : Dpq.DT0(lower_tail, log_p);
        }
        if (x < 0) {
            return Dpq.DT0(lower_tail, log_p);
        }
        if (Double.isInfinite(x)) {
            return Dpq.DT1(lower_tail, log_p);
        }

        x = Math.floor(x + 1e-7);
        return Beta.cdf(prob, size, x + 1, lower_tail, log_p);
    }

    public static double cdfMu(double x, double size, double mu)
    {
        return cdfMu(x, size, mu, true, false);
    }

    public static double cdfMu(double x, double size, double mu, boolean lower_tail)
    {
        return cdfMu(x, size, mu, lower_tail, false);
    }

    public static double cdfMu(double x, double size, double mu, boolean lower_tail, boolean log_p)
    {
        if (Double.isNaN(x) || Double.isNaN(size) || Double.isNaN(mu)) {
            return x + size + mu;
        }
        if (Double.isInfinite(mu)) {
            return Dpq.nanWarn();
        }
        if (size < 0 || mu < 0) {
            return Dpq.nanWarn();
        }
        if (size == 0) {
            return (x >= 0) ? Dpq.DT1(lower_tail, log_p) : Dpq.DT0(lower_tail, log_p);
        }

        if (x < 0) {
            return Dpq.DT0(lower_tail, log_p);
        }
        if (Double.isInfinite(x)) {
            return Dpq.DT1(lower_tail, log_p);
        }
        if (Double.isInfinite(size)) {
            return Poisson.cdf(x, mu, lower_tail, log_p);
        }

        x = Math.floor(x + 1e-7);

        {
            RefInt i_err = new RefInt(0);
            RefDouble w = new RefDouble(0.0);
            RefDouble wc = new RefDouble(0.0);
            Toms.betaRatio(size, x+1, size/(size+mu), mu/(size+mu), w, wc, i_err, log_p);
            if (i_err.val != 0) {
                System.out.printf("NBinomial.cdfMu() -> betaRatio() gave error code %d", i_err.val);
            }
            return lower_tail ? w.val : wc.val;
        }
    }

    public static double rand(double size, double prob)
    {
        if (Double.isInfinite(prob) || Double.isNaN(size) || size <= 0.0 || prob <= 0.0 || prob > 1.0) {
            return Dpq.nanWarn();
        }
        if (Double.isInfinite(size)) {
            size = Double.MAX_VALUE / 2.0;
        }
        return (prob == 1.0) ? 0.0 : Poisson.rand(Gamma.rand(size, (1.0 - prob) / prob));
    }

    public static double randMu(double size, double mu)
    {
        if (Double.isInfinite(mu) || Double.isNaN(size) || size <= 0.0 || mu < 0.0) {
            return Dpq.nanWarn();
        }
        if (Double.isInfinite(size)) {
            size = Double.MAX_VALUE / 2.0;
        }
        return (mu == 0.0) ? 0.0 : Poisson.rand(Gamma.rand(size, mu / size));
    }

    public static void main(String[] args)
    {
        System.out.println(pdf(9, 2, 0.23));
        System.out.println(pdfMu(9, 2, 5.6));
        System.out.println(pdf(9, 2, 5.6));
        System.out.println(cdfMu(90, 2, 5.6));
        System.out.println(cdfMu(90, 2, 5.6, false));

        for (int i = 0; i < 10; i++) {
            System.out.print(randMu(7, 0.8) + " ");
        }
        System.out.println();
    }
}
