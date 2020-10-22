package distributions;

import distributions.detail.*;

public class Beta {
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

        double x1 = 0.5 - x + 0.5;
        RefInt ierr = new RefInt(0);
        RefDouble w = new RefDouble(0.0);
        RefDouble wc = new RefDouble(0.0);

        Toms.bratio(a, b, x, x1, w, wc, ierr, log_p);
        if (ierr.val != 0 && ierr.val != 11 && ierr.val != 14) {
            System.out.printf("Beta.cdfRaw(%g, a = %g, b = %g, ..) -> bratio() gave error code %d", x, a, b, ierr.val);
        }
        return lower_tail ? w.val : wc.val;
    }
}
