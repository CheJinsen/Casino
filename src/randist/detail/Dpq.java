package randist.detail;

public class Dpq {
    public final static double M_LN2 = 0.693147180559945309417232121458;	/* ln(2) */
    public final static double DBL_EPSILON = 2.2204460492503131E-16;

    public static double nanWarn()
    {
        System.out.println("Argument out of domain.");
        return Double.NaN;
    }

    public static double D0(boolean log_p)
    {
        return log_p ? Double.NEGATIVE_INFINITY : 0.0;
    }

    public static double D1(boolean log_p)
    {
        return log_p ? 0.0 : 1.0;
    }

    public static double DT0(boolean lower_tail, boolean log_p)
    {
        return lower_tail ? D0(log_p) : D1(log_p);
    }

    public static double DT1(boolean lower_tail, boolean log_p)
    {
        return lower_tail ? D1(log_p) : D0(log_p);
    }

    public static double DHalf(boolean log_p)
    {
        return log_p ? Math.log(2.0) : 0.5;
    }

    public static double DLVal(double p, boolean lower_tail)
    {
        return lower_tail ? p : 1.0 - p;
    }

    public static double DCVal(double p, boolean lower_tail)
    {
        return lower_tail ? 1.0 - p : p;
    }

    public static double DVal(double x, boolean log_p)
    {
        return log_p ? Math.log(x) : x;
    }

    public static double DqIv(double p, boolean log_p)
    {
        return log_p ? Math.exp(p) : p;
    }

    public static double DExp(double x, boolean log_p)
    {
        return log_p ? x : Math.exp(x);
    }

    public static double DLog(double p, boolean log_p)
    {
        return log_p ? p : Math.log(p);
    }

    public static double DCLog(double p, boolean log_p)
    {
        return log_p ? Math.log1p(-p) : 1.0 - p;
    }

    public static double log1Exp(double x)
    {
        return x > Math.log(2.0) ? Math.log(-Math.expm1(x)) : Math.log1p(-Math.exp(x));
    }

    public static double DLExp(double x, boolean log_p)
    {
        return log_p ? log1Exp(x) : Math.log1p(-x);
    }

    public static double DTVal(double x, boolean lower_tail, boolean log_p)
    {
        return lower_tail ? DVal(x, log_p) : DCLog(x, log_p);
    }

    public static double DTCVal(double x, boolean lower_tail, boolean log_p)
    {
        return lower_tail ? DCLog(x, log_p) : DVal(x, log_p);
    }

    public static double DTqIv(double p, boolean lower_tail, boolean log_p)
    {
        return log_p ? (lower_tail ? Math.exp(p) : -Math.expm1(p))
                : DLVal(p, lower_tail);
    }

    public static double DTCIv(double p, boolean lower_tail, boolean log_p)
    {
        return log_p ? (lower_tail ? -Math.expm1(p) : Math.exp(p))
                : DCVal(p, lower_tail);
    }

    public static double DTExp(double x, boolean lower_tail, boolean log_p)
    {
        return DExp(DLVal(x, lower_tail), log_p);
    }

    public static double DTCExp(double x, boolean lower_tail, boolean log_p)
    {
        return DExp(DCVal(x, lower_tail), log_p);
    }

    public static double DT_log(double p, boolean lower_tail, boolean log_p)
    {
        return lower_tail ? DLog(p, log_p) : DLExp(p, log_p);
    }

    public static double DTCLog(double p, boolean lower_tail, boolean log_p)
    {
        return lower_tail ? DLExp(p, log_p) : DLog(p, log_p);
    }

    public static double DT_Log(double p, boolean lower_tail)
    {
        return lower_tail ? p : log1Exp(p);
    }

    public static double QP01Check(double p, boolean log_p)
    {
        if ((log_p && p > 0) || (!log_p && (p < 0 || p > 1)))
            return Double.NaN;
        return 0.0;
    }

    public static double QP01Boundaries(double p, double l, double r, boolean lower_tail, boolean log_p)
    {
        if (log_p) {
            if (p > 0)
                return Double.NaN;
            if (p == 0) {
                return lower_tail ? r : l;
            }
            if (p == Double.NEGATIVE_INFINITY) {
                return lower_tail ? l : r;
            }
        } else {
            if (p < 0 || p > 1)
                return Double.NaN;
            if (p == 0) {
                return lower_tail ? l : r;
            }
            if (p == 1) {
                return lower_tail ? r : l;
            }
        }
        return 0.0;
    }

    public static double pBounds01(double x, double xMin, double xMax, boolean lower_tail, boolean log_p)
    {
        if (x <= xMin)  return DT0(lower_tail, log_p);
        if (x >= xMax)  return DT1(lower_tail, log_p);
        return 0.0;
    }

    public static double pBoundsInf01(double x, boolean lower_tail, boolean log_p)
    {
        if (Double.isInfinite(x)) {
            if (x > 0)  return DT1(lower_tail, log_p);
            else        return DT0(lower_tail, log_p);
        }
        return 0.0;
    }

    public static double DFExp(double f, double x, boolean give_log)
    {
        if (give_log)   return -0.5 * Math.log(f) + x;
        else            return Math.exp(x) / Math.sqrt(f);
    }

    public static boolean negINonInt(double x)
    {
        return x < 0.0 || nonInt(x);
    }

    public static boolean nonInt(double x)
    {
        return Math.abs((x) - Math.round(x)) > 1e-7 * Math.max(1.0, Math.abs(x));
    }

    public static double nonIntCheck(double x, boolean log_p)
    {
        if (nonInt(x))
            System.out.printf("non-integer x = %f", x);
        return D0(log_p);
    }
}
