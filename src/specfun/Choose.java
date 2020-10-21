package specfun;

import randist.detail.Dpq;
import specfun.detail.Sign;
import java.math.BigDecimal;
import java.math.MathContext;

public class Choose {
    public static double choose(double n, double k)
    {
        double r, k0 = k;
        k = Math.round(k);
        final int k_small_max = 30;

        if (Double.isNaN(n) || Double.isNaN(k)) {
            return n + k;
        }

        if (Math.abs(k - k0) > 1e-7) {
            System.out.printf("'k' (%.2f) must be integer, rounded to %.0f\n", k0, k);
        }

        if (k < k_small_max) {
            if (n - k < k && n >= 0.0 && isInt(n)) {
                k = Math.round(n - k);
            }
            if (k < 0.0) {
                return 0.0;
            }
            if (k == 0.0) {
                return 1.0;
            }
            r = n;
            for (int j = 2; j <= k; j++) {
                r *= (n - j + 1.0) / j;
            }

            if (r >= Long.MAX_VALUE) {
                BigDecimal temp = BigDecimal.valueOf(r);
                return isInt(n) ? temp.round(MathContext.DECIMAL128).doubleValue() : r;
            }

            return isInt(n) ? Math.round(r) : r;
        }

        if (n < 0) {
            r = choose(-n + k - 1, k);
            if (isOdd(k)) {
                r = -r;
            }
            return r;
        } else if (isInt(n)) {
            n = Math.round(n);
            if (n < k)
                return 0.0;
            if (n - k < k_small_max)
                return choose(n, n - k);

            r = Math.exp(logFastChoose(n, k));
            if (r >= Long.MAX_VALUE) {
                BigDecimal temp = BigDecimal.valueOf(r);
                return temp.round(MathContext.DECIMAL128).doubleValue();
            }
            return Math.round(Math.exp(logFastChoose(n, k)));
        }

        if (n < k - 1.0) {
            Sign s_choose = new Sign(0);
            r = logFastChoose2(n, k, s_choose);
            return s_choose.val * Math.exp(r);
        }
        return Math.exp(logFastChoose(n, k));
    }

    public static double logFastChoose(double n, double k)
    {
        return -Math.log(n + 1.0) - Beta.logBeta(n - k + 1.0, k + 1.0);
    }

    private static boolean isOdd(double k)
    {
        return k != 2.0 * Math.floor(k / 2.0);
    }

    private static boolean isInt(double x)
    {
        return !Dpq.nonInt(x);
    }

    private static double logFastChoose2(double n, double k, Sign s_choose)
    {
        double r = Lgamma.lgammafn_sign(n - k + 1.0, s_choose);
        return Lgamma.lammafn(n + 1.0) - Lgamma.lammafn(k + 1.0) - r;
    }

    public static double logChoose(double n, double k)
    {
        double k0 = k;
        k = Math.round(k);

        if (Double.isNaN(n) || Double.isNaN(k)) {
            return n + k;
        }

        if (Math.abs(k - k0) > 1e-7) {
            System.out.printf("'k' (%.2f) must be integer, rounded to %.0f\n", k0, k);
        }
        if (k < 2.0) {
            if (k < 0.0)
                return Double.NEGATIVE_INFINITY;
            if (k == 0.0)
                return 0.0;
            return Math.log(Math.abs(n));
        }
        if (n < 0.0) {
            logChoose(-n + k - 1, k);
        } else if (isInt(n)) {
            n = Math.round(n);
            if (n < k)
                return Double.NEGATIVE_INFINITY;
            if (n - k < 2)
                return logChoose(n, n - k);
            return logFastChoose(n, k);
        }

        if (n < k - 1) {
            Sign s = new Sign(0);
            return logFastChoose2(n, k, s);
        }
        return logFastChoose(n, k);
    }

    public static void main(String[] args)
    {
        System.out.println(choose(50.321, 50));
        System.out.println(logChoose(0.8956, 50));
    }
}
