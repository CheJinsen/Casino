package specfun.detail;

import specfun.Lgamma;

public class Catherine {

    public static double bd0(double x, double np)
    {
        double ej, s, s1, v;
        int j;

        if(Double.isInfinite(x) || Double.isInfinite(np) || np == 0.0) {
            System.out.println("Argument out of domain in bd0.");
            return Double.NaN;
        }

        if (Math.abs(x-np) < 0.1*(x+np)) {
            v = (x-np)/(x+np);  // might underflow to 0
            s = (x-np)*v;/* s using v -- change by MM */
            if(Math.abs(s) < Double.MIN_VALUE) return s;
            ej = 2*x*v;
            v = v*v;
            for (j = 1; j < 1000; j++) { /* Taylor series; 1000: no infinite loop
						as |v| < .1,  v^2000 is "zero" */
                ej *= v;// = v^(2j+1)
                s1 = s+ej/((j<<1)+1);
                if (s1 == s) /* last term was effectively 0 */
                    return s1 ;
                s = s1;
            }
        }
        /* else:  | x - np |  is not too small */
        return(x*Math.log(x/np)+np-x);
    }

    public static double stirlerr(double n)
    {
        final double S0 = 0.083333333333333333333;       /* 1/12 */
        final double S1 = 0.00277777777777777777778;     /* 1/360 */
        final double S2 = 0.00079365079365079365079365;  /* 1/1260 */
        final double S3 = 0.000595238095238095238095238; /* 1/1680 */
        final double S4 = 0.0008417508417508417508417508;/* 1/1188 */

/*
  error for 0, 0.5, 1.0, 1.5, ..., 14.5, 15.0.
*/
        double[] sferr_halves = {
                0.0, /* n=0 - wrong, place holder only */
                0.1534264097200273452913848,  /* 0.5 */
                0.0810614667953272582196702,  /* 1.0 */
                0.0548141210519176538961390,  /* 1.5 */
                0.0413406959554092940938221,  /* 2.0 */
                0.03316287351993628748511048, /* 2.5 */
                0.02767792568499833914878929, /* 3.0 */
                0.02374616365629749597132920, /* 3.5 */
                0.02079067210376509311152277, /* 4.0 */
                0.01848845053267318523077934, /* 4.5 */
                0.01664469118982119216319487, /* 5.0 */
                0.01513497322191737887351255, /* 5.5 */
                0.01387612882307074799874573, /* 6.0 */
                0.01281046524292022692424986, /* 6.5 */
                0.01189670994589177009505572, /* 7.0 */
                0.01110455975820691732662991, /* 7.5 */
                0.010411265261972096497478567, /* 8.0 */
                0.009799416126158803298389475, /* 8.5 */
                0.009255462182712732917728637, /* 9.0 */
                0.008768700134139385462952823, /* 9.5 */
                0.008330563433362871256469318, /* 10.0 */
                0.007934114564314020547248100, /* 10.5 */
                0.007573675487951840794972024, /* 11.0 */
                0.007244554301320383179543912, /* 11.5 */
                0.006942840107209529865664152, /* 12.0 */
                0.006665247032707682442354394, /* 12.5 */
                0.006408994188004207068439631, /* 13.0 */
                0.006171712263039457647532867, /* 13.5 */
                0.005951370112758847735624416, /* 14.0 */
                0.005746216513010115682023589, /* 14.5 */
                0.005554733551962801371038690  /* 15.0 */
        };
        double nn;
        final double M_LN_SQRT_2PI = 0.918938533204672741780329736406;	// log(sqrt(2*pi))
        if (n <= 15.0) {
            nn = n + n;
            if (nn == (int)nn) return(sferr_halves[(int)nn]);
            return(Lgamma.lammafn(n + 1.) - (n + 0.5)*Math.log(n) + n - M_LN_SQRT_2PI);
        }

        nn = n*n;
        if (n>500) return((S0-S1/nn)/n);
        if (n> 80) return((S0-(S1-S2/nn)/nn)/n);
        if (n> 35) return((S0-(S1-(S2-S3/nn)/nn)/nn)/n);
        /* 15 < n <= 35 : */
        return((S0-(S1-(S2-(S3-S4/nn)/nn)/nn)/nn)/n);
    }
}
