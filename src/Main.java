import randist.Beta;
import specfun.Gamma;

public class Main {

    public static void main(String[] args) {
        System.out.println(Beta.cdf(0.975, 10, 2, true, false));
        System.out.println(Beta.cdf(0.005, 79, 50, true, false));
        System.out.println(Gamma.gammafn(171.0171));


        System.out.println("Well done...");
    }
}
