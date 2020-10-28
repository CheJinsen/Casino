import distributions.Beta;
import special_functions.Gamma;

public class Main
{
    public static void main(String[] args) {
        System.out.println(Beta.cdf(0.975, 10, 2, true, false));
        System.out.println(Beta.cdf(0.005, 79, 50, true, false));
        System.out.println(Gamma.gamma(0.710171));
        System.out.println(Gamma.gamma(17.123));

        System.out.println(0xffffffffL);
        System.out.println(0x7fffffffL);

        System.out.println("This foo user: ");
        System.out.println("Well done...");
    }
}
