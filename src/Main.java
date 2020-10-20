import randist.Beta;

public class Main {

    public static void main(String[] args) {
        System.out.println(Beta.cdf(0.975, 10, 2, true, false));
        System.out.println(Beta.cdf(0.005, 79, 50, true, false));
        System.out.println("Well done...");
    }
}
