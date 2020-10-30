import distributions.Beta;
import special_functions.Gamma;

public class Main
{
    public static void main(String[] args)
    {
        System.out.println("nothing happen.");
        label1:

        for (int i = 0; i < 50; i++) {
            System.out.print(" " + i);
            if (i == 25)
                break label1;
        }
        System.out.println();

        double a = 100.1;
        double b = 200;
        label2:
        if (a > 100)
            if (b == 200) {
                System.out.println("break to label2");
                break label2;
            }


        System.out.println("End for loop");

    }
}
