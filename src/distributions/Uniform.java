package distributions;

import java.util.Random;
import random.MersenneTwister;

public class Uniform {

    public static double rand()
    {
        Random r = new Random();
        MersenneTwister.initGenRand(r.nextLong());
        return MersenneTwister.genRandReal2();
    }

    public static void main(String[] args)
    {
        int length = 1000;
        double sum = 0.0;
        for (int i = 0; i < length; i++) {
            double temp = rand();
            System.out.println(temp);
            sum += temp;
        }
        System.out.println("Avg: " + sum / length);
    }
}
