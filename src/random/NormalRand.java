// Random variates from the STANDARD normal distribution  N(0,1).

package random;

import statistics.Stats;
import distributions.Normal;

public class NormalRand
{
    public static double rand()
    {
        final double BIG = 134217728;
        double u1 = UniformRand.rand();
        u1 = (int)(BIG * u1) + UniformRand.rand();
        return Normal.quantile(u1 / BIG, 0.0, 1.0, true, false);
    }

    public static void main(String[] args)
    {
        final int len = 1000;
        double[] arr = new double[len];

        for (int i = 0; i < len; i++) {
            arr[i] = rand();
        }
        System.out.println("mean: " + Stats.mean(arr));
        System.out.println("standard deviation: " + Stats.stdDev(arr));
    }
}
