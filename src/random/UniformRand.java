package random;

import java.util.Random;

public class UniformRand
{
    public static double rand()
    {
        Random r = new Random();
        MersenneTwister.initGenRand(r.nextLong());
        return MersenneTwister.genRandReal2();
    }
}
