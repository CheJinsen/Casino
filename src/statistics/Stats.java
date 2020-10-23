package statistics;

public class Stats
{
    public static double mean(double[] data)
    {
        nullWarn(data);
        double ret = 0.0;
        for (int i = 0; i < data.length; i++) {
            ret += (data[i] - ret) / (i + 1);
        }
        return ret;
    }

    public static double var(double[] data)
    {
        nullWarn(data);
        double ret = 0.0;
        double mean = mean(data);
        for (int i = 0; i < data.length; i++) {
            double delta = (data[i] - mean);
            ret += (delta * delta - ret) / (i + 1);
        }
        return ret * (double)data.length / (data.length - 1.0);
    }

    public static double stdDev(double[] data)
    {
        nullWarn(data);
        return Math.sqrt(var(data));
    }

    private static void nullWarn(Object data)
    {
        if (data == null) {
            throw new IllegalArgumentException("Argument is null.");
        }
    }

    public static void main(String[] args)
    {
        double[] arr = {1.23, 6, 7.98, 98.84, 34.86, 90.87, 34.567};
        System.out.println(mean(arr));
        System.out.println(var(arr));
        System.out.println(stdDev(arr));
    }
}
