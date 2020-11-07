import java.util.Vector;

public class Main
{
    public static void main(String[] args)
    {
        Vector<Integer> v = new Vector<>();

        v.setSize(10);
        for (int i = 0; i < 10; i++) {
            v.set(i, i * i);
        }
        System.out.println(v.capacity());
        System.out.println(v.get(9));

        Vector<Vector<Double>> vv = new Vector<>();
        vv.setSize(10);
        for (int i = 0; i < 10; i++) {
            Vector<Double> tmp = new Vector<>();
            tmp.setSize(10);
            vv.set(i, tmp);
        }

        for (int i = 0; i < 10; i++) {
            for (int j = 0; j < 10; j++) {
                vv.get(i).set(j, (double)i * j);
            }
        }

        System.out.println(vv.get(9).capacity());
        System.out.println(vv.get(9).get(1));
        System.out.println("Done...");
    }
}
