
public class Main {

    public static void main(String[] args) {
        long time = System.nanoTime();
        try {
            /* Args:
             * - psm file path (string)
             * - glycan database path (string)
             */
            Params analyzer = new Params(args[0].trim(), args[1].trim(), args[2].trim());
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        System.out.printf("Done in %.1f s.%n", (System.nanoTime() - time) * 1e-9f);
    }

}
