import java.util.ArrayList;
import java.util.HashMap;

public class Main {

    public static void main(String[] args) {
        long time = System.nanoTime();
        try {
            /* Args:
             * - psm file path (string)
             * - glycan database path (string)
             */
            Params params = new Params(args[1].trim(), args[2].trim());
            GlycanMatching matcher = new GlycanMatching(args[0].trim(), params);

            matcher.groupGlycosites();
            ArrayList<PSM> filledInPSMs = matcher.fillInGlycoSites();
            matcher.psmFile.printNewPSMs(filledInPSMs);

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        System.out.printf("Done in %.1f s.%n", (System.nanoTime() - time) * 1e-9f);
    }

}
