import java.util.ArrayList;
import java.util.HashMap;

public class Main {

    public static void main(String[] args) {
        long time = System.nanoTime();
        try {
            /* Args:
             * - psm file path (string)
             * - glycan database path (string, optional)
             * - glycan residue definitions path (string, optional)
             */
            String glycanResPath = args.length >= 3 ? args[2].trim() : "";
            String glycanDBpath = args.length >= 2 ? args[1].trim() : "";

            Params params = new Params(glycanDBpath, glycanResPath, 0.01);
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
