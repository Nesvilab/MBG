package com.mbg;

import java.util.ArrayList;
import java.util.Map;
import java.util.Objects;

import static java.lang.System.out;

public class MBG {

    private static final String name = "Glycoform Inference";
    private static final String version = "0.2.1";

    public static void main(String[] args) {
        long time = System.nanoTime();
        if (args.length == 0) {
            out.println("No parameters passed. Usage:");
            printHelp();
        }
        try {
            /* Args:
             *      --psm: psm file path (string)
             *      --glycandb: glycan database path (string, optional)
             *      --residuedb: residue definitions path (string, optional)
             *      --maxq: maximum glycan q-value to accept as an input glycoPSM (float, optional)
             *      --minpsms: minimum number of glycoPSMs for a given glycan to be used for inference
             */
            if (Objects.equals(args[0], "--match")) {
                out.printf("Running glycoform inference v%s\n", version);
                Params params = new Params(args);

                PSM_File psmFile = new PSM_File(params);
                ArrayList<PSM> newPSMs = new ArrayList<>();
                for (Map.Entry<String, ArrayList<PSM>> rawfileEntry: psmFile.rawfilePSMmap.entrySet()) {
                    // infer glycoforms within each raw file
                    GlycanMatching matcher = new GlycanMatching(params, psmFile, rawfileEntry.getValue());
                    matcher.groupGlycosites();
                    newPSMs.addAll(matcher.fillInGlycoSites(rawfileEntry.getKey()));
                }
                psmFile.printNewPSMs(newPSMs);

            } else {
                out.printf("%s %s\n", name, version);
                out.print("Usage:\n");
                out.printf("\tTo perform glycoform inference:\n" +
                        "\t\tjava -jar mbg-%s.jar --match [args]\n", version);
                out.println();
                printHelp();
                System.exit(0);
            }
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        out.printf("Done in %.1f s.%n", (System.nanoTime() - time) * 1e-9f);
    }

    private static void printHelp() {
        out.print("Args:\n" +
                " --match: main mode to perform glycoform inference\n" +
                "      --psm: psm file path (string)\n" +
                "      --residuedb: residue definitions path (string, optional)\n" +
                "      --maxq: maximum glycan q-value to accept as an input glycoPSM (float, optional)\n" +
                "      --minpsms: minimum number of glycoPSMs for a given glycan to be used for inference\n" +
                " --clean: alternate mode that removes inferred glycans from the PSM table if no intensity is found\n" +
                "      --psm: psm table to clean (string)\n");
    }

}
