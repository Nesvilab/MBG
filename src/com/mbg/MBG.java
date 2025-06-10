package com.mbg;

import java.util.ArrayList;
import java.util.Map;
import java.util.Objects;

import static java.lang.System.out;

/**
 * MBG (Match-between-glycans) - A tool for expanding glycopeptide identification
 * by generating and validating dummy glycopeptides with modified glycan structures.
 * 
 * This class serves as the main entry point for the MBG application, handling
 * command-line arguments and orchestrating the glycoform inference process.
 */
public class MBG {

    // Application metadata
    private static final String APPLICATION_NAME = "Glycoform Inference";
    private static final String APPLICATION_VERSION = "0.2.21";

    /**
     * Main entry point for the MBG application.
     * Processes command line arguments and initiates the glycoform inference workflow.
     * 
     * @param args Command line arguments:
     *             --match: Main mode for glycoform inference
     *             --psm: PSM file path
     *             --glycandb: Glycan database path (optional)
     *             --residuedb: Residue definitions path (optional)
     *             --maxq: Maximum glycan q-value threshold (optional)
     *             --minpsms: Minimum number of glycoPSMs for inference (optional)
     */
    public static void main(String[] args) {
        long startTime = System.nanoTime();
        
        if (args.length == 0) {
            out.println("No parameters passed. Usage:");
            printHelp();
        }
        
        try {
            if (Objects.equals(args[0], "--match")) {
                out.printf("Running glycoform inference v%s\n", APPLICATION_VERSION);
                Params params = new Params(args);

                // Initialize PSM file processing
                PSM_File psmFile = new PSM_File(params);
                out.printf("Reading PSM file %s\n", psmFile.file.toString());
                
                // Process PSMs for each raw file
                ArrayList<PSM> inferredPSMs = new ArrayList<>();
                for (Map.Entry<String, ArrayList<PSM>> rawFileEntry: psmFile.rawfilePSMmap.entrySet()) {
                    // Perform glycan matching and inference for each raw file
                    GlycanMatching matcher = new GlycanMatching(params, psmFile, rawFileEntry.getValue());
                    matcher.groupGlycosites();
                    inferredPSMs.addAll(matcher.fillInGlycoSites(rawFileEntry.getKey()));
                }
                
                // Output the results
                psmFile.printNewPSMs(inferredPSMs);

            } else {
                // Display usage information for invalid commands
                out.printf("%s %s\n", APPLICATION_NAME, APPLICATION_VERSION);
                out.print("Usage:\n");
                out.printf("\tTo perform glycoform inference:\n" +
                        "\t\tjava -jar mbg-%s.jar --match [args]\n", APPLICATION_VERSION);
                out.println();
                printHelp();
                System.exit(0);
            }
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        
        // Print execution time
        out.printf("Done in %.1f s.%n", (System.nanoTime() - startTime) * 1e-9f);
    }

    /**
     * Displays help information about command line arguments and usage.
     */
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
