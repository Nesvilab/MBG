package com.mbg;

import java.io.IOException;
import java.util.ArrayList;

public class TableCleaner {

    private Params params;
    private PSM_File psmFile;

    public TableCleaner(Params params, PSM_File psmFile) {
        this.params = params;
        this.psmFile = psmFile;
    }

    /**
     * Remove all inferred glycoform PSMs with no intensity found
     */
    public void cleanPSMtable() throws IOException {
        ArrayList<PSM> printPSMs = new ArrayList<>();
        for (PSM psm: psmFile.allPSMs) {
            if (!psm.inferredGlycoform) {
                printPSMs.add(psm);
            }
        }
        psmFile.printCleanedPSMs(printPSMs);
    }

}
