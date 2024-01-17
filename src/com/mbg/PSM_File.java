package com.mbg;

import java.io.*;
import java.util.*;

public class PSM_File {

    public String[] headers;
    public HashMap<String, ArrayList<PSM>> rawfilePSMmap;
    private final ArrayList<PSM> allPSMs;
    public File file;
    public HashMap<PSMattr, Integer> PSMcolumns;
    private static final HashMap<PSMattr, String> columnNames;
    public int addCols;

    static {
        columnNames = new HashMap<>();
        columnNames.put(PSMattr.RT, "Retention");
        columnNames.put(PSMattr.peptide, "Peptide");
        columnNames.put(PSMattr.protein, "Protein");
        columnNames.put(PSMattr.protein_start, "Protein Start");
        columnNames.put(PSMattr.assigned_mods, "Assigned Modifications");
        columnNames.put(PSMattr.glycan, "Total Glycan Composition");
        columnNames.put(PSMattr.spectrum, "Spectrum");
        columnNames.put(PSMattr.delta_mass, "Delta Mass");
        columnNames.put(PSMattr.mod_peptide, "Modified Peptide");
        columnNames.put(PSMattr.calc_pep_mass, "Calculated Peptide Mass");
        columnNames.put(PSMattr.calc_pep_mz, "Calculated M/Z");
        columnNames.put(PSMattr.charge, "Charge");
        columnNames.put(PSMattr.glyco_q_val, "Glycan q-value");
        columnNames.put(PSMattr.inferred_glycoform, "Inferred Glycoform?");
        columnNames.put(PSMattr.spectrum_file, "Spectrum File");
    }

    public PSM_File(Params params) throws Exception {
        BufferedReader in = new BufferedReader(new FileReader(params.psmPath.toFile()), 1 << 22);
        this.file = new File(String.valueOf(params.psmPath));
        this.headers = in.readLine().split("\t");
        this.PSMcolumns = initPSMcolumns();
        allPSMs = new ArrayList<>();

        // handle inserting new header(s) if needed
        if (PSMcolumns.getOrDefault(PSMattr.inferred_glycoform, -1) == -1) {
            addCols = 1;
            ArrayList<String> newHeaders = new ArrayList<>(Arrays.asList(this.headers));
            newHeaders.add(PSMcolumns.get(PSMattr.glyco_q_val) + 1, columnNames.get(PSMattr.inferred_glycoform));
            headers = newHeaders.toArray(new String[0]);
            this.PSMcolumns = initPSMcolumns();
        } else {
            addCols = 0;
        }

        // read PSMs into Map by rawfile
        rawfilePSMmap = new HashMap<>();
        String line;
        while((line = in.readLine()) != null) {
            if (line.length() > 0) {
                PSM psm = new PSM(line, PSMcolumns, params, false, addCols);
                allPSMs.add(psm);
                if (rawfilePSMmap.containsKey(psm.rawfile)) {
                    rawfilePSMmap.get(psm.rawfile).add(psm);
                } else {
                    ArrayList<PSM> psms = new ArrayList<>();
                    psms.add(psm);
                    rawfilePSMmap.put(psm.rawfile, psms);
                }
            }
        }
        in.close();
    }

    private HashMap<PSMattr, Integer> initPSMcolumns() {
        HashMap<PSMattr, Integer> psmColumns = new HashMap<>();
        for (Map.Entry<PSMattr, String> nameEntry: columnNames.entrySet()) {
            int index = getColumn(nameEntry.getValue());
            if (index == -1) {
                // todo: error checking
            } else {
                psmColumns.put(nameEntry.getKey(), index);
            }
        }
        return psmColumns;
    }

    public int getColumn(String head) {
        for(int i = 0; i < headers.length; i++)
            if(headers[i].equals(head))
                return i;
        return -1;
    }

    public void printNewPSMs(ArrayList<PSM> newPSMs) throws IOException {
        String tempName = this.file.getAbsolutePath() + ".mbg.tsv";
        PrintWriter out = new PrintWriter(new FileWriter(tempName));

        // write header
        out.println(String.join("\t", this.headers));

        // write original and new PSMs
        newPSMs.addAll(allPSMs);
        Collections.sort(newPSMs);
        for (PSM newPSM: newPSMs) {
            out.println(String.join("\t", newPSM.splits));
        }
        out.close();

        boolean success = this.file.delete();
        File newFileName = new File(tempName);
        success = newFileName.renameTo(this.file);
    }

    public enum PSMattr {
        RT,
        peptide,
        glycan,
        assigned_mods,
        protein,
        protein_start,
        spectrum,
        mod_peptide,
        delta_mass,
        calc_pep_mass,
        calc_pep_mz,
        charge,
        glyco_q_val,
        inferred_glycoform,
        spectrum_file
    }


}
