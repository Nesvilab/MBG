package com.mbg;

import java.io.*;
import java.util.*;

public class PSM_File {

    public String[] headers;
    public HashMap<String, ArrayList<PSM>> rawfilePSMmap;
    public final ArrayList<PSM> allPSMs;
    public File file;
    public HashMap<PSMattr, Integer> PSMcolumns;
    private static final HashMap<PSMattr, String> columnNames;
    public int addCols;
    public static int insertIndex;

    static {
        columnNames = new HashMap<>();
        columnNames.put(PSMattr.RT, "Retention");
        columnNames.put(PSMattr.peptide, "Peptide");
        columnNames.put(PSMattr.protein, "Protein");
        columnNames.put(PSMattr.IM, "Ion Mobility");
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
        columnNames.put(PSMattr.has_ms2_spec, "Has MS2 Spectrum?");
        columnNames.put(PSMattr.has_Y1, "Has Y1 ion?");
        columnNames.put(PSMattr.MBG_score, "MBG Score");
        columnNames.put(PSMattr.parent_glycan, "Parent Glycan");
        columnNames.put(PSMattr.added_residual, "Added residual");
        columnNames.put(PSMattr.spectrum_file, "Spectrum File");
        columnNames.put(PSMattr.intensity, "Intensity");
    }

    public PSM_File(Params params) throws Exception {
        BufferedReader in = new BufferedReader(new FileReader(params.psmPath.toFile()), 1 << 22);
        this.file = new File(String.valueOf(params.psmPath));
        this.headers = in.readLine().split("\t");
        this.PSMcolumns = initPSMcolumns();
        allPSMs = new ArrayList<>();

        // handle inserting new header(s) if needed
        if (PSMcolumns.getOrDefault(PSMattr.inferred_glycoform, -1) == -1) {
            addCols = 6;
            if (PSMcolumns.containsKey(PSMattr.glyco_q_val)) {
                insertIndex = PSMcolumns.get(PSMattr.glyco_q_val)+1;
            } else {
                insertIndex = PSMcolumns.get(PSMattr.glycan)+1;
            }
            ArrayList<String> newHeaders = new ArrayList<>(Arrays.asList(this.headers));
            newHeaders.add(insertIndex, columnNames.get(PSMattr.inferred_glycoform));
            newHeaders.add(insertIndex, columnNames.get(PSMattr.has_ms2_spec));
            newHeaders.add(insertIndex, columnNames.get(PSMattr.has_Y1));
            newHeaders.add(insertIndex, columnNames.get(PSMattr.MBG_score));
            newHeaders.add(insertIndex, columnNames.get(PSMattr.parent_glycan));
            newHeaders.add(insertIndex, columnNames.get(PSMattr.added_residual));
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
                PSM psm = new PSM(line, PSMcolumns, params, false, false, 0.0, "", "", addCols);
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

    public static void massWriter(List<Double> massesToWrite) throws IOException {
        // The file path to create or write to
        String filePath = "MBG_masses_list.txt";
        // Using try-with-resources to ensure the writer is closed automatically
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(filePath))) {
            // Writing masses to the file
            for (Double mass: massesToWrite) {
                writer.write(mass + "\n");
            }
            System.out.println("Glycan mass written successfully!");
        } catch (IOException e) {
            System.out.println("An error occurred: " + e.getMessage());
        }
        appendMassesFileContents("modmasses_ionquant.txt", "MBG_masses_list.txt");
    }

    private static void appendMassesFileContents(String fileAppendTo, String fileAppendFrom) throws IOException {
        File fileFrom = new File(fileAppendTo);
        if (!fileFrom.exists()) {
            return;
        }
        BufferedReader appendFromReader = new BufferedReader(new FileReader(fileAppendFrom));
        BufferedReader appendToReader = new BufferedReader(new FileReader(fileAppendTo));

        BufferedWriter writer = new BufferedWriter(new FileWriter(fileAppendTo, true));

        Set<Long> masses = new HashSet<>();
        // read all masses in the file to append to so we can avoid duplicates while appending
        String readline;
        while ((readline = appendToReader.readLine()) != null) {
            long mass = Math.round(Double.parseDouble(readline.trim()) * 100);
            masses.add(mass);
        }

        // read new masses and append, checking for duplicates
        while ((readline = appendFromReader.readLine()) != null) {
            long mass = Math.round(Double.parseDouble(readline.trim()) * 100);
            if (!masses.contains(mass)) {
                writer.write(readline + "\n");
            }
        }
        writer.close();
    }


    public enum PSMattr {
        RT,
        peptide,
        glycan,
        assigned_mods,
        protein,
        IM,
        protein_start,
        spectrum,
        mod_peptide,
        delta_mass,
        calc_pep_mass,
        calc_pep_mz,
        charge,
        glyco_q_val,
        inferred_glycoform,
        spectrum_file,
        has_ms2_spec,
        has_Y1,
        MBG_score,
        ionMobilityDiff,
        intensity,
        parent_glycan,
        added_residual
    }


}
