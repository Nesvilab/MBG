import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;

public class Params {

    public ArrayList<Glycan> glycanDatabase;
    private ArrayList<GlycanResidue> glycanResidues;
    public double maxGlycanQval;
    private static final String defaultResiduePath = "glycan_residues.txt";
    private static final String defaultGlycanPath = "glycan_database.txt";

    public Params(String glycanDBpath, String glycanResiduePath, double maxGlycanQval) throws Exception {
        glycanResidues = parseGlycoResiduesDB(glycanResiduePath);
        glycanDatabase = parseGlycanDB(glycanDBpath);
        this.maxGlycanQval = maxGlycanQval;
    }

    private ArrayList<Glycan> parseGlycanDB(String glycanDBpath) {
        ArrayList<Glycan> glycanDB = new ArrayList<>();
        try {
            BufferedReader in;
            if (glycanDBpath.equals("")) {
                // no glycan database provided - fall back to default glycan list
                in = new BufferedReader(new InputStreamReader(Params.class.getResourceAsStream(defaultGlycanPath)));
            } else {
                Path path = null;
                try {
                    path = Paths.get(glycanDBpath.replaceAll("['\"]", ""));
                } catch (Exception e) {
                    e.printStackTrace();
                    System.out.printf("Malformed glycan path string: [%s]%n", glycanDBpath);
                }
                if (path == null || !Files.exists(path)) {
                    System.out.printf("Glycan database file does not exist: [%s]%n", glycanDBpath);
                }
                in = new BufferedReader(new FileReader(path.toFile()));
            }

            HashMap<String, Boolean> glycansInDB = new HashMap<>();
            String readline;
            while ((readline = in.readLine()) != null) {
                Glycan glycan = new Glycan(readline, this);
                if (glycan.composition.size() == 0) {
                    continue;
                }
                String compStr = glycan.toCompString();
                // prevent duplicated glycans
                if (!glycansInDB.containsKey(compStr)) {
                    glycansInDB.put(compStr, true);
                    glycanDB.add(glycan);
                }
            }
        } catch (FileNotFoundException e) {
            System.out.printf("Glycan database file not found: [%s]%n", glycanDBpath);
        } catch (IOException e) {
            e.printStackTrace();
            System.out.println("IO Exception while reading database file");
        }
        return glycanDB;
    }

    private ArrayList<GlycanResidue> parseGlycoResiduesDB(String glycanResiduesPath) {
        ArrayList<GlycanResidue> residues = new ArrayList<>();
        BufferedReader in;
        try {
            // no glycan database provided - fall back to default glycan list in PeakAnnotator
            if (glycanResiduesPath.matches("")) {
                in = new BufferedReader(new InputStreamReader(Params.class.getResourceAsStream(defaultResiduePath)));
            } else {
                in = new BufferedReader(new FileReader(glycanResiduesPath));
            }
            String line;
            while ((line = in.readLine()) != null) {
                if (line.startsWith("#")) {
                    continue;
                }
                GlycanResidue residue = GlycanResidue.parseResidue(line);
                residues.add(residue);
            }
        } catch (IOException ex) {
            ex.printStackTrace();
            System.out.printf("Error parsing input monosaccharide table %s%n", glycanResiduesPath);
        }
        return residues;
    }

    /**
     * Search all GlycanResidue names and alternate names to try to match this input to one of them.
     * @param inputName string
     * @return
     */
    public GlycanResidue findResidueName(String inputName) {
        for (GlycanResidue residue: glycanResidues) {
            if (residue.name.equalsIgnoreCase(inputName)) {
                return residue;
            }
        }
        // no exact match: check alternate names
        for (GlycanResidue residue: glycanResidues) {
            for (String altName : residue.alternateNames) {
                if (altName.equalsIgnoreCase(inputName)) {
                    return residue;
                }
            }
        }
        // no match - return null and warn user
        return null;
    }


}
