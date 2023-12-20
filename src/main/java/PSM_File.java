import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

public class PSM_File {

    public String[] headers;
    public ArrayList<PSM> psms;
    public File file;
    public HashMap<PSMattr, Integer> PSMcolumns;
    private static final HashMap<PSMattr, String> columnNames;

    static {
        columnNames = new HashMap<>();
        columnNames.put(PSMattr.RT, "Retention");
        columnNames.put(PSMattr.peptide, "Peptide");
        columnNames.put(PSMattr.protein, "Protein");
        columnNames.put(PSMattr.protein_start, "Protein Start");
        columnNames.put(PSMattr.assigned_mods, "Assigned Modifications");
        columnNames.put(PSMattr.glycan, "Total Glycan Composition");
    }

    public PSM_File(String filename, Params params) throws Exception {
        BufferedReader in = new BufferedReader(new FileReader(filename), 1 << 22);
        this.file = new File(filename);
        this.headers = in.readLine().split("\t");
        this.PSMcolumns = initPSMcolumns();

        this.psms = new ArrayList<>();
        String line;
        while((line = in.readLine()) != null) {
            //if(cline.trim().length() > 0)
            if (line.length() > 0) {
                PSM psm = new PSM(line, PSMcolumns, params);
                // ignore non-glyco PSMs
                if (psm.glycanComp != null) {
                    this.psms.add(psm);
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

    public enum PSMattr {
        RT,
        peptide,
        glycan,
        assigned_mods,
        protein,
        protein_start
    }


}
