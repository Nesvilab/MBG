import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;

public class PSM_File {

    public String[] headers;
    public ArrayList<String> psms;
    public File file;

    public PSM_File(String filename) throws Exception {
        this(new File(filename));
    }

    public PSM_File(File f) throws Exception {
        BufferedReader in = new BufferedReader(new FileReader(f), 1 << 22);
        this.file = f;
        this.headers = in.readLine().split("\t");

        this.psms = new ArrayList<>();
        String line;
        while((line = in.readLine()) != null) {
            //if(cline.trim().length() > 0)
            if (line.length() > 0)
                this.psms.add(line);
        }
        in.close();
    }
}
