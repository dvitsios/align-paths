import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.Writer;
import java.net.MalformedURLException;
import java.net.URL;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Arrays;
import java.util.Date;
import java.util.StringTokenizer;

class GlobalInitializer{

    public static String baseDirectory;
    public static String genDataDirStr;
    public static String fastaDirStr;
    public static String xmlDirStr;
    public static String blastDbDirStr;
    public static String colMapsDirStr;
    public static String imgsDirForCurRun;
    
    

    public GlobalInitializer() throws MalformedURLException, IOException{

        try{

            baseDirectory = new File(".").getAbsolutePath();
            baseDirectory = baseDirectory.substring(0, baseDirectory.length() - 6);
            FileInputStream fstream = new FileInputStream(baseDirectory+"/input.txt");
            DataInputStream in = new DataInputStream(fstream);
            BufferedReader br = new BufferedReader(new InputStreamReader(in));

            System.out.println("Input parameters");

            
            //read pathway id
            br.readLine();
            MainClass.pathId = br.readLine();
            StringTokenizer stokenz = new StringTokenizer(MainClass.pathId);
            MainClass.pathId = stokenz.nextToken();

            
            System.out.println("> Pathway id: "+MainClass.pathId);

            //check if MCL is enabled
            br.readLine();
            MainClass.mclFlag = Boolean.parseBoolean(br.readLine());
            System.out.println("> MCL: "+MainClass.mclFlag);

            //read inflation parameter value, assign only one value
            br.readLine();
            MainClass.inflationParams[0] = Double.parseDouble(br.readLine());
            System.out.println("> Inflation parameter: "+MainClass.inflationParams[0]);

            //check if EM is enabled
            br.readLine();
            MainClass.emFlag = Boolean.parseBoolean(br.readLine());
            System.out.println("> EM: "+MainClass.emFlag);

            //read eValue value
            br.readLine();
            String eValStr = br.readLine();
            StringTokenizer st = new StringTokenizer(eValStr, ",");
            String eValBase = st.nextToken();
            String eValExp = st.nextToken();
            
            BlastXMLParser.expValueThreshold = Math.pow(Double.parseDouble(eValBase), Double.parseDouble(eValExp));
            System.out.println("> e-Value: "+BlastXMLParser.expValueThreshold.toString());

            //read number of genomes
            br.readLine();
            int numOfGenomes = Integer.parseInt(br.readLine());
            System.out.println("> Number of genomes: "+numOfGenomes);

            //read genomes names
            MainClass.orgsIds = new String[numOfGenomes];

            br.readLine();
            for(int i=0; i<numOfGenomes; i++){
                MainClass.orgsIds[i] = br.readLine();
                stokenz = new StringTokenizer(MainClass.orgsIds[i]);
                MainClass.orgsIds[i] = stokenz.nextToken();
                
                System.out.println(MainClass.orgsIds[i]);
            }

            in.close();

            Arrays.sort(MainClass.orgsIds);

        } catch (Exception e) {
                System.err.println("Error: " + e.getMessage());
        }

        //create directories
        genDataDirStr = baseDirectory+"/gen_data";
        boolean success = (new File(genDataDirStr)).mkdir();
        if (success) {
            System.out.println("Directory: " + genDataDirStr + " was created succesfully.");
        }
        fastaDirStr = genDataDirStr+"/FASTA_Files";
        success = (new File(fastaDirStr)).mkdir();
        if (success) {
            System.out.println("Directory: " + fastaDirStr + " was created succesfully.");
        }
        blastDbDirStr = genDataDirStr+"/Custom_BLAST_DBs";
        success = (new File(blastDbDirStr)).mkdir();
        if (success) {
            System.out.println("Directory: " + blastDbDirStr + " was created succesfully.");
        }
        xmlDirStr = genDataDirStr+"/Blast_Output_XML";
        success = (new File(xmlDirStr)).mkdir();
        if (success) {
            System.out.println("Directory: " + xmlDirStr + " was created succesfully.");
        }
        colMapsDirStr = baseDirectory+"/output";
        success = (new File(colMapsDirStr)).mkdir();
        if (success) {
            System.out.println("Directory: " + colMapsDirStr + " was created succesfully.");
        }
        String tmpGenomeDir = colMapsDirStr+"/"+(Integer.toString(MainClass.orgsIds.length))+"-genomes";
        success = (new File(tmpGenomeDir)).mkdir();
        if (success) {
            System.out.println("Directory: " + tmpGenomeDir + " was created succesfully.");
        }
        imgsDirForCurRun = tmpGenomeDir+"/";


        String imgsDirForCurRun_appendStr = "";
        DateFormat dateFormat = new SimpleDateFormat("dd/MM/yyyyy HH:mm:ss");
        Date date = new Date();
        imgsDirForCurRun_appendStr += dateFormat.format(date);

        imgsDirForCurRun_appendStr = imgsDirForCurRun_appendStr.replace("/", "-");
        imgsDirForCurRun_appendStr = imgsDirForCurRun_appendStr.replace(" ", "_");
        imgsDirForCurRun_appendStr = imgsDirForCurRun_appendStr.replaceFirst(":", "h");
        imgsDirForCurRun_appendStr = imgsDirForCurRun_appendStr.replaceFirst(":", "m");
        imgsDirForCurRun_appendStr += "s";
        imgsDirForCurRun += imgsDirForCurRun_appendStr;

        success = (new File(imgsDirForCurRun)).mkdir();
        if (success) {
            System.out.println("Directory: " + imgsDirForCurRun + " was created succesfully.\n");
        }

        //create a .txt file inside the 'imgsDirForCurRun' folder,
        //containing the genomes of the current test case.
        String genomesStringIdentifier = "";
        String[] genomesSortedList = new String[MainClass.orgsIds.length];
        for(int copyOrgsIds = 0; copyOrgsIds<MainClass.orgsIds.length; copyOrgsIds++)
            genomesSortedList[copyOrgsIds] = MainClass.orgsIds[copyOrgsIds];
        Arrays.sort(genomesSortedList);

        for(int m=0; m<MainClass.orgsIds.length; m++){
            genomesStringIdentifier += genomesSortedList[m]+"\n";
        }


        Writer writer = null;

        File file = new File(imgsDirForCurRun+"/genomes.txt");
        writer = new BufferedWriter(new FileWriter(file));
        writer.write(genomesStringIdentifier);

        if (writer != null)
            writer.close();


        boolean valid = false;
        while(!valid){
            try{
                //download NCBI_BlastOutput.dtd
                BufferedInputStream in1 = new BufferedInputStream(new URL("http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd").openStream());
                FileOutputStream fos1 = new FileOutputStream(genDataDirStr+"/Blast_Output_XML/NCBI_BlastOutput" + ".dtd");
                BufferedOutputStream bout1 = new BufferedOutputStream(fos1, 1024);
                byte data1[] = new byte[1024];
                int x1 = 0;
                while ((x1 = in1.read(data1, 0, 1024)) >= 0) {
                    bout1.write(data1, 0, x1);
                }
                bout1.close();
                in1.close();

                //download NCBI_BlastOutput.mod.dtd
                BufferedInputStream in2 = new BufferedInputStream(new URL("http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.mod.dtd").openStream());
                FileOutputStream fos2 = new FileOutputStream(genDataDirStr+"/Blast_Output_XML/NCBI_BlastOutput" + ".mod" + ".dtd");
                BufferedOutputStream bout2 = new BufferedOutputStream(fos2, 1024);
                byte data2[] = new byte[1024];
                int x2 = 0;
                while ((x2 = in2.read(data2, 0, 1024)) >= 0) {
                    bout2.write(data2, 0, x2);
                }
                bout2.close();
                in2.close();

                //download NCBI_Entity.mod.dtd
                BufferedInputStream in3 = new BufferedInputStream(new URL("http://www.ncbi.nlm.nih.gov/dtd/NCBI_Entity.mod.dtd").openStream());
                FileOutputStream fos3 = new FileOutputStream(genDataDirStr+"/Blast_Output_XML/NCBI_Entity" + ".mod" + ".dtd");
                BufferedOutputStream bout3 = new BufferedOutputStream(fos3, 1024);
                byte data3[] = new byte[1024];
                int x3 = 0;
                while ((x3 = in3.read(data3, 0, 1024)) >= 0) {
                    bout3.write(data3, 0, x3);
                }
                bout3.close();
                in3.close();

                valid = true;
            } catch(Exception e){
                    System.err.println("Error: " + e.getMessage());
            }
        }
    }
}