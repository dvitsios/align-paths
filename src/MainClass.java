/*
 *  Copyright (C) 2012 by Vitsios Dimitrios
 *
 *  Permission is hereby granted, free of charge, to any person obtaining a copy
 *  of this software and associated documentation files (the "Software"), to deal
 *  in the Software without restriction, including without limitation the rights
 *  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 *  copies of the Software, and to permit persons to whom the Software is
 *  furnished to do so, subject to the following conditions:
 *
 *  The above copyright notice and this permission notice shall be included in
 *  all copies or substantial portions of the Software.
 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 *  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 *  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 *  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 *  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 *  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 *  THE SOFTWARE.
 */

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.StringTokenizer;
import com.mathworks.toolbox.javabuilder.*;
import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.InputStreamReader;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import javax.xml.rpc.ServiceException;
import keggapi.KEGGLocator;
import keggapi.KEGGPortType;


class MainClass {


    public static String[] orgsIds;

  

    public static String pathId;

    public static double[] inflationParams = new double[1]; 

   
    public static FileWriter logFileStream;
    public static String outputLogStrMcl;
    public static String outputLogStrMclAppend;

    public static String outputLogStrEm;
    public static String outputLogStrEmAppend;
    
    public static int totalGenesNumber;
    public static volatile int createFastaThreadCnt;
    public static volatile int colorPathsThreadCnt;
    public static volatile int algoThreadsCnt;
    public static volatile int ecThreadCnt;
    public static volatile int enzThreadCnt;
    public static volatile int xmlParserThreadCnt;
    public static volatile int customDbThreadCnt;
    public static volatile int blastpThreadCnt;
    public static volatile int nestedGetEcsThreadsCnt;
    public static volatile int nestedGetColorMapsThreadCnt;
    public static volatile int nestedGetColorPathsThreadCntMcl;
    public static volatile int nestedGetColorPathsThreadCntEm;
    public static volatile int getBlackElsThreadMclCnt;

    public static boolean mclFlag = true;
    public static boolean emFlag = false;

    public static volatile List<String> blackECsListMcl = new ArrayList<String>();
    public static volatile List<Integer> blackClustersListMcl = new ArrayList<Integer>();
    public static volatile List<Integer> blackCurClustersListMcl = new ArrayList<Integer>();
    public static volatile List<String> blackCurGeneList = new ArrayList<String>();
    public static volatile List<String> blackExamGeneList = new ArrayList<String>();


    public static List<String> blackECsListEm = new ArrayList<String>();
    public static List<Integer> blackClustersListEm = new ArrayList<Integer>();
    public static List<Integer> blackCurClustersListEm = new ArrayList<Integer>();
    public static List<String> blackCurGeneListEm = new ArrayList<String>();
    public static List<String> blackExamGeneListEm = new ArrayList<String>();


    public static String clustersFromblackECsStr;
    public static String clustersFromblackECsStrEm;

    public static long elapsedTimeMillis;
    public static long startTime;


    // pr[orgsNum][GenesNumberInEachOrganism]: each line holds the genes strings for the corresponding organism (and for the specified path).
    public static volatile String[][] pr;
    // C[totalGenesNumber][totalGenesNumber]: contains the evalues calculated from blastp for every pair of genes.
    public static volatile double[][] C;
    // P[totalGenesNumber][orgsNum]: a value "1" in index (i,j) means that there is an homologous gene of gene "i" in organism "j". If there is not, the value is "0".
    public static volatile double[][] P;

    // ec[totalGenesNumber][]: contains the ec numbers corresponding to each of the genes involved in the examined pathway.
    public static volatile String[][][] ec;
    //enz[orgsNum][EnzymesNumberInEachOrganism]: contains the enzymes involved in each genome's specific pathway
    public static volatile String[][] enz;

    //contains the ec numbers corresponding to each gene
    public static String[][][] ecNumsList;

    //contains the EC numbers appeared in the examined pathway in alpharithmetical order.
    public static String[] ecsFinal;

    //contains the genes that correspond to each EC number appeared in the examined pathway.
    public static String[][] ec2GenesList;

    public static List<List<String>> clustersListMcl;
    public static List<List<String>> clustersListEm;

    public MainClass() throws IOException, MWException {

        

        clustersListMcl = new ArrayList<List<String>>();
        clustersListEm = new ArrayList<List<String>>();

        
    }



    public static void main(String[] args) throws Exception
    {

        
        startTime = System.currentTimeMillis();


        createFastaThreadCnt = 0;
        colorPathsThreadCnt = 0;
        algoThreadsCnt = 0;
        ecThreadCnt = 0;
        enzThreadCnt = 0;
        xmlParserThreadCnt = 0;
        customDbThreadCnt = 0;
        blastpThreadCnt = 0;
        nestedGetEcsThreadsCnt = 0;
        nestedGetColorMapsThreadCnt = 0;
        nestedGetColorPathsThreadCntMcl = 0;
        nestedGetColorPathsThreadCntEm = 0;
        getBlackElsThreadMclCnt = 0;
        outputLogStrMcl = "Log for MCL Algorithm.\n\n";
        outputLogStrEm = "Log for EM Algorithm.\n\n";
        outputLogStrMclAppend = "";
        outputLogStrEmAppend = "";
        clustersFromblackECsStr = "";
        clustersFromblackECsStrEm = "";
        int orgsNum;

        GlobalInitializer globInit = new GlobalInitializer();
        MainClass mainObject = new MainClass();

        orgsNum = MainClass.orgsIds.length;
        
        KeggClientInitializer keggClientObj = new KeggClientInitializer();


        //initialize pr[][]
        pr = new String[orgsIds.length][];
        System.out.println("--- Number of genes for pathway map: " + pathId + "  ---");
        {

            int parallelThreadsForGenesRetrieval = 10;

            if(MainClass.orgsIds.length>300)
                parallelThreadsForGenesRetrieval = 5;
            
            int iterations = MainClass.orgsIds.length/parallelThreadsForGenesRetrieval;
            int rem = MainClass.orgsIds.length%parallelThreadsForGenesRetrieval;

            for(int iter=0; iter<iterations; iter++){

                for(int i=0; i<parallelThreadsForGenesRetrieval; i++){

                    int idx = iter*parallelThreadsForGenesRetrieval + i;

                    Runnable r = new getGenesByGenomeThread(orgsIds[idx], idx);
                    Thread thr = new Thread(r);
                    thr.start();
                }

                while(true){
                    if(MainClass.createFastaThreadCnt > 0){
                        continue;
                    }
                    else
                        break;
                }

            }


            for(int i=0; i<rem; i++){
                
                int idx = iterations*parallelThreadsForGenesRetrieval + i;

                Runnable r = new getGenesByGenomeThread(orgsIds[idx], idx);
                Thread thr = new Thread(r);
                thr.start();
            }

            while(true){
                if(MainClass.createFastaThreadCnt > 0){
                    continue;
                }
                else
                    break;
            }


            
        }


        System.out.println("FASTA files were created successfully!");

        totalGenesNumber = 0;
        for (int i = 0; i < orgsIds.length; i++) {
            System.out.println(i);
            totalGenesNumber += pr[i].length;
        }
        System.out.println("Total number of genes: "+totalGenesNumber);

        //initialize C[][]
        C = new double[totalGenesNumber][totalGenesNumber];
        
        for (int i = 0; i < totalGenesNumber; i++) {
            for (int j = 0; j < totalGenesNumber; j++) {
                C[i][j] = 1000;
            }
        }

        //initialize P[][]
        P = new double[totalGenesNumber][orgsIds.length];
        for (int i = 0; i < totalGenesNumber; i++) {
            for (int j = 0; j < orgsIds.length; j++) {
                P[i][j] = 0;
            }
        }
        
        for (int i = 0; i<orgsIds.length; i++) {
            int offset = 0;
            for (int j = 0; j < i; j++) {
                offset += MainClass.pr[j].length;
            }
            for (int k = 0; k < MainClass.pr[i].length; k++) {
                P[k + offset][i] = 1;
            }
        }


        //initialize ec[][][]
        ec = new String[orgsIds.length][][];
        for(int i=0; i<orgsIds.length; i++)
            ec[i] = new String[pr[i].length][];

        //initialize enz[][]
        enz = new String[orgsIds.length][];


        //get enzymes by genome for the examined pathway
        {
            int parallelThreadsForECsRetrieval = 10;

            int iterations = MainClass.orgsIds.length/parallelThreadsForECsRetrieval;
            int rem = MainClass.orgsIds.length%parallelThreadsForECsRetrieval;

            for(int iter=0; iter<iterations; iter++){

                for(int i=0; i<parallelThreadsForECsRetrieval; i++){

                    int idx = iter*parallelThreadsForECsRetrieval + i;

                    Runnable ecR = new getECsByGenomeThread(idx);
                    Thread ecThr = new Thread(ecR);
                    ecThr.start();
                }

                while(true){
                    if(MainClass.ecThreadCnt > 0){
                        continue;
                    }
                    else
                        break;
                }

            }

            System.out.println("retrieve remaining ECs");

            for(int i=0; i<rem; i++){
                System.out.println("in rem ecs");
                int idx = iterations*parallelThreadsForECsRetrieval + i;

                Runnable ecR = new getECsByGenomeThread(idx);
                Thread ecThr = new Thread(ecR);
                ecThr.start();
            }
        

            while(true){
                if(MainClass.ecThreadCnt > 0){
                    continue;
                }
                else
                    break;
            }
        }
	
        System.out.println("********************************************\n********************************************\n");


        //get enzymes by pathway for each genome
        for(int idx=0; idx<MainClass.orgsIds.length; idx++){
            
            Runnable enzR = new getEnzsByPathwayThread(idx);
            Thread enzThr = new Thread(enzR);
            enzThr.start();
            Thread.sleep(1000);
           
        }

	
	while(true){
            if(MainClass.enzThreadCnt > 0){
                Thread.sleep(500);
                continue;
            }
            else
                break;
        }

        
        System.out.println("EC nums retrieved!");


        
        for (int i = 0; i < orgsIds.length; i++) {
            int baseIdx = i;
            String baseOrg = orgsIds[i];
            BlastDbWrapper custBlastDbObject = new BlastDbWrapper();
	    System.out.println("customDb: "+i);
            custBlastDbObject.createDb(baseOrg);
            for (int j = 0; j < orgsIds.length; j++) {
                System.out.println("balstp: "+i+", "+j);
                
                //include comparisons of type [queryOrg]2[baseOrg] where queryOrg=baseOrg.
                int queryIdx = j;
                String queryOrg = orgsIds[j];
                BlastpWrapper blastpObject = new BlastpWrapper();
                blastpObject.alignSequences(baseOrg, queryOrg);
                String xmlFileStr = GlobalInitializer.genDataDirStr+"/Blast_Output_XML/" + queryOrg + "2" + baseOrg + ".xml";
                BlastXMLParser parser = new BlastXMLParser();
                parser.parseXMLFile(xmlFileStr, baseOrg, queryOrg, baseIdx, queryIdx);


               
            }
        }
        

        
        //C matrix was created.
		
        /*
        for (int i=0; i<C.length; i++) {
           for (int j=0; j<C[i].length; j++) {
               System.out.print(C[i][j]+" ");
           }
           System.out.print("\n");
        }
         */

        
        //P matrix was created.


        if(mclFlag){

            System.out.println("MCL clustering started!");

            boolean success = (new File(GlobalInitializer.imgsDirForCurRun+"/MCL")).mkdir();
            if (success) {
                //System.out.println("Directory: " + (GlobalInitializer.imgsDirForCurRun+"/MCL") + " was created succesfully.\n");
            }

            Runnable r = new MclClusteringThread();
            Thread thr = new Thread(r);
            thr.start();
            Thread.sleep(2000);
            

            while(true){
            if(MainClass.algoThreadsCnt > 0){
                Thread.sleep(500);
                continue;
            }
            else
                break;
            }
            
            
            //System.out.println("MCL clustering completed!");
        } 

        if(emFlag){
            
            System.out.println("EM clustering started!");

            boolean success = (new File(GlobalInitializer.imgsDirForCurRun+"/EM")).mkdir();
            if (success) {
                //System.out.println("Directory: " + (GlobalInitializer.imgsDirForCurRun+"/EM") + " was created succesfully.\n");
            }
            Thread.sleep(1000);

            Runnable r = new EmClusteringThread();
            Thread thr = new Thread(r);
            thr.start();
            Thread.sleep(1000);
            
        }
        
        
	

        //create ecNumsList: for each gene we keep only the EC numbers that exit in the examined pathway.
        int totalEnzsFinal = 0;

        ecNumsList = new String[orgsNum][][];
        for(int i=0; i<ecNumsList.length; i++)
            ecNumsList[i] = new String[pr[i].length][];

        for(int i=0; i<ec.length; i++){
            for(int j=0; j<ec[i].length; j++){
                List<String> tmp = new ArrayList<String>();
                
                for(int k=0; k<ec[i][j].length; k++){
                    for(int m=0; m<enz[i].length; m++){
                        if(ec[i][j][k].equals(enz[i][m])){
                            tmp.add(enz[i][m]);
                            totalEnzsFinal++;
                            break;
                        }
                    }
                }

                ecNumsList[i][j] = new String[tmp.size()];
                String[] rowTmp = new String[tmp.size()];
                rowTmp = tmp.toArray(rowTmp);
                System.arraycopy(rowTmp, 0, ecNumsList[i][j], 0, rowTmp.length);
            }
        }

        
        
        
        /*
        //print (genes, ec numbers) pairs
        
        System.out.println("Genes \t| EC numbers");
        System.out.println("-------------------------");

        for(int i=0; i<ec.length; i++){
            for(int j=0; j<ec[i].length; j++){
                System.out.print(pr[i][j]+": ");
                for(int k=0; k<ec[i][j].length; k++)
                    System.out.print(ec[i][j][k]+" ");
                System.out.println("");
            }
        }


        System.out.println("\n");

        //print ec numbers per pathway.
        for(int i=0; i<enz.length; i++){
            System.out.println("EC numbers for pathway "+MainClass.orgsIds[i]+MainClass.pathId+":\n");
            for(int j=0; j<enz[i].length; j++){
                System.out.println(enz[i][j]);
            }
        }
        */


        ArrayList<String>[] enzList = new ArrayList[totalEnzsFinal];
        for(int i = 0; i < enzList.length; i++) {
            enzList[i] = new ArrayList<String>();
            enzList[i].add((String)"");
        }


        int ecOffset = 0;

        for(int i=0; i<ecNumsList.length; i++){
            for(int j=0; j<ecNumsList[i].length; j++){
                for(int n=0; n<ecNumsList[i][j].length; n++){
                    String curEc = ecNumsList[i][j][n];

                    int flag = 0;
                    for(int k=0; k<enzList.length; k++){
                        if(curEc.equals(enzList[k].get(0))){
                            enzList[k].add((String)MainClass.pr[i][j]);
                            flag = 1;
                            break;
                        }
                    }
                    if(flag == 0){
                        enzList[ecOffset].set(0, curEc);
                        enzList[ecOffset].add((String)MainClass.pr[i][j]);
                        ecOffset++;
                    }
                    else{
                        flag = 0;
                    }
                }
            }
        }


        
        int ec2GenesSize = 0;
        for(int i=0; i<enzList.length; i++){
            if(!(enzList[i].get(0).equals(""))){
                ec2GenesSize++;
            }
            else
                break;
        }

        ecsFinal = new String[ec2GenesSize];
        for(int i=0; i<ecsFinal.length; i++){
            ecsFinal[i] = (String)enzList[i].get(0);
        }

        Arrays.sort(ecsFinal);

        
        ec2GenesList = new String[ec2GenesSize][];

        for(int i=0; i<ec2GenesList.length; i++){
            String ecTmp = ecsFinal[i];

            for(int j=0; j<enzList.length; j++){
                if(ecTmp.equals(enzList[j].get(0))){
                    ec2GenesList[i] = new String[enzList[j].size() - 1];
                    for(int k=0; k<ec2GenesList[i].length; k++)
                        ec2GenesList[i][k] = (String)enzList[j].get(k+1);

                    break;
                }
            }
            
        }

        while(true){
            if(MainClass.algoThreadsCnt > 0){
                Thread.sleep(500);
                continue;
            }
            else
                break;
        }
        
        
        System.out.println("\nInitial Enzymes List (Genes -> EC numbers):");
        
        for(int i=0; i<ecNumsList.length; i++){
            for(int j=0; j<ecNumsList[i].length; j++){
                System.out.print(pr[i][j]+": ");
                for(int k=0; k<ecNumsList[i][j].length; k++)
                    System.out.print(ecNumsList[i][j][k]+" ");
                System.out.println("");
            }
        }

        /*
        System.out.println("\n\nEnzymes List -Unsorted- (EC numbers to Genes):");

        for(int i = 0; i < enzList.length; i++){
            for(int j=0; j<enzList[i].size(); j++){
                System.out.print(enzList[i].get(j)+" ");
            }
            System.out.println("");
        }
        System.out.println("\n");
        */



        String suf = "";

        if(mclFlag){
            System.out.println("printUniqueGenesListsMcl");
            MclClustering.printUniqueGenesLists();
            MainClass.outputLogStrMcl += MainClass.outputLogStrMclAppend;
	    System.out.println("getColoredMapsByClusterMcl");

            
            MclClustering.getColoredMapsByCluster();

            suf="I120";
	    System.out.println("getColoredPathwaysFunctionMcl");
            MclClustering.getColoredPathwaysFunction(suf);

            MclClustering.getPathwayWithBlackElems();
            MclClustering.getMultipleClustersFromECs();
             
        }
        if(emFlag){
            EmClustering.printUniqueGenesLists();
            MainClass.outputLogStrEm += MainClass.outputLogStrEmAppend;
            EmClustering.getColoredMapsByClusterEm();
            suf="";
            
            EmClustering.getColoredPathwaysFunctionEm(suf);
            EmClustering.getPathwayWithBlackElemsEm();
            
        }


        while(true){
            if(getBlackElsThreadMclCnt > 0){
                continue;
            }
            else
                break;
        }

/*
        if(mclFlag){

            PhylogeneticDistance phyloObj = new PhylogeneticDistance();
            phyloObj.getGenesListsForPhyloDist(suf);
            PhylogeneticDistance.createGenomesToClustersMatrix();
            PhylogeneticDistance.printGenomesToClustersMatrix();

            PhylogeneticDistance.createEcsToClustersMatrix();
            PhylogeneticDistance.printEcsToClustersMatrix();

        

            suf="I120";
	    System.out.println("getColoredMapsByClusterWithGroups - Mcl");
          

        }
*/

        
        //*** Terminating program and displaying output ***.

        finalizeAndPrintOutput();

    }

    







    // ********************************************

    // OTHER FUNCTIONS

    // ********************************************


    public static String getGeneNameByIdxInP(int geneIndex, int orgsNum) {
        int geneOffset = 0;
        int orgIdx = 0;
        int geneFinalIdx = geneIndex;
        for (int i = 0; i < orgsNum; i++) {
            geneOffset += pr[i].length;
            if (geneIndex < geneOffset) {
                orgIdx = i;
                break;
            } else {
                geneFinalIdx -= pr[i].length;
            }
        }

        return pr[orgIdx][geneFinalIdx];
    }

    public static String getOrganismNameByIdxInP(int geneIndex, int orgsNum) {

        int geneOffset = 0;
        int orgIdx = 0;
        for (int i = 0; i < orgsNum; i++) {
            geneOffset += pr[i].length;
            if (geneIndex < geneOffset) {
                orgIdx = i;
                break;
            }
        }
        return orgsIds[orgIdx];
    }

    public static int getIdxByNameInP(String geneStr) {
        int index=-1;
        int orgIdx=-1;
        int orgOffsetInP = 0;
        int geneOffsetInP = 0;

        StringTokenizer stokenz = new StringTokenizer(geneStr, ":");
        String orgStr = stokenz.nextToken();
        for(int i=0; i<orgsIds.length; i++){
            if(orgStr.equals(orgsIds[i])){
                orgIdx = i;
                break;
            }
        }

        for(int i=0; i<pr[orgIdx].length; i++){
            if(geneStr.equals(pr[orgIdx][i])){
                geneOffsetInP = i;
                break;
            }
        }

        for(int i=0; i<orgIdx; i++)
            orgOffsetInP += pr[i].length;

        index = orgOffsetInP + geneOffsetInP;

        return index;
    }


    


    


    


   
    public static void finalizeAndPrintOutput(){
        
    
        elapsedTimeMillis = System.currentTimeMillis() - startTime;
        //System.out.println("\n\nTotal time elapsed: " + (elapsedTimeMillis / 1000) + " sec");

        if(mclFlag){

            //print information about phylogenetic distances
            MainClass.outputLogStrMcl += PhylogeneticDistance.phyloDistancesOutputLog;

            //print information about genes from the same genome with same EC but in different clusters
            MainClass.outputLogStrMcl += clustersFromblackECsStr;

            

            MainClass.outputLogStrMcl += "\n\n>> Appendix.\n";
            MainClass.outputLogStrMcl += "> Clusters generated.\n";

            for(int i=0; i<clustersListMcl.size(); i++){
                MainClass.outputLogStrMcl += "\n\n- Cluster "+(i+1)+":\n";
                for(int j=0; j<clustersListMcl.get(i).size(); j++){
                    MainClass.outputLogStrMcl += clustersListMcl.get(i).get(j)+"  ";
                }

            }
            MainClass.outputLogStrMcl += "\n\n\n";

            MainClass.outputLogStrMcl += "> EC numbers -> Genes mapping:\n";


            for(int i = 0; i < ec2GenesList.length; i++){
                MainClass.outputLogStrMcl += ecsFinal[i]+" -> [";
                for(int j=0; j<ec2GenesList[i].length-1; j++){
                    MainClass.outputLogStrMcl += ec2GenesList[i][j]+", ";
                }
                MainClass.outputLogStrMcl += ec2GenesList[i][ec2GenesList[i].length-1]+"]\n";

            }
            MainClass.outputLogStrMcl += "\n";

            MainClass.outputLogStrMcl += "\nTotal time elapsed: " + (elapsedTimeMillis / 1000) + " sec\n";

            System.out.println(MainClass.outputLogStrMcl);

            try{
                String outputFileStr = GlobalInitializer.imgsDirForCurRun+"/MCL/outputLog.txt";
                logFileStream = new FileWriter(outputFileStr);
                BufferedWriter logBuf = new BufferedWriter(logFileStream);
                logBuf.write(MainClass.outputLogStrMcl);
                logBuf.close();

            } catch (Exception e){
                        System.err.println("Error: " + e.getMessage());
            }
        }
        
        System.out.println("\n\n");

        if(emFlag){


            MainClass.outputLogStrEm += clustersFromblackECsStrEm;

            MainClass.outputLogStrEm += "\n\n>> Appendix.\n";
            
            for(int i=0; i<clustersListEm.size(); i++){
                MainClass.outputLogStrEm += "\n\n- Cluster "+(i+1)+":\n";
                for(int j=0; j<clustersListEm.get(i).size(); j++){
                    MainClass.outputLogStrEm += clustersListEm.get(i).get(j)+"  ";
                }

            }
            MainClass.outputLogStrEm += "\n\n\n";


            MainClass.outputLogStrEm += "> EC numbers -> Genes mapping:\n";


            for(int i = 0; i < ec2GenesList.length; i++){
                MainClass.outputLogStrEm += ecsFinal[i]+" -> [";
                for(int j=0; j<ec2GenesList[i].length-1; j++){
                    MainClass.outputLogStrEm += ec2GenesList[i][j]+", ";
                }
                MainClass.outputLogStrEm += ec2GenesList[i][ec2GenesList[i].length-1]+"]\n";

            }
            MainClass.outputLogStrEm += "\n";

            MainClass.outputLogStrEm += "\nTotal time elapsed: " + (elapsedTimeMillis / 1000) + " sec\n";


            System.out.println(MainClass.outputLogStrEm);

            try{
                String outputFileStr = GlobalInitializer.imgsDirForCurRun+"/EM/outputLog.txt";
                logFileStream = new FileWriter(outputFileStr);
                BufferedWriter logBuf = new BufferedWriter(logFileStream);
                logBuf.write(MainClass.outputLogStrEm);
                logBuf.close();

            } catch (Exception e){
                        System.err.println("Error: " + e.getMessage());
            }
        }
    }


    public static void getOrganismsFromKeggByPathway(String refMapId) throws FileNotFoundException, IOException, ServiceException{


        KEGGLocator  locator;
        KEGGPortType serv;
        locator = new KEGGLocator();
        serv = locator.getKEGGPort();
        
        FileInputStream fstream = new FileInputStream("/Users/djifos/Desktop/kegg_org_defs.txt");
        DataInputStream in = new DataInputStream(fstream);
        BufferedReader br = new BufferedReader(new InputStreamReader(in));


        Writer writer = null;
        File file = new File("/Users/djifos/Desktop/genomes_map"+refMapId+".txt");
        writer = new BufferedWriter(new FileWriter(file));


        int cnt = 0;
        String strLine;
        while ((strLine = br.readLine()) != null){

            String orgId = strLine.substring(0,3);

            String query = "path:"+orgId+refMapId;
	    String[] retString = null;

                boolean valid = false;
                while(!valid){
                    try{

                        String[] genesByPathResults = serv.get_genes_by_pathway(query);
                        retString = genesByPathResults;

                        int genesNumber = genesByPathResults.length;

                        System.out.println("Line "+ (++cnt));

                        if(genesNumber>0){
                            writer.write(strLine+"\n");
                            System.out.println(strLine+" : PASS");
                        }
                        else
                            System.out.println(strLine + " : FAIL");

                        valid = true;
                    } catch(Exception e){
                        System.err.println("Error: " + e.getMessage());
                    }
                }

        }


        if (writer != null)
            writer.close();
    }


}






