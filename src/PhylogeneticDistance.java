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

import com.mathworks.toolbox.javabuilder.MWArray;
import com.mathworks.toolbox.javabuilder.MWCellArray;
import com.mathworks.toolbox.javabuilder.MWCharArray;
import com.mathworks.toolbox.javabuilder.MWClassID;
import com.mathworks.toolbox.javabuilder.MWException;
import com.mathworks.toolbox.javabuilder.MWNumericArray;
import com.mathworks.toolbox.javabuilder.services.ServiceException;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.Writer;
import java.net.URL;
import java.rmi.RemoteException;
import java.text.DateFormat;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import keggapi.KEGGLocator;
import keggapi.KEGGPortType;
import phylogeny_pkg.PhyloDistance;


public class PhylogeneticDistance{


    private static int outRowIndex = 0;
    
    public static String[][] dataForMatFun;
    public static String phyloDistancesOutputLog;

    
    public static ArrayList<double[]> distancesList;
    public static ArrayList<Integer> clustIdsList;
    public static ArrayList<Integer> ecIdsList;
    public static ArrayList<String> groupingsList;
    public static ArrayList<Integer> groupIdsList;
    public static ArrayList<Integer> groupCorrespondIdsList;
    public static ArrayList<Integer> genomSizePergroup;
    public static ArrayList<ArrayList<String>> genomeStringsPergroup;
    public static ArrayList<String> koIdsInlineGlobList;
    public static ArrayList<String> koIdsGlobalList;
    public static boolean[] groupsCntForEachEc;


    public static String pathForPhyloTreeImagesDir = GlobalInitializer.imgsDirForCurRun+"/MCL/phylo_trees";
    public static String pathForPhyloTreeGroupDir = "";
    public static String pathForPhyloTreeCustGenomGroupDir = "";
    public static String pathForPhyloTreesDir = "";

    public static String pathForKeggOrgDefs = GlobalInitializer.baseDirectory+"/lib/kegg_org_defs.txt";


    
    public static String perlScriptPath = GlobalInitializer.baseDirectory+"/lib/topd_v3.3.pl";



    private static int[][] genomesToClustersMatrix;
    private static double[][] ecsToClustersMatrix;

    
    public PhylogeneticDistance(){
        
    }



    public static String[] getKoByGene(String gene) throws ServiceException, RemoteException, javax.xml.rpc.ServiceException{

        String[] ko = null;
        boolean valid =false;

        while(!valid){
            try{
                KEGGLocator  locator;
                KEGGPortType serv;

                locator = new KEGGLocator();
                serv = locator.getKEGGPort();

                ko = serv.get_ko_by_gene(gene);
                
                valid = true;
            } catch(Exception e){
                  System.out.println("Error: "+e.getMessage());
            }
        }

        return ko;
    }


    public void getGenesListsForPhyloDist(String suf) throws FileNotFoundException, IOException, InterruptedException, MWException, ServiceException, RemoteException, javax.xml.rpc.ServiceException, Exception{

        phyloDistancesOutputLog = "";
        boolean foundAtLeastOnePureSet = false;
        

        for(int clustId=0; clustId<MainClass.clustersListMcl.size(); clustId++){

            distancesList = new ArrayList<double[]>();
            clustIdsList = new ArrayList<Integer>();
            ecIdsList = new ArrayList<Integer>();
            groupingsList = new ArrayList<String>();
            groupIdsList = new ArrayList<Integer>();
            groupCorrespondIdsList = new ArrayList<Integer>();

            genomSizePergroup = new ArrayList<Integer>();
            genomeStringsPergroup = new ArrayList<ArrayList<String>>();
            koIdsGlobalList = new ArrayList<String>();

            
            int extGenomSize = -1;
            boolean calcgroupEtc = false;
            boolean getColoredMapsWithGroupsFlag = false;


            int genomSize = 0;

            for(int ecId=0; ecId<MainClass.ec2GenesList.length; ecId++){

                

                ArrayList<ArrayList<String>> genesWithComEcAndCluster = new ArrayList<ArrayList<String>>();
                //each 'row' contains genes from a single genome having a common EC number and being in the same cluster
                ArrayList<ArrayList<String>> genesWithComClusterByGenome = new ArrayList<ArrayList<String>>();

                ArrayList<String> tmpGenomesList = new ArrayList<String>();

                for(int curOrg=0; curOrg<MainClass.orgsIds.length; curOrg++){

                    ArrayList<String> genesFromEachGenome = new ArrayList<String>();

                    for(int curGeneId=0; curGeneId<MainClass.ec2GenesList[ecId].length; curGeneId++){

                        String curGene = MainClass.ec2GenesList[ecId][curGeneId];

                            if(curGene.substring(0,3).equals(MainClass.orgsIds[curOrg]) && findClusterByGene(curGene)==clustId){
                                genesFromEachGenome.add(curGene);
                            }

                    }
                    if(genesFromEachGenome.size()>0){
                        genesWithComClusterByGenome.add(genesFromEachGenome);
                        tmpGenomesList.add(MainClass.orgsIds[curOrg]);
                    }
                }

                if(genesWithComClusterByGenome.size() == MainClass.orgsIds.length){
                    genomSize = MainClass.orgsIds.length;
                    break;
                }
                else if(genesWithComClusterByGenome.size() > genomSize)
                    genomSize = genesWithComClusterByGenome.size();
            }


        ArrayList<String> clusterGenomesList = new ArrayList<String>();

        if(genomSize>2){
            for(int ecId=0; ecId<MainClass.ec2GenesList.length; ecId++){      
                
                ArrayList<ArrayList<String>> genesWithComEcAndCluster = new ArrayList<ArrayList<String>>();
                //each 'row' contains genes from a single genome having a common EC number and being in the same cluster
                ArrayList<ArrayList<String>> genesWithComClusterByGenome = new ArrayList<ArrayList<String>>();

                ArrayList<String> tmpGenomesList = new ArrayList<String>();

                for(int curOrg=0; curOrg<MainClass.orgsIds.length; curOrg++){

                    ArrayList<String> genesFromEachGenome = new ArrayList<String>();
                    
                    for(int curGeneId=0; curGeneId<MainClass.ec2GenesList[ecId].length; curGeneId++){

                        String curGene = MainClass.ec2GenesList[ecId][curGeneId];
                        
                            if(curGene.substring(0,3).equals(MainClass.orgsIds[curOrg]) && findClusterByGene(curGene)==clustId){
                                genesFromEachGenome.add(curGene);
                            }

                    }
                    if(genesFromEachGenome.size()>0){
                        genesWithComClusterByGenome.add(genesFromEachGenome);
                        tmpGenomesList.add(MainClass.orgsIds[curOrg]);
                    }
                }

           if(genesWithComClusterByGenome.size() == genomSize){
               
                //retrieve all KO ids for the current EC
                ArrayList<String> koIdsList = new ArrayList<String>();

                ArrayList<ArrayList<String[]>> holdKoIdsList = new ArrayList<ArrayList<String[]>>();

                for(int koI=0; koI<genesWithComClusterByGenome.size(); koI++){
                    
                    ArrayList<String[]> tmpRowKoIdsList = new ArrayList<String[]>();
                    for(int koJ=0; koJ<genesWithComClusterByGenome.get(koI).size(); koJ++){
                        String[] koArrStr = getKoByGene(genesWithComClusterByGenome.get(koI).get(koJ));

                        tmpRowKoIdsList.add(koArrStr);

                            for(int koArrId=0; koArrId<koArrStr.length; koArrId++){
                                boolean addKoId = true;

                                for(int koIdsId=0; koIdsId<koIdsList.size(); koIdsId++){
                                    if((koIdsList.get(koIdsId)).equals(koArrStr[koArrId])){
                                        addKoId = false;
                                        break;
                                    }
                                }

                                if(addKoId)
                                    koIdsList.add(koArrStr[koArrId]);
                            }
                    }

                    holdKoIdsList.add(tmpRowKoIdsList);


                }

                System.out.println("KO ids in EC: "+MainClass.ecsFinal[ecId]+" :");
                for(int koId=0; koId<koIdsList.size(); koId++)
                    System.out.print(koIdsList.get(koId)+" ");

                koIdsInlineGlobList = new ArrayList<String>();

                for(int koId=0; koId<koIdsList.size(); koId++){

                    ArrayList<ArrayList<String>> genesWithComClusterEcKoIdByGenome = new ArrayList<ArrayList<String>>();

                    for(int koI=0; koI<genesWithComClusterByGenome.size(); koI++){
                        ArrayList<String> genesWithComKoByGenome = new ArrayList<String>();

                        for(int koJ=0; koJ<genesWithComClusterByGenome.get(koI).size(); koJ++){
                            String[] curKoIdArr = holdKoIdsList.get(koI).get(koJ);

                            for(int curId=0; curId<curKoIdArr.length; curId++){
                                if(curKoIdArr[curId].equals(koIdsList.get(koId)))
                                    genesWithComKoByGenome.add(genesWithComClusterByGenome.get(koI).get(koJ));
                            }
                        }

                        if(genesWithComKoByGenome.size()>0)
                                genesWithComClusterEcKoIdByGenome.add(genesWithComKoByGenome);
                    }
                   
                    
                    System.out.println("\ngenesWithComClusterEcKoIdByGenome: "+genesWithComClusterEcKoIdByGenome);
                    System.out.println("\ngenomSize: "+genomSize);

                    if(genesWithComClusterEcKoIdByGenome.size() == genomSize){

                        getColoredMapsWithGroupsFlag = true;

                        System.out.println("ready to calc phylogenetic distance");

                        calcgroupEtc = true;
                        extGenomSize = genomSize;
                        
                        pathForPhyloTreeGroupDir = pathForPhyloTreeImagesDir+"/"+genomSize+"_genomes";



                        String timestampForPhyloTrees = "";
                        DateFormat dateFormat = new SimpleDateFormat("dd/MM/yyyyy HH:mm:ss");
                        Date date = new Date();
                        timestampForPhyloTrees += dateFormat.format(date);

                        timestampForPhyloTrees = timestampForPhyloTrees.replace("/", "-");
                        timestampForPhyloTrees = timestampForPhyloTrees.replace(" ", " _ ");
                        timestampForPhyloTrees = timestampForPhyloTrees.replaceFirst(":", "h");
                        timestampForPhyloTrees = timestampForPhyloTrees.replaceFirst(":", "m");
                        timestampForPhyloTrees += "s";


                        String genomSublistStr = "";
                        for(int genSub=0; genSub<tmpGenomesList.size()-1; genSub++){
                            genomSublistStr += tmpGenomesList.get(genSub)+"_";
                        }
                        genomSublistStr += tmpGenomesList.get(tmpGenomesList.size()-1);
                        

                        pathForPhyloTreesDir = pathForPhyloTreeGroupDir+"/trees";


                        System.out.println("\nGenes with common cluster and EC by genome:");
                        for(int ki=0; ki<genesWithComClusterEcKoIdByGenome.size(); ki++){
                            for(int kj=0; kj<genesWithComClusterEcKoIdByGenome.get(ki).size(); kj++)
                                System.out.print(genesWithComClusterEcKoIdByGenome.get(ki).get(kj)+" ");
                            System.out.println("");
                        }

                        
                        //take all possible combinations between elements from all the rows of
                        //genesWithComClusterEcKoIdByGenome.
                        int combinedArraySize = 1;
                        for(int siz=0; siz<genesWithComClusterEcKoIdByGenome.size(); siz++)
                            combinedArraySize *= genesWithComClusterEcKoIdByGenome.get(siz).size();

                        System.out.println("Combined array size: "+combinedArraySize);

                        String[][] combinedArray = new String[combinedArraySize][genomSize];


                        for(int col=genomSize-1; col>=0; col--){

                            int iterations = 1;

                            if(col!=genomSize-1){
                                for(int iter=genomSize-1; iter>col; iter--){
                                    iterations *= genesWithComClusterEcKoIdByGenome.get(iter).size();
                                }
                            }

                            int cnt=0;

                            while(cnt<combinedArraySize){

                                for(int elemIter=0; elemIter<iterations; elemIter++){
                                    for(int inIter=0; inIter<genesWithComClusterEcKoIdByGenome.get(col).size(); inIter++){

                                        combinedArray[cnt][col] = genesWithComClusterEcKoIdByGenome.get(col).get(inIter);
                                        cnt++;
                                    }
                                }

                            }

                        }


                        /*
                         * Combined array:
                         * each row contains genes from all genomes
                         * having the same EC number and belonging to the same cluster
                         */


                        //Add the rows of the combined array to the general
                        //'genesWithComEcAndCluster' matrix which will be
                        //used for the 'pure' phylogenetic distance
                        //calculations with Matlab.
                        System.out.println("Combined array");
                        for(int cr=0; cr<combinedArray.length; cr++){

                            ArrayList<String> rowList = new ArrayList<String>();

                            for(int cc=0; cc<combinedArray[cr].length; cc++){
                                rowList.add(combinedArray[cr][cc]);
                                System.out.print(combinedArray[cr][cc]+" ");
                            }
                            genesWithComEcAndCluster.add(rowList);
                            System.out.println("");
                        }

                        for(int i=0; i<genesWithComEcAndCluster.size(); i++){
                            System.out.print("genesWithComEcAndCluster no. "+i+": ");
                            for(int j=0; j<genesWithComEcAndCluster.get(i).size();j++)
                                System.out.print(genesWithComEcAndCluster.get(i).get(j)+" ");
                            System.out.println("");
                        }

                        clustIdsList.add(clustId);
                        ecIdsList.add(ecId);
                        koIdsInlineGlobList.add(koIdsList.get(koId));
                        
                        System.out.println("->Calculating phylogenetic distance for cluster "+(clustId+1)+" and "+MainClass.ecsFinal[ecId]);
                        calcPhylogeneticDistance(genesWithComEcAndCluster, clustId, ecId, koId, genomSize, tmpGenomesList);

                        
                      }

                    }

                     for(int cpGlobList=0; cpGlobList<koIdsInlineGlobList.size(); cpGlobList++){
                         koIdsGlobalList.add(koIdsInlineGlobList.get(cpGlobList));
                     }

                  }

                      
                }

                System.out.println("ClustIds and ecIds");
                for(int cl=0; cl<clustIdsList.size(); cl++){
                    System.out.println("clustId: "+(clustIdsList.get(cl)+1)+" ec:"+MainClass.ecsFinal[ecIdsList.get(cl)]);
                }
                System.out.print("\n");


                System.out.println("$$$$$$$$$$$$$$$:");
                System.out.println("koIdsGlobalList:");
                for(int glId=0; glId<koIdsGlobalList.size(); glId++)
                    System.out.println(koIdsGlobalList.get(glId));


                
                //calc grouping
                for(int dist=0; dist<distancesList.size(); dist++){
                    System.out.println("calc group iteration:"+(dist+1));

                    int groupId = calcgrouping(clustIdsList.get(dist), ecIdsList.get(dist), dist, genomSize);
                    groupIdsList.add(groupId);


                    System.out.println("\n-----Groups List-----");
                    for(int g=0; g<groupingsList.size(); g++)
                        System.out.println(groupingsList.get(g));

                    System.out.println("\n");

                    System.out.println("/calc group iteration: "+(dist+1));
                }


                for(int rL=0; rL<groupingsList.size(); rL++){
                    phyloDistancesOutputLog += "*** GROUP: "+(rL+1)+" ***\n";

                    phyloDistancesOutputLog += groupingsList.get(rL).toString()+"\n";

                    for(int rId=0; rId<groupIdsList.size(); rId++){
                        if(groupIdsList.get(rId) == rL){
                            phyloDistancesOutputLog += "> Cluster: "+(clustIdsList.get(rId)+1)+", EC: "+MainClass.ecsFinal[ecIdsList.get(rId)]+"\n";
                        }
                    }
                    phyloDistancesOutputLog += "------------------------\n\n";
                }

                phyloDistancesOutputLog += "\n";

                phyloDistancesOutputLog += calcAverageDistancesByGroup(extGenomSize);

                phyloDistancesOutputLog += "\n";

                //add the corresponding group id to the tree images' name.
                renamePhyloTreeImages();


                groupsCntForEachEc = new boolean[ecIdsList.size()];
                for(int grCntId=0; grCntId<groupsCntForEachEc.length; grCntId++)
                    groupsCntForEachEc[grCntId] = true;

                for(int chkGroupsExt=0; chkGroupsExt<ecIdsList.size(); chkGroupsExt++){

                    for(int chkGroupsInt=0; chkGroupsInt<ecIdsList.size(); chkGroupsInt++){
                        
                        if(chkGroupsExt!=chkGroupsInt){
                            if((ecIdsList.get(chkGroupsExt)).equals(ecIdsList.get(chkGroupsInt))){
                                if(!(groupIdsList.get(chkGroupsExt)).equals(groupIdsList.get(chkGroupsInt))){
                                    groupsCntForEachEc[chkGroupsExt] = false;
                                }
                            }
                        }
                        
                    }
                }
                

                if(getColoredMapsWithGroupsFlag)
                    getColoredMapsByClusterWithGroups(clustId);

                getColoredPathwaysFunctionWithGroups(suf);
                    
            }            

            }

    }




    /*
     * Calculate phylogenetic distances using the Matlab-generated class PhyloDistance
     */
    public static void calcPhylogeneticDistance(ArrayList<ArrayList<String>> genesWithComEcAndCluster, int clustId, int ecId, int koId, int genomSize, ArrayList<String> tmpGenomesList){



        try{
            
            
            double[] averagePhyloDistances = new double[genomSize *(genomSize-1)/2];
            for(int i=0; i<averagePhyloDistances.length; i++)
                averagePhyloDistances[i] = 0.0;

            
            int geneGroupsForMatlab;

            for(geneGroupsForMatlab=1; geneGroupsForMatlab*genomSize<2000; geneGroupsForMatlab++){
            }

            int iterations = genesWithComEcAndCluster.size()/geneGroupsForMatlab;
            int rem = genesWithComEcAndCluster.size() - (iterations*geneGroupsForMatlab);

            System.out.println("Iterations: "+iterations);
            System.out.println("Remaining:"+rem);

            for(int iter=0; iter<iterations; iter++){

                MWArray cAr = new MWCellArray(geneGroupsForMatlab*genomSize, 2);
                int cArIndex = 1;

                for(int i=iter*geneGroupsForMatlab; i<((iter*geneGroupsForMatlab) + geneGroupsForMatlab); i++){

                    dataForMatFun = new String[genomSize][2];

                    for(int j=0; j<genesWithComEcAndCluster.get(i).size(); j++){

                        String gene = genesWithComEcAndCluster.get(i).get(j);
                        
                        
                        dataForMatFun[j][0] = gene.substring(0,3);
                        dataForMatFun[j][1] = getGeneSequenceFromFasta(gene);

                    }


                    for (int m=0; m<dataForMatFun.length; m++){
                        for(int n=0; n<dataForMatFun[m].length; n++)
                        {

                            int[] index = {cArIndex, n+1};
                            System.out.println("cArIndex: "+cArIndex);
                            cAr.set(index, dataForMatFun[m][n]);

                        }
                        cArIndex++;
                        System.out.println("cArIndex++: "+cArIndex);
                    }
                }


                System.out.println("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
                System.out.println("cAr set!");
                
                System.out.println("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");

                PhyloDistance ph = new PhyloDistance();

                Object[] o = ph.calc_phylogen_dist(1, cAr, genomSize);

                MWNumericArray x = (MWNumericArray)o[0];

                System.out.println("\nReturned distances matrix: ");
                System.out.print(x.toString());

                int[] dims = x.getDimensions();
                System.out.println("\ndims[0]: "+dims[0]+", dims[1]: "+dims[1]+"\n");


                //'distances' is a 1-by-(M*(M-1)/2) row vector corresponding to
                //the M*(M-1)/2 pairs of sequences in Seqs. The output 'distances' is arranged
                //in the order ((2,1),(3,1),..., (M,1),(3,2),...(M,2),...(M,M-1)).
                //This is the lower-left triangle of the full M-by-M distance matrix.
                //To get the distance between the Ith and the Jth sequences for I > J,
                //use the formula distances((J-1)*(M-J/2)+I-J).
                double[] distances = new double[dims[0] * dims[1]];


                int distIndex = 0;

                for(int distIdx=0; distIdx<dims[0]; distIdx++){
                    for(int distIdy=0; distIdy<dims[1]; distIdy++){
                        distances[distIndex++] = x.getDouble(distIdy+1);
                    }
                }

                PhyloDistance.disposeAllInstances();

                System.out.println("Iteration "+(iter+1)+"/"+(iterations+1)+", distances matrix:");
                for(int distId=0; distId<distances.length; distId++){
                    System.out.print(distances[distId]+" ");
                }

                for(int avgIdx=0; avgIdx<averagePhyloDistances.length; avgIdx++)
                    averagePhyloDistances[avgIdx] += distances[avgIdx];

            }


            //calc distances for the remaining genes
            MWArray cAr = new MWCellArray(rem*genomSize, 2);
            int cArIndex = 1;
            System.out.println("->Calc distances for the remaining genes:");
            for(int i=iterations*geneGroupsForMatlab; i<((iterations*geneGroupsForMatlab) + rem); i++){

                System.out.println("in rem, i:"+(i+1));

                dataForMatFun = new String[genomSize][2];

                for(int j=0; j<genesWithComEcAndCluster.get(i).size(); j++){

                    String gene = genesWithComEcAndCluster.get(i).get(j);

                    dataForMatFun[j][0] = gene.substring(0, 3);
                    dataForMatFun[j][1] = getGeneSequenceFromFasta(gene);

                }


                for (int m=0; m<dataForMatFun.length; m++){
                    for(int n=0; n<dataForMatFun[m].length; n++)
                    {

                        int[] index = {cArIndex, n+1};
                        System.out.println("cArIndex: "+cArIndex);
                        cAr.set(index, dataForMatFun[m][n]);

                    }
                    cArIndex++;
                    System.out.println("cArIndex++: "+cArIndex);
                }

                System.out.println("debug test");
            }

            System.out.println("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
            System.out.println("cAr set!");
            
            System.out.println("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");

            PhyloDistance ph = new PhyloDistance();

            Object[] o = ph.calc_phylogen_dist(1, cAr, genomSize);

            MWNumericArray x = (MWNumericArray)o[0];

            System.out.println("\nReturned distances matrix: ");
            System.out.print(x.toString());

            int[] dims = x.getDimensions();
            System.out.println("\ndims[0]: "+dims[0]+", dims[1]: "+dims[1]+"\n");


            //'distances' is a 1-by-(M*(M-1)/2) row vector corresponding to
            //the M*(M-1)/2 pairs of sequences in Seqs. The output 'distances' is arranged
            //in the order ((2,1),(3,1),..., (M,1),(3,2),...(M,2),...(M,M-1)).
            //This is the lower-left triangle of the full M-by-M distance matrix.
            //To get the distance between the Ith and the Jth sequences for I > J,
            //use the formula distances((J-1)*(M-J/2)+I-J).
            double[] distances = new double[dims[0] * dims[1]];


            int distIndex = 0;

            for(int distIdx=0; distIdx<dims[0]; distIdx++){
                for(int distIdy=0; distIdy<dims[1]; distIdy++){
                    distances[distIndex++] = x.getDouble(distIdy+1);
                }
            }

            PhyloDistance.disposeAllInstances();

                System.out.println("Remaining genes, distances matrix:");
                for(int distId=0; distId<distances.length; distId++){
                    System.out.print(distances[distId]+" ");
                }

                for(int avgIdx=0; avgIdx<averagePhyloDistances.length; avgIdx++)
                    averagePhyloDistances[avgIdx] += distances[avgIdx];

                
            //get final distances matrix
            for(int avgIdx=0; avgIdx<averagePhyloDistances.length; avgIdx++)
                    averagePhyloDistances[avgIdx] /= (iterations+1);

            
            System.out.println("before distance add");
            distancesList.add(averagePhyloDistances);
            System.out.println("after distance add");

            
            //create phylogenetic tree
            createPhyloTree(averagePhyloDistances, clustId, ecId, koId, "", genomSize, tmpGenomesList);

            //print information about the calculated phylogenetic distances
            printPhyloDistances(averagePhyloDistances, clustId, ecId, genomSize);
            

        }  catch (Exception e){
            System.err.println("Error: " + e.getMessage());
        }
        finally{

        }
    }


    private static String getTreeImageNameForCalcgroup(int clustId, int ecId, int koId, String extraDefStr){

        if(extraDefStr.equals("")){
            String trImgName = pathForPhyloTreeGroupDir+"/PhylTree.Cluster"+(clustId+1)+".";

            String ec = MainClass.ecsFinal[ecId];
            ec = ec.replace( ':', '_' );
            trImgName += ec;

            String ko = koIdsGlobalList.get(koId);
            ko = ko.replace( "ko:", "koId_" );
            trImgName += "."+ko;

            trImgName += ".png";

            return trImgName;
        } else{
            String trImgName = pathForPhyloTreeGroupDir+"/"+extraDefStr+".png";
            return trImgName;
        }


    }

    private static String getTreeFileNameForCalcgroup(int clustId, int ecId, int koId, String extraDefStr){

        if(extraDefStr.equals("")){
            String trFileName =  pathForPhyloTreesDir+"/PhylTree.Cluster"+(clustId+1)+".";

            String ec = MainClass.ecsFinal[ecId];
            ec = ec.replace( ':', '_' );
            trFileName += ec;

            String ko = koIdsGlobalList.get(koId);
            ko = ko.replace( "ko:", "koId_" );
            trFileName += "."+ko;

            trFileName += ".txt";

            return trFileName;
        } else{
            String trFileName = pathForPhyloTreesDir+"/"+extraDefStr+".txt";
            return trFileName;
        }
    }

    
    private static String getTreeImageName(int clustId, int ecId, int koId, String extraDefStr){

        if(extraDefStr.equals("")){
            String trImgName = pathForPhyloTreeGroupDir+"/PhylTree.Cluster"+(clustId+1)+".";

            String ec = MainClass.ecsFinal[ecId];
            ec = ec.replace( ':', '_' );
            trImgName += ec;

            String ko = koIdsInlineGlobList.get(koId);
            ko = ko.replace( "ko:", "koId_" );
            trImgName += "."+ko;

            trImgName += ".png";
            
            return trImgName;
        } else{
            String trImgName = pathForPhyloTreeGroupDir+"/"+extraDefStr+".png";
            return trImgName;
        }

        
    }


    private static String getTreeFileName(int clustId, int ecId, int koId, String extraDefStr){

        if(extraDefStr.equals("")){
            String trFileName =  pathForPhyloTreesDir+"/PhylTree.Cluster"+(clustId+1)+".";
            
            String ec = MainClass.ecsFinal[ecId];
            ec = ec.replace( ':', '_' );
            trFileName += ec;

            String ko = koIdsInlineGlobList.get(koId);
            ko = ko.replace( "ko:", "koId_" );
            trFileName += "."+ko;

            trFileName += ".txt";

            return trFileName;
        } else{
            String trFileName = pathForPhyloTreesDir+"/"+extraDefStr+".txt";
            return trFileName;
        }
    }


    private static void createPhyloTree(double[] averagePhyloDistances, int clustId, int ecId, int koId, String extraDefStr, int genomSize, ArrayList<String> tmpGenomesList) throws MWException, IOException{


        
        String pathForPhyloTreeImage = getTreeImageName(clustId, ecId, koId, extraDefStr);
        String pathForPhyloTreeFile = getTreeFileName(clustId, ecId, koId, extraDefStr);

        System.out.println("pathForPhyloTreeImage: "+pathForPhyloTreeImage);
        System.out.println("pathForPhyloTreeFile: "+pathForPhyloTreeFile);


        boolean successImg = (new File(pathForPhyloTreeImagesDir)).mkdir();
        if (successImg) {
            System.out.println("Directory: " + pathForPhyloTreeImagesDir + " was created succesfully.");
        }

        boolean successGroup = (new File(pathForPhyloTreeGroupDir)).mkdir();
        if (successGroup) {
            System.out.println("Directory: " + pathForPhyloTreeGroupDir + " was created succesfully.");
        }


        
        //create a .txt file inside the 'pathForPhyloTreeCustGenomGroupDir' folder,
        //containing the genomes of the current test case.
        String genomesStringIdentifier = "";
        ArrayList<String> genomesSortedList = new ArrayList<String>();
        for(int copyList=0; copyList<tmpGenomesList.size(); copyList++)
            genomesSortedList.add(tmpGenomesList.get(copyList));
        Collections.sort(genomesSortedList);


        Writer writer = null;

        File file = new File(pathForPhyloTreeGroupDir+"/genomes_included.txt");
        writer = new BufferedWriter(new FileWriter(file));

        for(int m=0; m<genomesSortedList.size(); m++){

            String orgDef = "";

            FileInputStream fstream = new FileInputStream(pathForKeggOrgDefs);
            DataInputStream in = new DataInputStream(fstream);
            BufferedReader br = new BufferedReader(new InputStreamReader(in));

            String strLine;
            while ((strLine = br.readLine()) != null){
                if(strLine.startsWith(genomesSortedList.get(m))){
                    orgDef = strLine.substring(5,strLine.length());
                    break;
                }
            }

            writer.write(genomesSortedList.get(m)+" "+orgDef+"\n");
        }


        if (writer != null)
            writer.close();
        
        

        ArrayList<String> exludedGenomesList = new ArrayList<String>();

        for(int searchMainOrgs=0; searchMainOrgs<MainClass.orgsIds.length; searchMainOrgs++){

            String curMainOrg = MainClass.orgsIds[searchMainOrgs];
            boolean included = false;

            for(int searchTmpGenomes=0; searchTmpGenomes<tmpGenomesList.size(); searchTmpGenomes++){

                String curTmpOrg = tmpGenomesList.get(searchTmpGenomes);
                if(curMainOrg.equals(curTmpOrg)){
                    included = true;
                    break;
                }
            }

            if(!included)
                exludedGenomesList.add(MainClass.orgsIds[searchMainOrgs]);
        }


        Writer writerEx = null;
        File fileEx = new File(pathForPhyloTreeGroupDir+"/genomes_excluded.txt");
        writerEx = new BufferedWriter(new FileWriter(fileEx));

        for(int exclGenome=0; exclGenome<exludedGenomesList.size(); exclGenome++){
            String orgDef = "";

            FileInputStream fstream = new FileInputStream(pathForKeggOrgDefs);
            DataInputStream in = new DataInputStream(fstream);
            BufferedReader br = new BufferedReader(new InputStreamReader(in));

            String strLine;
            while ((strLine = br.readLine()) != null){
                if(strLine.startsWith(exludedGenomesList.get(exclGenome))){
                    orgDef = strLine.substring(5,strLine.length());
                    break;
                }
            }

            writerEx.write(exludedGenomesList.get(exclGenome)+" "+orgDef+"\n");

        }

        if (writerEx != null)
            writerEx.close();



        boolean successTrees = (new File(pathForPhyloTreesDir)).mkdir();
        if (successTrees) {
            System.out.println("Directory: " + pathForPhyloTreesDir + " was created succesfully.");
        }

        dataForMatFun = new String[genomSize][2];

        for(int genomId=0; genomId<genomSize; genomId++){


            String orgId = tmpGenomesList.get(genomId);
            String orgDef = "";

            FileInputStream fstream = new FileInputStream(pathForKeggOrgDefs);
            DataInputStream in = new DataInputStream(fstream);
            BufferedReader br = new BufferedReader(new InputStreamReader(in));

            String strLine;
            while ((strLine = br.readLine()) != null)   {
                if(strLine.startsWith(orgId)){
                    orgDef = strLine.substring(5,strLine.length());
                    break;
                }
            }

            in.close();

            System.out.println("org definition: "+orgDef);

            dataForMatFun[genomId][0] = orgDef+" ("+orgId+")";
            dataForMatFun[genomId][1] = "";

        }

        MWArray cAr = new MWCellArray(genomSize, 2);

        for (int m=0; m<dataForMatFun.length; m++){
            for(int n=0; n<dataForMatFun[m].length; n++)
            {

                int[] index = {m+1, n+1};
                cAr.set(index, dataForMatFun[m][n]);

            }
        }


        
        MWNumericArray distancesForMatlab = new MWNumericArray(averagePhyloDistances, MWClassID.DOUBLE);
        MWCharArray pathForMatlab = new MWCharArray(pathForPhyloTreeImage);


        PhyloDistance ph = new PhyloDistance();

        Object[] o = ph.create_phylogenetic_tree(1, cAr, distancesForMatlab, pathForMatlab, pathForPhyloTreeFile);


        MWCharArray x = (MWCharArray)o[0];

        System.out.println("\nPath where the tree image is stored: ");
        System.out.println(x.toString());

        PhyloDistance.disposeAllInstances();

        System.out.println("Phylogenetic tree image for cluster "+(clustId+1)+" and EC "+ecId+" was created succesfully!");

    }


    public static void printPhyloDistances(double[] averagePhyloDistances, int clustId, int ecId, int genomSize){


        DecimalFormat twoDForm = new DecimalFormat("#.###");

        System.out.println("\n\n------------------------->\n");
        System.out.println("Average phylogenetic distances:");
        for(int i=0; i<averagePhyloDistances.length; i++)
            System.out.print(averagePhyloDistances[i]+" ");
        System.out.println("");
        System.out.println("\n<-------------------------\n\n");



        phyloDistancesOutputLog += "\nCluster: "+(clustId+1)+", EC: "+MainClass.ecsFinal[ecId]+"\n";

        //write distances output in matrix display format to log file
        phyloDistancesOutputLog += "\n- Average phylogenetic distances:\n\n";
        phyloDistancesOutputLog += "Genomes  |   ";

        for(int k=0; k<genomSize; k++)
            phyloDistancesOutputLog +=  MainClass.orgsIds[k]+"   |   ";

        phyloDistancesOutputLog += "\n";

        for(int k=0; k<(genomSize*10+10); k++){
            phyloDistancesOutputLog += "-";
        }


        int distRegCounter = 0;

        for(int i=0; i<genomSize; i++){
            phyloDistancesOutputLog += "\n   "+MainClass.orgsIds[i]+"   |  ";

            //between the Ith and the Jth sequences for I > J,
            //use the formula distances((J-1)*(M-J/2)+I-J).
            int j=0;
            while(j<i){
                double dist = Double.valueOf(twoDForm.format(averagePhyloDistances[distRegCounter++]));

                //count number of decimal places
                String tmp = Double.toString(dist);
                String[] res = tmp.split("\\.");
                phyloDistancesOutputLog += dist;
                if(res[1].length() == 4)
                    phyloDistancesOutputLog += " |  ";
                else if(res[1].length() == 3)
                        phyloDistancesOutputLog += "  |  ";
                else if(res[1].length() == 2)
                        phyloDistancesOutputLog += "   |  ";
                else
                    phyloDistancesOutputLog += "    |  ";

                j++;
            }
        }
        phyloDistancesOutputLog += "\n\n";

        System.out.println(phyloDistancesOutputLog);

    }

    public static int calcgrouping(int clustId, int ecId, int distId, int genomSize) throws FileNotFoundException, IOException, InterruptedException{


            int indexIngroupList = -1;
            boolean newgroupMat = true;

            String curTreeWrappedString = "";
            String curTreePath = getTreeFileNameForCalcgroup(clustId, ecId, distId, "");

            System.out.println("curTreePath: "+curTreePath);

            FileInputStream fstream = new FileInputStream(curTreePath);
            DataInputStream in = new DataInputStream(fstream);
            BufferedReader br = new BufferedReader(new InputStreamReader(in));

            String strLine;
            
            while ((strLine = br.readLine()) != null)   {
                curTreeWrappedString += strLine;
            }
            
            in.close();


            for(int group=0; group<groupingsList.size(); group++){

                String baseTreeWrappedString = "";
                String baseTreePath = getTreeFileNameForCalcgroup(clustIdsList.get(groupCorrespondIdsList.get(group)), ecIdsList.get(groupCorrespondIdsList.get(group)), groupCorrespondIdsList.get(group), "");

                System.out.println("baseTreePath: "+baseTreePath);

                FileInputStream fBasestream = new FileInputStream(baseTreePath);
                DataInputStream inBase = new DataInputStream(fBasestream);
                BufferedReader brBase = new BufferedReader(new InputStreamReader(inBase));

                String baseStrLine;
                
                while ((baseStrLine = brBase.readLine()) != null)   {
                    baseTreeWrappedString += baseStrLine;
                }
                inBase.close();

                System.out.println("curTreeWrappedString: "+curTreeWrappedString);
                System.out.println("baseTreeWrappedString: "+baseTreeWrappedString);

                if(curTreeWrappedString.equals(baseTreeWrappedString)){

                    System.out.println("Not new group. cur: "+curTreeWrappedString+", base: "+baseTreeWrappedString);
                    indexIngroupList = group;
                    newgroupMat = false;
                    break;
                }
            }

            

            if(newgroupMat){
                groupingsList.add(curTreeWrappedString);
                int genomSizeForCurgroup;

                String simplifStr = curTreeWrappedString;

                simplifStr = simplifStr.replace("(","");
                simplifStr = simplifStr.replace(")","");
                simplifStr = simplifStr.replace(",","");
                simplifStr = simplifStr.replace(";","");

                genomSizeForCurgroup = simplifStr.length()/genomSize;

                genomSizePergroup.add(genomSizeForCurgroup);
                ArrayList<String> genomesForCurgroup = new ArrayList<String>();
                
                for(int strIndex=0; strIndex<genomSizeForCurgroup; strIndex++){
                    genomesForCurgroup.add(simplifStr.substring(strIndex*3, (strIndex*3)+3));
                }


                genomeStringsPergroup.add(genomesForCurgroup);

                indexIngroupList = groupingsList.size()-1;
                groupCorrespondIdsList.add(distId);
                System.out.println("New group. cur: "+curTreeWrappedString+" groupIndex: "+indexIngroupList);

            }
           
            return indexIngroupList;
    }

    
    private static String calcAverageDistancesByGroup(int genomSize) throws MWException{

        DecimalFormat twoDForm = new DecimalFormat("#.###");
        String distByGroupStr = "";

        for(int rId=0; rId<groupingsList.size(); rId++){

            int listCnt = 0;
            distByGroupStr += "GROUP "+(rId+1)+":\n";
            
            int M = genomSize;
            double[] avgDistbyGroup = new double[M*(M-1)/2];
            for(int avgId=0; avgId<avgDistbyGroup.length; avgId++)
                avgDistbyGroup[avgId] = 0;
            
            for(int d=0; d<distancesList.size(); d++){
                if(groupIdsList.get(d) == rId){
                    listCnt++;
                    for(int i=0; i<distancesList.get(d).length; i++)
                        avgDistbyGroup[i] += distancesList.get(d)[i];
                }
            }

            for(int avgId=0; avgId<avgDistbyGroup.length; avgId++)
                avgDistbyGroup[avgId] /= listCnt;

            distByGroupStr += "\n- Average phylogenetic distances:\n\n";
            distByGroupStr += "Genomes  |   ";

            for(int k=0; k<genomSize; k++)
                distByGroupStr +=  MainClass.orgsIds[k]+"   |   ";

            distByGroupStr += "\n";

            for(int k=0; k<(genomSize*10+10); k++){
                distByGroupStr += "-";
            }


            int distRegCounter = 0;

            for(int i=0; i<genomSize; i++){
                distByGroupStr += "\n   "+MainClass.orgsIds[i]+"   |  ";

                //between the Ith and the Jth sequences for I > J,
                //use the formula distances((J-1)*(M-J/2)+I-J).
                int j=0;
                while(j<i){
                    double dist = Double.valueOf(twoDForm.format(avgDistbyGroup[distRegCounter++]));

                    //count number of decimal places
                    String tmp = Double.toString(dist);
                    String[] res = tmp.split("\\.");
                    distByGroupStr += dist;
                    if(res[1].length() == 4)
                        distByGroupStr += " |  ";
                    else if(res[1].length() == 3)
                            distByGroupStr += "  |  ";
                    else if(res[1].length() == 2)
                            distByGroupStr += "   |  ";
                    else
                        distByGroupStr += "    |  ";

                    j++;
                }
            }
            distByGroupStr += "\n\n";

            String extraDefStr = "Group_"+(rId+1)+"_avg_dist";
            
        }

        return distByGroupStr;
    }


    public static void renamePhyloTreeImages(){

        for(int dist=0; dist<distancesList.size(); dist++){

            String oldImageName = getTreeImageNameForCalcgroup(clustIdsList.get(dist), ecIdsList.get(dist), dist, "");


            String newImageName = pathForPhyloTreeGroupDir+"/Group_"+(groupIdsList.get(dist)+1)+".Tree.Cl_"+(clustIdsList.get(dist)+1)+".";
            String ec = MainClass.ecsFinal[ecIdsList.get(dist)];
            ec = ec.replace( ':', '_' );
            newImageName += ec;

            String ko = koIdsGlobalList.get(dist);
            ko = ko.replace( "ko:", "koId_" );
            newImageName += "."+ko;

            newImageName += ".png";

            
            File file = new File(oldImageName);
            File file2 = new File(newImageName);

            boolean success = file.renameTo(file2);
            if (success) {
                System.out.println("Image "+oldImageName+" was succesfully renamed to "+newImageName);
            }

        }
    }




    public static void getColoredMapsByClusterWithGroups(int clustId) throws RemoteException, IOException, ServiceException, javax.xml.rpc.ServiceException{

            try{
                boolean success = (new File((GlobalInitializer.imgsDirForCurRun)+"/MCL/images_with_groups").mkdir());
                if (success) {
                    System.out.println("Directory: " + (GlobalInitializer.imgsDirForCurRun)+"/MCL/images_with_groups" + " was created succesfully.\n");
                }

            } catch (Exception e){
                System.err.println("Error: " + e.getMessage());
            }

            String whiteFgColor = "#ffffff";

            int baseFgColor = 10000000;
            int fgColorOffset = 2000000;


            String[] bgColors = {"#0000ff", "#ff0000", "#f0f000", "#00f0f0", "#f000f0", "#0f0f00", "000f0f", "0f000f", "fff000", "0fff00", "00fff0", "000fff"};

            
            int i = clustId;
            
                System.out.println("test simul colored maps");
                try {
                    KEGGLocator  locator;
                    KEGGPortType serv;
                    locator = new KEGGLocator();
                    serv = locator.getKEGGPort();

                    ArrayList<String> obj_list = new ArrayList<String>();
                    ArrayList<String> fg_list = new ArrayList<String>();
                    ArrayList<String> bg_list = new ArrayList<String>();


                    int objIndex = 0;
                    for (int j = 0; j < MainClass.clustersListMcl.get(i).size(); j++) {
                        String curGene = MainClass.clustersListMcl.get(i).get(j);
                        int genomeIndex = -1;
                        String genomeString = curGene.substring(0,3);

                        for (int iter1 = 0; iter1 < MainClass.orgsIds.length; iter1++) {
                            if(genomeString.equals(MainClass.orgsIds[iter1]))
                                genomeIndex = iter1;
                        }
                        
                        for (int iter2 = 0; iter2 < MainClass.pr[genomeIndex].length; iter2++) {
                            if (((String) MainClass.clustersListMcl.get(i).get(j)).equals(MainClass.pr[genomeIndex][iter2])) {
                                if (MainClass.ecNumsList[genomeIndex][iter2].length > 0) {

                                    for (int k = 0; k < MainClass.ecNumsList[genomeIndex][iter2].length; k++) {
                                        int blackFlag = 0;
                                        boolean addEc = true;

                                        for(int n=0; n<obj_list.size(); n++){
                                            if(obj_list.get(n).equals(MainClass.ecNumsList[genomeIndex][iter2][k])){
                                                addEc = false;
                                                break;
                                            }
                                        }
                                        if(addEc){
                                            obj_list.add(MainClass.ecNumsList[genomeIndex][iter2][k]);
                                            String curEC = MainClass.ecNumsList[genomeIndex][iter2][k];
                                            int ecIdxExam = -1;

                                            for (int t1 = 0; t1 < MainClass.ecsFinal.length; t1++) {
                                                if (curEC.equals(MainClass.ecsFinal[t1])) {
                                                    ecIdxExam = t1;
                                                    break;
                                                }
                                            }
                                            for (int t1 = 0; t1 < MainClass.ec2GenesList[ecIdxExam].length; t1++) {
                                                String examGene = MainClass.ec2GenesList[ecIdxExam][t1];

                                                for (int ii = 0; ii < MainClass.clustersListMcl.size(); ii++) {
                                                    if (ii != i) {
                                                        for (int jj = 0; jj < MainClass.clustersListMcl.get(ii).size(); jj++) {

                                                            if ((examGene.equals(MainClass.clustersListMcl.get(ii).get(jj))) && (!examGene.equals(curGene))) {
                                                                blackFlag = 1;
                                                                MainClass.blackECsListMcl.add(MainClass.ecsFinal[ecIdxExam]);
                                                                MainClass.blackCurClustersListMcl.add(i);
                                                                MainClass.blackClustersListMcl.add(ii);
                                                                MainClass.blackCurGeneList.add(curGene);
                                                                MainClass.blackExamGeneList.add(examGene);

                                                                break;
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                            if (blackFlag == 1) {
                                                blackFlag = 0;
                                                bg_list.add("#0f0f00");
                                            } else {
                                                bg_list.add(bgColors[i]);
                                            }

                                            int groupIndex = -1;
                                            boolean groupParticipant = false;

                                            for(int findDistId = 0; findDistId<distancesList.size(); findDistId++){

                                                if(groupsCntForEachEc[findDistId]){
                                                    if((MainClass.ecNumsList[genomeIndex][iter2][k]).equals(MainClass.ecsFinal[ecIdsList.get(findDistId)])){

                                                        int cntCurEc = 0;
                                                        for(int chkEcId=0; chkEcId<ecIdsList.size(); chkEcId++){
                                                            if(ecIdsList.get(findDistId).equals(ecIdsList.get(chkEcId)))
                                                                cntCurEc++;
                                                        }
                                                        groupIndex = findDistId;
                                                        groupParticipant = true;
                                                        break;
                                                    }
                                                }
                                            }

                                            if(groupParticipant){
                                                int finalFgColor = baseFgColor - (fgColorOffset*groupIdsList.get(groupIndex));
                                                String finalFgColorStr = "#"+Integer.toHexString(finalFgColor);

                                                fg_list.add(finalFgColorStr);
                                            } else{
                                                fg_list.add(whiteFgColor);
                                            }

                                            objIndex++;
                                        }
                                    }
                                } else {

                                        String[] res = serv.get_ko_by_gene(MainClass.clustersListMcl.get(i).get(j));
                                        for (int k = 0; k < res.length; k++) {
                                            obj_list.add(res[k]);
                                            bg_list.add(bgColors[i]);

                                            int groupIndex = -1;
                                            boolean groupParticipant = false;

                                            for(int findDistId = 0; findDistId<distancesList.size(); findDistId++){
                                                if(groupsCntForEachEc[findDistId]){
                                                    if((res[k]).equals(koIdsGlobalList.get(findDistId))){
                                                        groupIndex = findDistId;
                                                        groupParticipant = true;
                                                        break;
                                                    }
                                                }
                                            }

                                            if(groupParticipant){
                                                int finalFgColor = baseFgColor - (fgColorOffset*groupIdsList.get(groupIndex));
                                                String finalFgColorStr = "#"+Integer.toHexString(finalFgColor);

                                                fg_list.add(finalFgColorStr);
                                            } else{
                                                fg_list.add(whiteFgColor);
                                            }
                                            
                                            objIndex++;
                                        }

                                }
                            }
                        }
                    }

                    System.out.println("obj_list size (with groups): "+ obj_list.size());
                    System.out.println("fg_list size (with groups): "+ fg_list.size());
                    System.out.println("bg_list size (with groups): "+ bg_list.size());



                    String[] obj_arr = new String[obj_list.size()];
                    String[] fg_arr = new String[fg_list.size()];
                    String[] bg_arr = new String[bg_list.size()];

                    obj_list.toArray(obj_arr);
                    fg_list.toArray(fg_arr);
                    bg_list.toArray(bg_arr);

                    System.out.println("Path String (with groups): "+"map"+MainClass.pathId);

                    System.out.println("obj_arr (with groups):");
                    for(int itest=0; itest<obj_arr.length; itest++){
                        System.out.print(obj_arr[itest]+" ");
                    }

                    System.out.println("\nfg_arr (with groups):");
                    for(int itest=0; itest<fg_arr.length; itest++){
                        System.out.print(fg_arr[itest]+" ");
                    }

                    System.out.println("\nbg_arr (with groups):");
                    for(int itest=0; itest<bg_arr.length; itest++){
                        System.out.print(bg_arr[itest]+" ");
                    }


                    boolean valid = false;
                    while(!valid){
                        try{
                            System.out.println("downloading map img (with groups), cluster: "+(i+1));
                            String pathStr = "map"+MainClass.pathId;
                            String url;
                            url = serv.color_pathway_by_objects(pathStr, obj_arr, fg_arr, bg_arr);
                            System.out.println(url);
                            String imFileStr = GlobalInitializer.imgsDirForCurRun + "/MCL/images_with_groups/cluster_" + (i + 1) + "_map.png";
                            saveImage(url, imFileStr);

                            valid = true;
                        } catch (Exception e){
                            System.err.println("Error: "+e.getMessage());
                        }
                    }
                    System.out.println("/downloading map img (with groups), cluster: "+(i+1));
                } catch (IOException ex) {
                    Logger.getLogger(MclClustering.class.getName()).log(Level.SEVERE, null, ex);
                }

    }




    public static void getColoredPathwaysFunctionWithGroups(String suf) throws Exception{



            boolean success = (new File((GlobalInitializer.imgsDirForCurRun)+"/MCL/images_with_groups").mkdir());
            if (success) {
                System.out.println("Directory: " + (GlobalInitializer.imgsDirForCurRun)+"/MCL/images_with_groups" + " was created succesfully.\n");
            }


             class nestedGetColorPathsThread implements Runnable{

                String genomeStr;
                String suf;

                
                ArrayList<String> obj_list = new ArrayList<String>();
                ArrayList<String> fg_list = new ArrayList<String>();
                ArrayList<String> bg_list = new ArrayList<String>();
                private String[] bgColors = {"#0000ff", "#ff0000", "#f0f000", "#00f0f0", "#f000f0", "#0f0f00", "000f0f", "0f000f", "fff000", "0fff00", "00fff0", "000fff"};
                private String whiteFgColor = "#ffffff";

                int baseFgColor = 10000000;
                int fgColorOffset = 2000000;

                nestedGetColorPathsThread(String genomeStr, String suf){
                    this.genomeStr = genomeStr;
                    this.suf = suf;
                }

                public void run(){
                    MainClass.nestedGetColorPathsThreadCntMcl++;
                    System.out.println("test started");

                    try{
                        int objIndex = 0;
                        int clusterOffset = 0;
                        

                        for(int i=0; i<MainClass.clustersListMcl.size(); i++){
                            for(int j=0; j<MainClass.clustersListMcl.get(i).size(); j++){

                                String curGene = MainClass.clustersListMcl.get(i).get(j);



                                boolean addEc = false;
                                if(curGene.substring(0,3).equals(this.genomeStr)){
                                    addEc = true;
                                }



                                int blackFlag = 0;

                                if(addEc){
                                    System.out.println("in genome: "+genomeStr+", inside addEc");

                                    obj_list.add(curGene);

                                    int genomeIndex = -1;
                                    String genomeString = curGene.substring(0,3);
                                    for (int iter1 = 0; iter1 < MainClass.orgsIds.length; iter1++) {
                                        if(genomeString.equals(MainClass.orgsIds[iter1]))
                                            genomeIndex = iter1;
                                    }

                               
                                    for(int iter2=0; iter2<MainClass.pr[genomeIndex].length; iter2++){
                                        if( ((String)MainClass.clustersListMcl.get(i).get(j)).equals(MainClass.pr[genomeIndex][iter2]) ){

                                            if(MainClass.ecNumsList[genomeIndex][iter2].length>0){

                                                for(int k=0; k<MainClass.ecNumsList[genomeIndex][iter2].length; k++){

                                                    String curEC = MainClass.ecNumsList[genomeIndex][iter2][k];

                                                    int ecIdxExam = -1;
                                                    for(int t1=0; t1<MainClass.ecsFinal.length; t1++){
                                                        if(curEC.equals(MainClass.ecsFinal[t1])){
                                                            ecIdxExam = t1;
                                                            break;
                                                        }
                                                    }

                                                    for(int t1=0; t1<MainClass.ec2GenesList[ecIdxExam].length; t1++){
                                                        String examGene = MainClass.ec2GenesList[ecIdxExam][t1];

                                                        for(int ii=0; ii<MainClass.clustersListMcl.size(); ii++)
                                                            if(ii!=i){
                                                                for(int jj=0; jj<MainClass.clustersListMcl.get(ii).size(); jj++){
                                                                    if( (examGene.equals(MainClass.clustersListMcl.get(ii).get(jj))) && (!examGene.equals(curGene))){
                                                                        if( (examGene.substring(0,3)).equals(curGene.substring(0,3))){
                                                                            blackFlag = 1;
                                                                            break;
                                                                        }
                                                                    }

                                                                }
                                                            }
                                                        if(blackFlag == 1)
                                                            break;
                                                    }

                                                    if(blackFlag == 1)
                                                            break;


                                                    int groupIndex = -1;
                                                    boolean groupParticipant = false;

                                                    for(int findDistId = 0; findDistId<distancesList.size(); findDistId++){

                                                        if(groupsCntForEachEc[findDistId]){
                                                            if((MainClass.ecNumsList[genomeIndex][iter2][k]).equals(MainClass.ecsFinal[ecIdsList.get(findDistId)])){

                                                                int cntCurEc = 0;
                                                                for(int chkEcId=0; chkEcId<ecIdsList.size(); chkEcId++){
                                                                    if(ecIdsList.get(findDistId).equals(ecIdsList.get(chkEcId)))
                                                                        cntCurEc++;
                                                                }
                                                                groupIndex = findDistId;
                                                                groupParticipant = true;
                                                                break;
                                                            }
                                                        }
                                                    }

                                                    if(groupParticipant){
                                                        int finalFgColor = baseFgColor - (fgColorOffset*groupIdsList.get(groupIndex));
                                                        String finalFgColorStr = "#"+Integer.toHexString(finalFgColor);

                                                        fg_list.add(finalFgColorStr);
                                                    } else{
                                                        fg_list.add(whiteFgColor);
                                                    }
                                                    
                                                }


                                            }

                                        }
                                        if(blackFlag == 1)
                                            break;

                                    }
                                   
                                if(blackFlag == 1){
                                    blackFlag = 0;
                                    bg_list.add("#0f0f00");
                                }
                                else{
                                    bg_list.add(bgColors[i]);
                                }

                                objIndex++;

                                }
                                
                            }

                            clusterOffset += MainClass.clustersListMcl.get(i).size();
                        }


                        System.out.println("obj_list with groups size: "+ obj_list.size());
                        System.out.println("fg_list with groups size: "+ fg_list.size());
                        System.out.println("bg_list with groups size: "+ bg_list.size());

                        String[] obj_arr = new String[obj_list.size()];
                        String[] fg_arr = new String[fg_list.size()];
                        String[] bg_arr = new String[bg_list.size()];

                        obj_list.toArray(obj_arr);
                        fg_list.toArray(fg_arr);
                        bg_list.toArray(bg_arr);

                        boolean valid = false;
                        while(!valid){
                            try{
                                System.out.println("inColoringMaps with groups, genome: "+genomeStr);
                                KEGGLocator  locator;
                                KEGGPortType serv;
                                locator = new KEGGLocator();
                                serv = locator.getKEGGPort();
                                String pathStr = "path:"+genomeStr+MainClass.pathId;
                                String url = serv.color_pathway_by_objects(pathStr, obj_arr, fg_arr, bg_arr);
                                System.out.println(url);
                                String imFileStr = GlobalInitializer.imgsDirForCurRun+"/MCL/images_with_groups/"+genomeStr+"Path"+suf+".png";
                                saveImage(url, imFileStr);

                                valid = true;
                                System.out.println("/inColoringMaps with groups, genome: "+genomeStr);
                            } catch(Exception e){
                                    System.out.println("Error: "+e.getMessage());
                            }
                        }
                    } catch(Exception ex){
                        Logger.getLogger(MclClustering.class.getName()).log(Level.SEVERE, null, ex);
                    }

                    MainClass.nestedGetColorPathsThreadCntMcl--;
                 }
            }

            int parallelThreadsForPathsRetrieval = 10;

            int iterations = MainClass.orgsIds.length/parallelThreadsForPathsRetrieval;
            int rem = MainClass.orgsIds.length%parallelThreadsForPathsRetrieval;

            for(int iter=0; iter<iterations; iter++){

                for(int i=0; i<parallelThreadsForPathsRetrieval; i++){

                    int idx = iter*parallelThreadsForPathsRetrieval + i;

                        try{
                            Runnable ecR = new nestedGetColorPathsThread(MainClass.orgsIds[idx], suf);
                            Thread ecThr = new Thread(ecR);
                            ecThr.start();

                        } catch(Exception e){

                                System.err.println("Error: " + e.getMessage());

                        }
                }


                while(true){
                    if(MainClass.nestedGetColorPathsThreadCntMcl > 0){
                        continue;
                    }
                    else
                        break;
                }

            }


            for(int i=0; i<rem; i++){
                int idx = iterations*parallelThreadsForPathsRetrieval + i;

                    try{
                        Runnable ecR = new nestedGetColorPathsThread(MainClass.orgsIds[idx], suf);
                        Thread ecThr = new Thread(ecR);
                        ecThr.start();
                    } catch(Exception e){

                            System.err.println("Error: " + e.getMessage());
                    }
            }

            while(true){
                if(MainClass.nestedGetColorPathsThreadCntMcl > 0){
                    continue;
                }
                else
                    break;
            }

        }


    public static void createEcsToClustersMatrix(){
        
        ecsToClustersMatrix = new double[MainClass.ecsFinal.length][MclClustering.clustersList.size()];
        for(int i=0; i<ecsToClustersMatrix.length; i++)
           for(int j=0; j<ecsToClustersMatrix[i].length; j++)
                ecsToClustersMatrix[i][j] = 0;
        

        for(int i=0; i < MainClass.ec2GenesList.length; i++){

            String curEcNumber = MainClass.ecsFinal[i];
            for(int j=0; j<MainClass.ec2GenesList[i].length; j++){
                String curGene = MainClass.ec2GenesList[i][j];

                int curCluster = findClusterByGene(curGene);

                ecsToClustersMatrix[i][curCluster] += 1;
            }
        }


        for(int i=0; i<ecsToClustersMatrix.length; i++){
            int sumOfGeneHits = 0;
            for(int j=0; j<ecsToClustersMatrix[i].length; j++){
                sumOfGeneHits += ecsToClustersMatrix[i][j];
            }

            for(int j=0; j<ecsToClustersMatrix[i].length; j++){
                ecsToClustersMatrix[i][j] /= sumOfGeneHits;
            }
        }
    }


    public static void printEcsToClustersMatrix(){


        DecimalFormat twoDForm = new DecimalFormat("#.##");
        String ecsToClustersMatrixOutputStr = "";

        ecsToClustersMatrixOutputStr += "EC Numbers  |   ";

        for(int k=0; k<MclClustering.clustersList.size(); k++)
            ecsToClustersMatrixOutputStr +=  "Cluster "+(k+1)+"   |   ";

        ecsToClustersMatrixOutputStr += "\n";

        for(int k=0; k<(MclClustering.clustersList.size()*16+13); k++){
            ecsToClustersMatrixOutputStr += "-";
        }

        for(int i=0; i<MainClass.ecsFinal.length; i++){

           
            ecsToClustersMatrixOutputStr += "\n"+MainClass.ecsFinal[i];
            
            for(int spaceCnt=0; spaceCnt<(12 - MainClass.ecsFinal[i].length()); spaceCnt++){
                ecsToClustersMatrixOutputStr+=" ";
            }
            ecsToClustersMatrixOutputStr += "|       ";
            

            
            for(int j=0; j<MclClustering.clustersList.size(); j++){

                double dist = Double.valueOf(twoDForm.format(ecsToClustersMatrix[i][j]));
                String distStr = Double.toString(dist);

                if(distStr.length() == 3)
                    ecsToClustersMatrixOutputStr += dist+"     |       ";
                else
                    ecsToClustersMatrixOutputStr += dist+"    |       ";
            }
        }

        System.out.println("\n<><><><><><><><><><><><><><><>\n");
        System.out.println(ecsToClustersMatrixOutputStr);
        System.out.println("\n<><><><><><><><><><><><><><><>\n");


        phyloDistancesOutputLog += "\n- (EC numbers - Clusters) Matrix:\n";
        phyloDistancesOutputLog += ecsToClustersMatrixOutputStr;
        phyloDistancesOutputLog += "\n";

    }


    public static void createGenomesToClustersMatrix(){
        genomesToClustersMatrix = new int[MainClass.orgsIds.length][MclClustering.clustersList.size()];

        for(int curOrgId=0; curOrgId<MainClass.orgsIds.length; curOrgId++){
            String curOrgStr = MainClass.orgsIds[curOrgId];

            for(int clustId=0; clustId<MclClustering.clustersList.size(); clustId++){
                
                boolean clustIncluded = false;
                
                for(int curGeneInClust=0; curGeneInClust<MclClustering.clustersList.get(clustId).size(); curGeneInClust++){
                    String curGeneStr = MclClustering.clustersList.get(clustId).get(curGeneInClust);
                    if(curOrgStr.equals(curGeneStr.substring(0, 3))){
                        clustIncluded = true;
                        break;
                    }
                }

                if(clustIncluded)
                    genomesToClustersMatrix[curOrgId][clustId] = 1;
                else
                    genomesToClustersMatrix[curOrgId][clustId] = 0;
            }
        }



    }


    public static void printGenomesToClustersMatrix(){

        String genomesToClustersMatrixOutputStr = "";

        genomesToClustersMatrixOutputStr += "Genomes  |   ";

        for(int k=0; k<MclClustering.clustersList.size(); k++)
            genomesToClustersMatrixOutputStr +=  "Cluster "+(k+1)+"   |   ";

        genomesToClustersMatrixOutputStr += "\n";

        for(int k=0; k<(MclClustering.clustersList.size()*16+10); k++){
            genomesToClustersMatrixOutputStr += "-";
        }


        
        for(int i=0; i<MainClass.orgsIds.length; i++){
            genomesToClustersMatrixOutputStr += "\n   "+MainClass.orgsIds[i]+"   |       ";

            
            for(int j=0; j<MclClustering.clustersList.size(); j++){
                int dist = genomesToClustersMatrix[i][j];

                genomesToClustersMatrixOutputStr += dist+"       |       ";
                
            }
        }

        System.out.println("\n<><><><><><><><><><><><><><><>\n");
        System.out.println(genomesToClustersMatrixOutputStr);
        System.out.println("\n<><><><><><><><><><><><><><><>\n");


        phyloDistancesOutputLog += "\n- (Genomes - Clusters) Matrix:\n";
        phyloDistancesOutputLog += genomesToClustersMatrixOutputStr;
        phyloDistancesOutputLog += "\n";
        
    }


    public static void saveImage(String imageUrl, String destinationFile) throws IOException {
                URL url = new URL(imageUrl);
                InputStream is = url.openStream();
                OutputStream os = new FileOutputStream(destinationFile);

                byte[] b = new byte[2048];
                int length;

                while ((length = is.read(b)) != -1) {
                        os.write(b, 0, length);
                }

                is.close();
                os.close();
        }



    private static String callPerlScriptForgrouping(String inputFile) throws IOException, InterruptedException{

                List<String> args = new ArrayList<String>();
                ProcessBuilder builder;

                System.out.println("->calling perl script...");
                
                args.add("perl");
                args.add(perlScriptPath);
                args.add("-f");
                args.add(inputFile);
                args.add("-r");
                args.add("random");
                args.add("-n");
                args.add("100");
                args.add("-m");
                args.add("split");


                builder = new ProcessBuilder(args);
                Process p = builder.start();

                InputStream is = p.getInputStream();
                InputStreamReader isr = new InputStreamReader(is);
                BufferedReader br = new BufferedReader(isr);
                String line;

                String outStr = "";

                while ((line = br.readLine()) != null) {
                        outStr += line+"\n";
                }
                
                System.out.println(outStr.substring(1,2));
                

                System.out.println("->perl script before waitFor()...");
                p.waitFor();
                System.out.println("->perl script finished!");

                String perlScriptOutput = GlobalInitializer.baseDirectory+"/src/topd_result.txt";
                
                return perlScriptOutput;
                
    }


    public static Double[][] sort2Darray(Double[][] twoDarray){

            Arrays.sort(twoDarray, new Comparator<Double[]>() {
            @Override
            public int compare(final Double[] entry1, final Double[] entry2) {
                final Double distance1 = entry1[1];
                final Double distance2 = entry2[1];
                return distance1.compareTo(distance2);
            }
            });

            return twoDarray;
    }
    

    public static String getGeneSequenceFromFasta(String gene) throws FileNotFoundException, IOException{

            String sequence = "";
            Integer curCluster = findClusterByGene(gene);
            int geneStrLength = gene.length();

            FileInputStream fstream = new FileInputStream(GlobalInitializer.imgsDirForCurRun+"/MCL/fasta/cluster_"+curCluster.toString()+"_FastaFile.txt");
            DataInputStream in = new DataInputStream(fstream);
            BufferedReader br = new BufferedReader(new InputStreamReader(in));

            String strLine;

            while ((strLine = br.readLine()) != null){

                if(strLine.substring(0,1).equals(">")){

                    if(gene.equals(strLine.substring(1, geneStrLength+1))){

                        String subseq = "";

                        while(true){
                            subseq = br.readLine();
                            if((subseq != null) && (!subseq.substring(0,1).equals(">"))){
                                sequence += subseq;
                            }
                            else
                                break;
                        }
                        break;
                    }
                }
            }

            return sequence;

    }


    public static boolean checkMatricesEquality(int[] mat1, int[] mat2){

        boolean equal = true;
        if(mat1.length == mat2.length){
            for(int i=0; i<mat1.length; i++){
                if(mat1[i] != mat2[i]){
                    equal = false;
                    break;
                }
            }
        } else{
            equal = false;
        }

        return equal;
    }

    public static String findEcByGene(String gene){

        int ecId = -1;

        for(int i=0; i<MainClass.ec2GenesList.length; i++){
            for(int j=0; j<MainClass.ec2GenesList[i].length; j++){
                if((MainClass.ec2GenesList[i][j]).equals(gene)){
                    ecId = i;
                    break;
                }
            }
        }

        return MainClass.ecsFinal[ecId];
    }


    public static int findClusterByGene(String gene){

        int clustId = -1;

        for(int i=0; i<MainClass.clustersListMcl.size(); i++){
            for(int j=0; j<MainClass.clustersListMcl.get(i).size(); j++){
                if(gene.equals(MainClass.clustersListMcl.get(i).get(j))){
                    clustId = i;
                    break;
                }
            }
        }

        
        return clustId;
    }

    
}




