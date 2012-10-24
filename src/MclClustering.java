import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.net.URL;
import java.rmi.RemoteException;
import java.text.DecimalFormat;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.xml.rpc.ServiceException;
import keggapi.KEGGLocator;
import keggapi.KEGGPortType;
import net.sf.javaml.distance.JaccardIndexDistance;




class MclClustering{

	
	private static String  mclDirStr;
        public static List<List<String>> clustersList;
	private static int[][] genesByGenomeInClusters;

        private static double[][] clustersDistValidMatMcl;
        private static double[][] clustersDistValidMatMcl_Max;
        private static double[][] clustersDistValidMatMcl_Min;
        private static double[][] clustersDistValidMatMcl_Std;

        private static double[][] clustersHomogenValidMatMcl;
        private static double[][] clustersHomogenValidMatMcl_Max;
        private static double[][] clustersHomogenValidMatMcl_Min;
        private static double[][] clustersHomogenValidMatMcl_Std;

        private DecimalFormat twoDForm;
	// *** PARAMETERS ***
	//define distance-similarity metric
	private static JaccardIndexDistance jacDist = new JaccardIndexDistance();
	
	
	
	MclClustering(){
		
		System.out.println("\n\n");
                clustersList = new ArrayList<List<String>>();
                
                twoDForm = new DecimalFormat("#.###");

	}
	
	public static void createAbcFile(int orgsNum, int totalGenesNumber){
		
		//create MCL data directory
		mclDirStr = GlobalInitializer.genDataDirStr+"/MCL_Files";
		boolean success = (new File(mclDirStr)).mkdir();
		if (success) {
			System.out.println("Directory: "+mclDirStr+" was created succesfully.");
		}
		
		
		try{
                    int cnt = 0;
                    float tmplength;
			FileWriter abcFstream = new FileWriter(GlobalInitializer.genDataDirStr+"/MCL_Files/graph");
			BufferedWriter abcOut = new BufferedWriter(abcFstream);
			
			for (int io=0; io<totalGenesNumber; io++) {
				//DenseInstance geneA = new DenseInstance(MainClass.P[io]);
				String geneAStr = MainClass.getGeneNameByIdxInP(io, orgsNum);
				for (int ii=io; ii<totalGenesNumber; ii++) {
                                        cnt++;
					//DenseInstance geneB = new DenseInstance(MainClass.P[ii]);
					
					// *** Using distance-metric PARAMETER
					
                                        //double edgeLength = jacDist.measure(geneA, geneB);
					String geneBStr = MainClass.getGeneNameByIdxInP(ii, orgsNum);
                                        float edgeLength = calcJaccardDist(MainClass.P[io], MainClass.P[ii]);
                                        abcOut.write(geneAStr+" "+geneBStr+" "+edgeLength+"\n");
				}
			}
			
			abcOut.close();
		}catch (Exception e){
			System.err.println("Error: " + e.getMessage());
		}
		
	}
	
	public void mcxload(double inflArg){

            //mcxload
            try {
                List<String> args = new ArrayList<String>();
                ProcessBuilder builder;

                System.out.println("mcxload");
                String infl = Double.toString(inflArg);
                String suf = Integer.toString((int)(10*inflArg));
                args.add("mcxload");
                args.add("-abc");
                args.add(GlobalInitializer.genDataDirStr+"/MCL_Files/graph");
                args.add("--stream-mirror");
                args.add("-write-tab");
                args.add("data.tab");
                args.add("-o");
                
                args.add(GlobalInitializer.genDataDirStr+"/MCL_Files/data.mci");
                
                builder = new ProcessBuilder(args);
                Process p = builder.start();
                p.waitFor();
                System.out.println("/mcxload");
            } catch (IOException e) {
                    System.err.println("Error: " + e.getMessage());
            } catch(InterruptedException ie){
                    System.err.println("Error: " + ie.getMessage());
            }

        }

        public void mcl(double inflArg){

            //mcl
            try {
                List<String> args = new ArrayList<String>();
                ProcessBuilder builder;

                System.out.println("mcl");
                String infl = Double.toString(inflArg);
                String suf = Integer.toString((int)(10*inflArg));
                args.add("mcl");
                args.add(GlobalInitializer.genDataDirStr+"/MCL_Files/data.mci");
                args.add("-I");
                args.add(infl);
                args.add("-o");
                args.add(GlobalInitializer.genDataDirStr+"/MCL_Files/out.data.mci.I"+suf);
                
                builder = new ProcessBuilder(args);
                Process p = builder.start();
                p.waitFor();
                System.out.println("/mcl");
            } catch (IOException e) {
                    System.err.println("Error: " + e.getMessage());
            } catch(InterruptedException ie){
                    System.err.println("Error: " + ie.getMessage());
            }

        }

        public void mcxdump(double inflArg){
            //mcxdump
            try {
                List<String> args = new ArrayList<String>();
                ProcessBuilder builder;

                System.out.println("mcxdump");
                String infl = Double.toString(inflArg);
                String suf = Integer.toString((int)(10*inflArg));
                args.add("mcxdump");
                args.add("-icl");
                args.add(GlobalInitializer.genDataDirStr+"/MCL_Files/out.data.mci.I"+suf);
                args.add("-tabr");
                args.add("data.tab");
                args.add("-o");
                args.add(GlobalInitializer.genDataDirStr+"/MCL_Files/dump.data.mci.I"+suf);
                //args.add("-append-log");
                //args.add("y");
                builder = new ProcessBuilder(args);
                Process p = builder.start();
                p.waitFor();
                System.out.println("/mcxdump");
            } catch (IOException e) {
                    System.err.println("Error: " + e.getMessage());
            } catch(InterruptedException ie){
                    System.err.println("Error: " + ie.getMessage());
            }
        }

	public void cluster(double inflArg){

            mcxload(inflArg);
            mcl(inflArg);
            mcxdump(inflArg);
	}
	
	public void printOutput(double inflArg){
		
		System.out.println("\n\n*** MCL Output, Inflation: "+inflArg+" ***\n\n");
		try{
			String suf = Integer.toString((int)(10*inflArg));
			FileInputStream fstreamMcl = new FileInputStream(GlobalInitializer.genDataDirStr+"/MCL_Files/dump.data.mci.I"+suf);
			DataInputStream inMcl = new DataInputStream(fstreamMcl);
			BufferedReader brMcl = new BufferedReader(new InputStreamReader(inMcl));
			String strLine;
			int printFlag = 0;
			while ((strLine = brMcl.readLine()) != null){
				if (strLine.equals("(mclruninfo")) {
					printFlag = 1;
				}
				if (printFlag == 1)
					System.out.println (strLine);
			}
			inMcl.close();
		} catch (Exception e){
			System.err.println("Error: " + e.getMessage());
		}
		System.out.println("\n\n*** End of MCL Output, Inflation: "+inflArg+" ***\n\n");
		
	}

        public List<List<String>> getClusters(double inflArg){
		try{
                    System.out.println("in getClusters");
			String suf = Integer.toString((int)(10*inflArg));
			FileInputStream fstreamMcl = new FileInputStream(GlobalInitializer.genDataDirStr+"/MCL_Files/dump.data.mci.I"+suf);
			DataInputStream inMcl = new DataInputStream(fstreamMcl);
			BufferedReader brMcl = new BufferedReader(new InputStreamReader(inMcl));

			int i = 0;
			int j = 0;
			int numOfClusters = 0;
			
			int[] numOfElsInClusters = new int[MainClass.orgsIds.length];
			for (i=0; i<numOfElsInClusters.length; i++) {
				numOfElsInClusters[i] = 0;
			}
                        
                        
			while(true){
                            String line = brMcl.readLine();
                            if(line != null){
                                numOfClusters++;
                                List<String> cluster = new ArrayList<String>();
                                
                                StringTokenizer stokenz = new StringTokenizer(line, "\t");
                                
                                while (stokenz.hasMoreTokens()) {
                                    String tmpGenStr = stokenz.nextToken();
                                    int geneIdx = MainClass.getIdxByNameInP(tmpGenStr);
                                    //System.out.print(tmpGenStr+" ");
                                    //for(int k=0; k<MainClass.P[geneIdx].length; k++)
                                        //System.out.print(MainClass.P[geneIdx][k]+" ");
                                    cluster.add(tmpGenStr);
                                   
                                }
                                
                                clustersList.add(cluster);
                                //System.out.println("Num of elements in cluster "+numOfClusters+": "+cluster.size());
                               
                            }
                            else
                                break;
                        }

                        brMcl.close();
			inMcl.close();
                        fstreamMcl.close();

                        //count genes per genome in each cluster
                        genesByGenomeInClusters = new int[clustersList.size()][MainClass.orgsIds.length];
                        for(i=0; i<clustersList.size(); i++)
                            for(j=0; j<MainClass.orgsIds.length; j++)
                                genesByGenomeInClusters[i][j] = 0;
                        
                        for(j=0; j<clustersList.size(); j++){

                            ArrayList curList = (ArrayList)clustersList.get(j);
                            
                            String queryGenome = null;

                            for(int k=clustersList.get(j).size()-1; k>=0; k--){
                                queryGenome = ((String)curList.get(k)).substring(0,3);
                                genesByGenomeInClusters[j][findGenomeIdByGene(queryGenome)]++;
                            }

                        }

                        //Write to output log
                        MainClass.outputLogStrMcl += "> Pathway: map"+MainClass.pathId+".\n\n";
                        MainClass.outputLogStrMcl += "> Genomes examined: ";
                        for(i=0; i<MainClass.orgsIds.length-1; i++)
                            MainClass.outputLogStrMcl += MainClass.orgsIds[i]+", ";
                        MainClass.outputLogStrMcl += MainClass.orgsIds[i]+".\n\n";

                        MainClass.outputLogStrMcl += "> Total number of genes: ";
                        MainClass.outputLogStrMcl += MainClass.totalGenesNumber+".\n\n";

                        MainClass.outputLogStrMcl += "> Number of clusters: ";
                        MainClass.outputLogStrMcl += clustersList.size()+".\n";

                        MainClass.outputLogStrMcl += "\n";
                        
                        MainClass.outputLogStrMcl += "> Number of genes in each cluster:\n";
                        for( i=0; i<clustersList.size(); i++){
                            MainClass.outputLogStrMcl += "- Cluster "+(i+1)+": ";
                            MainClass.outputLogStrMcl += clustersList.get(i).size()+" [";
                            int k;
                            for(k=0; k<MainClass.orgsIds.length-1; k++)
                                MainClass.outputLogStrMcl += MainClass.orgsIds[k]+"("+genesByGenomeInClusters[i][k]+"), ";
                            MainClass.outputLogStrMcl += MainClass.orgsIds[k]+"("+genesByGenomeInClusters[i][k]+")].\n";
                        }
                        MainClass.outputLogStrMcl += "\n";

                        

                    System.out.println("end of getClusters");
		} catch (Exception e){
                        System.err.println("Error: " + e.getMessage());
		}

                return clustersList;

	}

        public void getUniqueGenesByGenome(){

            

            int[] foundUniqeGenesInCluster = new int[clustersList.size()];
            String[] genomeInClusterWithUniqeGenes = new String[clustersList.size()];


            for(int i=0; i<foundUniqeGenesInCluster.length; i++){
                foundUniqeGenesInCluster[i] = 1;
                genomeInClusterWithUniqeGenes[i] = null;
            }


            for(int j=0; j<clustersList.size(); j++){

                ArrayList curList = (ArrayList)clustersList.get(j);
                String baseGenome = ((String)curList.get(0)).substring(0,3);
                genomeInClusterWithUniqeGenes[j] = baseGenome;
                String queryGenome = null;

                for(int k=clustersList.get(j).size()-1; k>=0; k--){
                    queryGenome = ((String)curList.get(k)).substring(0,3);
                    if(baseGenome.equals(queryGenome)){
                        continue;
                    }
                    else{
                        foundUniqeGenesInCluster[j] = 0;
                        break;
                    }
                }

            }

        }
        


        public void validateClusters(double inflArg){
		try{

                    MainClass.outputLogStrMclAppend += "> Clustering validation:\n\n";
                    
                    clustersDistValidMatMcl = new double[clustersList.size()][clustersList.size()];
                    clustersDistValidMatMcl_Max = new double[clustersList.size()][clustersList.size()];
                    clustersDistValidMatMcl_Min = new double[clustersList.size()][clustersList.size()];
                    clustersDistValidMatMcl_Std = new double[clustersList.size()][clustersList.size()];



                    clustersHomogenValidMatMcl = new double[clustersList.size()][clustersList.size()];
                    

                    for(int i=0; i<clustersList.size(); i++){
                        for(int j=0; j<clustersList.size(); j++){
                            clustersDistValidMatMcl_Max[i][j] = 0.0;
                            clustersDistValidMatMcl_Min[i][j] = 10000.0;
                            clustersDistValidMatMcl_Std[i][j] = 0.0;

                        }
                    }

                    //internal clustering validation
                    for(int iter=0; iter<clustersList.size(); iter++){
                        float avgDistInCluster = 0;
                        int homologiesTotalNumber = 0;
                        double homologiesMean = 0;

                        ArrayList curList = (ArrayList)clustersList.get(iter);
                        int curListSize = curList.size();

                        for(int i=0; i<curListSize; i++){

                            int iGeneIdx = MainClass.getIdxByNameInP((String)curList.get(i));
                            
                            for(int j=0; j<curListSize; j++){
                                int jGeneIdx = MainClass.getIdxByNameInP((String)curList.get(j));
				
				// *** Using distance-metric PARAMETER
				float edgeLength = calcJaccardDist(MainClass.P[iGeneIdx], MainClass.P[jGeneIdx]);
                                avgDistInCluster += edgeLength;

                                //update max
                                if(edgeLength > clustersDistValidMatMcl_Max[iter][iter])
                                    clustersDistValidMatMcl_Max[iter][iter] = edgeLength;

                                //update min
                                if(edgeLength<clustersDistValidMatMcl_Min[iter][iter])
                                    clustersDistValidMatMcl_Min[iter][iter] = edgeLength;


                                if(MainClass.C[iGeneIdx][jGeneIdx] <= BlastXMLParser.expValueThreshold){
                                    homologiesTotalNumber++;
                                }
                            }
			}
                        avgDistInCluster = avgDistInCluster/(curListSize*curListSize);
                        homologiesMean = (double)homologiesTotalNumber/(curListSize*curListSize);


                        //find std
                        for(int i=0; i<curListSize; i++){

                            int iGeneIdx = MainClass.getIdxByNameInP((String)curList.get(i));
                            
                            for(int j=0; j<curListSize; j++){
                                int jGeneIdx = MainClass.getIdxByNameInP((String)curList.get(j));
				
				// *** Using distance-metric PARAMETER
				float edgeLength = calcJaccardDist(MainClass.P[iGeneIdx], MainClass.P[jGeneIdx]);

                                clustersDistValidMatMcl_Std[iter][iter] += Math.pow(edgeLength-avgDistInCluster,2);

                            }
                        }

                        clustersDistValidMatMcl_Std[iter][iter] = Math.sqrt(clustersDistValidMatMcl_Std[iter][iter]/(curListSize*curListSize - 1));


                        clustersDistValidMatMcl[iter][iter] = avgDistInCluster;
                        clustersHomogenValidMatMcl[iter][iter] = homologiesMean;

                    }


                    //external clustering validation
                    if(clustersList.size()>1){
                        
                        for(int ii=0; ii<clustersList.size(); ii++){
                            for(int ij=0; ij<clustersList.size(); ij++){
                                int homologiesTotalNumber = 0;
                                double homologiesMean = 0.0;

                                if(ii!=ij){
                                    ArrayList curList0 = (ArrayList)clustersList.get(ii);
                                    int curListSize0 = curList0.size();
                                    ArrayList curList1 = (ArrayList)clustersList.get(ij);
                                    int curListSize1 = curList1.size();
                                    float avgDistInCluster = 0;
                                    for(int i=0; i<curListSize0; i++){

                                        int iGeneIdx = MainClass.getIdxByNameInP((String)curList0.get(i));
                                        
                                        for(int j=0; j<curListSize1; j++){
                                            int jGeneIdx = MainClass.getIdxByNameInP((String)curList1.get(j));
                                            
                                            // *** Using distance-metric PARAMETER
                                            float edgeLength = calcJaccardDist(MainClass.P[iGeneIdx], MainClass.P[jGeneIdx]);
                                            avgDistInCluster += edgeLength;


                                            //update max
                                            if(edgeLength>clustersDistValidMatMcl_Max[ii][ij])
                                                clustersDistValidMatMcl_Max[ii][ij] = edgeLength;

                                            //update min
                                            if(edgeLength<clustersDistValidMatMcl_Min[ii][ij])
                                                clustersDistValidMatMcl_Min[ii][ij] = edgeLength;


                                            if(MainClass.C[iGeneIdx][jGeneIdx] <= BlastXMLParser.expValueThreshold){
                                                homologiesTotalNumber++;
                                            }
                                        }
                                    }
                                    avgDistInCluster = avgDistInCluster/(curListSize0*curListSize1);
                                    homologiesMean = (double)homologiesTotalNumber/(curListSize0*curListSize1);

                                    //find std
                                    for(int i=0; i<curListSize0; i++){

                                        int iGeneIdx = MainClass.getIdxByNameInP((String)curList0.get(i));
                                        
                                        for(int j=0; j<curListSize1; j++){
                                            int jGeneIdx = MainClass.getIdxByNameInP((String)curList1.get(j));
                                            
                                            // *** Using distance-metric PARAMETER
                                            float edgeLength = calcJaccardDist(MainClass.P[iGeneIdx], MainClass.P[jGeneIdx]);

                                            clustersDistValidMatMcl_Std[ii][ij] += Math.pow(edgeLength-avgDistInCluster,2);

                                        }
                                    }

                                    clustersDistValidMatMcl_Std[ii][ij] = Math.sqrt(clustersDistValidMatMcl_Std[ii][ij]/(curListSize0*curListSize1 - 1));



                                    clustersDistValidMatMcl[ii][ij] = avgDistInCluster;
                                    clustersHomogenValidMatMcl[ii][ij] = homologiesMean;

                                }
                            }
                        }
                    }


                   MainClass.outputLogStrMclAppend += "- Average ~Similarity~ between clusters:\n\n";
                   MainClass.outputLogStrMclAppend += "Clusters |    ";

                   for(int k=0; k<clustersList.size(); k++)
                       MainClass.outputLogStrMclAppend += (k+1)+"    |    ";
                   MainClass.outputLogStrMclAppend += "\n";
                   for(int k=0; k<(clustersList.size()*10+10); k++){
                       MainClass.outputLogStrMclAppend += "-";
                   }

                   for(int i=0; i<clustersList.size(); i++){
                       MainClass.outputLogStrMclAppend += "\n    "+(i+1)+"    |  ";
                       for(int j=0; j<clustersList.size(); j++){
                           double dist = Double.valueOf(twoDForm.format(clustersDistValidMatMcl[i][j]));
                           //count number of decimal places
                           String tmp = Double.toString(dist);
                           String[] res = tmp.split("\\.");
                           MainClass.outputLogStrMclAppend += dist;
                           if(res[1].length() == 4)
                               MainClass.outputLogStrMclAppend +=" |  ";
                           else if(res[1].length() == 3)
                                   MainClass.outputLogStrMclAppend +="  |  ";
                           else if(res[1].length() == 2)
                                   MainClass.outputLogStrMclAppend +="   |  ";
                           else
                               MainClass.outputLogStrMclAppend +="    |  ";
                       }
                   }
                   MainClass.outputLogStrMclAppend += "\n\n";


                   MainClass.outputLogStrMclAppend += "- Maximum ~Similarity~ between clusters:\n\n";
                   MainClass.outputLogStrMclAppend += "Clusters |    ";

                   for(int k=0; k<clustersList.size(); k++)
                       MainClass.outputLogStrMclAppend += (k+1)+"    |    ";
                   MainClass.outputLogStrMclAppend += "\n";
                   for(int k=0; k<(clustersList.size()*10+10); k++){
                       MainClass.outputLogStrMclAppend += "-";
                   }

                   for(int i=0; i<clustersList.size(); i++){
                       MainClass.outputLogStrMclAppend += "\n    "+(i+1)+"    |  ";
                       for(int j=0; j<clustersList.size(); j++){
                           double dist = Double.valueOf(twoDForm.format(clustersDistValidMatMcl_Max[i][j]));

                           //count number of decimal places
                           String tmp = Double.toString(dist);
                           String[] res = tmp.split("\\.");
                           MainClass.outputLogStrMclAppend += dist;
                           if(res[1].length() == 4)
                               MainClass.outputLogStrMclAppend +=" |  ";
                           else if(res[1].length() == 3)
                                   MainClass.outputLogStrMclAppend +="  |  ";
                           else if(res[1].length() == 2)
                                   MainClass.outputLogStrMclAppend +="   |  ";
                           else
                               MainClass.outputLogStrMclAppend +="    |  ";
                       }
                   }
                   MainClass.outputLogStrMclAppend += "\n\n";


                   MainClass.outputLogStrMclAppend += "- Minimum ~Similarity~ between clusters:\n\n";
                   MainClass.outputLogStrMclAppend += "Clusters |    ";

                   for(int k=0; k<clustersList.size(); k++)
                       MainClass.outputLogStrMclAppend += (k+1)+"    |    ";
                   MainClass.outputLogStrMclAppend += "\n";
                   for(int k=0; k<(clustersList.size()*10+10); k++){
                       MainClass.outputLogStrMclAppend += "-";
                   }

                   for(int i=0; i<clustersList.size(); i++){
                       MainClass.outputLogStrMclAppend += "\n    "+(i+1)+"    |  ";
                       for(int j=0; j<clustersList.size(); j++){
                           double dist = Double.valueOf(twoDForm.format(clustersDistValidMatMcl_Min[i][j]));

                           //count number of decimal places
                           String tmp = Double.toString(dist);
                           String[] res = tmp.split("\\.");
                           MainClass.outputLogStrMclAppend += dist;
                           if(res[1].length() == 4)
                               MainClass.outputLogStrMclAppend +=" |  ";
                           else if(res[1].length() == 3)
                                   MainClass.outputLogStrMclAppend +="  |  ";
                           else if(res[1].length() == 2)
                                   MainClass.outputLogStrMclAppend +="   |  ";
                           else
                               MainClass.outputLogStrMclAppend +="    |  ";
                       }
                   }
                   MainClass.outputLogStrMclAppend += "\n\n";




                   MainClass.outputLogStrMclAppend += "- Standard deviation of ~Similarity~ between clusters:\n\n";
                   MainClass.outputLogStrMclAppend += "Clusters |    ";

                   for(int k=0; k<clustersList.size(); k++)
                       MainClass.outputLogStrMclAppend += (k+1)+"    |    ";
                   MainClass.outputLogStrMclAppend += "\n";
                   for(int k=0; k<(clustersList.size()*10+10); k++){
                       MainClass.outputLogStrMclAppend += "-";
                   }

                   for(int i=0; i<clustersList.size(); i++){
                       MainClass.outputLogStrMclAppend += "\n    "+(i+1)+"    |  ";
                       for(int j=0; j<clustersList.size(); j++){
                           double dist = Double.valueOf(twoDForm.format(clustersDistValidMatMcl_Std[i][j]));

                           //count number of decimal places
                           String tmp = Double.toString(dist);
                           String[] res = tmp.split("\\.");
                           MainClass.outputLogStrMclAppend += dist;
                           if(res[1].length() == 4)
                               MainClass.outputLogStrMclAppend +=" |  ";
                           else if(res[1].length() == 3)
                                   MainClass.outputLogStrMclAppend +="  |  ";
                           else if(res[1].length() == 2)
                                   MainClass.outputLogStrMclAppend +="   |  ";
                           else
                               MainClass.outputLogStrMclAppend +="    |  ";
                       }
                   }
                   MainClass.outputLogStrMclAppend += "\n\n**************************************\n\n";






                   MainClass.outputLogStrMclAppend += "- ~Homologies/Gene~ between clusters:\n\n";
                   MainClass.outputLogStrMclAppend += "Clusters |    ";

                   for(int k=0; k<clustersList.size(); k++)
                       MainClass.outputLogStrMclAppend += (k+1)+"    |    ";
                   MainClass.outputLogStrMclAppend += "\n";
                   for(int k=0; k<(clustersList.size()*10+10); k++){
                       MainClass.outputLogStrMclAppend += "-";
                   }

                   for(int i=0; i<clustersList.size(); i++){
                       MainClass.outputLogStrMclAppend += "\n    "+(i+1)+"    |  ";
                       for(int j=0; j<clustersList.size(); j++){
                           double homog = Double.valueOf(twoDForm.format(clustersHomogenValidMatMcl[i][j]));
                           MainClass.outputLogStrMclAppend += homog;

                           //count number of decimal places
                           String tmp = Double.toString(homog);
                           String[] res = tmp.split("\\.");

                           if(res[1].length() == 4)
                                   MainClass.outputLogStrMclAppend +=" |  ";
                           else if(res[1].length() == 3)
                                   MainClass.outputLogStrMclAppend +="  |  ";
                           else if(res[1].length() == 2)
                               MainClass.outputLogStrMclAppend +="   |  ";
                           else
                               MainClass.outputLogStrMclAppend +="    |  ";
                       }
                   }

                   MainClass.outputLogStrMclAppend += "\n\n";

                } catch(Exception e){
                    System.err.println("Error: " + e.getMessage());
		}

	}



        public void getColoredPathwaysByGenome(double inflArg){

            boolean success = (new File((GlobalInitializer.imgsDirForCurRun)+"/MCL/images").mkdir());
            if (success) {
                System.out.println("Directory: " + (GlobalInitializer.imgsDirForCurRun)+"/MCL/images" + " was created succesfully.\n");
            }

            String suf = "I"+Integer.toString((int)(10*inflArg));
            String algo = "MCL";
            MainClass.colorPathsThreadCnt = 0;
            try{

                for(int iter=0; iter<MainClass.orgsIds.length; iter++){
                    
                    String genomeStr = MainClass.orgsIds[iter];
                    Runnable r = new getColoredPathwaysThread(algo, genomeStr, suf, iter, (ArrayList<List<String>>) clustersList);
                    Thread thr = new Thread(r);
                    thr.start();
                    
                }

                while(true){
                    if(MainClass.colorPathsThreadCnt > 0){
                        
                        continue;
                    }
                    else
                        break;
                }
                    
            } catch(Exception e){
                    System.err.println("Error: " + e.getMessage());
            }
        }


        public static void printUniqueGenesLists(){

        int[] foundUniqeGenesInCluster = new int[MainClass.clustersListMcl.size()];
        String[] genomeInClusterWithUniqeGenes = new String[MainClass.clustersListMcl.size()];


            for(int i=0; i<foundUniqeGenesInCluster.length; i++){
                foundUniqeGenesInCluster[i] = 1;
                genomeInClusterWithUniqeGenes[i] = null;
            }


            for(int j=0; j<MainClass.clustersListMcl.size(); j++){

                ArrayList curList = (ArrayList)MainClass.clustersListMcl.get(j);
                String baseGenome = ((String)curList.get(0)).substring(0,3);
                genomeInClusterWithUniqeGenes[j] = baseGenome;
                String queryGenome = null;

                for(int k=MainClass.clustersListMcl.get(j).size()-1; k>=0; k--){
                    queryGenome = ((String)curList.get(k)).substring(0,3);
                    if(baseGenome.equals(queryGenome)){
                        continue;
                    }
                    else{
                        foundUniqeGenesInCluster[j] = 0;
                        break;
                    }
                }

            }

            try{
                int foundAtLeastOne = 0;
                MainClass.outputLogStrMcl += "> List of clusters with genes belonging only to a single genome:";
                for(int i=0; i<MainClass.clustersListMcl.size(); i++){
                    if(foundUniqeGenesInCluster[i] == 1){
                        foundAtLeastOne = 1;
                        MainClass.outputLogStrMcl += "\n- Cluster "+(i+1)+":\n";
                        MainClass.outputLogStrMcl += "# Genome: "+genomeInClusterWithUniqeGenes[i]+"\n";
                        MainClass.outputLogStrMcl += "# Genes: \n";



                        for(int j=0; j<MainClass.clustersListMcl.get(i).size(); j++){
                            MainClass.outputLogStrMcl += MainClass.clustersListMcl.get(i).get(j)+" (";
                            
                            //find gene in pr[][] matrix
                            for(int iter1=0; iter1<MainClass.pr.length; iter1++){
                            for(int iter2=0; iter2<MainClass.pr[iter1].length; iter2++){
                                
                                    if( ((String)MainClass.clustersListMcl.get(i).get(j)).equals(MainClass.pr[iter1][iter2]) ){
                                        for(int k=0; k<MainClass.ecNumsList[iter1][iter2].length; k++){
                                            MainClass.outputLogStrMcl += MainClass.ecNumsList[iter1][iter2][k];
                                            if(k<MainClass.ecNumsList[iter1][iter2].length-1)
                                                MainClass.outputLogStrMcl += " ";
                                        }
                                    }
                                }
                            }

                            MainClass.outputLogStrMcl += ")\n";
                        }


                        }
                    }

                MainClass.outputLogStrMcl += "\n";

                if(foundAtLeastOne == 0){
                    MainClass.outputLogStrMcl += "empty";
                    MainClass.outputLogStrMcl += "\n\n";
                }


            } catch(Exception e){

                System.err.println("Error: " + e.getMessage());
            }
        }


        public static void getColoredMapsByCluster() throws RemoteException, IOException, ServiceException{

            try{
                boolean success = (new File((GlobalInitializer.imgsDirForCurRun)+"/MCL/images").mkdir());
                if (success) {
                    System.out.println("Directory: " + (GlobalInitializer.imgsDirForCurRun)+"/MCL/images" + " was created succesfully.\n");
                }

            } catch (Exception e){
                System.err.println("Error: " + e.getMessage());
            }


          
            String fgColor = "#ffffff";

            String[] bgColors = {"#0000ff", "#ff0000", "#f0f000", "#00f0f0", "#f000f0", "#0f0f00", "000f0f", "0f000f", "fff000", "0fff00", "00fff0", "000fff"};
            
            
                for(int i=0; i<MainClass.clustersListMcl.size(); i++){
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
                                                    fg_list.add(fgColor);
                                                    objIndex++;
                                                }
                                            }
                                        } else {

                                                String[] res = serv.get_ko_by_gene(MainClass.clustersListMcl.get(i).get(j));
                                                for (int k = 0; k < res.length; k++) {
                                                    obj_list.add(res[k]);
                                                    fg_list.add(fgColor);
                                                    bg_list.add(bgColors[i]);
                                                    objIndex++;
                                                }

                                        }
                                    }
                                }
                            
                        }

                        System.out.println("obj_list size: "+ obj_list.size());
                        System.out.println("fg_list size: "+ fg_list.size());
                        System.out.println("bg_list size: "+ bg_list.size());

                        

                        String[] obj_arr = new String[obj_list.size()];
                        String[] fg_arr = new String[fg_list.size()];
                        String[] bg_arr = new String[bg_list.size()];

                        obj_list.toArray(obj_arr);
                        fg_list.toArray(fg_arr);
                        bg_list.toArray(bg_arr);

                        System.out.println("Path String: "+"map"+MainClass.pathId);

                        System.out.println("obj_arr:");
                        for(int itest=0; itest<obj_arr.length; itest++){
                            System.out.print(obj_arr[itest]+" ");
                        }

                        System.out.println("\nfg_arr:");
                        for(int itest=0; itest<fg_arr.length; itest++){
                            System.out.print(fg_arr[itest]+" ");
                        }

                        System.out.println("\nbg_arr:");
                        for(int itest=0; itest<bg_arr.length; itest++){
                            System.out.print(bg_arr[itest]+" ");
                        }
                        

                        boolean valid = false;
                        while(!valid){
                            try{
                                System.out.println("downloading map img, cluster: "+(i+1));
                                String pathStr = "map"+MainClass.pathId;
                                String url;
                                url = serv.color_pathway_by_objects(pathStr, obj_arr, fg_arr, bg_arr);
                                System.out.println(url);
                                String imFileStr = GlobalInitializer.imgsDirForCurRun + "/MCL/images/cluster_" + (i + 1) + "_map.png";
                                saveImage(url, imFileStr);

                                valid = true;
                            } catch (Exception e){
                                System.err.println("Error: "+e.getMessage());
                            }
                        }
                        System.out.println("/downloading map img, cluster: "+(i+1));
                    } catch (IOException ex) {
                        Logger.getLogger(MclClustering.class.getName()).log(Level.SEVERE, null, ex);
                    }
            }
                  
        }


        public void createFastaForClusters() throws Exception{

            System.out.println("in createFastaForClusters");
            boolean success = (new File((GlobalInitializer.imgsDirForCurRun)+"/MCL/fasta").mkdir());
            if (success) {
                System.out.println("Directory: " + (GlobalInitializer.imgsDirForCurRun)+"/MCL/fasta" + " was created succesfully.\n");
            }

            FastaCreator clFasta = new FastaCreator();

            for(int i=0; i<clustersList.size(); i++){
                String[] clusterString = new String[clustersList.get(i).size()];
                for(int j=0; j<clustersList.get(i).size(); j++){
                    clusterString[j] = (String)clustersList.get(i).get(j);
                }
                String fastaString = clFasta.getFastaString(clusterString);
                String fastaId = GlobalInitializer.imgsDirForCurRun+"/MCL/fasta/cluster_"+i+"_FastaFile.txt";
                clFasta.createFastaFile(fastaId, fastaString);
            }

            System.out.println("end of createFastaForClusters");
        }

        
        public static void getColoredPathwaysFunction(String suf) throws Exception{

        

            boolean success = (new File((GlobalInitializer.imgsDirForCurRun)+"/MCL/images").mkdir());
            if (success) {
                System.out.println("Directory: " + (GlobalInitializer.imgsDirForCurRun)+"/MCL/images" + " was created succesfully.\n");
            }


             class nestedGetColorPathsThread implements Runnable{

                String genomeStr;
                String suf;

                private String fgColor = "#ffffff";
                ArrayList<String> obj_list = new ArrayList<String>();
                ArrayList<String> fg_list = new ArrayList<String>();
                ArrayList<String> bg_list = new ArrayList<String>();
                private String[] bgColors = {"#0000ff", "#ff0000", "#f0f000", "#00f0f0", "#0f0f00", "000f0f"};

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
                                fg_list.add(fgColor);
                                objIndex++;

                                }
                                
                            }

                            clusterOffset += MainClass.clustersListMcl.get(i).size();
                        }


                        System.out.println("obj_list size: "+ obj_list.size());
                        System.out.println("fg_list size: "+ fg_list.size());
                        System.out.println("bg_list size: "+ bg_list.size());

                        String[] obj_arr = new String[obj_list.size()];
                        String[] fg_arr = new String[fg_list.size()];
                        String[] bg_arr = new String[bg_list.size()];

                        obj_list.toArray(obj_arr);
                        fg_list.toArray(fg_arr);
                        bg_list.toArray(bg_arr);

                        boolean valid = false;
                        while(!valid){
                            try{
                                System.out.println("inColoringMaps, genome: "+genomeStr);
                                KEGGLocator  locator;
                                KEGGPortType serv;
                                locator = new KEGGLocator();
                                serv = locator.getKEGGPort();
                                String pathStr = "path:"+genomeStr+MainClass.pathId;
                                String url = serv.color_pathway_by_objects(pathStr, obj_arr, fg_arr, bg_arr);
                                System.out.println(url);
                                String imFileStr = GlobalInitializer.imgsDirForCurRun+"/MCL/images/"+genomeStr+"Path"+suf+".png";
                                saveImage(url, imFileStr);

                                valid = true;
                                System.out.println("/inColoringMaps, genome: "+genomeStr);
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

        public static void getPathwayWithBlackElems() throws RemoteException, IOException, ServiceException{

            

            KEGGLocator  locator;
            KEGGPortType serv;

            locator = new KEGGLocator();
            serv = locator.getKEGGPort();


            HashSet hashSetBlackECsListMcl = new HashSet(MainClass.blackECsListMcl);
            ArrayList<String> newBlackECsListMcl = new ArrayList<String>(hashSetBlackECsListMcl);


            int size = newBlackECsListMcl.size();



            for(int i=0; i<size; i++)
                System.out.println(newBlackECsListMcl.get(i));
            
            String fgColor = "#ffffff";
            String[] obj_list = new String[size];
            String[] fg_list = new String[size];
            String[] bg_list = new String[size];
            for (int i = 0; i < size; i++) {
                obj_list[i] = newBlackECsListMcl.get(i);
                bg_list[i] = "#0f0f00";
                fg_list[i] = fgColor;
            }
            String pathStr = "map"+MainClass.pathId;
            boolean valid =false;
            while(!valid){
                try{
                    String url = serv.color_pathway_by_objects(pathStr, obj_list, fg_list, bg_list);
                    System.out.println(url);
                    String imFileStr = GlobalInitializer.imgsDirForCurRun + "/MCL/images/blackElements_map.png";
                    saveImage(url, imFileStr);

                    valid = true;
                } catch(Exception e){
                                    System.out.println("Error: "+e.getMessage());
                }
            }
        }


        public static void getMultipleClustersFromECs(){

            List<String> sortedBlackECsListMcl = new ArrayList<String>();
            HashSet hs = new HashSet();
            hs.addAll(MainClass.blackECsListMcl);
            sortedBlackECsListMcl.clear();
            sortedBlackECsListMcl.addAll(hs);
            //Collections.copy(sortedBlackECsListMcl, MainClass.blackECsListMcl);
            Collections.sort(sortedBlackECsListMcl);


            List<List<Integer>> sortedBlackClustersFromECsListMcl = new ArrayList<List<Integer>>();


            for(int iter=0; iter<sortedBlackECsListMcl.size(); iter++){
                
                List<Integer> insertedClusters = new ArrayList<Integer>();
                insertedClusters.clear();

                for(int curEcId=0; curEcId<MainClass.blackECsListMcl.size(); curEcId++){

                    if(sortedBlackECsListMcl.get(iter).equals(MainClass.blackECsListMcl.get(curEcId))){

                        //add cluster from 'cur' list
                        boolean insertNewCluster = true;
                        for(int k=0; k<insertedClusters.size(); k++){
                            if(MainClass.blackCurClustersListMcl.get(curEcId)==insertedClusters.get(k)){
                                insertNewCluster = false;
                                break;
                            }
                        }
                        if(insertNewCluster)
                            insertedClusters.add(MainClass.blackCurClustersListMcl.get(curEcId));

                        //add cluster from 'exam' list
                        insertNewCluster = true;
                        for(int k=0; k<insertedClusters.size(); k++){
                            if(MainClass.blackClustersListMcl.get(curEcId)==insertedClusters.get(k)){
                                insertNewCluster = false;
                                break;
                            }
                        }
                        if(insertNewCluster)
                            insertedClusters.add(MainClass.blackClustersListMcl.get(curEcId));

                    }

                }

                sortedBlackClustersFromECsListMcl.add(insertedClusters);
            }


            List<List<String[]>> sortedBlackGenesListMcl = new ArrayList<List<String[]>>();
            sortedBlackGenesListMcl.clear();

            for(int iter=0; iter<sortedBlackECsListMcl.size(); iter++){

                List<String[]> intermedGenesListToInsert = new ArrayList<String[]>();

                for(int clustId=0; clustId<sortedBlackClustersFromECsListMcl.get(iter).size(); clustId++){

                    int curCluster = sortedBlackClustersFromECsListMcl.get(iter).get(clustId);

                    List<String> genesToInsert = new ArrayList<String>();

                    for(int curEcId=0; curEcId<MainClass.blackECsListMcl.size(); curEcId++){
                        if(sortedBlackECsListMcl.get(iter).equals(MainClass.blackECsListMcl.get(curEcId))){

                            if(MainClass.blackCurClustersListMcl.get(curEcId)==curCluster){

                                boolean insertGeneFlag = true;
                                for(int k=0; k<genesToInsert.size(); k++){
                                    if(genesToInsert.get(k).equals(MainClass.blackCurGeneList.get(curEcId))){
                                        insertGeneFlag = false;
                                        break;
                                    }
                                }

                                if(insertGeneFlag)
                                    genesToInsert.add(MainClass.blackCurGeneList.get(curEcId));
                            }

                            if(MainClass.blackClustersListMcl.get(curEcId)==curCluster){

                                boolean insertGeneFlag = true;
                                for(int k=0; k<genesToInsert.size(); k++){
                                    if(genesToInsert.get(k).equals(MainClass.blackExamGeneList.get(curEcId))){
                                        insertGeneFlag = false;
                                        break;
                                    }
                                }

                                if(insertGeneFlag)
                                    genesToInsert.add(MainClass.blackExamGeneList.get(curEcId));
                            }
                        }
                    }

                    String[] tmpGenesToInsert = new String[genesToInsert.size()];
                    tmpGenesToInsert = genesToInsert.toArray(tmpGenesToInsert);

                    intermedGenesListToInsert.add(tmpGenesToInsert);

                }

                    sortedBlackGenesListMcl.add(intermedGenesListToInsert);
            }

            

            MainClass.outputLogStrMcl += "{Total number of 'BLACK' elements: "+sortedBlackECsListMcl.size()+"}\n"+
                                    "- Clusters/Genes extracted from 'BLACK' elements:\n";
                                        
            for(int i=0; i<sortedBlackECsListMcl.size(); i++){
                MainClass.outputLogStrMcl += "> "+sortedBlackECsListMcl.get(i)+": ";
                for(int j=0; j<sortedBlackClustersFromECsListMcl.get(i).size(); j++){
                    MainClass.outputLogStrMcl += "\nCluster "+(sortedBlackClustersFromECsListMcl.get(i).get(j)+1)+" -> [";
                    for(int k=0; k<sortedBlackGenesListMcl.get(i).get(j).length; k++)
                        MainClass.outputLogStrMcl +=sortedBlackGenesListMcl.get(i).get(j)[k]+" ";
                    MainClass.outputLogStrMcl +="]";
                }

                MainClass.outputLogStrMcl += "\n";
            }

        }


        


        public static float calcJaccardDist(double[] a, double[] b){

           float distance = -1;
           float m00, m01, m10, m11;
           m00=0;
           m01=0;
           m10=0;
           m11=0;


           for(int i=0; i<a.length; i++){
                if(a[i]==0.0){
                    if(b[i]==1.0)
                        m01++;
                    else
                        m00++;
                }
                else if(a[i] == 1.0){
                    if(b[i]==0.0)
                        m10++;
                    else
                        m11++;
                }
                else{
                    System.out.println("An error occured in function calcJaccardDist()");
                }
           }

           distance = (m11+m00)/(m01 + m10 + m11);

           return distance;
        }



        public int findGenomeIdByGene(String geneSubStr){

            int genomeId = -1;
            for(int i=0; i<MainClass.orgsIds.length; i++){
                if(geneSubStr.equals(MainClass.orgsIds[i])){
                    genomeId =  i;
                    break;
                }
            }

            return genomeId;
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
}
