import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.PrintStream;
import java.net.URL;
import java.rmi.RemoteException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.xml.rpc.ServiceException;
import keggapi.KEGGLocator;
import keggapi.KEGGPortType;
import weka.clusterers.EM;
import weka.core.Instances;


class EmClustering {

        
	public Instances graphData;
        public static int numOfClusters;
        private static String[][] clustersArray;
        private static List<List<String>> clustersList;
	private static int[][] genesByGenomeInClusters;
        private static double[][] clustersHomogenValidMatEm;
        private DecimalFormat twoDForm;
       
	public EmClustering(){
            clustersList = new ArrayList<List<String>>();
            twoDForm = new DecimalFormat("#.000");
        }

	public void createArff(int orgsNum, String[] orgsIds, int totalGenesNumber){
		//create an arff file for the graph representation
		try{
			//create EM files directory
			String emDirStr = GlobalInitializer.genDataDirStr+"/EM_Files";
			boolean success = (new File(emDirStr)).mkdir();
			if (success) {
				System.out.println("Directory: "+emDirStr+" was created succesfully.");
			}
			
			FileWriter arffFstream = new FileWriter(GlobalInitializer.genDataDirStr+"/EM_Files/graph.arff");
			BufferedWriter arffOut = new BufferedWriter(arffFstream);
			arffOut.write("@RELATION P_matrix_data\n\n");
			for (int i=0; i<orgsNum; i++) {
				arffOut.write("@ATTRIBUTE "+orgsIds[i]+"_homology"+"\tNUMERIC\n");
			}
                        
			
			arffOut.write("\n@DATA\n");
			for (int i=0; i<totalGenesNumber; i++) {
				String geneStr = MainClass.getGeneNameByIdxInP(i, orgsNum);
				for (int j=0; j<orgsNum; j++) {
					arffOut.write(Double.toString(MainClass.P[i][j]));
                                        if(j<orgsNum-1)
                                            arffOut.write(",");
				}
				
				if(i != (totalGenesNumber-1))
					arffOut.write("\r");
			}
			arffOut.close();
		} catch (Exception e){
			System.err.println("Error: " + e.getMessage());
		}
		
		try{
			BufferedReader arffReader = new BufferedReader(new FileReader(GlobalInitializer.genDataDirStr+"/EM_Files/graph.arff"));
			graphData = new Instances(arffReader);
			arffReader.close();
		} catch (Exception e){
			System.err.println("Error: " + e.getMessage());
		}
		
	}	
		
	public void cluster(){

                PrintStream orgStream = null;
		PrintStream fileStream 	= null;
                String emOutput = null;

		try{

                    // Saving the orginal stream
                    orgStream = System.out;
                    fileStream = new PrintStream(new FileOutputStream(GlobalInitializer.genDataDirStr+"/EM_Files/clustersTemp.txt"));
                    // Redirecting console output to file
                    System.setOut(fileStream);
                    // Redirecting runtime exceptions to file
                    System.setErr(fileStream);

                    String[] options = {"-I", "100", "-t", GlobalInitializer.genDataDirStr+"/EM_Files/graph.arff", "-T", GlobalInitializer.genDataDirStr+"/EM_Files/graph.arff", "-p", "0"};
                    EM.main(options);

                    
                    //repeat algorithm execution to print standard output to console
                    EM emCluster = new EM();

                    emCluster.setOptions(options);
                    emCluster.buildClusterer(graphData);
                    emOutput = emCluster.toString();
                    System.out.println(emOutput);
                    
                    throw new Exception("Test Exception");
                    
                }
                catch (FileNotFoundException fnfEx)
		{
			System.out.println("Error in IO Redirection");
		}
		catch (Exception ex)
		{
			System.out.println("Redirecting output & exceptions to file");
		}
		finally
		{
			//Restoring back to console
			System.setOut(orgStream);
                        System.setErr(orgStream);
                        
		}
		
	}

        public void getNumOfClusters(){

            numOfClusters = 0;
           
            try{
                    FileInputStream fstreamEm = new FileInputStream(GlobalInitializer.genDataDirStr+"/EM_Files/clustersTemp.txt");
                    DataInputStream inEm = new DataInputStream(fstreamEm);
                    BufferedReader brEm = new BufferedReader(new InputStreamReader(inEm));


                    int i,j;
                    i = j = 0;

                    for(int geneIdx=0; geneIdx<MainClass.totalGenesNumber; geneIdx++){

                        String line = brEm.readLine();

                    }

                    String clustLine = null;
                    for(i=0; i<6; i++){
                        clustLine = brEm.readLine();
                    }
                    String clStr = clustLine.replace("Number of clusters selected by cross validation: ", "");
                    numOfClusters = Integer.parseInt(clStr);
                    
                    brEm.close();
                    inEm.close();
                    fstreamEm.close();


		} catch (Exception e){
                        System.err.println("Error: " + e.getMessage());
		}

        }

        public List<List<String>> getClusters(){

            
            clustersArray = new String[numOfClusters][MainClass.totalGenesNumber];
            for(int i=0; i<clustersArray.length; i++){
                        for (int j=0; j<clustersArray[i].length; j++){
                            clustersArray[i][j] = "";
                        }
            }
            
            try{
                    FileInputStream fstreamEm = new FileInputStream(GlobalInitializer.genDataDirStr+"/EM_Files/clustersTemp.txt");
                    DataInputStream inEm = new DataInputStream(fstreamEm);
                    BufferedReader brEm = new BufferedReader(new InputStreamReader(inEm));

                    
                    int[] indexInCluster = new int[numOfClusters];
                    for (int i=0; i<indexInCluster.length; i++) {
                            indexInCluster[i] = -1;
                    }

                    int i,j;
                    i = j = 0;

                    for(int geneIdx=0; geneIdx<MainClass.totalGenesNumber; geneIdx++){

                        String line = brEm.readLine();

                        StringTokenizer stokenz = new StringTokenizer(line, " ");

                        String nextToken = stokenz.nextToken();
                        int clusterIdx = Integer.parseInt(stokenz.nextToken());
                        String tmpGeneStr = MainClass.getGeneNameByIdxInP(geneIdx, MainClass.orgsIds.length);

                        indexInCluster[clusterIdx]++;
                        clustersArray[clusterIdx][indexInCluster[clusterIdx]] = tmpGeneStr;
                        
                    }

                    
                    brEm.close();
                    inEm.close();
                    fstreamEm.close();

                    
                    List<List<String>> clustersListTemp = new ArrayList<List<String>>();

                    for(i=0; i<clustersArray.length; i++){
                        List<String> cluster = new ArrayList<String>();
                        for (j=0; j<clustersArray[i].length; j++){
                            if (!clustersArray[i][j].equals(""))
                                cluster.add(clustersArray[i][j]);
                        }
                        clustersListTemp.add(cluster);
                    }

                    for(i=0; i<clustersListTemp.size(); i++){
                        if(!clustersListTemp.get(i).isEmpty()){
                            clustersList.add(clustersListTemp.get(i));
                        }
                    }

                    
                    
                    
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
                    MainClass.outputLogStrEm += "> Pathway: map"+MainClass.pathId+".\n\n";
                    MainClass.outputLogStrEm += "> Genomes examined: ";
                    for(i=0; i<MainClass.orgsIds.length-1; i++)
                        MainClass.outputLogStrEm += MainClass.orgsIds[i]+", ";
                    MainClass.outputLogStrEm += MainClass.orgsIds[i]+".\n\n";

                    MainClass.outputLogStrEm += "> Total number of genes: ";
                    MainClass.outputLogStrEm += MainClass.totalGenesNumber+".\n\n";

                    MainClass.outputLogStrEm += "> Number of clusters: ";
                    MainClass.outputLogStrEm += clustersList.size()+".\n";

                    MainClass.outputLogStrEm += "\n";

                    MainClass.outputLogStrEm += "> Number of genes in each cluster:\n";
                    for( i=0; i<clustersList.size(); i++){
                        MainClass.outputLogStrEm += "- Cluster "+(i+1)+": ";
                        MainClass.outputLogStrEm += clustersList.get(i).size()+" [";
                        int k;
                        for(k=0; k<MainClass.orgsIds.length-1; k++)
                            MainClass.outputLogStrEm += MainClass.orgsIds[k]+"("+genesByGenomeInClusters[i][k]+"), ";
                        MainClass.outputLogStrEm += MainClass.orgsIds[k]+"("+genesByGenomeInClusters[i][k]+")].\n";
                    }
                    MainClass.outputLogStrEm += "\n";


                    

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


            try{
                int foundAtLeastOne = 0;
                MainClass.outputLogStrEm += "> List of clusters with genes belonging only to a single genome:";
                for(int i=0; i<clustersList.size(); i++){
                    if(foundUniqeGenesInCluster[i] == 1){
                        foundAtLeastOne = 1;
                        MainClass.outputLogStrEm += "\n- Cluster "+(i+1)+":\n";
                        MainClass.outputLogStrEm += "# Genome: "+genomeInClusterWithUniqeGenes[i]+"\n";
                        MainClass.outputLogStrEm += "# Genes: \n";

                        for(int j=0; j<clustersList.get(i).size(); j++)
                            MainClass.outputLogStrEm += clustersList.get(i).get(j)+"\n";
                    }
                }
                MainClass.outputLogStrEm += "\n";
                if(foundAtLeastOne == 0){
                    MainClass.outputLogStrEm += "empty";
                    MainClass.outputLogStrEm += "\n\n";
                }




            } catch(Exception e){
                System.err.println("Error: " + e.getMessage());
            }


        }

        public void validateClusters(){

                clustersHomogenValidMatEm = new double[clustersList.size()][clustersList.size()];

                try{

                    MainClass.outputLogStrEmAppend += "> Clustering validation:\n\n";
                    
                    //internal clustering validation
                    for(int iter=0; iter<clustersList.size(); iter++){

                        int homologiesTotalNumber = 0;
                        double homologiesMean = 0;

                        ArrayList curList = (ArrayList)clustersList.get(iter);
                        int curListSize = curList.size();

                        for(int i=0; i<curListSize; i++){

                            int iGeneIdx = MainClass.getIdxByNameInP((String)curList.get(i));

                            for(int j=0; j<curListSize; j++){
                                int jGeneIdx = MainClass.getIdxByNameInP((String)curList.get(j));

                                if(MainClass.C[iGeneIdx][jGeneIdx] <= BlastXMLParser.expValueThreshold){
                                    homologiesTotalNumber++;
                                }
                            }
			}
                        homologiesMean = (double)homologiesTotalNumber/(curListSize*curListSize);

                        clustersHomogenValidMatEm[iter][iter] = homologiesMean;
                        
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
                                    for(int i=0; i<curListSize0; i++){

                                        int iGeneIdx = MainClass.getIdxByNameInP((String)curList0.get(i));
                                        
                                        for(int j=0; j<curListSize1; j++){
                                            int jGeneIdx = MainClass.getIdxByNameInP((String)curList1.get(j));
                                            
                                            if(MainClass.C[iGeneIdx][jGeneIdx] <= BlastXMLParser.expValueThreshold){
                                                homologiesTotalNumber++;
                                            }
                                        }
                                    }
                                    homologiesMean = (double)homologiesTotalNumber/(curListSize0*curListSize1);

                                    clustersHomogenValidMatEm[ii][ij] = homologiesMean;

                                }
                            }
                        }
                    }

                    MainClass.outputLogStrEmAppend += "- ~Homologies/Gene~ between clusters:\n\n";
                    MainClass.outputLogStrEmAppend += "Clusters |    ";

                    for(int k=0; k<clustersList.size(); k++)
                        MainClass.outputLogStrEmAppend += (k+1)+"    |    ";
                    MainClass.outputLogStrEmAppend += "\n";
                    for(int k=0; k<(clustersList.size()*10+10); k++){
                        MainClass.outputLogStrEmAppend += "-";
                    }

                    for(int i=0; i<clustersList.size(); i++){
                        MainClass.outputLogStrEmAppend += "\n    "+(i+1)+"    |  ";
                        for(int j=0; j<clustersList.size(); j++){
                            double homog = Double.valueOf(twoDForm.format(clustersHomogenValidMatEm[i][j]));
                            MainClass.outputLogStrEmAppend += homog;

                            //count number of decimal places
                            String tmp = Double.toString(homog);
                            String[] res = tmp.split("\\.");

                            if(res[1].length() == 4)
                                    MainClass.outputLogStrEmAppend +=" |  ";
                            else if(res[1].length() == 3)
                                    MainClass.outputLogStrEmAppend +="  |  ";
                            else if(res[1].length() == 2)
                                MainClass.outputLogStrEmAppend +="   |  ";
                            else
                                MainClass.outputLogStrEmAppend +="    |  ";
                        }
                    }

                    MainClass.outputLogStrEmAppend += "\n\n";




                } catch(Exception e){
                    System.err.println("Error: " + e.getMessage());
		}

	}

        public void getColoredPathwaysByGenome(){

            boolean success = (new File((GlobalInitializer.imgsDirForCurRun)+"/EM/images").mkdir());
            if (success) {
                System.out.println("Directory: " + (GlobalInitializer.imgsDirForCurRun)+"/EM/images" + " was created succesfully.\n");
            }

            String suf = "";
            String algo = "EM";
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

        int[] foundUniqeGenesInCluster = new int[MainClass.clustersListEm.size()];
        String[] genomeInClusterWithUniqeGenes = new String[MainClass.clustersListEm.size()];


            for(int i=0; i<foundUniqeGenesInCluster.length; i++){
                foundUniqeGenesInCluster[i] = 1;
                genomeInClusterWithUniqeGenes[i] = null;
            }


            for(int j=0; j<MainClass.clustersListEm.size(); j++){

                ArrayList curList = (ArrayList)MainClass.clustersListEm.get(j);
                String baseGenome = ((String)curList.get(0)).substring(0,3);
                genomeInClusterWithUniqeGenes[j] = baseGenome;
                String queryGenome = null;

                for(int k=MainClass.clustersListEm.get(j).size()-1; k>=0; k--){
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
                MainClass.outputLogStrEm += "> List of clusters with genes belonging only to a single genome:";
                for(int i=0; i<MainClass.clustersListEm.size(); i++){
                    if(foundUniqeGenesInCluster[i] == 1){
                        foundAtLeastOne = 1;
                        MainClass.outputLogStrEm += "\n- Cluster "+(i+1)+":\n";
                        MainClass.outputLogStrEm += "# Genome: "+genomeInClusterWithUniqeGenes[i]+"\n";
                        MainClass.outputLogStrEm += "# Genes: \n";



                        for(int j=0; j<MainClass.clustersListEm.get(i).size(); j++){
                            MainClass.outputLogStrEm += MainClass.clustersListEm.get(i).get(j)+" (";
                            
                            //find gene in pr[][] matrix
                            for(int iter1=0; iter1<MainClass.pr.length; iter1++){
                            for(int iter2=0; iter2<MainClass.pr[iter1].length; iter2++){
                                
                                    if( ((String)MainClass.clustersListEm.get(i).get(j)).equals(MainClass.pr[iter1][iter2]) ){
                                        for(int k=0; k<MainClass.ecNumsList[iter1][iter2].length; k++){
                                            MainClass.outputLogStrEm += MainClass.ecNumsList[iter1][iter2][k];
                                            if(k<MainClass.ecNumsList[iter1][iter2].length-1)
                                                MainClass.outputLogStrEm += " ";
                                        }
                                    }
                                }
                            }

                            MainClass.outputLogStrEm += ")\n";
                        }


                        }
                    }

                MainClass.outputLogStrEm += "\n";

                if(foundAtLeastOne == 0){
                    MainClass.outputLogStrEm += "empty";
                    MainClass.outputLogStrEm += "\n\n";
                }


            } catch(Exception e){

                System.err.println("Error: " + e.getMessage());
            }
        }

        public static void getColoredMapsByClusterEm() throws RemoteException, IOException, ServiceException{

            try{
                 boolean success = (new File((GlobalInitializer.imgsDirForCurRun)+"/EM/images").mkdir());
                if (success) {
                    System.out.println("Directory: " + (GlobalInitializer.imgsDirForCurRun)+"/EM/images" + " was created succesfully.\n");
                }

            } catch (Exception e){
                System.err.println("Error: " + e.getMessage());
            }

			
			String fgColor = "#ffffff";
			
            String[] bgColors = {"#0000ff", "#ff0000", "#f0f000", "#00f0f0", "000f0f","0f0f0f", "000fff"};
			
            
			for(int i=0; i<MainClass.clustersListEm.size(); i++){
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
					for (int j = 0; j < MainClass.clustersListEm.get(i).size(); j++) {
						String curGene = MainClass.clustersListEm.get(i).get(j);
						int genomeIndex = -1;
						String genomeString = curGene.substring(0,3);
						for (int iter1 = 0; iter1 < MainClass.orgsIds.length; iter1++) {
							if(genomeString.equals(MainClass.orgsIds[iter1]))
								genomeIndex = iter1;
						}
						for (int iter2 = 0; iter2 < MainClass.pr[genomeIndex].length; iter2++) {
							if (((String) MainClass.clustersListEm.get(i).get(j)).equals(MainClass.pr[genomeIndex][iter2])) {
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
												for (int ii = 0; ii < MainClass.clustersListEm.size(); ii++) {
													if (ii != i) {
														for (int jj = 0; jj < MainClass.clustersListEm.get(ii).size(); jj++) {
															if ((examGene.equals(MainClass.clustersListEm.get(ii).get(jj))) && (!examGene.equals(curGene))) {
																blackFlag = 1;
																
																MainClass.blackECsListEm.add(MainClass.ecsFinal[ecIdxExam]);
																MainClass.blackCurClustersListEm.add(i);
																MainClass.blackClustersListEm.add(ii);
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
									
									String[] res = serv.get_ko_by_gene(MainClass.clustersListEm.get(i).get(j));
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
					
					boolean valid = false;
					while(!valid){
						try{
							System.out.println("downloading map img, cluster: "+(i+1));
							String pathStr = "map"+MainClass.pathId;
							String url;
							url = serv.color_pathway_by_objects(pathStr, obj_arr, fg_arr, bg_arr);
							System.out.println(url);
							String imFileStr = GlobalInitializer.imgsDirForCurRun + "/EM/images/cluster_" + (i + 1) + "_map.png";
							saveImage(url, imFileStr);
							
							valid = true;
						} catch (Exception e){
							System.err.println("Error: "+e.getMessage());
						}
					}
					System.out.println("/downloading map img, cluster: "+(i+1));
				} catch (IOException ex) {
					Logger.getLogger(EmClustering.class.getName()).log(Level.SEVERE, null, ex);
				}
            }


        }


        public static void getColoredPathwaysFunctionEm(String suf){

            boolean success = (new File((GlobalInitializer.imgsDirForCurRun)+"/EM/images").mkdir());
                if (success) {
                    System.out.println("Directory: " + (GlobalInitializer.imgsDirForCurRun)+"/EM/images" + " was created succesfully.\n");
                }

			
	 class nestedGetColorPathsThread implements Runnable{
				
                String genomeStr;
                String suf;
				
                private String fgColor = "#ffffff";
                ArrayList<String> obj_list = new ArrayList<String>();
                ArrayList<String> fg_list = new ArrayList<String>();
                ArrayList<String> bg_list = new ArrayList<String>();
                private String[] bgColors = {"#0000ff", "#ff0000", "#f0f000", "#00f0f0", "#f000f0", "#0f0f00", "000f0f", "0f000f", "fff000", "0fff00", "00fff0", "000fff"};

                
                nestedGetColorPathsThread(String genomeStr, String suf){
                    this.genomeStr = genomeStr;
                    this.suf = suf;
                }
                
                public void run(){
                    MainClass.nestedGetColorPathsThreadCntEm++;
                    System.out.println("test started");
                    
                    try{
                        int objIndex = 0;
                        int clusterOffset = 0;
                        
						
                        for(int i=0; i<MainClass.clustersListEm.size(); i++){
                            for(int j=0; j<MainClass.clustersListEm.get(i).size(); j++){
								
                                String curGene = MainClass.clustersListEm.get(i).get(j);
								
                                
								
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
                                        if( ((String)MainClass.clustersListEm.get(i).get(j)).equals(MainClass.pr[genomeIndex][iter2]) ){
											
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
														
                                                        for(int ii=0; ii<MainClass.clustersListEm.size(); ii++)
                                                            if(ii!=i){
                                                                for(int jj=0; jj<MainClass.clustersListEm.get(ii).size(); jj++){
                                                                    if( (examGene.equals(MainClass.clustersListEm.get(ii).get(jj))) && (!examGene.equals(curGene))){
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
							
                            clusterOffset += MainClass.clustersListEm.get(i).size();
                           
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
                                String imFileStr = GlobalInitializer.imgsDirForCurRun+"/EM/images/"+genomeStr+"Path"+suf+".png";
                                saveImage(url, imFileStr);
								
                                valid = true;
                                System.out.println("/inColoringMaps, genome: "+genomeStr);
                            } catch(Exception e){
								System.out.println("Error: "+e.getMessage());
                            }
                        }
                    } catch(Exception ex){
                        Logger.getLogger(EmClustering.class.getName()).log(Level.SEVERE, null, ex);
                    }
					
                    MainClass.nestedGetColorPathsThreadCntEm--;
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
                    if(MainClass.nestedGetColorPathsThreadCntEm > 0){
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
                if(MainClass.nestedGetColorPathsThreadCntEm > 0){
                    continue;
                }
                else
                    break;
            }
			


        }

        public static void getPathwayWithBlackElemsEm() throws RemoteException, IOException{

            int size = MainClass.blackECsListEm.size();
            String fgColor = "#ffffff";

            String[] obj_list = new String[size];
            String[] fg_list = new String[size];
            String[] bg_list = new String[size];

            for(int i=0; i<size; i++){
                obj_list[i] = MainClass.blackECsListEm.get(i);
                bg_list[i] = "#0f0f00";
                fg_list[i] = fgColor;
            }

            String pathStr = "map00010";
            String url = KeggClientInitializer.serv.color_pathway_by_objects(pathStr, obj_list, fg_list, bg_list);
            System.out.println(url);
            String imFileStr = GlobalInitializer.imgsDirForCurRun+"/EM/images/blackElements_map.png";
            saveImage(url, imFileStr);

        }


        public static void getMultipleClustersFromECsEm(){

            


            List<String> checkedECs = new ArrayList<String>();

            MainClass.clustersFromblackECsStrEm +="\n\n>Clusters for black colored EC numbers:\n";

            for(int i=0; i<MainClass.blackECsListEm.size(); i++){
                
                int flagChecked = 0;

                for(int m=0; m<checkedECs.size(); m++){
                    if(MainClass.blackECsListEm.get(i).equals(checkedECs.get(m))){
                        flagChecked = 1;
                        break;
                    }
                }

                if(flagChecked == 1){
                    flagChecked = 0;
                    continue;
                }
                else{
                    List<String> curGenesChecked = new ArrayList<String>();
                    List<String> examGenesChecked = new ArrayList<String>();


                    MainClass.clustersFromblackECsStrEm += "\n"+MainClass.blackECsListEm.get(i)+" :: \n";

                    curGenesChecked.add(MainClass.blackCurGeneListEm.get(i));
                    examGenesChecked.add(MainClass.blackExamGeneListEm.get(i));

                    MainClass.clustersFromblackECsStrEm += "- ( "+MainClass.blackCurGeneListEm.get(i)+" -> "+MainClass.blackCurClustersListEm.get(i)+", "+MainClass.blackExamGeneListEm.get(i)+" -> "+MainClass.blackClustersListEm.get(i)+" )\n";

                    checkedECs.add(MainClass.blackECsListEm.get(i));


                    for(int j=0; j<MainClass.blackECsListEm.size(); j++){
                        if( (i!=j) && ((MainClass.blackECsListEm.get(i)).equals(MainClass.blackECsListEm.get(j)) ) ){

                            int pairOfGenesCheckedFlag = 0;

                            for(int n=0; n<curGenesChecked.size(); n++){
                                if((MainClass.blackCurGeneListEm.get(j)).equals(curGenesChecked.get(n)) && (MainClass.blackExamGeneListEm.get(j)).equals(examGenesChecked.get(n))){
                                    pairOfGenesCheckedFlag = 1;
                                    break;
                                }
                                else if((MainClass.blackCurGeneListEm.get(j)).equals(examGenesChecked.get(n)) && (MainClass.blackExamGeneListEm.get(j)).equals(curGenesChecked.get(n))){
                                    pairOfGenesCheckedFlag = 1;
                                    break;
                                }
                            }

                            if(pairOfGenesCheckedFlag!=1){
                                curGenesChecked.add(MainClass.blackCurGeneListEm.get(j));
                                examGenesChecked.add(MainClass.blackExamGeneListEm.get(j));
                                MainClass.clustersFromblackECsStrEm += "- ( "+MainClass.blackCurGeneListEm.get(j)+" -> "+MainClass.blackCurClustersListEm.get(j)+", "+MainClass.blackExamGeneListEm.get(j)+" -> "+MainClass.blackClustersListEm.get(j)+" )\n";
                            }
                        }
                    }
                }
            }

        }


        public void createFastaForClusters() throws Exception{

            boolean success = (new File((GlobalInitializer.imgsDirForCurRun)+"/EM/fasta").mkdir());
            if (success) {
                System.out.println("Directory: " + (GlobalInitializer.imgsDirForCurRun)+"/EM/fasta" + " was created succesfully.\n");
            }

            FastaCreator clFasta = new FastaCreator();

            for(int i=0; i<clustersList.size(); i++){
                String[] clusterString = new String[clustersList.get(i).size()];
                for(int j=0; j<clustersList.get(i).size(); j++){
                    clusterString[j] = (String)clustersList.get(i).get(j);
                }
                String fastaString = clFasta.getFastaString(clusterString);
                String fastaId = GlobalInitializer.imgsDirForCurRun+"/EM/fasta/cluster_"+i+"_FastaFile.txt";
                clFasta.createFastaFile(fastaId, fastaString);
            }
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