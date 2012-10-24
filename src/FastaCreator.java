import java.io.BufferedWriter;
import java.io.FileWriter;
import javax.xml.rpc.ServiceException;
import keggapi.KEGGLocator;
import keggapi.KEGGPortType;


class FastaCreator {
	
	private KEGGLocator  locator;
	private KEGGPortType serv;

	
	public FastaCreator() throws ServiceException{
		locator = new KEGGLocator();
		serv = locator.getKEGGPort();
	}

        public String[] getGenesByGenome(String orgId, String pathId) throws Exception{

                String query = "path:"+orgId+pathId;
		String[] retString = null;

                boolean valid = false;
                while(!valid){
                    try{

                        String[] genesByPathResults = serv.get_genes_by_pathway(query);
                        
                        retString = genesByPathResults;

                        int genesNumber = genesByPathResults.length;
                        
                        System.out.println(orgId+": "+genesNumber+" genes\n");
                        valid = true;
                    } catch(Exception e){
                        System.err.println("Error: " + e.getMessage());
                    }
                }

                return retString;
        }


	public String getFastaString(String[] genesQuery) throws Exception{
		
		String fastaString = null;
		String getFastaString = "-f -n a";
                int genesNumber = genesQuery.length;

                boolean valid = false;
                while(!valid){
                    try {
			//construct query string as an argument for bget()
			if (genesNumber <= 100) {
				for (int i=0; i<genesNumber; i++) 
					getFastaString += " "+genesQuery[i];
			
				//retrieve amino acid sequence in FASTA format
				fastaString = serv.bget(getFastaString);
			}
			else {
				fastaString ="";
				int iter = (genesNumber / 100);
				int rem = (genesNumber % 100);
			
				for (int i=0; i<iter; i++) {
					for (int j=0; j<100; j++) {
						getFastaString += " "+genesQuery[j + (100*i)];
					}
					fastaString += serv.bget(getFastaString);
					getFastaString = "-f -n a";
				}	
				for (int i=0; i<rem; i++) {
					getFastaString += " "+genesQuery[i +(100*iter)];
				}
				fastaString += serv.bget(getFastaString);
			}
                        valid = true;
                    } catch(Exception e){
                        System.err.println("Error: " + e.getMessage());
                    }
                }

                return fastaString;
		
	}

        public void createFastaFile(String fastaFilePath, String fastaString){
            //create FASTA-formatted txt files
            try{
                    FileWriter fstream = new FileWriter(fastaFilePath);
                    BufferedWriter out = new BufferedWriter(fstream);
                    out.write(fastaString);
                    out.close();
            }catch (Exception e){
                    System.err.println("Error: " + e.getMessage());
            }
        }
	
}