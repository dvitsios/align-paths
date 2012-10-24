import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

class BlastpWrapper {
	
	private List<String> args;
	private ProcessBuilder builder;
	
	BlastpWrapper(){
		args = new ArrayList<String>();
	}
	
	public void alignSequences(String baseOrg, String queryOrg){
		args.add("blastp"); 
		args.add("-db");
		args.add(GlobalInitializer.genDataDirStr+"/Custom_BLAST_DBs/"+baseOrg+"DB.db");
		args.add("-query");
		args.add(GlobalInitializer.genDataDirStr+"/FASTA_Files/"+queryOrg+"FastaFile.txt");
		args.add("-outfmt");
		args.add("5");
		args.add("-out");
		args.add(GlobalInitializer.genDataDirStr+"/Blast_Output_XML/"+queryOrg+"2"+ baseOrg+".xml");

		builder = new ProcessBuilder(args);
		try {
			Process p = builder.start();
			p.waitFor();
			
		} catch (IOException e) {
			System.err.println("Error: " + e.getMessage());
		}
		catch(InterruptedException ie){
			System.out.println("Blastp<"+baseOrg+", "+queryOrg+"> execution thread was intrrupted by another thread while it was waiting!");
		}
	}
}
