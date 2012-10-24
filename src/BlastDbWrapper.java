import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

class BlastDbWrapper {
	
	private List<String> args;
	private ProcessBuilder builder;
	
	BlastDbWrapper(){
		args = new ArrayList<String>();
	}
	
	public void createDb(String org){
		args.add("makeblastdb"); 
		args.add("-in");
		args.add(GlobalInitializer.genDataDirStr+"/FASTA_Files/"+org+"FastaFile.txt");
		args.add("-dbtype");
		args.add("prot");
		args.add("-out");
		args.add(GlobalInitializer.genDataDirStr+"/Custom_BLAST_DBs/"+org+"DB.db");
		builder = new ProcessBuilder(args);
		try {
			Process p = builder.start();
			p.waitFor();
		} catch (IOException e) {
			System.err.println("Error: " + e.getMessage());
		} catch(InterruptedException ie){
			System.out.println("Custom database creation thread was intrrupted by another thread while it was waiting!");
		}
	}
}