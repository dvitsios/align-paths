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
