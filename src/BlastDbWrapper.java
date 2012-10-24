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