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

public class getGenesByGenomeThread implements Runnable {

    private String genomeStr;
    private int idx;

    public getGenesByGenomeThread(String genomeStr, int idx){
        this.genomeStr = genomeStr;
        this.idx = idx;
    }
    
    public void run() {
    
        MainClass.createFastaThreadCnt++;
        
        try {

            FastaCreator myFastaCreator = new FastaCreator();

            String[] tempPr = myFastaCreator.getGenesByGenome(genomeStr, MainClass.pathId);
            String fastaString = myFastaCreator.getFastaString(tempPr);
            myFastaCreator.createFastaFile(GlobalInitializer.genDataDirStr+"/FASTA_Files/"+genomeStr+"FastaFile.txt", fastaString);
            MainClass.pr[idx] = new String[tempPr.length];
            System.arraycopy(tempPr, 0, MainClass.pr[idx], 0, tempPr.length);

            
        } catch(Exception e){
                System.err.println("Error: " + e.getMessage());
        }

        MainClass.createFastaThreadCnt--;
    }
};