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