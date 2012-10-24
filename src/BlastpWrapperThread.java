class BlastpWrapperThread implements Runnable{

    private String baseOrg;
    private String queryOrg;

    BlastpWrapperThread(String baseOrg, String queryOrg){

        this.baseOrg = baseOrg;
        this.queryOrg = queryOrg;
        System.out.println("blastp: "+baseOrg+" "+queryOrg);
    }

    public void run(){

        MainClass.blastpThreadCnt++;

        BlastpWrapper blastpObject = new BlastpWrapper();
        blastpObject.alignSequences(baseOrg, queryOrg);

        MainClass.blastpThreadCnt--;
    }
}