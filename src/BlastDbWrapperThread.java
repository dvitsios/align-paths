class BlastDbWrapperThread implements Runnable{

    private String baseOrg;

    BlastDbWrapperThread(String baseOrg){

        this.baseOrg = baseOrg;    
        System.out.println("customDb: "+baseOrg);
    }

    public void run(){

        MainClass.customDbThreadCnt++;

        BlastDbWrapper custBlastDbObject = new BlastDbWrapper();
        custBlastDbObject.createDb(baseOrg);

        MainClass.customDbThreadCnt--;
    }
}