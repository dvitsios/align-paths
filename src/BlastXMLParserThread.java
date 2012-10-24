class BlastXMLParserThread implements Runnable{

    private String xmlFile;
    private String baseOrg;
    private String queryOrg;
    private int baseIdx;
    private int queryIdx;

    BlastXMLParserThread(String baseOrg, String queryOrg, int baseIdx, int queryIdx){

        this.baseOrg = baseOrg;
        this.queryOrg = queryOrg;
        this.baseIdx = baseIdx;
        this.queryIdx = queryIdx;

        xmlFile = GlobalInitializer.genDataDirStr+"/Blast_Output_XML/" + queryOrg + "2" + baseOrg + ".xml";
                    
    }

    public void run(){

        MainClass.xmlParserThreadCnt++;

        BlastXMLParser parser = new BlastXMLParser();
        parser.parseXMLFile(xmlFile, baseOrg, queryOrg, baseIdx, queryIdx);
    
        MainClass.xmlParserThreadCnt--;
    }
}