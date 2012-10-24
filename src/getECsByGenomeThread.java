import java.rmi.RemoteException;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.xml.rpc.ServiceException;
import keggapi.KEGGLocator;
import keggapi.KEGGPortType;

public class getECsByGenomeThread implements Runnable{

    private static KEGGLocator  locator;
    private static KEGGPortType serv;

    private String[] prTmp;
    private int orgIdx;
    private int simulNestedThreads;
   


    public getECsByGenomeThread(int orgIdx) throws ServiceException{

        prTmp = new String[MainClass.pr[orgIdx].length];
        System.arraycopy(MainClass.pr[orgIdx], 0, prTmp, 0, MainClass.pr[orgIdx].length);
        this.orgIdx = orgIdx;
        locator = new KEGGLocator();
        serv = locator.getKEGGPort();

        simulNestedThreads = 1;

    }

    public void run() {




        MainClass.ecThreadCnt++;

        for(int index=0; index<prTmp.length; index++){
            boolean valid = false;
            while(!valid){
                try {
                        String[] tmp = serv.get_enzymes_by_gene(prTmp[index]);

                        System.out.println(MainClass.orgsIds[orgIdx]+" "+prTmp[index]);

                        if(tmp.length>0){
                                MainClass.ec[orgIdx][index] = new String[tmp.length];
                                System.arraycopy(tmp, 0, MainClass.ec[orgIdx][index], 0, tmp.length);
                        }
                        else{
                            MainClass.ec[orgIdx][index] = new String[1];
                            for(int t=0; t<MainClass.ec[orgIdx][index].length; t++)
                                MainClass.ec[orgIdx][index][t] = "";
                        }

                        valid = true;
                    } catch (RemoteException ex) {
                        Logger.getLogger(getECsByGenomeThread.class.getName()).log(Level.SEVERE, null, ex);
                }
            }
        }

        
        MainClass.ecThreadCnt--;
    }

    
}
