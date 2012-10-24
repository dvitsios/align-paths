import java.rmi.RemoteException;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.xml.rpc.ServiceException;
import keggapi.KEGGLocator;
import keggapi.KEGGPortType;

class getEnzsByPathwayThread implements Runnable{

            private static KEGGLocator  locator;
            private static KEGGPortType serv;
            private int orgIdx;

            public getEnzsByPathwayThread(int orgIdx) throws ServiceException{

                locator = new KEGGLocator();
                serv = locator.getKEGGPort();
                this.orgIdx = orgIdx;
            }

            public void run() {

                MainClass.enzThreadCnt++;
                boolean valid = false;
                while(!valid){
                    try {
                        String[] tmp = serv.get_enzymes_by_pathway("path:"+MainClass.orgsIds[orgIdx]+MainClass.pathId);
                        MainClass.enz[orgIdx] = new String[tmp.length];
                        System.arraycopy(tmp, 0, MainClass.enz[orgIdx], 0, tmp.length);

                        valid = true;
                    } catch (RemoteException ex) {
                        Logger.getLogger(MainClass.class.getName()).log(Level.SEVERE, null, ex);
                    }
                }



                MainClass.enzThreadCnt--;
            }

        }