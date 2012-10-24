
import java.util.logging.Level;
import java.util.logging.Logger;

public class EmClusteringThread implements Runnable{

    public void run(){
        try {
            MainClass.algoThreadsCnt++;
            EmClustering emCl = new EmClustering();
            emCl.createArff(MainClass.orgsIds.length, MainClass.orgsIds, MainClass.totalGenesNumber);
            emCl.cluster();
            emCl.getNumOfClusters();
            MainClass.clustersListEm = emCl.getClusters();
         
            emCl.validateClusters();
            
            emCl.createFastaForClusters();
            MainClass.algoThreadsCnt--;
            
        } catch (Exception ex) {
            Logger.getLogger(EmClusteringThread.class.getName()).log(Level.SEVERE, null, ex);
        }
        }
}
