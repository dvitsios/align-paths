import java.util.logging.Level;
import java.util.logging.Logger;

public class MclClusteringThread implements Runnable{

     public void run() {

        MainClass.algoThreadsCnt++;
        //MCL clustering
        MclClustering.createAbcFile(MainClass.orgsIds.length, MainClass.totalGenesNumber);
        for (int i = 0; i < MainClass.inflationParams.length; i++) {
            try {
                MclClustering mclCl = new MclClustering();
                mclCl.cluster(MainClass.inflationParams[i]);
                
                MainClass.clustersListMcl = mclCl.getClusters(MainClass.inflationParams[i]);
                
                mclCl.validateClusters(MainClass.inflationParams[i]);
                
                mclCl.createFastaForClusters();
            }
            
            catch (Exception ex) {
                Logger.getLogger(MclClusteringThread.class.getName()).log(Level.SEVERE, null, ex);
            }

        }
       
        MainClass.algoThreadsCnt--;
     }

}
