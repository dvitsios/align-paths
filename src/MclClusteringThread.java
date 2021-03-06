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
