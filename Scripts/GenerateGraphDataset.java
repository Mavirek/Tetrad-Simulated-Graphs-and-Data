package edu.pitt.csb.mgm;

import edu.cmu.tetrad.algcomparison.graph.RandomForward;
import edu.cmu.tetrad.algcomparison.graph.RandomGraph;
import edu.cmu.tetrad.algcomparison.independence.FisherZ;
import edu.cmu.tetrad.algcomparison.simulation.ConditionalGaussianSimulation;
import edu.cmu.tetrad.algcomparison.simulation.LeeHastieSimulation;
import edu.cmu.tetrad.algcomparison.simulation.STARS;
import edu.cmu.tetrad.algcomparison.simulation.SemSimulation;
import edu.cmu.tetrad.algcomparison.statistic.*;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DataUtils;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.search.*;
import edu.cmu.tetrad.util.Parameters;
import edu.pitt.csb.stability.DataGraphSearch;
import edu.pitt.csb.stability.SearchWrappers;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.io.*;
import java.io.PrintStream;


/**
 * Created by Allen Poon on 6/20/2017.
 */
public class GenerateGraphDataset {
     public static void main(String[] args) throws Exception
     {
        int numRuns = 20;
        int[] numVars = new int[]{50,100};
        int[] sampleSizes = new int[]{300};//,500};//{100,1000,3000,5000};//{100, 300, 500};
         NormalDistribution nd = new NormalDistribution(3,1);
         Parameters parameters = new Parameters();
         parameters.set("numRuns", numRuns);
        parameters.set("numLatents",0);
        parameters.set("avgDegree", nd.sample());
        parameters.set("maxDegree",100);
        parameters.set("maxIndegree",100);
        parameters.set("maxOutdegree",100);
        parameters.set("standardize",true);


        File dir = new File("mgmTrueGraphDataSets");
        dir.mkdirs();
        //True Graph and DataSet Generating
        //C = Continuous, D = Discrete
         for(int n=0; n<numVars.length; n++) {
             parameters.set("numMeasures", numVars[n]);
             for (int ss = 0; ss < sampleSizes.length; ss++) {
                 parameters.set("sampleSize", sampleSizes[ss]);
                 RandomGraph g = new RandomForward();
                 //Continuous Graph/DataSet for algorithms 1-5
//                 SemSimulation s = new SemSimulation(g);
//                 s.createData(parameters);

                 //Mixed Graph/DataSet for algorithms 1-5
                 ConditionalGaussianSimulation cgs = new ConditionalGaussianSimulation(g);
                 LeeHastieSimulation lhs = new LeeHastieSimulation(g);
                 parameters.set("percentDiscrete", 50);
                 parameters.set("numCategories", 4);
                 parameters.set("differentGraphs",true);
                 lhs.createData(parameters);

                 for (int i = 0; i < numRuns; i++) {

//                     //Continuous Graph/DataSet
//                     File graphOut = new File(dir, "Graph_" + numVars[n] + "_" + (i + 1 + sampleSizes[ss]) + "_C.txt");
//                     DataSet d = (DataSet) s.getDataModel(i);
//                     GraphUtils.saveGraph(s.getTrueGraph(i), graphOut, false);
//                     PrintStream dataOut = new PrintStream(new File(dir, "DataSet_" + numVars[n] + "_" + parameters.get("sampleSize") + "_" + (i + 1 + sampleSizes[ss]) + "_CGS.txt"));
//                     dataOut.println(d.toString());
//                     dataOut.close();


//                     //Mixed Graph/DataSet
                     File graphOut2 = new File(dir, "Graph_" + numVars[n] + "_" + (i + 1 + sampleSizes[ss]) + "_M.txt");
                     DataSet d2 = (DataSet) lhs.getDataModel(i);
                     GraphUtils.saveGraph(lhs.getTrueGraph(i), graphOut2, false);
                     PrintStream dataOut2 = new PrintStream(new File(dir, "DataSet_" + numVars[n] + "_" + parameters.get("sampleSize") + "_" + (i + 1 + sampleSizes[ss]) + "_M.txt"));
                     dataOut2.println(d2.toString());
                     dataOut2.close();

                 }
             }
         }
    }
}
