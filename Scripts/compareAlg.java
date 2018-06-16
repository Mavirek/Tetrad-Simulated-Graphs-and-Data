package edu.cmu.tetrad.algcomparison.simulation;

import com.sun.org.apache.xerces.internal.impl.xpath.XPath;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DataUtils;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.search.*;
import edu.pitt.csb.mgm.*;
import edu.pitt.csb.stability.DataGraphSearch;
import edu.pitt.csb.stability.SearchWrappers;

import java.io.File;
import java.io.PrintStream;

/**
 * Created by Allen Poon on 6/23/2017.
 */

/**
 * Generates estimated graphs for the listed algorithms based on the generated data sets
 */
public class compareAlg {
    public static void main(String[] args) throws Exception {
        int numRuns = 20;
        int[] varSize = {50, 100};
        double[] alphaValues = {.0001, .001, .01, .03, .05, .08, .1};
        double[] penaltyDiscounts = {.5, 1, 2, 4, 8, 10, 20};
        int[] sampleSizes = {300,500};//,500};//{3000,5000};

        long startTime, endTime, totalTime = 0;

        //double[] lambdaValues = {.1, .2, .3, .35, .4, .5, .6}; //lambda values for STEPS for MGM
        double[] lambdaValues = new double[30];
        for (int i = 0; i < lambdaValues.length; i++) {
            lambdaValues[i] = 0.05 + i * (0.9) / 30;
        }
        //File starValues = new File("STAR Alpha+Penalty Values.txt");
        //PrintStream starOut = new PrintStream(starValues);

        File dir = new File("mgm Algorithm Comparison");
        dir.mkdirs();

        File curveDir = new File("mgm Curve Data");
        curveDir.mkdirs();

        File timeResults = new File(dir, "Algorithm Timed Results.txt");
        PrintStream timeOut = new PrintStream(timeResults);
        timeOut.println("GraphID, DataSetID, Algorithm, Elapsed Time");


        for (int n = 0; n < varSize.length; n++) {
            for (int ss = 0; ss < sampleSizes.length; ss++) {
                String dataSetID;
                DataSet data;
                DataGraphSearch c;
                STARS stars;
                double alphaToUse;
                File algResults;
                File curveResults;

                for (int x = 0; x < numRuns; x++) {
                    dataSetID = "DataSet_" + varSize[n] + "_" + (x + sampleSizes[ss] + 1);
                    data = MixedUtils.loadDataSet("mgmTrueGraphDataSets", "DataSet_" + varSize[n] + "_" + sampleSizes[ss] + "_" + (x + sampleSizes[ss] + 1) + "_M.txt");
                    //System.out.println("is mixed? "+data.isMixed());
                    //System.out.println("CPC " + dataSetID);
                    //CPC Stable
                    STEPS steps = new STEPS(data, lambdaValues, 0.01, 20);
                    Graph temp = steps.runSteps();

                    //starOut.println("CPC\t"+sampleSizes[ss]+"\t"+alphaToUse);
                    //System.out.println("CPC Alpha = "+alphaToUse);
                    //double CPCStableAlpha = alphaToUse;

//                    c = new SearchWrappers.CpcStableWrapper("Fisher Z", .05);
//                    stars = new STARS(c, data,false);
//                    alphaToUse = stars.runStars(alphaValues, dataSetID);
//                    CpcStable cpcs = new CpcStable(new IndTestFisherZ(data, alphaToUse));
//
//                    algResults = new File(dir, "CPC_" + dataSetID + ".txt");
//                    startTime = System.nanoTime();
//                    Graph CPCGraph = cpcs.search();
//                    endTime = System.nanoTime();
//                    GraphUtils.saveGraph(CPCGraph, algResults, false);
//                    totalTime = endTime - startTime;
//                    timeOut.println("Graph_" + (x + 1 + sampleSizes[ss]) + ", " + dataSetID + ", CPC Stable, " + (totalTime / 1000000000.0) + " s");

//                    for (int a = 0; a < alphaValues.length; a++) {
//                        //curveResults = new File(curveDir, "CPC_" + dataSetID + "_" + alphaValues[a] + ".txt");
//                        cpcs = new CpcStable(new IndTestFisherZ(data, alphaValues[a]));
//                        CPCGraph = cpcs.search();
//                        //GraphUtils.saveGraph(CPCGraph, curveResults, false);
//                    }

                    //mgm CPC Stable
                    System.out.println("mgm CPC Stable "+dataSetID);
                    c = new SearchWrappers.CpcStableWrapper("Multinomial AJ",.05);
                    stars = new STARS(c, data,false);
                    alphaToUse = stars.runStars(alphaValues,dataSetID);
                    algResults = new File(dir, "mgm-CPCStable_" + dataSetID + ".txt");
                    CpcStable cpcs = new CpcStable(new IndTestMultinomialAJ(data,alphaToUse));
                    cpcs.setInitialGraph(temp);
                    startTime = System.nanoTime();
                    Graph CPCGraph = cpcs.search();
                    GraphUtils.saveGraph(CPCGraph, algResults, false);
                    endTime = System.nanoTime();
                    totalTime = endTime - startTime;
                    timeOut.println("Graph_" + (x + 1 + sampleSizes[ss]) + ", " + dataSetID + ", CPC Stable, " + (totalTime / 1000000000.0) + " s");

                    for (int a = 0; a < alphaValues.length; a++) {
                        curveResults = new File(curveDir, "mgm-CPCStable_" + dataSetID + "_" + alphaValues[a] + ".txt");
                        cpcs = new CpcStable(new IndTestMultinomialAJ(data, alphaValues[a]));
                        cpcs.setInitialGraph(temp);
                        CPCGraph = cpcs.search();
                        GraphUtils.saveGraph(CPCGraph,curveResults,false);
                    }

                    //timeOut.println("\n\n");

                    //PC Max

//                    dataSetID = "DataSet_" +varSize[n]+"_" + (x + sampleSizes[ss] + 1);
//                    data = MixedUtils.loadDataSet("mgmTrueGraphDataSets", "DataSet_"+varSize[n]+"_" + sampleSizes[ss] + "_" + (x + sampleSizes[ss] + 1) + "_M.txt");
//                    System.out.println("PC Max " + dataSetID);
//
//                    //starOut.println("PCMax\t"+sampleSizes[ss]+"\t"+alphaToUse);
//                    //double PCMaxAlphaToUse = alphaToUse;
//                    //System.out.println("PCMax Alpha = "+alphaToUse);
//                    c = new SearchWrappers.PcMaxWrapper("Fisher Z", .05);
//                    stars = new STARS(c, data,false);
//                    alphaToUse = stars.runStars(alphaValues, dataSetID);
//                    PcMax pcm = new PcMax(new IndTestFisherZ(data, alphaToUse));
//
//                    algResults = new File(dir, "PCMax_" + dataSetID + ".txt");
//                    startTime = System.nanoTime();
//                    Graph PCMaxGraph = pcm.search();
//                    endTime = System.nanoTime();
//                    GraphUtils.saveGraph(PCMaxGraph, algResults, false);
//                    totalTime = endTime - startTime;
//                    timeOut.println("Graph_" + (x + 1 + sampleSizes[ss]) + ", " + dataSetID + ", PC Max, " + (totalTime / 1000000000.0) + " s");


//                    for (int a = 0; a < alphaValues.length; a++) {
//                        curveResults = new File(curveDir, "PCMax_" + dataSetID + "_" + alphaValues[a] + ".txt");
//                        pcm = new PcMax(new IndTestFisherZ(data, alphaValues[a]));
//                        PCMaxGraph = pcm.search();
//                        GraphUtils.saveGraph(PCMaxGraph, curveResults, false);
//                    }
//
                    //mgm-PC Max
                    System.out.println("mgm PC Max "+dataSetID);
                    c = new SearchWrappers.PcMaxWrapper("Multinomial AJ",.05);
                    stars = new STARS(c, data,false);
                    alphaToUse = stars.runStars(alphaValues,dataSetID);
                    algResults = new File(dir,"mgm-PCMax_"+dataSetID+".txt");
                    PcMax pcm = new PcMax(new IndTestMultinomialAJ(data, alphaToUse));
                    pcm.setInitialGraph(temp);
                    startTime = System.nanoTime();
                    Graph PCMaxGraph = pcm.search();
                    GraphUtils.saveGraph(PCMaxGraph,algResults,false);
                    endTime = System.nanoTime();
                    totalTime = endTime - startTime;
                    timeOut.println("Graph_"+(x+1 + sampleSizes[ss])+", " +dataSetID+", mgm-PC Max, "+(totalTime/1000000000.0)+" s");

                    for (int a = 0; a < alphaValues.length; a++) {
                        curveResults = new File(curveDir, "mgm-PCMax_" + dataSetID + "_" + alphaValues[a] + ".txt");
                        pcm = new PcMax(new IndTestMultinomialAJ(data, alphaValues[a]));
                        pcm.setInitialGraph(temp);
                        PCMaxGraph = pcm.search();
                        GraphUtils.saveGraph(PCMaxGraph,curveResults,false);
                    }

                    //PC Stable
//                    System.out.println("PC Stable " + dataSetID);
//                    c = new SearchWrappers.PcStableWrapper("Fisher Z", .05);
//                    stars = new STARS(c, data,false);
//                    alphaToUse = stars.runStars(alphaValues, dataSetID);
//                    PcStable pcs = new PcStable(new IndTestFisherZ(data, alphaToUse));
//
//                    algResults = new File(dir, "PCStable_" + dataSetID + ".txt");
//                    startTime = System.nanoTime();
//                    Graph PCStableGraph = pcs.search();
//                    endTime = System.nanoTime();
//                    GraphUtils.saveGraph(PCStableGraph, algResults, false);
//                    totalTime = endTime - startTime;
//                    timeOut.println("Graph_" + (x + 1 + sampleSizes[ss]) + ", " + dataSetID + ", PC Stable, " + (totalTime / 1000000000.0) + " s");

//                    for (int a = 0; a < alphaValues.length; a++) {
//                        curveResults = new File(curveDir, "PCStable_" + dataSetID + "_" + alphaValues[a] + ".txt");
//                        pcs = new PcStable(new IndTestFisherZ(data, alphaValues[a]));
//                        PCStableGraph = pcs.search();
//                        GraphUtils.saveGraph(PCStableGraph, curveResults, false);
//                    }

                    //mgm-PC Stable
                        System.out.println("mgm PC Stable "+dataSetID);
                        c = new SearchWrappers.PcStableWrapper("Multinomial AJ",.05);
                        stars = new STARS(c, data,false);
                        alphaToUse = stars.runStars(alphaValues,dataSetID);
                       algResults = new File(dir,"mgm-PCStable_"+dataSetID+".txt");
                        //double[] stepsResults = steps.runSteps();
                       PcStable pcs = new PcStable(new IndTestMultinomialAJ(data, alphaToUse));
                       pcs.setInitialGraph(temp);
                       startTime = System.nanoTime();
                       Graph PCStableGraph = pcs.search();
                        GraphUtils.saveGraph(PCStableGraph,algResults,false);
                        endTime = System.nanoTime();
                       totalTime = endTime - startTime;
                       timeOut.println("Graph_"+(x+1 + sampleSizes[ss])+", " +dataSetID+", mgm-PC Stable, "+(totalTime/1000000000.0)+" s");

                    for (int a = 0; a < alphaValues.length; a++) {
                        curveResults = new File(curveDir, "mgm-PCStable_" + dataSetID + "_" + alphaValues[a] + ".txt");
                        pcs = new PcStable(new IndTestMultinomialAJ(data, alphaValues[a]));
                        pcs.setInitialGraph(temp);
                        PCStableGraph = pcs.search();
                        GraphUtils.saveGraph(PCStableGraph,curveResults,false);
                    }


                    //timeOut.println("\n\n");

                    //FGES
//
//                    dataSetID = "DataSet_" +varSize[n]+"_" + (x + sampleSizes[ss] + 1);
//                    data = MixedUtils.loadData("TrueGraphDataSets", "DataSet_"+varSize[n]+"_" + sampleSizes[ss] + "_" + (x + sampleSizes[ss] + 1) + "_C.txt");
//                    System.out.println("FGES " + dataSetID);
//                    c = new SearchWrappers.FgesWrapper(1);
//                    stars = new STARS(c, data,true);
//                    double penaltyToUse = stars.runStars(penaltyDiscounts, dataSetID);
//                    //starOut.println("FGES\t"+sampleSizes[ss]+"\t"+penaltyToUse);
//                    Fges f = new Fges(new SemBicScore(DataUtils.getCovMatrix(data), penaltyToUse));
//                    algResults = new File(dir, "FGES_" + dataSetID + ".txt");
//                    startTime = System.nanoTime();
//                    Graph FGESGraph = f.search();
//                    endTime = System.nanoTime();
//                    GraphUtils.saveGraph(FGESGraph, algResults, false);
//                    totalTime = endTime - startTime;
//                    timeOut.println("Graph_" + (x + 1 + sampleSizes[ss]) + ", " + dataSetID + ", FGES, " + (totalTime / 1000000000.0) + " s");

//                    for (int a = 0; a < penaltyDiscounts.length; a++) {
//                        curveResults = new File(curveDir, "FGES_" + dataSetID + "_" + penaltyDiscounts[a] + ".txt");
//                        f = new Fges(new SemBicScore(DataUtils.getCovMatrix(data), penaltyDiscounts[a]));
//                        FGESGraph = f.search();
//                        GraphUtils.saveGraph(FGESGraph, curveResults, false);
//                    }
//                }
                    //timeOut.println("\n\n");

                }
            }
            timeOut.close();
            //starOut.close();
        }
    }
}