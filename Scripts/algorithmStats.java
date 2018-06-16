package edu.pitt.csb.mgm;

import edu.cmu.tetrad.algcomparison.statistic.*;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.graph.Edge;
import edu.cmu.tetrad.graph.Endpoint;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.util.StatUtils;
import jdk.net.SocketFlow;

import java.io.File;
import java.io.PrintStream;
import java.util.List;

/**
 * Created by Allen Poon on 7/7/2017.
 */
public class algorithmStats {
    public static void main(String[] args) throws Exception
    {

        int numRuns = 20;
        int[] varSize = {100, 50};
        double[] alphaValues = {.0001, .001, .01, .03, .05, .08, .1};
        double[] penaltyDiscounts = {.5, 1, 2, 4, 8, 10, 20};
        int[] sampleSizes = {100,1000,3000,5000,10000};

        double[] lambdaValues = {.1, .2, .3, .35, .4, .5, .6}; //lambda values for STEPS for MGM

        String[] algorithms = new String[]{"CPC","PCMax","PCStable","FGES"};

        //adjacency stats
        double[][][][][][] adjStatTable = new double[varSize.length][algorithms.length][sampleSizes.length][3][alphaValues.length][numRuns]; //3 = prec+rec+fpr
        //alphaValues.length == penaltyDiscounts.length -> interchangeable
        //arrowhead stats
        double[][][][][][] ahStatTable = new double[varSize.length][algorithms.length][sampleSizes.length][3][alphaValues.length][numRuns]; //3 = prec+rec+fpr

        //average stats (prec/rec/fpr)
        double[][][][][] avgAdjTable = new double[varSize.length][algorithms.length][sampleSizes.length][3][alphaValues.length];
        double[][][][][] avgAhTable = new double[varSize.length][algorithms.length][sampleSizes.length][3][alphaValues.length];

        double[][][][][][] avgStarsTable = new double[varSize.length][algorithms.length][sampleSizes.length][2][3][numRuns]; //2 = adj+ah, 3 = prec+rec+fpr
        double[][][][][] avgStarsValues = new double[varSize.length][algorithms.length][sampleSizes.length][2][3]; //2 = adj+ah, 3 = prec+rec+fpr

        //avg stats (SD)
        double[][][][][] adjSDTable = new double[varSize.length][algorithms.length][sampleSizes.length][3][alphaValues.length];
        double[][][][][] ahSDTable = new double[varSize.length][algorithms.length][sampleSizes.length][3][alphaValues.length];

        Statistics stats = new Statistics();
        //stats.add(new AdjacencyPrecision());
        //stats.add(new AdjacencyRecall());
        //stats.add(new ArrowheadPrecision());
        //stats.add(new ArrowheadRecall());
        stats.add(new SHD());

        List<Statistic> statList = stats.getStatistics();



        System.out.println("Comparing True Graphs and Estimated Graphs");
        for(int n=0; n<varSize.length; n++) {
            File graphComparison = new File("Cont True+Est Graph Comparison - Size "+varSize[n]+".txt");
            PrintStream compareOut = new PrintStream(graphComparison);
            compareOut.println("Algorithm \t TrueGraphID \t Sample Size \t AP \t AR \t AHP \t AHR \t SHD \t adj F1 \t ah F1");

            //Avg adj and arrowhead stats
            File avgStats = new File("Cont Average Adj and AH - Size " + varSize[n] + ".txt");
            PrintStream avgOut = new PrintStream(avgStats);
            avgOut.println("Algorithm, Sample Size, Alpha/Penalty, \t AP, \t AR, \t AHP, \t AHR, \t SHD \t adj F1 \t ah F1");

            File sderrStats = new File("Cont SDErr Adj and AH - Size "+varSize[n]+".txt");
            PrintStream sderr = new PrintStream(sderrStats);
            sderr.println("Algorithm, Sample Size, Alpha/Penalty, \t apSD, \t arSD, \t ahpSD, \t ahrSD");

            File allAlphas = new File("Cont All Alphas Adj and AH - Size " + varSize[n] + ".txt");
            PrintStream allOut = new PrintStream(allAlphas);
            allOut.println("Algorithm \t TrueGraphID \t Sample Size \t Alpha/Penalty \t AP \t AR \t AHP \t AHR,\t SHD \t adj F1 \t ah F1");

            boolean isFGES = false;

            for (int i = 0; i < algorithms.length; i++) {
                if(i==3)
                {
                    isFGES = true;
                }
                System.out.println(algorithms[i]);
                for (int x = 0; x < sampleSizes.length; x++) {
                    double bestadjHM = 0.0;
                    double bestahHM = 0.0;
                    double bestadjAlpha = 0.0;
                    double bestahAlpha = 0.0;
                    double fgesAdjHM = 0.0;
                    double fgesAhHM = 0.0;
                    double fgesAdjPenalty = 0.0;
                    double fgesAhPenalty = 0.0;
                    String adjTrueGraphID = "";
                    String ahTrueGraphID = "";
                    String fgesAdjID = "";
                    String fgesAhID = "";
                    for (int y = 0; y < numRuns; y++) {
                        String trueGraphID = "Graph_"+varSize[n]+"_" + (y + 1 + sampleSizes[x]) + "_C";
                        Graph trueGraph = GraphUtils.loadGraphTxt(new File("TrueGraphDataSets", trueGraphID + ".txt"));
                        DataSet data = MixedUtils.loadData("TrueGraphDataSets","DataSet_"+varSize[n]+"_" + sampleSizes[x] + "_" + (x + sampleSizes[x] + 1) + "_C.txt");
                        //STARS continuous
                        Graph estGraph = GraphUtils.loadGraphTxt(new File("Algorithm Comparison", algorithms[i] + "_DataSet_" +varSize[n]+"_" +(sampleSizes[x] + y + 1) + ".txt"));
                        compareOut.print(algorithms[i] + "\t" + trueGraphID + "\t" + sampleSizes[x]+"\t");

                        double[] adjacencyResults = apPrecRec(estGraph, trueGraph, varSize[n]);
                        compareOut.print(adjacencyResults[0] + "\t" + adjacencyResults[1]+"\t");
                        avgStarsTable[n][i][x][0][0][y] = adjacencyResults[0];
                        avgStarsTable[n][i][x][0][1][y] = adjacencyResults[1];
                        avgStarsTable[n][i][x][0][2][y] = adjacencyResults[2];

                        double[] arrowheadResults = ahPrecRec(estGraph, trueGraph, varSize[n]);
                        compareOut.print(arrowheadResults[0] + "\t" + arrowheadResults[1] + "\t");
                        avgStarsTable[n][i][x][1][0][y] = arrowheadResults[0];
                        avgStarsTable[n][i][x][1][1][y] = arrowheadResults[1];
                        avgStarsTable[n][i][x][1][2][y] = arrowheadResults[2];

//                        for (Statistic stat : statList) {
//                            //System.out.println(stat.getAbbreviation() + ":" + stat.getValue(trueGraph,estGraph));
//                            double val = stat.getValue(trueGraph, estGraph);
//                            compareOut.print(val + "\t");
//                        }

                        //SHD
                        double[] shdValues = MixedUtils.structuralHamming(trueGraph,estGraph,data);
                        compareOut.print(shdValues[0]+"\t");//+" "+shdValues[1]+" "+shdValues[2]+" Sum: "+shdValues[3]);
                        //compareOut.println(sampleSizes[x]);

                        //F1
                        double starsAdjFOne = 2*((avgStarsTable[n][i][x][0][0][y]*avgStarsTable[n][i][x][0][1][y])/(avgStarsTable[n][i][x][0][0][y]+avgStarsTable[n][i][x][0][1][y]));
                        double starsAhFOne = 2*((avgStarsTable[n][i][x][1][0][y]*avgStarsTable[n][i][x][1][1][y])/(avgStarsTable[n][i][x][1][0][y]+avgStarsTable[n][i][x][1][1][y]));
                        compareOut.println(starsAdjFOne+"\t"+starsAhFOne);

                        //All alpha values continuous
                        //CPC, PCMax, and PCStable
                        if (i < algorithms.length-1) {
                            for (int a = 0; a < alphaValues.length; a++) {
                                Graph alphaEstGraph = GraphUtils.loadGraphTxt(new File("Curve Data", algorithms[i] + "_DataSet_" +varSize[n]+"_" + (sampleSizes[x] + y + 1) + "_" + alphaValues[a] + ".txt"));
                                allOut.print(algorithms[i]+"\t"+trueGraphID+"\t"+sampleSizes[x]+"\t"+alphaValues[a]+"\t");
                                //Adjacency Prec/Rec/fpr
                                double[] rocAdj = apPrecRec(alphaEstGraph, trueGraph, varSize[n]);
                                allOut.print(rocAdj[0]+ "\t"+rocAdj[1]+"\t");
                                adjStatTable[n][i][x][0][a][y] = rocAdj[0];
                                adjStatTable[n][i][x][1][a][y] = rocAdj[1];
                                //adjStatTable[n][i][x][2][a][y] = rocAdj[2];

                                //Arrowhead Prec/Rec/fpr
                                double[] rocAH = ahPrecRec(alphaEstGraph, trueGraph, varSize[n]);
                                allOut.print(rocAH[0]+"\t"+rocAH[1]+"\t");
                                ahStatTable[n][i][x][0][a][y] = rocAH[0];
                                ahStatTable[n][i][x][1][a][y] = rocAH[1];
                                ahStatTable[n][i][x][2][a][y] = rocAH[2];

                                //SHD
                                double[] alphasSHD = MixedUtils.structuralHamming(trueGraph,alphaEstGraph,data);
                                adjStatTable[n][i][x][2][a][y] = alphasSHD[0];
                                allOut.print(alphasSHD[0]+"\t");

                                //F1
                                double alphasAdjFOne = 2*((rocAdj[0]*rocAdj[1])/(rocAdj[0]+rocAdj[1]));
                                double alphasAhFOne = 2*((rocAH[0]*rocAH[1])/(rocAH[0]+rocAH[1]));
                                if(alphasAdjFOne>bestadjHM)
                                {
                                    bestadjHM = alphasAdjFOne;
                                    adjTrueGraphID = algorithms[i]+ " "+trueGraphID;
                                    bestadjAlpha = alphaValues[a];
                                }
                                if(alphasAhFOne>bestahHM)
                                {
                                    bestahHM = alphasAhFOne;
                                    ahTrueGraphID = algorithms[i]+" "+trueGraphID;
                                    bestahAlpha = alphaValues[a];
                                }
                                allOut.println(alphasAdjFOne+"\t"+alphasAhFOne);

                            }
                        }
                        //FGES
                        else {
                            for (int p = 0; p < penaltyDiscounts.length; p++) {
                                Graph penaltyEstGraph = GraphUtils.loadGraphTxt(new File("Curve Data", "FGES_DataSet_" +varSize[n]+"_" +(sampleSizes[x] + y + 1) + "_" + penaltyDiscounts[p] + ".txt"));
                                allOut.print(algorithms[i]+"\t"+trueGraphID+"\t"+sampleSizes[x]+"\t"+penaltyDiscounts[p]+"\t");

                                //Adjacency Prec/Rec/fpr
                                double[] rocAdj = apPrecRec(penaltyEstGraph, trueGraph,varSize[n]);
                                allOut.print(rocAdj[0]+ "\t"+rocAdj[1]+"\t");
                                adjStatTable[n][i][x][0][p][y] = rocAdj[0];
                                adjStatTable[n][i][x][1][p][y] = rocAdj[1];
                                //adjStatTable[n][i][x][2][p][y] = rocAdj[2];


                                //Arrowhead Prec/Rec/fpr
                                double[] rocAH = ahPrecRec(penaltyEstGraph, trueGraph, varSize[n]);
                                allOut.print(rocAH[0]+"\t"+rocAH[1]+"\t");
                                ahStatTable[n][i][x][0][p][y] = rocAH[0];
                                ahStatTable[n][i][x][1][p][y] = rocAH[1];
                                ahStatTable[n][i][x][2][p][y] = rocAH[2];

                                //SHD
                                double[] alphasSHD = MixedUtils.structuralHamming(trueGraph,penaltyEstGraph,data);
                                adjStatTable[n][i][x][2][p][y] = alphasSHD[0];
                                allOut.print(alphasSHD[0]+"\t");

                                //F1
                                double alphasAdjFOne = 2*((rocAdj[0]*rocAdj[1])/(rocAdj[0]+rocAdj[1]));
                                double alphasAhFOne = 2*((rocAH[0]*rocAH[1])/(rocAH[0]+rocAH[1]));
                                if(alphasAdjFOne>fgesAdjHM)
                                {
                                    fgesAdjHM = alphasAdjFOne;
                                    fgesAdjID = algorithms[i]+ " "+trueGraphID;
                                    fgesAdjPenalty = penaltyDiscounts[p];
                                }
                                if(alphasAhFOne>fgesAhHM)
                                {
                                    fgesAhHM = alphasAhFOne;
                                    fgesAhID = algorithms[i]+" "+trueGraphID;
                                    fgesAhPenalty = penaltyDiscounts[p];
                                }
                                allOut.println(alphasAdjFOne+"\t"+alphasAhFOne);
                            }
                        }
                    }

                    //average alpha/penalty discount values after 20 runs
                    for (int a = 0; a < alphaValues.length; a++) {
                        //averaging adj prec and rec
                        avgAdjTable[n][i][x][0][a] = StatUtils.mean(adjStatTable[n][i][x][0][a]);
                        avgAdjTable[n][i][x][1][a] = StatUtils.mean(adjStatTable[n][i][x][1][a]);
                        //averaging SHD
                        avgAdjTable[n][i][x][2][a] = StatUtils.mean(adjStatTable[n][i][x][2][a]);
                        //SD err adj prec/rec
                        adjSDTable[n][i][x][0][a] = StatUtils.sd(adjStatTable[n][i][x][0][a])/Math.sqrt(20);
                        adjSDTable[n][i][x][1][a] = StatUtils.sd(adjStatTable[n][i][x][1][a])/Math.sqrt(20);



                        //averaging ah prec and rec
                        avgAhTable[n][i][x][0][a] = StatUtils.mean(ahStatTable[n][i][x][0][a]);
                        avgAhTable[n][i][x][1][a] = StatUtils.mean(ahStatTable[n][i][x][1][a]);
                        //avgAhTable[n][i][x][2][a] = StatUtils.mean(ahStatTable[n][i][x][2][a]);
                        //SD err ah prec/rec
                        ahSDTable[n][i][x][0][a] = StatUtils.sd(ahStatTable[n][i][x][0][a])/Math.sqrt(20);
                        ahSDTable[n][i][x][1][a] = StatUtils.sd(ahStatTable[n][i][x][1][a])/Math.sqrt(20);


                        if (i < algorithms.length-1) //CPC, PCMax, PCStable
                        {
                            avgOut.print(algorithms[i] + "\t" + sampleSizes[x] + "\t" + alphaValues[a] + "\t" + avgAdjTable[n][i][x][0][a] + "\t" + avgAdjTable[n][i][x][1][a] + "\t" + avgAhTable[n][i][x][0][a] + "\t" + avgAhTable[n][i][x][1][a] + "\t" + avgAdjTable[n][i][x][2][a]);
                            double adjHMean = 2*((avgAdjTable[n][i][x][0][a]*avgAdjTable[n][i][x][1][a])/(avgAdjTable[n][i][x][0][a]+avgAdjTable[n][i][x][1][a]));
                            double ahHMean = 2*((avgAhTable[n][i][x][0][a]*avgAhTable[n][i][x][1][a])/(avgAhTable[n][i][x][0][a]+avgAhTable[n][i][x][1][a]));
                            avgOut.println("\t"+adjHMean+"\t"+ahHMean);
                            sderr.println(algorithms[i] + "\t" + sampleSizes[x] + "\t" + alphaValues[a] + "\t" + adjSDTable[n][i][x][0][a] + "\t" + adjSDTable[n][i][x][1][a] + "\t" + ahSDTable[n][i][x][0][a] + "\t" + ahSDTable[n][i][x][1][a]);

                        }
                        else //FGES
                        {
                            avgOut.print(algorithms[i] + "\t" + sampleSizes[x] + "\t" + penaltyDiscounts[a] + "\t" + avgAdjTable[n][i][x][0][a] + "\t" + avgAdjTable[n][i][x][1][a] + "\t" + avgAhTable[n][i][x][0][a] + "\t" + avgAhTable[n][i][x][1][a] + "\t" + avgAdjTable[n][i][x][2][a]);
                            double adjHMean = 2*((avgAdjTable[n][i][x][0][a]*avgAdjTable[n][i][x][1][a])/(avgAdjTable[n][i][x][0][a]+avgAdjTable[n][i][x][1][a]));
                            double ahHMean = 2*((avgAhTable[n][i][x][0][a]*avgAhTable[n][i][x][1][a])/(avgAhTable[n][i][x][0][a]+avgAhTable[n][i][x][1][a]));
                            avgOut.println("\t"+adjHMean+"\t"+ahHMean);
                            sderr.println(algorithms[i] + "\t" + sampleSizes[x] + "\t" + penaltyDiscounts[a] + "\t" + adjSDTable[n][i][x][0][a] + "\t" + adjSDTable[n][i][x][1][a] + "\t" + ahSDTable[n][i][x][0][a] + "\t" + ahSDTable[n][i][x][1][a]);
                        }
                    }
                    if(!isFGES)
                    {
                        System.out.println("Best F1 for Adj in: "+adjTrueGraphID+"\tF1: "+bestadjHM+"\tAlpha: "+bestadjAlpha);
                        System.out.println("Best F1 for AH in: "+ahTrueGraphID+"\tF1: "+bestahHM+"\tAlpha: "+bestahAlpha);
                    }
                    else
                    {
                        System.out.println("Best F1 for FGES Adj in: "+fgesAdjID+"\tF1: "+fgesAdjHM+"\tPenalty: "+fgesAdjPenalty);
                        System.out.println("Best F1 for FGES AH in: "+fgesAhID+"\tF1: "+fgesAhHM+"\tPenalty: "+fgesAhPenalty);
                    }
                }
            }

        }

        for(int n=0; n<varSize.length; n++) {
            File avgStarStats = new File("Average STARS Adj and AH - Size " + varSize[n] + ".txt");
            PrintStream starOut = new PrintStream(avgStarStats);
/*
            File avgStarsSD = new File("STARS SD Adj and AH - Size "+varSize[n]+".txt");
            PrintStream sdOut = new PrintStream(avgStarsSD);
            sdOut.println("apSD,\tarSD,\tahpSD,\tahrSD");
*/
            File avgStarsSDErr = new File("SDErr STARS Adj and AH - Size "+varSize[n]+".txt");
            PrintStream sdErrOut = new PrintStream(avgStarsSDErr);
            sdErrOut.println("apSDErr,\tarSDErr,\tahpSDErr,\tahrSDErr");

            for (int i = 0; i < algorithms.length; i++) //algorithm
            {
                for (int ss = 0; ss < sampleSizes.length; ss++) //sample size
                {
                    starOut.print(algorithms[i] + "\t" + sampleSizes[ss] + "\t");
                    for (int a = 0; a < 2; a++) //adj/ah
                    {
                        for (int p = 0; p < 3; p++) //prec/rec/fpr
                        {
                            avgStarsValues[n][i][ss][a][p] = StatUtils.mean(avgStarsTable[n][i][ss][a][p]);
                            starOut.print(avgStarsValues[n][i][ss][a][p] + "\t");
                            if(p<2) //exclude fpr
                            {
                                //sdOut.print(currSD + "\t");
                                sdErrOut.print(StatUtils.sd(avgStarsTable[n][i][ss][a][p])/Math.sqrt(sampleSizes[ss])+"\t");

                            }
                        }
                    }
                    starOut.println();
                    sdErrOut.println();
                }
            }
            starOut.close();
            sdErrOut.close();
        }
    }


    public static double[] apPrecRec(Graph est, Graph truth, int numVars)
    {
        double tp = 0;
        double fp = 0;
        double fn = 0;
        for(Edge e:est.getEdges())
        {
            /*
            e.getNode1()
                    e.getNode2() instanceof
                    */
            if(truth.getEdge(truth.getNode(e.getNode1().getName()),truth.getNode(e.getNode2().getName()))!=null)
            {
                tp++;
            }
            else
            {
                fp++;
            }
        }
        for(Edge e:truth.getEdges())
        {
            if(est.getEdge(est.getNode(e.getNode1().getName()),est.getNode(e.getNode2().getName()))==null)
            {
                fn++;
            }
        }
        //TN = N(N-1)/2 - TP-FP-FN
        //N(N-1)/2 = total number of edges

        double tn = ((numVars*(numVars-1))/2)-tp-fp-fn;

        double [] result = new double[3];
        result[0] = tp/(tp+fp); //adjacency precision
        result[1] = tp/(tp+fn);//adjacency recall or tpr
        result[2] = fp/(fp+tn); //false pos rate
        return result;
    }

    public static double [] ahPrecRec(Graph est, Graph truth, int numVars)
    {
        double arrowsTp = 0;
        double arrowsFp = 0;
        double   arrowsFn = 0;

        // Get edges from the true Graph to compute TruePositives, TrueNegatives and FalseNeagtives
        //    System.out.println(this.truth.getEdges());

        for (Edge edge : truth.getEdges()) {

            List<Edge> edges1 = est.getEdges(est.getNode(edge.getNode1().getName()), est.getNode(edge.getNode2().getName()));
            Edge edge1;

            if (edges1.size() == 1) {
                edge1 = edges1.get(0);
            } else {
                edge1 = est.getDirectedEdge(est.getNode(edge.getNode1().getName()), est.getNode(edge.getNode2().getName()));
            }

            //      System.out.println(edge1 + "(est)");

            Endpoint e1Est = null;
            Endpoint e2Est = null;

            if (edge1 != null) {
                e1Est = edge1.getProximalEndpoint(est.getNode(edge.getNode1().getName()));
                e2Est = edge1.getProximalEndpoint(est.getNode(edge.getNode2().getName()));
            }
            //      System.out.println(e1Est);
            //      System.out.println(e2Est);

            List<Edge> edges2 = truth.getEdges(edge.getNode1(), edge.getNode2());
            Edge edge2;

            if (edges2.size() == 1) {
                edge2 = edges2.get(0);
            } else {
                edge2 = truth.getDirectedEdge(edge.getNode1(), edge.getNode2());
            }

            //       System.out.println(edge2 + "(truth)");

            Endpoint e1True = null;
            Endpoint e2True = null;

            if (edge2 != null) {
                e1True = edge2.getProximalEndpoint(edge.getNode1());
                e2True = edge2.getProximalEndpoint(edge.getNode2());
            }
            //       System.out.println(e1True);
            //       System.out.println(e2True);


            if (e1True == Endpoint.ARROW && e1Est != Endpoint.ARROW) {
                arrowsFn++;
            }

            if (e2True == Endpoint.ARROW && e2Est != Endpoint.ARROW) {
                arrowsFn++;
            }

            if (e1True == Endpoint.ARROW && e1Est == Endpoint.ARROW) {
                arrowsTp++;
            }

            if (e2True == Endpoint.ARROW && e2Est == Endpoint.ARROW) {
                arrowsTp++;
            }



        }
        // Get edges from the estimated graph to compute only FalsePositives
        // System.out.println(this.est.getEdges());

        for (Edge edge : est.getEdges()) {

            List<Edge> edges1 = est.getEdges(edge.getNode1(), edge.getNode2());
            Edge edge1;

            if (edges1.size() == 1) {
                edge1 = edges1.get(0);
            } else {
                edge1 = est.getDirectedEdge(edge.getNode1(), edge.getNode2());
            }
            //      System.out.println(edge1 + "(est)");

            Endpoint e1Est = null;
            Endpoint e2Est = null;

            if (edge1 != null) {
                e1Est = edge1.getProximalEndpoint(edge.getNode1());
                e2Est = edge1.getProximalEndpoint(edge.getNode2());
            }
            //       System.out.println(e1Est);
            //       System.out.println(e2Est);

            List<Edge> edges2 = truth.getEdges(truth.getNode(edge.getNode1().getName()), truth.getNode(edge.getNode2().getName()));
            Edge edge2;

            if (edges2.size() == 1) {
                edge2 = edges2.get(0);
            } else {
                edge2 = truth.getDirectedEdge(truth.getNode(edge.getNode1().getName()), truth.getNode(edge.getNode2().getName()));
            }

            //          System.out.println(edge2 + "(truth)");

            Endpoint e1True = null;
            Endpoint e2True = null;

            if (edge2 != null) {
                e1True = edge2.getProximalEndpoint(truth.getNode(edge.getNode1().getName()));
                e2True = edge2.getProximalEndpoint(truth.getNode(edge.getNode2().getName()));
            }
            //          System.out.println(e1True);
            //          System.out.println(e2True);


            if (e1Est == Endpoint.ARROW && e1True != Endpoint.ARROW) {
                arrowsFp++;
            }

            if (e2Est == Endpoint.ARROW && e2True != Endpoint.ARROW) {
                arrowsFp++;
            }


        }
        //TN = N(N-1) - TP-FP-FN
        double arrowsTn = (numVars*(numVars-1))-arrowsTp-arrowsFp-arrowsFn;
        double [] result = new double[3];
        result[0] = arrowsTp/(arrowsTp+arrowsFp); //prec
        result[1] = arrowsTp/(arrowsTp+arrowsFn); //recall
        result[2] = arrowsFp/(arrowsFp+arrowsTn); //false pos rate
        return result;
    }
}
