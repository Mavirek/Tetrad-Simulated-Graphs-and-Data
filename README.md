# Tetrad-Simulated-Graphs-and-Data
Contains true graphs, estimated graphs, and datasets (continuous, discrete, mixed)

File Name Scheme:

- True Graphs:

Graph_<var-Size>_<sampleSize+GraphNum>_<Continuous/Discrete/Mixed>.txt

ex. Graph_100_10020_C.txt

	* Variable Size = 100
	* Sample Size = 10000
	* Graph Number = 20
	* Data Type = Continuous

- Data Sets:

(I apologize for not keeping it consistent with the graph naming scheme)

DataSet_<var-Size>_<sampleSize>_<sampleSize+IterationNum>_<Continuous/Discrete/Mixed>.txt

ex. DataSet_100_500_518_M.txt
	
	* Variable Size = 100
	* Sample Size = 500
	* Iteration Num = 18
	* Data Type = Mixed
	
- Estimated Graphs:

<Algorithm>_DataSet_<var-Size>_<sampleSize+IterationNum>.txt

ex. CPC_DataSet_50_107.txt

	* Variable Size = 50
	* Sample Size = 100
	* Iteration Num = 7
	* Data Type = Continuous (Depends on which folder the estimated graph is from. Different experiments call for different data types.)
	
	
Folders and their contents:

* Algorithm Comparison
	- Continuous data, estimated graphs from CPC, PC, PCStable, and FGES. All with STARS and STEPS for FGES
	
* CPC Independence Tests
	- Mixed data, estimated graphs from FisherZ, Chi Square, Multinomial AJ, Conditional Gaussian LRT. All with STARS
	
	
	