# multi-omics-heterogeneity-analysis
MOHA algorithm published by John Graf and Maria Zavodszky, GE Global Research Center.

Please see the [published paper](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0188878) for details.


## Example Commands
First change the working directory to the location of the jar file that has 
a subdirectory called “example”.  The next line computes thresholds from a 
cell stat file with the biomarker measures called quant_006.csv    It will 
create the threshold filed called quant_006.csv.thresholds.txt  
```
java -jar MOHAtool.jar -computeThresholds=example/quant_006.csv -biomarkerMetricTag=_Cell_Mean
```
First change the working directory to the location of the jar file that has 
a subdirectory called “example”.  The next line computes thresholds from a 
cell stat file with the biomarker measures called quant_006.csv    It will 
create the threshold filed called quant_006.csv.thresholds.txt  
```
java -jar MOHAtool.jar -computeCellStates=example/quant_006.csv -thresholdFile=example/quant_006.csv.thresholds.txt
```
The final step is to compute the heterogeneity metrics  You need to give it the input 
file  “quant_006.csv.MarkerStates.txt” which will output the metrics to the 
file out_moha.txt.
```
java -jar MOHAtool.jar -computeHeterogeneity=example/quant_006.csv.MarkerStates.txt -outputFile=example/out_moha_006.txt -append=false
```

## Notes
It would be preferable to have enough build files stored here that MOHAtool.jar
can be easily rebuilt, rather than saving it in GitHub.  (to-do)