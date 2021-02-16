# Channel-Genes
R Code for differential mRNA Expression analysis. Looks for channel protein changes in the human genome. Could do other organism with small changes.


# Usage

Simply load in your results from DESeq2 (or other DE analysis). Best is if your rownames have been converted to HGNC symbols, otherwise basic modifications to the code need to be made. Feedback is greatly appreciated. There are many neuron channels which is only amplified by all the different isotypes that exist as well. Keep an eye out for false results by checking gene descriptions in the dataframe output.

3 User inputs are required. They are clearly notated in the code.

# Version

### Dependencies

dplyr_1.0.4  
biomaRt_2.42.1  
hash_2.2.6.1  

### My R Session

R version 3.6.1 (2019-07-05)  
Platform: x86_64-conda_cos6-linux-gnu (64-bit)  
Running under: Ubuntu 20.04.1 LTS  
