# DILS_3pop_SNP   
Version 1.0.0 of DILS_3pop_SNP used in Leroy et. __al.__ (2021).  This is a modified version of DILS:  
1. taking fasta file as input (ex: **Morpho_Peru_dils.fa.** in **/data**)  
2. investigate models of divergence with/without gene flow for 3 species/populations.  
3. designed for RAD datasets, by conditionning the mutations on the observed number of SNPs and not on the mutation rate.   
To adapt to your project, the best is to set a single directory where *all* analysis (different projects for different species, etc ...) will be conducted.   
Let's name this directory *projectpath* 
The path to such *projectpath* has to be specified in the Snakefile (named here *Snakefile_3pop*).  
   
## in DILS_3pop.sh  
line 11  
*binpath=*'' # full link to the bin directory of the DILSmcsnp depository  
  
## in Snakefile_3pop  
lines 44 and 47  
*binpath=*'' # full path to the bin directory of the DILSmcsnp depository  
*projectpath=* '' # full path to the directory where all DILS analysis will be performed, i.e, where yaml files have to be placed  
  
## in yaml file  
*infile:* # full path to the fasta file, placed within the projectpath  
*config_yaml:* # full path to the current yaml file
