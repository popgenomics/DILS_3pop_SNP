#!/usr/bin/bash
## launch DILS for 2 populations
## the provided argument is for --configfile, expecting the yaml file
module load java-jdk/8.0.112
module load pypy/2.7-5.10.0
module load snakemake/5.3.0
module load r/3.6.3
module load python/2.7
binpath="/shared/ifbstor1/home/croux/scratch/morpho/DILS_3pop_SNP/bin"
snakemake --snakefile ${binpath}/Snakefile_3pop -p -j 500 --configfile ${1} --cluster-config ${binpath}/cluster_3pop.json --cluster "sbatch --nodes={cluster.node} --ntasks={cluster.n} --cpus-per-task={cluster.cpusPerTask} --time={cluster.time}"

