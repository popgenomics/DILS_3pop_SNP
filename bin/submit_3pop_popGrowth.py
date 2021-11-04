#!/usr/bin/python
import sys
import os

# example:
# ./submit_3pop.py ab AbC 1M 2N 100 /home/croux/Programmes/ABC_4pop/3pop

if len(sys.argv) != 6:
	print('./submit_3pop_popGrowth.py ab_AbC_1M_2N 100 /shared/home/croux/scratch/morpho/bin_snp /shared/mfs/data/home/croux/scratch/morpho/bin/dev/config.yaml gof_RF/ABC')
	sys.exit()

model = sys.argv[1]

#AB_model = sys.argv[1] # ab = AM model between sp.A and sp.B; AB = SC model between sp.A and sp.B
#ABC_model = sys.argv[2] # abc = no ongoing migration between A<->C and B<->C, AbC ongoing migration only between A<->C (same for aBC); abc no ongoing migration
#M_model = sys.argv[3] # 1M or 2M
#N_model = sys.argv[4] # 1N or 2N
nMultilocusSim = int(sys.argv[2]) # number of multilocus simulations
binpath = sys.argv[3] # /shared/home/croux/scratch/morpho/bin_snp


AB_model = model.split('_')[0]
ABC_model = model.split('_')[1]
M_model = model.split('_')[2]
N_model = model.split('_')[3]

config_yaml = sys.argv[4]
analysis = sys.argv[5]

infile = open('bpfile', 'r')
tmp = infile.readline()
tmp = infile.readline().strip().split('\t')
nLoci = len(tmp)

#commande = 'python2 {0}/priorgen_3pop_popGrowth.py {1} {2} {3} {4} {5} {7} {8} | {0}/msnsam tbs {6} -t tbs -r tbs tbs -I 3 tbs tbs tbs 0 -n 1 tbs -en tbs 1 tbs -n 2 tbs -en tbs 2 tbs -n 3 tbs -en tbs 3 tbs -m 1 2 tbs -m 2 1 tbs -m 1 3 tbs -m 3 1 tbs -m 2 3 tbs -m 3 2 tbs -em tbs 1 2 tbs -em tbs 2 1 tbs -em tbs 1 3 tbs -em tbs 3 1 tbs -em tbs 2 3 tbs -em tbs 3 2 tbs -ej tbs 2 1 -en tbs 1 tbs -ej tbs 3 1 -eN tbs tbs | pypy {0}/mscalc_3pop.py'.format(binpath, AB_model, ABC_model, M_model, N_model, nMultilocusSim, nMultilocusSim*nLoci, config_yaml, analysis)
commande = 'python2 {0}/priorgen_3pop_popGrowth.py {1} {2} {3} {4} {5} {7} {8} | java -jar {0}/msms.jar tbs {6} -s tbs -r tbs tbs -I 3 tbs tbs tbs 0 -n 1 tbs -en tbs 1 tbs -n 2 tbs -en tbs 2 tbs -n 3 tbs -en tbs 3 tbs -m 1 2 tbs -m 2 1 tbs -m 1 3 tbs -m 3 1 tbs -m 2 3 tbs -m 3 2 tbs -em tbs 1 2 tbs -em tbs 2 1 tbs -em tbs 1 3 tbs -em tbs 3 1 tbs -em tbs 2 3 tbs -em tbs 3 2 tbs -ej tbs 2 1 -en tbs 1 tbs -ej tbs 3 1 -eN tbs tbs | pypy {0}/mscalc_3pop.py'.format(binpath, AB_model, ABC_model, M_model, N_model, nMultilocusSim, nMultilocusSim*nLoci, config_yaml, analysis)
#print(commande)
os.system(commande)

