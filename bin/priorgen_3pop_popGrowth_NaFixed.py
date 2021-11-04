#!/usr/bin/python
# #!/home/roux/python/Python-2.7.14/python
# -*- coding: utf-8 -*-
epsilon = 0.0001
import sys
from numpy.random import uniform
from numpy.random import binomial
from numpy.random import beta
from numpy.random import seed

if sys.argv[2] not in ['abc', 'AbC', 'aBC', 'ABC']:
	sys.exit('\n\tmodel {0} is not in the list of allowed models [abc, AbC, aBC, ABC]\n'.format(sys.argv[2]))

#seed(100)
# assumed topology :  ( (A; B); C )
help = "\n\t\033[1;31;40mTakes one indicator of migration within the focal group, another one for the outgroup, homo/hetero Ne, homo/hetero N.m and a number of multilocus simulations as arguments:\033[0m\n\n\t"
help += "\033[1;32;40m#Command line for ms: \033[0m\n\tjava -jar msms.jar tbs 10000 -s tbs -r tbs tbs -I 3 tbs tbs tbs 0 -n 1 tbs -en tbs 1 tbs -n 2 tbs -en tbs 2 tbs -n 3 tbs -en tbs 3 tbs -m 1 2 tbs -m 2 1 tbs -m 1 3 tbs -m 3 1 tbs -m 2 3 tbs -m 3 2 tbs -em tbs 1 2 tbs -em tbs 2 1 tbs -em tbs 1 3 tbs -em tbs 3 1 tbs -em tbs 2 3 tbs -em tbs 3 2 tbs -ej tbs 2 1 -en tbs 1 tbs -ej tbs 3 1 -eN tbs tbs\n"

help += "\t\033[1;32;40mExample: ./priorgen_3pop.py ab/AB abc/ABC 1M/2M 1N/2N 1000 /shared/mfs/data/home/croux/scratch/morpho/bin/dev/config.yaml gof_RF\033[0m\n"

if len(sys.argv) != 8:
	print(help)
	sys.exit()

def randomBeta(estimated_value, nMultilocus):
	a = 50.0
	b = 50.0
	scalar_tmp = beta(a=a, b=b, size=nMultilocus)
	scalar = [ i/(a/(a+b)) for i in scalar_tmp ]

	res = [ estimated_value * i for i in scalar ]
	return(res)

analysis = sys.argv[7]
# Configuration of the prior distribution
nMultilocus = int(sys.argv[5])

# prior from the yaml
N_bound = [0, 10]
T_bound = [0, 10]
M_bound = [0.4, 50]
shape_bound = [0.01, 20]

infile = open(sys.argv[6], 'r')
for line in infile:
	line = line.strip()
	line = line.split(':')
	if(line[0]=='Nref'):
		Nref = float(line[1])
	if(line[0]=='Nmin'):
		N_bound[0] = float(line[1])
	if(line[0]=='Nmax'):
		N_bound[1] = float(line[1])
	if(line[0]=='Tmin'):
		T_bound[0] = float(line[1])
	if(line[0]=='Tmax'):
		T_bound[1] = float(line[1])
	if(line[0]=='Mmin'):
		M_bound[0] = float(line[1])
	if(line[0]=='Mmax'):
		M_bound[1] = float(line[1])
	if(line[0]=='shapeMin'):
		shape_bound[0] = float(line[1])
	if(line[0]=='shapeMax'):
		shape_bound[1] = float(line[1])
infile.close()	

# read bpfile
infile = open("bpfile", "r")
tmp = infile.readline()
L = [ float(i) for i in infile.readline().strip().split("\t") ]
nsamA = [ int(i) for i in infile.readline().strip().split("\t") ]
nsamB = [ int(i) for i in infile.readline().strip().split("\t") ]
nsamC = [ int(i) for i in infile.readline().strip().split("\t") ]
nSNPs = [ int(i) for i in infile.readline().strip().split("\t") ]
rho = [ float(i) for i in infile.readline().strip().split("\t") ]
infile.close()

# number of loci
nLoci = len(L)

# sum of nsamA + nsamB
nsam_tot = [ nsamA[i] + nsamB[i] + nsamC[i] for i in range(nLoci) ]

if analysis == 'ABC':
	# param multilocus: values that will be printed in priorfile.txt
	## N = N_pop_i / Nref
	N1 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
#	N1 = [1] * nMultilocus
	
	N2 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	N3 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)

	N1a = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	N2a = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	N3a = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)

	Na_12 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus) # ancestor between N1 and N2
#	Na = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus) # N1N2-N3N4 
	Na = [1] * nMultilocus
	
	## factor of local reduction in Ne. Model of "background selection"
	shape_N_a = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
	shape_N_b = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)

	## times
	Tsplit = uniform(low = T_bound[0], high = T_bound[1], size = nMultilocus) # split ((A,B), C)
	Tsplit_12 = [ uniform(low = 0, high = Tsplit[i], size = 1)[0] for i in range(nMultilocus) ] # split (A, B)

	Tdem_1 = [ uniform(low = 0, high = Tsplit_12[i], size = 1)[0] for i in range(nMultilocus) ]
	Tdem_2 = [ uniform(low = 0, high = Tsplit_12[i], size = 1)[0] for i in range(nMultilocus) ]
	Tdem_3 = [ uniform(low = 0, high = Tsplit[i], size = 1)[0] for i in range(nMultilocus) ]

	if sys.argv[1] == 'AB': # if ongoing migration between A and B
		shape_T_a = 1
		shape_T_b = 5
	else:
		shape_T_a = 5
		shape_T_b = 1
	Tsmall_12 = [ beta(shape_T_a, shape_T_b, size = 1)[0]*((Tsplit_12[i]-T_bound[0])+T_bound[0]) for i in range(nMultilocus) ] # secondary contact or end of migration between (A,B) and C

	if sys.argv[2] == 'ABC': # if gene flow between (A;B) and C
		shape_T_a = 1
		shape_T_b = 5
	else:
		shape_T_a = 5
		shape_T_b = 1
	Tsmall_13 = [ beta(shape_T_a, shape_T_b, size = 1)[0]*((Tsplit[i]-T_bound[0])+T_bound[0]) for i in range(nMultilocus) ] # secondary contact or end of migration between (A,B) and C


	## Miration rates
	### priorgen AB AC BC  
	### between A and B
	M12 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
	M21 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)

	### between (A,B) and C
	if sys.argv[2] in ['ABC', 'AbC']:
		M13 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
		M31 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
	else:
		M13 = [0] * nMultilocus
		M31 = [0] * nMultilocus

	if sys.argv[2] in ['ABC', 'aBC']:
		M23 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
		M32 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
	else:
		M23 = [0] * nMultilocus
		M32 = [0] * nMultilocus

	### present and past migration rates
	if sys.argv[1] == 'AB':
		M12_current = M12 # uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
		M21_current = M21 # uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
		M12_past = [0] * nMultilocus
		M21_past = [0] * nMultilocus
	else:
		M12_current = [0] * nMultilocus
		M21_current = [0] * nMultilocus
		M12_past = M12 # uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
		M21_past = M21 # uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)

	if sys.argv[2] != 'abc':
		M13_current = M13 # uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
		M31_current = M31 # uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
		M23_current = M23 # uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
		M32_current = M32 # uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
		M13_past = [0] * nMultilocus
		M31_past = [0] * nMultilocus
		M23_past = [0] * nMultilocus
		M32_past = [0] * nMultilocus
	else:
		M13_current = [0] * nMultilocus
		M31_current = [0] * nMultilocus
		M23_current = [0] * nMultilocus
		M32_current = [0] * nMultilocus
		M13_past = M13 # uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
		M31_past = M31 # uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
		M23_past = M23 # uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
		M32_past = M32 # uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)

	## factor of local reduction in Me. One Beta distribution for each of the migration rates 
	shape_M12_a = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
	shape_M12_b = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
	shape_M21_a = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
	shape_M21_b = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)

	shape_M13_a = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus) # same Beta distribution for the migration A<->C than for B<->C
	shape_M13_b = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
	shape_M31_a = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus) # same Beta distribution for the migration A<->C than for B<->C
	shape_M31_b = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)

if analysis == 'gof_RF':
	param_RF = {}
	posteriorFile = open('../../modelComp/posterior_RF.txt', 'r')
	line = posteriorFile.readline()
	for line in posteriorFile:
		line = line.strip().split('\t')
		param_RF[line[0]] = [float(line[1])]
	posteriorFile.close()
	
	N1  = param_RF['N1']*nMultilocus
	N2  = param_RF['N2']*nMultilocus
	N3  = param_RF['N3']*nMultilocus
	N1a  = param_RF['N1a']*nMultilocus
	N2a  = param_RF['N2a']*nMultilocus
	N3a  = param_RF['N3a']*nMultilocus
	Na_12  = param_RF['Na_12']*nMultilocus
	Na  = param_RF['Na']*nMultilocus
	shape_N_a  = param_RF['shape_N_a']*nMultilocus
	shape_N_b  = param_RF['shape_N_b']*nMultilocus
	Tdem_1  = param_RF['Tdem_1']*nMultilocus
	Tdem_2  = param_RF['Tdem_2']*nMultilocus
	Tdem_3  = param_RF['Tdem_3']*nMultilocus
	Tsplit_12  = param_RF['Tsplit_12']*nMultilocus
	Tsplit  = param_RF['Tsplit']*nMultilocus
	Tsmall_12  = param_RF['Tsmall_AB']*nMultilocus
	Tsmall_13  = param_RF['Tsmall_ABC']*nMultilocus
	M12  = param_RF['M12']*nMultilocus
	shape_M12_a  = param_RF['shape_M12_a']*nMultilocus
	shape_M12_b  = param_RF['shape_M12_b']*nMultilocus
	M21  = param_RF['M21']*nMultilocus
	shape_M21_a  = param_RF['shape_M21_a']*nMultilocus
	shape_M21_b  = param_RF['shape_M21_b']*nMultilocus
	M13  = param_RF['M13']*nMultilocus
	M31  = param_RF['M31']*nMultilocus
	M23  = param_RF['M23']*nMultilocus
	M32  = param_RF['M32']*nMultilocus
	shape_M13_a  = param_RF['shape_M13_a']*nMultilocus
	shape_M13_b  = param_RF['shape_M13_b']*nMultilocus
	shape_M31_a  = param_RF['shape_M13_a']*nMultilocus
	shape_M31_b  = param_RF['shape_M13_b']*nMultilocus
	shape_M23_a  = param_RF['shape_M13_a']*nMultilocus
	shape_M23_b  = param_RF['shape_M13_b']*nMultilocus

	### present and past migration rates
	if sys.argv[1] == 'AB':
		M12_current = M12 # uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
		M21_current = M21 # uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
		M12_past = [0] * nMultilocus
		M21_past = [0] * nMultilocus
	else:
		M12_current = [0] * nMultilocus
		M21_current = [0] * nMultilocus
		M12_past = M12 # uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
		M21_past = M21 # uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)

	if sys.argv[2] != 'abc':
		M13_current = M13 # uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
		M31_current = M31 # uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
		M23_current = M23 # uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
		M32_current = M32 # uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
		M13_past = [0] * nMultilocus
		M31_past = [0] * nMultilocus
		M23_past = [0] * nMultilocus
		M32_past = [0] * nMultilocus
	else:
		M13_current = [0] * nMultilocus
		M31_current = [0] * nMultilocus
		M23_current = [0] * nMultilocus
		M32_current = [0] * nMultilocus
		M13_past = M13 # uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
		M31_past = M31 # uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
		M23_past = M23 # uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
		M32_past = M32 # uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)

if analysis == 'optimization1':
	param_RF = {}
	posteriorFile = open('../../modelComp/posterior_RF.txt', 'r')
	line = posteriorFile.readline()
	for line in posteriorFile:
		line = line.strip().split('\t')
		param_RF[line[0]] = float(line[1])
	posteriorFile.close()
	
	N1  = randomBeta(param_RF['N1'], nMultilocus)
	N2  = randomBeta(param_RF['N2'], nMultilocus)
	N3  = randomBeta(param_RF['N3'], nMultilocus)
	N1a  = randomBeta(param_RF['N1a'], nMultilocus)
	N2a  = randomBeta(param_RF['N2a'], nMultilocus)
	N3a  = randomBeta(param_RF['N3a'], nMultilocus)
	Na_12  = randomBeta(param_RF['Na_12'], nMultilocus)
	Na  = randomBeta(param_RF['Na'], nMultilocus)
	shape_N_a  = randomBeta(param_RF['shape_N_a'], nMultilocus)
	shape_N_b  = randomBeta(param_RF['shape_N_b'], nMultilocus)
	Tdem_1  = randomBeta(param_RF['Tdem_1'], nMultilocus)
	Tdem_2  = randomBeta(param_RF['Tdem_2'], nMultilocus)
	Tdem_3  = randomBeta(param_RF['Tdem_3'], nMultilocus)
	Tsplit_12  = randomBeta(param_RF['Tsplit_12'], nMultilocus)
	Tsplit  = randomBeta(param_RF['Tsplit'], nMultilocus)
	Tsmall_12  = randomBeta(param_RF['Tsmall_AB'], nMultilocus)
	Tsmall_12 = [ Tsmall_12[i] if Tsmall_12[i]<Tsplit_12[i] else Tsplit_12[i]-epsilon for i in range(nMultilocus) ] # to keep Tsmall_12 < Tsplit_12

	Tsmall_13  = randomBeta(param_RF['Tsmall_ABC'], nMultilocus)
	M12  = randomBeta(param_RF['M12'], nMultilocus)
	shape_M12_a  = randomBeta(param_RF['shape_M12_a'], nMultilocus)
	shape_M12_b  = randomBeta(param_RF['shape_M12_b'], nMultilocus)
	M21  = randomBeta(param_RF['M21'], nMultilocus)
	shape_M21_a  = randomBeta(param_RF['shape_M21_a'], nMultilocus)
	shape_M21_b  = randomBeta(param_RF['shape_M21_b'], nMultilocus)
	M13  = randomBeta(param_RF['M13'], nMultilocus)
	M31  = randomBeta(param_RF['M31'], nMultilocus)
	M23  = randomBeta(param_RF['M23'], nMultilocus)
	M32  = randomBeta(param_RF['M32'], nMultilocus)
	shape_M13_a  = randomBeta(param_RF['shape_M13_a'], nMultilocus)
	shape_M13_b  = randomBeta(param_RF['shape_M13_b'], nMultilocus)
	shape_M31_a  = randomBeta(param_RF['shape_M13_a'], nMultilocus)
	shape_M31_b  = randomBeta(param_RF['shape_M13_b'], nMultilocus)
	shape_M23_a  = randomBeta(param_RF['shape_M13_a'], nMultilocus)
	shape_M23_b  = randomBeta(param_RF['shape_M13_b'], nMultilocus)

	### present and past migration rates
	if sys.argv[1] == 'AB':
		M12_current = M12 # uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
		M21_current = M21 # uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
		M12_past = [0] * nMultilocus
		M21_past = [0] * nMultilocus
	else:
		M12_current = [0] * nMultilocus
		M21_current = [0] * nMultilocus
		M12_past = M12 # uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
		M21_past = M21 # uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)

	if sys.argv[2] != 'abc':
		M13_current = M13 # uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
		M31_current = M31 # uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
		M23_current = M23 # uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
		M32_current = M32 # uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
		M13_past = [0] * nMultilocus
		M31_past = [0] * nMultilocus
		M23_past = [0] * nMultilocus
		M32_past = [0] * nMultilocus
	else:
		M13_current = [0] * nMultilocus
		M31_current = [0] * nMultilocus
		M23_current = [0] * nMultilocus
		M32_current = [0] * nMultilocus
		M13_past = M13 # uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
		M31_past = M31 # uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
		M23_past = M23 # uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
		M32_past = M32 # uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)

# Three pop version
# param monolocus: values that will be read by ms
priorfile = "N1\tN2\tN3\tN1a\tN2a\tN3a\tNa_12\tNa\tshape_N_a\tshape_N_b\tTdem_1\tTdem_2\tTdem_3\tTsplit_12\tTsplit\tTsmall_AB\tTsmall_ABC\tM12\tshape_M12_a\tshape_M12_b\tM21\tshape_M21_a\tshape_M21_b\tM13\tM31\tM23\tM32\tshape_M13_a\tshape_M13_b\n"

for sim in range(nMultilocus):
	#priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}\t{14:.5f}\t{15:.5f}\t{16:.5f}\t{17:.5f}\t{18:.5f}\t{19:.5f}\t{20:.5f}\t{21:.5f}\t{22:.5f}\t{23:.5f}\t{24:.5f}\t{25:.5f}\t{26:.5f}\t{27:.5f}\t{28:.5f}\t{29:.5f}\t{30:.5f}\t{31:.5f}\t{32:.5f}\t{33:.5f}\t{34:.5f}\t{35:.5f}\t{36:.5f}\t{37:.5f}\t{38:.5f}\t{39:.5f}\t{40:.5f}\n".format(N1[sim], N2[sim], N3[sim], N4[sim], Na_12[sim], Na_34[sim], Na[sim], shape_N_a[sim], shape_N_b[sim], Tsplit_12[sim], Tsplit_34[sim], Tsplit[sim], Tsc_13[sim], Tsc_24[sim], M12[sim], M21[sim], M13[sim], shape_M13_a[sim], shape_M13_b[sim], M31[sim], shape_M31_a[sim], shape_M31_b[sim], M14[sim], shape_M14_a[sim], shape_M14_b[sim], M41[sim], shape_M41_a[sim], shape_M41_b[sim], M23[sim], shape_M23_a[sim], shape_M23_b[sim], M32[sim], shape_M32_a[sim], shape_M32_b[sim], M24[sim], shape_M24_a[sim], shape_M24_b[sim], M42[sim], shape_M42_a[sim], shape_M42_b[sim], M34[sim], M43[sim])
	priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}\t{14:.5f}\t{15:.5f}\t{16:.5f}\t{17:.5f}\t{18:.5f}\t{19:.5f}\t{20:.5f}\t{21:.5f}\t{22:.5f}\t{23:.5f}\t{24:.5f}\t{25:.5f}\t{26:.5f}\t{27:.5f}\t{28:.5f}\n".format(N1[sim], N2[sim], N3[sim], N1a[sim], N2a[sim], N3a[sim], Na_12[sim], Na[sim], shape_N_a[sim], shape_N_b[sim], Tdem_1[sim], Tdem_2[sim], Tdem_3[sim], Tsplit_12[sim], Tsplit[sim], Tsmall_12[sim], Tsmall_13[sim], M12[sim], shape_M12_a[sim], shape_M12_b[sim], M21[sim], shape_M21_a[sim], shape_M21_b[sim], M13[sim], M31[sim], M23[sim], M32[sim], shape_M13_a[sim], shape_M13_b[sim])
	# vectors of size 'nLoci' containing parameters
	## effective sizes
	if "1N" in sys.argv[4]:
		N1_vec = [N1[sim]] * nLoci
		N2_vec = [N2[sim]] * nLoci
		N3_vec = [N3[sim]] * nLoci
		N1a_vec = [N1a[sim]] * nLoci
		N2a_vec = [N2a[sim]] * nLoci
		N3a_vec = [N3a[sim]] * nLoci
		Na_12_vec = [Na_12[sim]] * nLoci
		Na_vec = [Na[sim]] * nLoci
	
	if "2N" in sys.argv[4]:
		scalar_N = beta(shape_N_a[sim], shape_N_b[sim], size=nLoci)
		scalar_N = [ i/(shape_N_a[sim]/(shape_N_a[sim]+shape_N_b[sim])) for i in scalar_N ]
		N1_vec = [ N1[sim]*i for i in scalar_N ]
		N2_vec = [ N2[sim]*i for i in scalar_N ]
		N3_vec = [ N3[sim]*i for i in scalar_N ]
		N1a_vec = [ N1a[sim]*i for i in scalar_N ]
		N2a_vec = [ N2a[sim]*i for i in scalar_N ]
		N3a_vec = [ N3a[sim]*i for i in scalar_N ]
		Na_12_vec = [ Na_12[sim]*i for i in scalar_N ]
		Na_vec = [ Na[sim]*i for i in scalar_N ]
	
	## migration rates
	
	if "1M" in sys.argv[3]:
		M12_current_vec = [M12_current[sim]] * nLoci
		M21_current_vec = [M21_current[sim]] * nLoci

		M13_current_vec = [M13_current[sim]] * nLoci
		M31_current_vec = [M31_current[sim]] * nLoci

		M23_current_vec = [M23_current[sim]] * nLoci
		M32_current_vec = [M32_current[sim]] * nLoci
	
		M12_past_vec = [M12_past[sim]] * nLoci
		M21_past_vec = [M21_past[sim]] * nLoci

		M13_past_vec = [M13_past[sim]] * nLoci
		M31_past_vec = [M31_past[sim]] * nLoci

		M23_past_vec = [M23_past[sim]] * nLoci
		M32_past_vec = [M32_past[sim]] * nLoci
	
	if "2M" in sys.argv[3]:	
		scalar_M12 = beta(shape_M12_a[sim], shape_M12_b[sim], size = nLoci)
		scalar_M12 = [ i/(shape_M12_a[sim]/(shape_M12_a[sim]+shape_M12_b[sim])) for i in scalar_M12 ]
		
		scalar_M21 = beta(shape_M21_a[sim], shape_M21_b[sim], size = nLoci)
		scalar_M21 = [ i/(shape_M21_a[sim]/(shape_M21_a[sim]+shape_M21_b[sim])) for i in scalar_M21 ]
		
		scalar_M13 = beta(shape_M13_a[sim], shape_M13_b[sim], size = nLoci)
		scalar_M13 = [ i/(shape_M13_a[sim]/(shape_M13_a[sim]+shape_M13_b[sim])) for i in scalar_M13 ]
		
		scalar_M31 = beta(shape_M31_a[sim], shape_M31_b[sim], size = nLoci)
		scalar_M31 = [ i/(shape_M31_a[sim]/(shape_M31_a[sim]+shape_M31_b[sim])) for i in scalar_M31 ]
		
		scalar_M23 = scalar_M13 # strong assumption here -> to test
		scalar_M32 = scalar_M31 # strong assumption here -> to test
		
		M12_current_vec = [ M12_current[sim] * i for i in scalar_M12 ]
		M21_current_vec = [ M21_current[sim] * i for i in scalar_M21 ]
	
		M13_current_vec = [ M13_current[sim] * i for i in scalar_M13 ]
		M31_current_vec = [ M31_current[sim] * i for i in scalar_M31 ]
		
		M23_current_vec = [ M23_current[sim] * i for i in scalar_M13 ]
		M32_current_vec = [ M32_current[sim] * i for i in scalar_M31 ]
	
		M12_past_vec = [ M12_past[sim] * i for i in scalar_M12 ]
		M21_past_vec = [ M21_past[sim] * i for i in scalar_M21 ]
	
		M13_past_vec = [ M13_past[sim] * i for i in scalar_M13 ]
		M31_past_vec = [ M31_past[sim] * i for i in scalar_M31 ]
		
		M23_past_vec = [ M23_past[sim] * i for i in scalar_M13 ]
		M32_past_vec = [ M32_past[sim] * i for i in scalar_M31 ]
	for locus in range(nLoci):
		# java -jar msms.jar tbs 10000 -s tbs -r tbs tbs -I 3 tbs tbs tbs 0 -n 1 tbs -en tbs 1 tbs -n 2 tbs -en tbs 2 tbs -n 3 tbs -en tbs 3 tbs -m 1 2 tbs -m 2 1 tbs -m 1 3 tbs -m 3 1 tbs -m 2 3 tbs -m 3 2 tbs -em tbs 1 2 tbs -em tbs 2 1 tbs -em tbs 1 3 tbs -em tbs 3 1 tbs -em tbs 2 3 tbs -em tbs 3 2 tbs -ej tbs 2 1 -en tbs 1 tbs -ej tbs 3 1 -eN tbs tbs
		print('{0}\t{1}\t{2:.5f}\t{3}\t{4}\t{5}\t{6}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}\t{14:.5f}\t{15:.5f}\t{16:.5f}\t{17:.5f}\t{18:.5f}\t{19:.5f}\t{20:.5f}\t{21:.5f}\t{22:.5f}\t{23:.5f}\t{24:.5f}\t{25:.5f}\t{26:.5f}\t{27:.5f}\t{28:.5f}\t{29:.5f}\t{30:.5f}\t{31:.5f}\t{32:.5f}\t{33:.5f}\t{34:.5f}\t{35:.5f}\t{36:.5f}\t{37:.5f}\t{38:.5f}\t{39:.5f}'.format(nsam_tot[locus], int(nSNPs[locus]), rho[locus], int(L[locus]), nsamA[locus], nsamB[locus], nsamC[locus], N1_vec[locus], Tdem_1[sim], N1a_vec[locus], N2_vec[locus], Tdem_2[sim], N2a_vec[locus], N3_vec[locus], Tdem_3[sim], N3a_vec[locus], M12_current_vec[locus], M21_current_vec[locus], M13_current_vec[locus], M31_current_vec[locus], M23_current_vec[locus], M32_current_vec[locus], Tsmall_12[sim], M12_past_vec[locus], Tsmall_12[sim], M21_past_vec[locus], Tsmall_13[sim], M13_past_vec[locus], Tsmall_13[sim], M31_past_vec[locus], Tsmall_13[sim], M23_past_vec[locus], Tsmall_13[sim], M32_past_vec[locus], Tsplit_12[sim], Tsplit_12[sim], Na_12_vec[locus], Tsplit[sim], Tsplit[sim], Na_vec[locus]))
outfile = open("priorfile.txt", "w")
outfile.write(priorfile)
outfile.close()


