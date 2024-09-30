import msprime as msp
import numpy as np
import scipy.stats as st
from plinkio import plinkfile as pl
import os
import glob
import random


# Source pops related by (S1,(S2,S3))
# T is target pop
# T is initially a sister to S2, but then recieves contribution from S1 and S3
# with no admixture, (S1, ((S2, T), S3))
# O are nonadmixed groups that split off of each branch
# NB: O split off more anciently than the most ancient event on their branch (either node, sample, or admixture)
# nS vector of sample sizes in each S, i.e. (nS1, nS2, nS3)
# nT sample size of target population
# nO sample size of unadmixed groups (assumed to be the same in all of them)
# tS vector ages of samples in each S, i.e. (tS1, tS2, tS3)
# ??: does it matter if tSi is more ancient than the admixture from that branch?
# tT sampling time for target population
# NB: unadmixed groups are assumed to all be modern samples
# t2T time of T,S2 split
# t23 time of S2,S3 split
# t123 time of S1, (S2, S3) split
# ta vector of admixture times from S1 and S3 to T, i.e. (ta1T, ta3T)
# f vector of admixture fractions from S1 and S3 to T, i.e. (f1T, f3T)
# tsO vector of split times of the O as a fraction of the branch they're on
# tsOi = 0 -> as recent as possible, tsOi = 1 -> as ancient as possible
# N vector of effective sizes in order of population numbers
# fs1s2 admixture proportion between S1 and S2, happens as recently as possible
# fo1o2 admixture proportion between O1 and O2, happens as recently as possible
# SIMULATES DIPLOIDS
def sim_admix(nS=[10,10,10],nT=10,nO=10,tS=[0,0,0],tT=0,t2T=100,t23=200,t123=300,ta=[50,50],f=[.1,.1],tsO=[0.1,0.1,0.1,0.1], N = [10000]*8,mu=1.25e-8,r=1e-8,length=1000000,num_rep=1,coverage=10, err = None, seed=None, fs1s2 = None, fo1o2 = None):
	#generate samples
	samples = set_up_pops(tT,tS,nT,nS,nO)
	#generate demography
	demography = set_up_demography(t2T, t23, t123, ta, tsO, f)
	#add admixture between sources if necessary
	if fs1s2 is not None: demography = add_source_admixture(demography, fs1s2, tS)
	#add admixture between outgroups if necessary
	if fo1o2 is not None: demography = add_outgroup_admixture(demography, fo1o2)
	#make some population configurations
	pops = [msp.PopulationConfiguration(initial_size = n) for n in N]
	#simulate some data
	sims = msp.simulate(samples=samples,Ne=N[0],population_configurations=pops, demographic_events = demography, mutation_rate = mu, length=length, num_replicates = num_rep, recombination_rate = r, random_seed = seed)

	#return whole simulations as an object

	return sims
	


def set_up_pops(tT, tS, nT, nS, nO):
	#target
	samples = [msp.Sample(population=0, time = tT)]*(2*nT) #pop T
	#sources
	samples.extend([msp.Sample(population=1, time = tS[0])]*(2*nS[0])) #pop S1
	samples.extend([msp.Sample(population=2, time = tS[1])]*(2*nS[1])) #pop S2
	samples.extend([msp.Sample(population=3, time = tS[2])]*(2*nS[2])) #pop S3
	#outgroups
	samples.extend([msp.Sample(population=4, time = 0)]*(2*nO)) #pop O1 (outgroup of S1)
	samples.extend([msp.Sample(population=5, time = 0)]*(2*nO)) #pop O2 (outgroup of (T, S2))
	samples.extend([msp.Sample(population=6, time = 0)]*(2*nO)) #pop  O3 (outgroup of S3)
	samples.extend([msp.Sample(population=7, time = 0)]*(2*nO)) #pop O4 (outgroup of ((T,S2),S3))
	return samples

def set_up_demography(t2T, t23, t123, ta, tsO, f):
	#divergence of source populations
	source_divergence = [msp.MassMigration(time=t2T,source=2,destination=0,proportion=1), # (T, S2)
				msp.MassMigration(time=t23,source=3,destination=0,proportion=1), # ((T, S2), S3)
				msp.MassMigration(time=t123,source=1,destination=0,proportion=1)] # (S1, ((T, S2), S3))
	#divergence of outgroup populations
	o = np.array(tsO)*np.array([t123-ta[0],t23-t2T,t23-ta[1],t123-t23]) #time along branch
	o += np.array([ta[0],t2T,ta[1],t23]) #starting point
	outgroup_divergence = [msp.MassMigration(time=o[0], source = 4, destination = 1, proportion = 1), # merge of O1
				msp.MassMigration(time=o[1], source = 5, destination = 0, proportion = 1), # merge O2
				msp.MassMigration(time=o[2], source = 6, destination = 3, proportion = 1), # merge O3
				msp.MassMigration(time=o[3], source = 7, destination = 0, proportion = 1)] #merge O4
	#admixture times
	#NB: fraction from S2 to T is implicitly 1 - f1T - f2T I think
	admix = [msp.MassMigration(time = ta[0], source = 0, destination = 1, proportion = f[0]), #fraction from S1 to T
			msp.MassMigration(time = ta[1], source = 0, destination = 3, proportion = f[1])] #fraction from S3 to T 
	#combine and sort the demography
	demography = source_divergence + outgroup_divergence + admix
	return sorted(demography, key = lambda x: x.time)

def add_source_admixture(demography, fs1s2, tS):
	ts1s2 = min(tS[:2])
	source_admix = msp.MassMigration(time = ts1s2, source = 1, destination = 2, proportion = fs1s2)
	demography.append(source_admix)
	return sorted(demography, key = lambda x: x.time)

def add_outgroup_admixture(demography, fo1o2):
	to1o2 = 0
	outgroup_admix = msp.MassMigration(time = to1o2, source = 4, destination = 5, proportion = fo1o2)
	demography.append(outgroup_admix)
	return sorted(demography, key = lambda x: x.time)


def draw_reads(GT, coverage, err):
	num_reads = st.poisson.rvs(mu=coverage,size=len(GT))
	p_der = GT/2.*(1-err)+(1-GT/2.)*err
	derived_reads = st.binom.rvs(num_reads,p_der)
	return (derived_reads,num_reads-derived_reads)

def generate_plink_bed(names,outname,pos,alleles,GT=None,reads=None, r = 1e-8, duplicate = False):
	if GT is None and reads is None:
		raise StandardError("Need to provide either GT or reads")
	elif GT is not None and reads is not None:
		raise StandardError("Must provide exactly one of GT or reads")
	elif alleles is None:
		raise StandardError("Alleles must be specified")
	elif GT is not None:
		generate_plink_bed_GT(names,GT,pos,alleles,outname, r = r, duplicate = duplicate)
	elif reads is not None:
		generate_plink_bed_reads(names,reads,pos,alleles,outname, r = r, duplicate = duplicate)

def generate_plink_bed_GT(names,GT,pos,alleles,outname , r = 1e-8, duplicate = False):
	pos_included = set()
	#create sample objects
	samples = [pl.Sample(0,name,0,0,0,0) for name in names]
	outfile = pl.create(outname,samples)
	#write the file
	for i, site in enumerate(GT):
		cur_label = "chr%i_%i"%(pos[i][0],pos[i][1]+1)
		if cur_label not in pos_included:
			cur_ind = 1
			pos_included.add(cur_label)
		elif cur_label in pos_included and duplicate:
			cur_label += "_%i"%cur_ind
			cur_ind += 1
		elif cur_label in pos_included and not duplicate:
			continue
		cur_locus = pl.Locus(pos[i][0],cur_label,r*pos[i][1],pos[i][1]+1,alleles[i][0],alleles[i][1])
		outfile.write_row(cur_locus,site)
	outfile.close()

def generate_plink_bed_reads(names,reads,pos,alleles,outname, r = 1e-8, duplicate = False):
	pos_included = set()
	#create sample objects
	samples = [pl.Sample(fid=0,iid=name,father_iid=0,mother_iid=0,sex=0,affection=5,phenotype=float(name.split("_")[1])) for name in names]
	outfile = pl.create(outname,samples)
	#write the file
	for i, site in enumerate(reads):
		cur_label = "chr%i_%i"%(pos[i][0],pos[i][1]+1)
		if cur_label not in pos_included:
			cur_ind = 1
			pos_included.add(cur_label)
		elif cur_label in pos_included and duplicate:
			cur_label += "_%i"%cur_ind
			cur_ind += 1
		elif cur_label in pos_included and not duplicate:
			continue
		cur_locus = pl.Locus("chr_%s" %pos[i][0],cur_label,r*pos[i][1],pos[i][1]+1,alleles[i][0],alleles[i][1])
		cur_reads = sample_random_read(site)
		outfile.write_row(cur_locus,cur_reads)
	outfile.close()

	
def sample_random_read(reads):
	total = reads[0]+reads[1]
	total[total==0] = -1 #hack to avoid divide by zero
	prob_derived = (reads[0]+0.0)/total
	random_read = 2*st.binom.rvs(1,prob_derived)
	random_read[total==-1] = 3 #set missing ones to actually be missing!
	return random_read	

def generate_beagle_GL(names,reads,err,pos,alleles,outname,duplicate=False, precision=3):
	pos_included = set()
	outfile = open("%s.gl"%outname,"w")
	outfile.write("marker\tallele1\tallele2\t")
	header = '\t'.join(['\t'.join([name]*3) for name in names])
	outfile.write("%s\n"%header)
	for i, read in enumerate(reads):
		cur_label = "chr%i_%i"%(pos[i][0],pos[i][1]+1)
		if cur_label not in pos_included:
			cur_ind = 1
			pos_included.add(cur_label)
		elif cur_label in pos_included and duplicate:
			cur_label += "_%i"%cur_ind
			cur_ind += 1
		elif cur_label in pos_included and not duplicate:
			continue
		GL = map(lambda x: str(round(x, precision)),compute_GL(read,err))
		outfile.write("%s\t%s\t%s\t"%(cur_label,alleles[i][0],alleles[i][1]))
		outfile.write("%s\n"%'\t'.join(GL))
	outfile.close()		

def compute_GL(reads,err):
	total = np.sum(reads,axis=0)
	GT = np.array([0,1,2])
	p_err = GT/2.*(1-err)+(1-GT/2.)*err
	GL0 = st.binom.pmf(reads[0],total,p_err[0])
	GL1 = st.binom.pmf(reads[0],total,p_err[1])
	GL2 = st.binom.pmf(reads[0],total,p_err[2])
	GL = np.vstack((GL0,GL1,GL2))
	GL /= np.sum(GL, axis=0)
	GL = GL.flatten("F")
	return GL

def randomDNA(length,gc_content=0.41):
	at_content = 1.0 - gc_content
	AT = int(round(100 * at_content))
	GC = int(round(100 * gc_content))
	return ''.join(random.choice("A"*AT+"C"*GC+"G"*GC+"T"*AT) for _ in xrange(length))

def mutate(base, no_transition = True):
	#assume all mutations are equally likely and give the option to exclude transition mutations
	nucs = ["A","C","G","T"]
	nucs.remove(base)
	if no_transition:
		if base == "A":
			nucs.remove("G")
		elif base == "C":
			nucs.remove("T")
		elif base == "G":
			nucs.remove("A")
		elif base == "T":
			nucs.remove("C")
	return random.choice(nucs)
	

