# -*- coding: utf-8 -*-
################################################################################
##  *
##  *  Function: single-species miRNA identification [Run miRDeep2 for multiple samples]
##  *  Writer:Siyuan Feng
##  *  Mail: siyuanfeng.bioinfo@gmail.com
##  *  Version:19.11.5
##  *
################################################################################

import os
from optparse import OptionParser

### help and usage ### 
usage = "usage: %prog [options] args"
description = "single-species miRNA identification."
version = '%prog 19.11.5'
parser = OptionParser(usage=usage,version=version, description = description)
parser.add_option("--species",
                    action="store",
                    dest="species",
                    help="name of the selected species, for example, 'pig'",
                    metavar = 'NAME')         
parser.add_option("--species_short",
                    action="store",
                    dest="species_short",
                    help="shortname of the selected species, for example, 'ssc'",
                    metavar = 'NAME')
parser.add_option("--queue",
                    action="store",
                    dest="queue",
                    help="name of the queue where you submit your task on the server ",
                    metavar = 'NAME')                 
parser.add_option("--genome",
                    action="store",
                    dest="genome",
                    help="path of the fasta file of selected genome\nPlease note, the identifier line must be like '>chr_1'. Delete the space-delimited information after the chromosome ID",
                    metavar = 'PATH')  
parser.add_option("--bowtie_index",
                    action="store",
                    dest="bowtie_index",
                    help="path of the bowtie index built for a selected genome",
                    metavar = 'PATH') 
parser.add_option("--work_dir",
                    action="store",
                    dest="work_dir",
                    help="path of your work direction",
                    metavar = 'PATH')
parser.add_option("--sample_list",
                    action="store",
                    dest="sample_list",
                    help="path of a tab-delimited list of sample name (first column) and sample path (second column)",
                    metavar = 'PATH')                         
parser.add_option("--mirbase_mature",
                    action="store",
                    dest="mirbase_mature",
                    help="path of the fasta file of miRBase mature sequence database (mature_format.fa) including all speceis (latest version: miRBase release 22). Please don't change this file (unless miRBase is updated).". ,
                    metavar = 'PATH')  
parser.add_option("--mirbase_hairpin",
                    action="store",
                    dest="mirbase_hairpin",
                    help="path of the fasta file of miRBase hairpin sequence database (hairpin_format.fa) including all speceis (latest version: miRBase release 22). Please don't change this file (unless miRBase is updated).",
                    metavar = 'PATH')                      
parser.add_option("--mirbase_mature_ref",
                    action="store",
                    dest="mirbase_mature_ref",
                    help="path of the fasta file of miRBase mature sequence database (mammal_mature.fa) including reference speceis (including the selected species). This file is used to idendify conserve microRNA",
                    metavar = 'PATH')                      
parser.add_option("--template_mapping",
                    action="store",
                    dest="template_mapping",
                    help="path of template perl script (template_mapping.pl) for mapping. Please don't change this file.",
                    metavar = 'PATH')                     
parser.add_option("--extract_mirna",
                    action="store",
                    dest="extract_mirna",
                    help="path of the perl script (extract_mirna.pl) to extract sequence from miRBase database. Please don't change this file.",
                    metavar = 'PATH')                     
                                                          
(options,args) = parser.parse_args()
#################

### parameters ###
species = options.species
species_short = options.species_short
queue = options.queue
genome = options.genome
bowtie_index = options.bowtie_index
work_dir = options.work_dir
sample_list = options.sample_list
mirbase_hairpin = options.mirbase_hairpin
mirbase_mature = options.mirbase_mature
mirbase_mature_ref = options.mirbase_mature_ref
template_mapping = options.template_mapping
extract_mirna = options.extract_mirna
#################

### function to extract microRNA from miRBase databse ###
def mirbase(spec, short, dir):
	hairpin = mirbase_hairpin
	mature = mirbase_mature
	x = open(mirbase_mature_ref, 'r').readlines()
	for line in x:
		if line.startswith('>' + short):
			x = x[0:x.index(line)] + x[x.index(line) + 2:]
	open(mirbase_mature_ref, 'r').close()

	f = open('{0}/{1}/mirbase/others_mature.fa'.format(dir, spec), 'w')
	for line in x:
		f.write(line)
	f.close()

	os.system('perl {2} {0} {1}'.format(hairpin, mature, extract_mirna) + ''.join(
		' {0}'.format(short)))
	os.system('mv {0}/{2}_hairpin.fa  {0}/{1}/mirbase/{2}_hairpin.fa'.format(dir, spec, short))
	os.system('mv {0}/{2}_mature.fa  {0}/{1}/mirbase/{2}_mature.fa'.format(dir, spec, short))
	return None
#################

### function to systematically identify 3' adaptor of small RNA-Seq data ###
def find_adaptor(fastq):
	f=open(fastq)
	
	##parse fastq file
	
	reads=[]
	x=f.readline()
	while x!='':
	    reads.append(f.readline().strip())
	    f.readline()
	    f.readline()
	    x=f.readline()
	
	##find k_mers
	
	k_mers={}
	k_pos={}
	for read in reads[:10000]:
	    for pos in range(0,len(read)):
	        k_mers[read[pos:pos+6]]=k_mers.get(read[pos:pos+6],0)+1
	        k_pos[read[pos:pos+6]]=k_pos.get(read[pos:pos+6],0)+pos
	
	##convert dictionary to list and sort
	
	mylist=zip(k_mers.values(),k_mers.keys())
	mylist.sort(reverse=True)
	
	
	
	##print most abundant k_mers
	adapter=""
	tmp_pos=9999999999999999
	top=mylist[0][0]
	temp_i=-1
	for i in mylist[:50]:
	    if k_pos[i[1]]<tmp_pos and i[0]*1.5 >top:
	        tmp_pos=k_pos[i[1]]
	        adapter=i[1]
	        temp_i=i
	
	if temp_i[0]*2>10000:
	    return adapter
	else:
	    return 'NNNNNN'
#################

os.makedirs('{0}/{1}/mirbase'.format(work_dir, species))
mirbase(species, species_short, work_dir)
########

for line in open(sample_list, 'r'):
	line = line.strip().split('\t')
	sample = line[0]
	rawdata = line[1]
	adaptor = find_adaptor(rawdata)
	myqsub = '{0}.pbs'.format(sample)
	f = open(myqsub, 'w')
	print >> f, '''#!/bin/sh
#PBS -N {0}
#PBS -l nodes=1:ppn=4
#PBS -l mem=2gb
#PBS -q {6}
echo start at time `date +%F'  '%H:%M`
#mapping
sample={0}
species={1}
name={2}
dir={3}

cp {7} $dir/$species/$sample\_mapping.pl
sed -i "s?dir_bowtie_index?{8}?g" $dir/$species/$sample\_mapping.pl
sed -i "s?NNNNNN?{9}?g" $dir/$species/$sample\_mapping.pl

cd $dir/$species
perl $sample\_mapping.pl -i {4} -a $sample

#expression
mkdir $dir/$species/$sample/identification/
cd $dir/$species/$sample/identification/
reads=$dir/$species/$sample/$sample\.fa
genome={5}
arf=$dir/$species/$sample/$sample\_vs_genome.arf
this_species_mature=$dir/$species/mirbase/$name\_mature.fa
ohthers_mature=$dir/$species/mirbase/others_mature.fa
this_species_hairpin=$dir/$species/mirbase/$name\_hairpin.fa
miRDeep2.pl $reads $genome $arf $this_species_mature $ohthers_mature  $this_species_hairpin -d 2> report.log
echo finish at time `date +%F'  '%H:%M`'''.format(sample, species, species_short, work_dir, rawdata, genome, queue, template_mapping, bowtie_index, adaptor)
	f.close()
	os.system('qsub {0}'.format(myqsub))
open(sample_list, 'r').close()



