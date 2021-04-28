# -*- coding: utf-8 -*-
################################################################################
##  *
##  *  Function: single-species miRNA profiling [Extract microRNA expression for multiple samples]
##  *  Writer:Siyuan Feng
##  *  Mail: siyuanfeng.bioinfo@gmail.com
##  *  Version:19.11.5
##  *
################################################################################

import csv
import os
import sys
import math
import Levenshtein
import types
import itertools
from optparse import OptionParser

### help and usage ### 
usage = "usage: %prog [options] args"
description = "single-species miRNA profiling based on output of identification.py"
version = '%prog 19.11.5'
parser = OptionParser(usage=usage,version=version, description = description)
parser.add_option("--expr_cutoff",
                    action="store",
                    dest="expr_cutoff",
                    help="Cutoff of readcount. Cutoff = 1 means we discard microRNA with readcount < 1 in all samples. Default = 3",
                    metavar = 'CUTOFF',
                    type = 'int',
                    default = 3) 
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
parser.add_option("--work_dir",
                    action="store",
                    dest="work_dir",
                    help="path of your work direction",
                    metavar = 'PATH')
parser.add_option("--work_dir_identification",
                    action="store",
                    dest="work_dir_identification",
                    help="path of your work direction of identification (first step)",
                    metavar = 'PATH')                    
parser.add_option("--sample_list",
                    action="store",
                    dest="sample_list",
                    help="path of a tab-delimited list of sample name (first column) and sample path (second column)",
                    metavar = 'PATH')                                                                                   
(options,args) = parser.parse_args()
#############################

# def a function for getting the middle part of a string
def GetMiddleStr(content, startStr, endStr):
	startIndex = content.index(startStr)
	if startIndex >= 0:
		startIndex += len(startStr)
	endIndex = content.index(endStr, startIndex)
	return content[startIndex:endIndex]
#############################

### parameters ###
species = options.species
species_short = options.species_short
expr_cutoff = options.expr_cutoff
work_dir = options.work_dir
work_dir_identification = options.work_dir_identification
sample_list = options.sample_list
###################

#######################part 1 to extract expression of each sample within one species##############
# parent_path = '{0}/{1}'.format(species_path, mismatch)
# os.mkdir(parent_path)
# name_exchange = '/lustre/fengsiyuan/data/01_longkeren/all_177samples/graph_analysis/quantification/ortholog/test1/name_exchange.csv'
dir_map = '{0}/each_mappingrate'.format(work_dir)
dir_result = '{0}/all_results'.format(work_dir)
os.mkdir('{0}'.format(dir_map))
os.mkdir('{0}'.format(dir_result))
#############################

### prepare for file reading ###
# must put the reference species at the first place of list "shortname" ###

# varname = []
corres = {}
allmir_species_known = {}
allmir_species_new = {}
allmir_collapsed = {}

corres.setdefault(species_short, [])
allmir_species_known.setdefault(species_short, {})
allmir_species_new.setdefault(species_short, {})

# # samplename exchange
# name_relation = {}
# with open(name_exchange, 'r') as f:
# 	for line in f.readlines():
# 		name = line.strip().split(',')
# 		name_relation[name[0]] = name[1] + '_' + name[3] + name[2]
# f.close()

# trans = {}
# for key, value in name_relation.items():
# 	trans[value] = key

### construct the known mirna set for each species, according to mirbase mature mirna
spec_no_mirbase = []
mirbase = open('{0}/{1}/mirbase/{2}_hairpin.fa'.format(work_dir_identification, species, species_short), 'r')
real_spec = species_short
for mir in mirbase:
	if mir.startswith('>'):
		mirname = '%s-%s-%s' % (real_spec, mir.rstrip().split('-')[1], mir.rstrip().split('-')[2])
		allmir_species_known[real_spec][mirname] = {}  # allmir_species_known={'bta':{'bta-mir-199a':{},bta-let-7a:{}},'gga':{{}},,,,}

# input all expression results
for line in open(sample_list, 'r'):
	line = line.strip().split('\t')
	line = line[0]
	new_name = line
	# varname.append(line)
	corres[species_short].append(new_name)
	os.system('cat {0}/{3}/{1}/identification/result*.csv > {2}/{4}.csv '.format(work_dir_identification, line, dir_result, species, new_name))
	l = {}
	f = open("{0}/{1}.csv".format(dir_result, new_name), 'r')
	x = f.readlines()
	f.close()

	# get the sum of cleandata readcount for each sample
	os.system('cat {0}/{1}.e* > {2}/{3}.txt'.format(work_dir_identification, line, dir_map, new_name))
	f = open('{0}/{1}.txt'.format(dir_map, new_name), 'r')
	a = f.readlines()
	info = a[a.index('Mapping statistics\n') + 3]
	total = float(info.split()[2])
	f.close()

	# read known mirna in dict
	for i in x[x.index('mature miRBase miRNAs detected by miRDeep2\n') + 2:x.index('#miRBase miRNAs not detected by miRDeep2\n') - 1]:
		mir = i.split()[11].split('-')
		unique_number = i.split()[0].split('_')[1]
		mir_id = '%s-%s-%s-%s-%s' % (new_name, species_short, mir[1].replace('miR', 'mir'), mir[2], unique_number)  # cattle_lung1-bta-mir-199a-22755
		mir_id_short = '%s-%s-%s' % (species_short, mir[1].replace('miR', 'mir'), mir[2])
		if mir_id_short in allmir_species_known[species_short]:
			pass
		else:
			for each in allmir_species_known[species_short]:
				if mir_id_short in each:
					mir_id_short = each

		sequence = i.split()[17]
		coordin = i.split()[18]
		total_expression = i.split()[6]
		mature_expression = i.split()[7]
		star_expression = i.split()[9]
		mature_expression_normalized = str(math.log(int(i.split()[7]) * 1000000 / total + 1, 2))
		star_expression_normalized = str(math.log(int(i.split()[9]) * 1000000 / total + 1, 2))
		mature_sequence = i.split()[15]
		star_sequence = i.split()[16]

		if sequence.index(mature_sequence) < sequence.index(star_sequence):
			allmir_species_known[species_short][mir_id_short][mir_id] = '%s#%s#%s#%s#%s#%s#%s#%s#%s' % (sequence, coordin, mature_expression, mature_expression_normalized, star_expression, star_expression_normalized, mature_sequence, star_sequence, total_expression)
		else:
			allmir_species_known[species_short][mir_id_short][mir_id] = '%s#%s#%s#%s#%s#%s#%s#%s#%s' % (sequence, coordin, star_expression, star_expression_normalized, mature_expression, mature_expression_normalized, star_sequence, mature_sequence, total_expression)


	# read conserved and novel mirna in dict
	identifier = 1
	for i in x[x.index('novel miRNAs predicted by miRDeep2\n') + 2:x.index('mature miRBase miRNAs detected by miRDeep2\n') - 3]:
		if i.split()[5]!='-' or i.split()[10] == 'no':
			continue
		if i.split()[12] == '-':
			mir_id = species_short + '-' + new_name + '-' + 'mir-' + str(identifier) + '-novel'  # bta-C1_L-mir-1-novel
			identifier += 1
		else:
			if i.split()[12].split('-')[1] == 'let':
				mir_id = species_short + '-' + new_name + '-' + 'let-' + i.split()[12].split('-')[2] + '-conserve'
			else:
				mir_id = species_short + '-' + new_name + '-' + 'mir-' + i.split()[12].split('-')[
					2] + '-conserve'  # bta-C1_L-mir-18a-conserve

		sequence = i.split()[17]
		coordin = i.split()[18]
		mature_expression = i.split()[7]
		star_expression = i.split()[9]
		mature_expression_normalized = str(math.log(int(i.split()[7]) * 1000000 / total + 1, 2))
		star_expression_normalized = str(math.log(int(i.split()[9]) * 1000000 / total + 1, 2))
		mature_sequence = i.split()[15]
		star_sequence = i.split()[16]

		if sequence.index(mature_sequence) < sequence.index(star_sequence):
			l[mir_id] = '%s#%s#%s#%s#%s#%s#%s#%s' % (
			sequence, coordin, mature_expression, mature_expression_normalized, star_expression,
			star_expression_normalized, mature_sequence, star_sequence)
		else:
			l[mir_id] = '%s#%s#%s#%s#%s#%s#%s#%s' % (
			sequence, coordin, star_expression, star_expression_normalized, mature_expression,
			mature_expression_normalized, star_sequence, mature_sequence)

		# put annotions from all samples into one dict:'l'
	allmir_species_new[species_short].update(l)
	
# collapse mirs with very close position

allmir_collapsed.setdefault(species_short, {})
j = {}
k = {}
for key in allmir_species_new[species_short]:
	j[key] = allmir_species_new[species_short][key]
for mir1 in j:
	for mir2 in j:
		coordin1 = j[mir1].split('#')[1].split(':')
		coordin2 = j[mir2].split('#')[1].split(':')
		sequence1_5p = j[mir1].split('#')[6]
		sequence2_5p = j[mir2].split('#')[6]
		sequence1_3p = j[mir1].split('#')[7]
		sequence2_3p = j[mir2].split('#')[7]

		if coordin1[0] == coordin2[0] and coordin1[2] == coordin2[2] and abs(int(coordin1[1].split('..')[0]) - int(coordin2[1].split('..')[0])) <= 40 and abs(int(coordin1[1].split('..')[1]) - int(coordin2[1].split('..')[1])) <= 40:
			j[mir2] = j[mir1]
		elif Levenshtein.ratio(sequence1_5p, sequence2_5p) >= 0.8 or Levenshtein.ratio(sequence1_3p,sequence2_3p) >= 0.8 or Levenshtein.ratio(sequence1_5p, sequence2_3p) >= 0.8 or Levenshtein.ratio(sequence1_3p, sequence2_5p) >= 0.8:
			j[mir2] = j[mir1]

for key, value in j.items():
	k[value] = key

allmir_collapsed[species_short].update(k)  # set of all mirrna after collapsing

# extract
miR = {}
miR.setdefault(species_short, {})
for mir in allmir_species_known[species_short]:
	if allmir_species_known[species_short][mir] != {}:
		miR[species_short][mir] = []


for key in allmir_collapsed[species_short]:
	name = allmir_collapsed[species_short][key]
	miR[name[0:3]][name] = []
###miR={'bta':{'bta-mir-7i':[[sample1],[sample2]]}}


for x in miR[species_short]:
	# novel and conserve mirna
	if x.endswith('novel') or x.endswith('conserve'):
		for sample in corres[species_short]:  # samples of one species
			thissample = []
			for mir in allmir_species_new[species_short]:
				if sample in mir:
					items1 = allmir_species_new[species_short][x].split('#')
					coordin1 = items1[1].split(':')
					items2 = allmir_species_new[species_short][mir].split('#')
					coordin2 = items2[1].split(':')
					# for sequence comparing, to reduce false negative
					sequence1_5p = items1[6]
					sequence2_5p = items2[6]
					sequence1_3p = items1[7]
					sequence2_3p = items2[7]
					precursor2 = items2[0]

					if coordin1[0] == coordin2[0] and coordin1[2] == coordin2[2] and abs(
									int(coordin1[1].split('..')[0]) - int(
									coordin2[1].split('..')[0])) <= 40 and abs(
									int(coordin1[1].split('..')[1]) - int(coordin2[1].split('..')[1])) <= 40:
						thissample = [items2[0], items2[1], items2[2], items2[3], items2[4], items2[5], items2[6],
									  items2[7]]
						already_extracted = mir
						break
					elif Levenshtein.ratio(sequence1_5p, sequence2_5p) >= 0.8 or Levenshtein.ratio(sequence1_3p,sequence2_3p) >= 0.8 or Levenshtein.ratio(sequence1_5p, sequence2_3p) >= 0.8 or Levenshtein.ratio(sequence1_3p, sequence2_5p) >= 0.8:
						thissample = [items2[0], items2[1], items2[2], items2[3], items2[4], items2[5], items2[6],
									  items2[7]]
						already_extracted = mir
					# make sure whether the mature/star_sequence is equal to 5p/3p-sequence,if not,just exchange them
			if thissample != []:
				if x != already_extracted and not already_extracted in miR[species_short].keys():
					del allmir_species_new[species_short][already_extracted]
			else:
				thissample = [items1[0], items1[1], '0', '0', '0', '0', items1[6], items1[7]]
			miR[species_short][x].append(thissample)


		# known mirna
	else:
		for sample in corres[species_short]:
			thissample = []
			optional_mir = {}
			for mir in allmir_species_known[species_short][x]:
				if sample in mir:
					optional_mir[mir] = int(allmir_species_known[species_short][x][mir].split('#')[8])
				else:
					continue
			if optional_mir == {}:
				thissample = ['null'] * 2 + ['0'] * 4 + ['null'] * 2
			else:
				selected_mir = max(optional_mir, key=optional_mir.get)
				for ele in allmir_species_known[species_short][x][selected_mir].split('#')[:-1]:
					thissample.append(ele)
			miR[species_short][x].append(thissample)

# write all infomation into one csv
varname_inorder = []
mir_selected = {}
varname_inorder = corres[species_short]
mir_selected.setdefault(species_short, [])
headers1 = ['miR'] + ['precursor_sequence'] + [''] * (len(varname_inorder) - 1) + ['precursor_coordinate'] + [''] * (len(varname_inorder) - 1) + [
			   'mature_expression'] + [''] * (len(varname_inorder) - 1) + ['mature_expression_normalized'] + [ ''] * ( len(varname_inorder) - 1) + [
			   'mature_sequence'] + [''] * (len(varname_inorder) - 1)
headers2 = ['miR'] + varname_inorder * 5
pre_num = {}
with open('{0}/expression_{1}_{2}.csv'.format(work_dir, expr_cutoff, species_short), 'wb') as f:
	f_csv = csv.writer(f)
	f_csv.writerow(headers1)
	f_csv.writerow(headers2)
	newname_conserve = {'-'.join(['conserve'] + [key.split('-')[0], key.split('-')[2], key.split('-')[3]]): 0 for key in miR[species_short] if 'conserve' in key}
	id_novel = 1
	for key in miR[species_short]:
		pre_num[key] = []
		elements = key.split('-')
		if 'conserve' in key:
			elements = key.split('-')
			newname = '-'.join(['conserve'] + [elements[0], elements[2], elements[3]])
			if newname_conserve[newname] == 0:
				newname_conserve[newname] += 1
				pass
			else:
				newname1 = newname + '-' + str(newname_conserve[newname])
				newname_conserve[newname] += 1
				newname = newname1

		elif 'novel' in key:
			elements = key.split('-')
			newname = '-'.join(['novel'] + [elements[0], elements[2], str(id_novel)])
			id_novel += 1
		else:
			newname = key


		row1 = [newname.replace('mir', 'miR') + '-5p']
		row2 = [newname.replace('mir', 'miR') + '-3p']
		for x1 in [0, 1, 2, 3, 6]:
			for lists in miR[species_short][key]:
				row1 += [lists[x1]]
		for x2 in [0, 1, 4, 5, 7]:
			for lists in miR[species_short][key]:
				row2 += [lists[x2]]

		expression1 = []
		expression2 = []
		for i in row1[len(varname_inorder) * 2 + 1:len(varname_inorder) * 3 + 1]:
			expression1.append(int(i))
		if max(expression1) >= expr_cutoff:
			f_csv.writerow(row1)
			mir_selected[species_short].append(key)
			pre_num[key].append('5p')

		for i in row2[len(varname_inorder) * 2 + 1:len(varname_inorder) * 3 + 1]:
			expression2.append(int(i))
		if max(expression2) >= expr_cutoff:
			f_csv.writerow(row2)
			mir_selected[species_short].append(key)
			pre_num[key].append('3p')

mir_selected[species_short] = list(set(mir_selected[species_short]))
f.close()
######count the number of nove_mir,conserve_mir,mir-5p,mir-3p
novel_pre = 0
conserve_pre = 0
known_pre = 0
novel_5p = 0
novel_3p = 0
conserve_5p = 0
conserve_3p = 0
known_5p = 0
known_3p = 0
novel_both = 0
conserve_both = 0
known_both = 0
for i in pre_num:
	if 'novel' in i:
		if len(pre_num[i]) == 2:
			novel_both += 1
		elif len(pre_num[i]) == 1:
			globals()['novel_%s' % pre_num[i][0]] += 1
		else:
			continue
		novel_pre += 1
	else:
		if 'conserve' in i:
			if len(pre_num[i]) == 2:
				conserve_both += 1
			elif len(pre_num[i]) == 1:
				globals()['conserve_%s' % pre_num[i][0]] += 1
			else:
				continue
			conserve_pre += 1

		else:
			if len(pre_num[i]) == 2:
				known_both += 1
			elif len(pre_num[i]) == 1:
				globals()['known_%s' % pre_num[i][0]] += 1
			else:
				continue
			known_pre += 1

with open('{0}/{1}_mirna_count_stats.csv'.format(work_dir, species_short), 'wb') as f:
	f_csv = csv.writer(f)
	head = ['type', 'Pre-miRNAs', 'miRNA-5p', 'miRNA-3p', 'Both', 'Mature miRNAs']
	f_csv.writerow(head)
	f_csv.writerow(['known', known_pre, known_5p, known_3p, known_both, known_5p + known_3p + known_both * 2])
	f_csv.writerow(['novel', novel_pre, novel_5p, novel_3p, novel_both, novel_5p + novel_3p + novel_both * 2])
	f_csv.writerow(['conserve', conserve_pre, conserve_5p, conserve_3p, conserve_both,
					conserve_5p + conserve_3p + conserve_both * 2])
f.close()

