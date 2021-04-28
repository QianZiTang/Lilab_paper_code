
import sys
import re
from optparse import OptionParser
import subprocess
import os
import time
import math
import numpy as np
import scipy.stats

class generate:

	def __init__ (self,path='/Lustre02/data/hic/forJinplot/crossspecies/',tissue='SAT'):

		self.path = '/Lustre02/data/hic/forJinplot/crossspecies/'
		## self.input1 = open('/Lustre02/data/hic/forJinOrth/mRNAOrth/nochicken/all_families_gene_one2one.10.txt','r')
		## self.input2 = open('/Lustre02/data/hic/forJinOrth/quant_summary/sus_scrofa.Tau.7tissues.xls','r')
		self.input3 = open('/Lustre01/tangqianzi/forJinplot/forJin_addHiC/crossspecies/exprs_sele/merge_codingtable_forcluster.final.tpm.nochicken.nomouse.normed.filtered.xls','r')
		self.input4 = open(path+'/formatch_results.nochicken.nomouse.correct.filtered.xls','r')
		self.input5 = open('/Lustre02/data/hic/forJinOrth/exprs_trees/divergence_time.txt','r')

		self.output11 = open(path+'/housekeeping_forplot11.'+tissue+'.all.txt','w')
		self.output12 = open(path+'/housekeeping_forplot12.'+tissue+'.all.txt','w')
		self.output21 = open(path+'/housekeeping_forplot21.'+tissue+'.all.txt','w')
		self.output22 = open(path+'/housekeeping_forplot22.'+tissue+'.all.txt','w')
		self.tissue = tissue

	def generate (self):

	        genes1=[]
	        genes2=[]
	        i=0
	        for line in self.input4:
	               line=line.rstrip()
	               parts=line.split()
	               if i!=0:
	                       if int(parts[1])==1:
	                               ## print line
	                               genes1.append(parts[0])
	                       else:
	                               genes2.append(parts[0])


	               i+=1


                i=0
                for line in self.input3:
                       line=line.rstrip()
                       parts=line.split()
                       if i==0:
                            ##allspecies=[]
                            index2species={}
                            allresults1={}
                            allresults2={}

                            spe2samples={}

                            for j in range(1,len(parts)):
                                    ## myt=parts[j].split('-')[0]
                                    myspe=parts[j].split('_')[-1]
                                    if myspe not in spe2samples:
                                            spe2samples[myspe]=[]
                                    spe2samples[myspe].append(parts[j])

                                    index2species[str(j)]=parts[j]

                                    allresults1[parts[j]]=[]
                                    allresults2[parts[j]]=[]

                       else:
                            if parts[0] in genes1:
                                    for j in range(1,len(parts)):
                                            if str(j) in index2species:
                                                    myname=index2species[str(j)]
                                                    allresults1[myname].append(float(parts[j]))

                            elif parts[0] in genes2:
                                    for j in range(1,len(parts)):
                                            if str(j) in index2species:
                                                    myname=index2species[str(j)]
                                                    allresults2[myname].append(float(parts[j]))

                       i+=1

                ## print

                allpairs={}
                i=0
                for line in self.input5:
                       line=line.rstrip()
                       parts=line.split('\t')
                       if i==0:
                            index2name={}
                            allnames=[]
                            for j in range(0,len(parts)):
                                    index2name[str(j+1)]=parts[j]
                                    if parts[j]!='chicken' and parts[j]!='guinea_pig' and parts[j]!='mouse':
                                            allnames.append(parts[j])

                       else:
                            for j in range(1,len(parts)):
                                    myname1=parts[0]
                                    myname2=index2name[str(j)]
                                    if myname1!='chicken' and myname1!='guinea_pig' and myname2!='chicken' and myname2!='guinea_pig' and myname1!='mouse' and myname2!='mouse':
                                            allpairs[myname1+'-'+myname2]=float(parts[j])

                       i+=1


                final1={}
                final2={}


                finalall=[]
                ##finalall2=[]

                for i in range(0,len(allnames)):
                       for j in range(0,len(allnames)):
                            if j>i:
                                   myname1=allnames[i]
                                   myname2=allnames[j]

                                   allcors1=[]
                                   allcors2=[]

                                   mydist=str(allpairs[myname1+'-'+myname2])

                                   for sample1 in spe2samples[myname1]:
                                        for sample2 in spe2samples[myname2]:

                                                mycor1=float(scipy.stats.spearmanr(allresults1[sample1],allresults1[sample2])[0])
                                                mycor2=float(scipy.stats.spearmanr(allresults2[sample1],allresults2[sample2])[0])

                                                allcors1.append(mycor1)
                                                allcors2.append(mycor2)

                                                finalall.append([float(mydist),mycor1,mycor2])
                                                ##finalall2.append([float(mydist),mycor2])

                                   ## print>>self.output11,mydist+'\t'+str(np.mean(allcors1))+'\t'+myname1+'-'+myname2
                                   ## print>>self.output12,mydist+'\t'+str(np.mean(allcors2))+'\t'+myname1+'-'+myname2

                                   if mydist not in final1:
                                        final1[mydist]=[]

                                   final1[mydist].append(float(mycor1))

                                   if mydist not in final2:
                                        final2[mydist]=[]

                                   final2[mydist].append(float(mycor2))



                finalall.sort()

                for parts in finalall:
                        print>>self.output11,str(parts[0])+'\t'+str(parts[1])
                        print>>self.output12,str(parts[0])+'\t'+str(parts[2])



                theresults1=[]
                theresults2=[]
                for mydist in final1:
                        finalval1=np.mean(final1[mydist])
                        finalval2=np.mean(final2[mydist])

                        theresults1.append([float(mydist),str(finalval1)])
                        theresults2.append([float(mydist),str(finalval2)])

                theresults1.sort()
                theresults2.sort()

                for parts in theresults1:
                        print>>self.output21,str(parts[0])+'\t'+parts[1]


                for parts in theresults2:
                        print>>self.output22,str(parts[0])+'\t'+parts[1]


		self.input3.close()
		self.input4.close()
		self.input5.close()

		self.output11.close()
                self.output12.close()
                self.output21.close()
                self.output22.close()


def main():

	usage = "usage: %prog [options] <pathandfiles>"
	description = "Generate jobs."

	optparser = OptionParser(version="%prog 0.1",description=description,usage=usage,add_help_option=False)
	optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")

	(options,pathandfiles) = optparser.parse_args()

	generate(tissue=pathandfiles[0]).generate()


if __name__ == '__main__':

	try:
		main()
	except KeyboardInterrupt:
		sys.stderr.write("User interrupt me! ;-) See you!\n")
		sys.exit(0)
