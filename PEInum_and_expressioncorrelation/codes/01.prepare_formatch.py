
import sys
import re
from optparse import OptionParser
import subprocess
import os
import time
import math
import numpy as np

#median of mean

class generate:

	def __init__ (self):

		self.path = '/Lustre02/data/hic/forJinplot/crossspecies/'

	        self.input01 = open('/Lustre01/tangqianzi/forJinplot/forJin_addHiC/crossspecies/exprs_sele/orth_genes.nochicken.nomouse.txt','r')
		self.input02 = open('/Lustre02/data/hic/forJinOrth/quant_summary/sus_scrofa.Tau.7tissues.xls','r')

		self.input1 = open('/Lustre02/data/hic/forJinplot/crossspecies/merge_Enum.nochicken.nomouse.filtered.xls','r')
		self.input2 = open('/Lustre01/tangqianzi/forJinplot/forJin_addHiC/crossspecies/sample_infor.rest.txt','r')
		self.input3 = open('/Lustre01/tangqianzi/forJinplot/forJin_addHiC/crossspecies/exprs_sele/merge_codingtable_forcluster.final.tpm.nochicken.nomouse.normed.filtered.xls','r')
		self.output = open(self.path+'/formatch.nochicken.nomouse.correct.filtered.xls','w')

	def generate (self):

	        pig2human={}
	        for line in self.input01:
	               line=line.rstrip()
	               parts=line.split('\t')
	               mypig=parts[2].split('_')[0]
	               myhuman=parts[1].split('_')[0]

	               pig2human[mypig]=myhuman


                mygenes=[]
                count0=0
                mygenes2=[]
                mygenes3=[]
	        for line in self.input02:
	               line=line.rstrip()
	               parts=line.split('\t')
	               if str(parts[1])!='NA' and float(parts[1])<0.1:
	                       if parts[0] in pig2human:
	                               mygenes.append(pig2human[parts[0]])
	                       count0+=1

	               elif str(parts[1])!='NA' and float(parts[1])<0.75 and float(parts[1])>0.1:
	                       if parts[0] in pig2human:
	                               mygenes2.append(pig2human[parts[0]])

	               elif str(parts[1])!='NA':
	                       if parts[0] in pig2human:
	                               mygenes3.append(pig2human[parts[0]])


	        print count0





	        species2sample={}
	        for line in self.input2:
	               line=line.rstrip()
	               parts=line.split()
	               if parts[0]=='chicken' or parts[0]=='mouse':
	                       continue

	               if parts[0] not in species2sample:
	                       species2sample[parts[0]]=[]

	               species2sample[parts[0]].append(parts[1]+'_'+parts[3])


	        keepdata=[]
	        i=0
	        for line in self.input3:
	               line=line.rstrip()
	               parts=line.split()
	               keepdata.append(parts)
	               if i==0:
	                       species2sample2={}
	                       for  j in range(1,len(parts)):
	                               myspe=parts[j].split('_')[-1]
	                               if myspe=='mouse':
	                                       continue

	                               if myspe not in species2sample2:
	                                       species2sample2[myspe]=[]
	                               species2sample2[myspe].append(parts[j])

	               i+=1

                print species2sample2

                mygenes_E=[]
                mygenes2_E=[]
                mygenes3_E=[]

                gene2index={}
	        i=0
	        mymax=0
	        for line in self.input1:
	               line=line.rstrip()
	               parts=line.split('\t')
	               if i==0:
	                       species2index={}

	                       for myspe in species2sample:
	                               species2index[myspe]=[]
	                               for mysample in species2sample[myspe]:
	                                       for j in range(1,len(parts)):
	                                               if mysample==parts[j]:
	                                                       species2index[myspe].append(j)

	                       print species2index



	               else:
	                       allvals=[]
                               for myspe in species2index:
                                    each=[]
                                    for j in species2index[myspe]:
                                            each.append(float(parts[j]))

                                    myval=np.mean(each)
                                    myval=int(round(myval,0))

                                    allvals.append(myval)

                               mymedian=np.median(allvals)
                               ## mymedian=np.mean(allvals)

                               if mymedian>=5:
                                    sig=1
                               elif mymedian<=1:
                                    sig=0
                               else:
                                    sig='NA'

                               gene2index[parts[0]]=str(sig)

                               mymax=max(mymax,mymedian)

                               if parts[0] in mygenes:
                                    mygenes_E.append(mymedian)

                               if parts[0] in mygenes2:
                                    mygenes2_E.append(mymedian)

                               if parts[0] in mygenes3:
                                    mygenes3_E.append(mymedian)



	               i+=1


                print "mymax",mymax

                print np.mean(mygenes_E),"housekeepings"
                print np.mean(mygenes2_E),"others"
                print np.mean(mygenes3_E),"specific"


	        i=0
	        count1=0
	        count2=0
	        print>>self.output,'treat'+'\t'+'mymeanexprs'
	        ## for line in self.input3:
	        for parts in keepdata:
	               ##line=line.rstrip()
	               ##parts=line.split('\t')
	               if i==0:
	                       species2index2={}
	                       for myspe in species2sample2:
	                               species2index2[myspe]=[]
	                               for mysample in species2sample2[myspe]:
	                                       for j in range(1,len(parts)):
	                                               if mysample==parts[j]:
	                                                       species2index2[myspe].append(j)

	                       print species2index2



	               elif i!=0:
	                       alleach=[]
	                       for myspe in species2index2:
	                               each=[]
	                               for j in species2index2[myspe]:
	                                       each.append(float(parts[j]))

	                               mymean=np.mean(each)
	                               alleach.append(mymean)

	                       myallmean=np.mean(alleach)

	                       if parts[0] in gene2index:
	                               sig=gene2index[parts[0]]
	                               ## if sig!='NA' and float(myallmean)>1:
	                               if sig!='NA':
	                                       print>>self.output,parts[0]+'\t'+str(sig)+'\t'+str(myallmean)

	               i+=1

                print count1
                print count2

		self.input1.close()
		self.input2.close()
                self.output.close()

		self.input01.close()
		self.input02.close()


def main():

	usage = "usage: %prog [options] <pathandfiles>"
	description = "Generate jobs."

	optparser = OptionParser(version="%prog 0.1",description=description,usage=usage,add_help_option=False)
	optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")

	(options,pathandfiles) = optparser.parse_args()

	generate().generate()


if __name__ == '__main__':

	try:
		main()
	except KeyboardInterrupt:
		sys.stderr.write("User interrupt me! ;-) See you!\n")
		sys.exit(0)
