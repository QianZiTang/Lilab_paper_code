
import sys
import re
from optparse import OptionParser
import subprocess
import os
import time
import numpy as np
import math
import random

#scipy.stats.pearsonr
#scipy.stats.spearmanr

class generate:

	def __init__ (self,path='F:'):

                self.input1 = open(path+'/3cluster-UMAP-distance.txt','r')
                self.input2 = open(path+'/3cluster-group.txt','r')
                self.output = open(path+'/clusterI-rep-spots187.txt','w')

	def generate (self):

                index2class={}
                i=0
	        for line in self.input2:
	               line=line.rstrip()
	               parts=line.split()
	               if i!=0:
	                       index2class[parts[0]]=parts[1]

	               i+=1

                i=0
                alldata={}
                allx=[]
                ally=[]
                for line in self.input1:
                       line=line.rstrip()
                       parts=line.split('\t')
                       if i!=0:
                            mytype=index2class[str(i)]
                            if mytype=='0':
                                    alldata[parts[0]]=[float(parts[1]),float(parts[2])]
                                    allx.append(float(parts[1]))
                                    ally.append(float(parts[2]))


                       i+=1


                mycenter_x=np.mean(allx)
                mycenter_y=np.mean(ally)

                alldists=[]
                for myname in alldata:
                        myx=alldata[myname][0]
                        myy=alldata[myname][1]

                        thedist=math.sqrt((myx-mycenter_x)**2+(myy-mycenter_y)**2)

                        alldists.append([thedist,myname])

                alldists.sort()

                startpoint_name=alldists[0][1]

                keeplist=[]
                keeplist.append(startpoint_name)
                for m in range(0,187):
                        eachdists=[]
                        for myname1 in alldata:
                                if myname1 not in keeplist:
                                        thedist2=[]
                                        myx1=alldata[myname1][0]
                                        myy1=alldata[myname1][1]
                                        for myname2 in keeplist:
                                                myx2=alldata[myname2][0]
                                                myy2=alldata[myname2][1]

                                                thedist=math.sqrt((myx1-myx2)**2+(myy1-myy2)**2)
                                                thedist2.append(thedist)

                                        meandist2=np.mean(thedist2)
                                        eachdists.append([meandist2,myname1])

                        eachdists.sort()
                        myaddname=eachdists[0][1]

                        keeplist.append(myaddname)


                print len(keeplist)

                for c in keeplist:
                        print>>self.output,c

                self.input1.close()
                self.input2.close()
                self.output.close()


def main():

	usage = "usage: %prog [options] <pathandfiles>"
	description = "Generate jobs."

	optparser = OptionParser(version="%prog 0.1",description=description,usage=usage,add_help_option=False)
	optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")

	(options,pathandfiles) = optparser.parse_args()

	#generate().generate()
	generate().generate()


if __name__ == '__main__':

	try:
		main()
	except KeyboardInterrupt:
		sys.stderr.write("User interrupt me! ;-) See you!\n")
		sys.exit(0)
