#!/usr/bin/env python
# Programmer : zhuxp
# Date: 
# Last-modified: 12-25-2015, 22:16:09 EST
from __future__ import print_function
VERSION="0.1"
import os,sys,argparse
from bam2x.Annotation import BED12
from bam2x import TableIO,Tools
from bam2x import IO
import time
import itertools
import multiprocessing as mp
from collections import defaultdict
def init():
    return [];
def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency : bam2x')
    p.add_argument('-v','--version',action='version',version='%(prog)s '+VERSION)
    p.add_argument('-i','--input',dest="input",default="stdin",type=str,help="input file DEFAULT: STDIN")
    p.add_argument('-b','--db',dest="db",type=str,help="input file db bed file (RFP gene coordinates file)")
    p.add_argument('-o','--output',dest="output",type=str,default="stdout",help="output file DEFAULT: STDOUT")
    #p.add_argument('-n','--num_cpus',dest="num_cpus",type=int,default=4,help="number of cpus DEFAULT: %(default)i")
    if len(sys.argv)==1:
        print(p.print_help(),file=sys.stderr)
        exit(0)
    return p.parse_args()
def Main():
    '''
    IO TEMPLATE
    '''
    global args,out
    args=ParseArg()
    fin=IO.fopen(args.input,"r")
    out=IO.fopen(args.output,"w")
    G2RFP=defaultdict(init);
    '''
    END OF IO TEMPLATE 
    '''
    m=[0.0 for i in xrange(200)] # 100 to 200
    for b in TableIO.parse(IO.fopen(args.db,"r"),"bed6"):
        G2RFP[b.chr].append(b)
    total_reads=0;
    for i0,i in enumerate(TableIO.parse(fin,"bed6")):
        spectral=[0 for j in xrange(200)] # 100 to 200
        for j in G2RFP[i.chr]:
            dis=i.start-j.start
            if(dis >=-100 and dis<100): # 50 to 100
                spectral[dis+100]+=j.score # 50 to 100
        spectral=norm(spectral)
        total_reads+=i.score
        m=[a*i.score+b for a,b in itertools.izip(spectral,m)]
        if i0%100==0:
            print("{} processed\r".format(i0),file=sys.stderr)
    print("pos\tvalue",file=out);
    for i,x in enumerate(m):
        print("{}\t{}".format(i,float(x)/total_reads),file=out)



def norm(spectrum):
    s=sum(spectrum);
    if s==0:
        s=1
    return [float(i)/s for i in spectrum]
    
if __name__=="__main__":
    Main()






