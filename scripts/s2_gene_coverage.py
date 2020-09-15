import numpy as np
from collections import Counter
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import pickle
#from os import *


NP29="all_NP29.bed.NP29.bed"
ref="317.111_no_dup.gff.NP29.bed"

def process_bed(infile, RPKM=True):
    '''this function takes the bed file, and counts the occurrence per gene of a hit then return a dictionary per gene'''
    bed = [line.strip().split('\t')[8] for line in open(infile).readlines()]
    bed_dict = dict(Counter(bed))
    bed_dict_correct=sum(bed_dict.values())
    for key in bed_dict:
        if RPKM == True:
            temp = (bed_dict[key]/float(bed_dict_correct))*1000
        else:
            temp=(bed_dict[key])
        bed_dict[key] = temp
    return bed_dict

def compare(host, host_name, full_genome, pd_hosts):
    pd_hosts.append(host_name)
    for rec in full_genome:
        if rec in host.keys():
                h1 = host[rec]
                print "YES"
        else:
            h1=0
        pd_hosts[rec][host_name]=h1
    return pd_hosts

def iterate_through_samples(sample_list, pd_hosts):
    '''give list of samples, make pd_hosts file combining all'''
    for sample in sample_list:
        temp=process_bed(sample, RPKM=False)
        pd_hosts=compare(temp, sample, full_genome, pd_hosts)
    return pd_hosts


sample_list=[line.strip() for line in os.popen('ls | grep NP29 | grep Barley').readlines()]
sample_list.append('NP29')
pd_hosts = pd.DataFrame(columns=reference.keys(), index=sample_list)
np = process_bed(NP29, RPKM=False)
pd_hosts = compare(np, "NP29", full_genome, pd_hosts

'''Process samples together as 'all'''
pd_hosts=iterate_through_samples(sample_list, pd_hosts)





plt.subplot(2,2,1)
x=pd_hosts.ix["NP29"]
y=pd_hosts.ix["peas"]
plt.title("Peas vs NP29")
plt.xlim(0,max(x))
plt.ylim(0,max(y))
plt.scatter(x, y, s=25, c="green", marker='o', alpha=0.5)

plt.subplot(2,2,2)
x=pd_hosts.ix["NP29"]
y=pd_hosts.ix["barley"]
plt.title("Barley vs NP29")
plt.xlim(0,max(x))
plt.ylim(0,max(y))
plt.scatter(x, y, s=25, c="green", marker='o', alpha=0.5)

plt.subplot(2,2,3)
x=pd_hosts.ix["NP29"]
y=pd_hosts.ix["col"]
plt.title("A. thaliana vs NP29")
plt.xlim(0,max(x))
plt.ylim(0,max(y))
plt.scatter(x, y, s=25, c="green", marker='o', alpha=0.5)

plt.subplot(2,2,4)
x=pd_hosts.ix["NP29"]
y=pd_hosts.ix["canola"]
plt.title("Canola vs NP29")
plt.xlim(0,max(x))
plt.ylim(0,max(y))
plt.scatter(x, y, s=25, c="green", marker='o', alpha=0.5)
plt.tight_layout()
plt.savefig("/ebio/abt6_projects9/Pseudomonas_diversity/Tnseq/processed_reads/hiseq0081/all_hosts.pdf")



x=[numpy.log10(rec) for rec in ((pd_hosts.ix["canola"]+0.000000001)/(pd_hosts.ix["NP29"]+0.000000001))]
y=[numpy.log10(rec) for rec in ((pd_hosts.ix["col"]+0.000000001)/(pd_hosts.ix["NP29"]+0.000000001))]
plt.xlim(0,max(x))
plt.ylim(0,max(y))
plt.scatter(x, y, s=25, c=pd_hosts.ix["NP29"], marker='o', alpha=0.5)
plt.tight_layout()
plt.savefig("/ebio/abt6_projects9/Pseudomonas_diversity/Tnseq/processed_reads/hiseq0081/try.pdf")


cols = pd_hosts.columns
bt = pd_hosts.apply(lambda x: x > 0)
non_zero=bt.apply(lambda x: list(cols[x.values]), axis=1)

together=[]
track=[]
for line in non_zero.keys():
   for rec in non_zero[line]:
       if rec not in together:
           together.append(rec)
   print len(together)
   track.append([line, len(together)])


for rec in pd_hosts.index:
    temp=[samp for samp in pd_hosts.ix[rec] if samp!=0]
    for thing in temp:
        if thing in present:
            pass
        else:
            present.append(thing)
    print len(present)
