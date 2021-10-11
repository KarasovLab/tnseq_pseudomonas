import numpy as np
from collections import Counter
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import pickle
import os

colum="all_Col-0_8d.bed.NP29.bed"
peas="all_Peas_8d.bed.NP29.bed"
NP29="all_NP29.bed.NP29.bed"
canola="all_Canola_8d.bed.NP29.bed"
barley="all_Barley_8d.bed.NP29.bed"
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
    for rec in full_genome:
        if rec in host.keys():
                h1 = host[rec]
                print "YES"
        else:
            h1=0
        pd_hosts[rec][host_name]=h1
    return pd_hosts


def combine_compare(host_list, full_genome):
    pd_hosts = pd.DataFrame(columns=reference.keys(), index=['col', 'canola', 'barley', 'peas', 'NP29'])

def fisher_test(host1, host2, pd_hosts):
    tot_1=sum(pd_hosts.ix[host1])
    tot_2=sum(pd_hosts.ix[host2])
    gene_dict={}
    for gene in pd_hosts.columns:
        g1=pd_hosts[gene][host1]
        g2=pd_hosts[gene][host2]
        oddratio, pvalue=scipy.stats.fisher_exact([[tot_1, g1], [tot_2, g2]])
        gene_dict[gene]=[oddratio, pvalue]
    return gene_dict

def pca_bed_files():
    '''this function is meant to perform'''
    cmd='ls | grep NP29.bed | grep -v all' 
    bed_list=os.popen(cmd).read().split()
    pd_within=pd.DataFrame(columns=reference.keys(), index=[rec.split("_tnseq")[0] for rec in bed_list])
    for samp in bed_list:
        samp_name=samp.split("_tnseq")[0]
        host=process_bed(samp, RPKM=False)
        pd_within=compare(host, samp_name, full_genome, pd_within)
    return pd_within


reference = process_bed(ref)
full_genome = reference
pd_all=pca_bed_files()

#write pd_all to file
pickle.dump(pd_all, open('/ebio/abt6_projects9/Pseudomonas_diversity/Tnseq/processed_reads/hiseq0081/pd_all.cpk', 'w'))







'''Process samples together as 'all'''
np = process_bed(NP29, RPKM=False)
pea = process_bed(peas, RPKM=False)
can = process_bed(canola, RPKM=False)
bar = process_bed(barley, RPKM=False)
col = process_bed(colum, RPKM=False)
pd_hosts = pd.DataFrame(columns=reference.keys(), index=['col', 'canola', 'barley', 'peas', 'NP29'])
pd_hosts = compare(np, "NP29", full_genome, pd_hosts)
pd_hosts = compare(pea, 'peas', full_genome, pd_hosts)
pd_hosts = compare(bar, 'barley', full_genome, pd_hosts)
pd_hosts = compare(can, 'canola', full_genome, pd_hosts)
pd_hosts = compare(col, 'col', full_genome, pd_hosts)


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
