# -*- coding: utf-8 -*-


import numpy  as np
import sys
import  math
import scipy.stats
from collections import OrderedDict

args = sys.argv

#k=7

inFile = open(args[1])
miR_file= open(args[2])
k =int(args[3])

#miR_file= open("../top50_ginger.txt")
out = open (str(k) +"mers.txt","w")
hits_file = open (str(k) +"mers_hits.txt","w")


def fastaToDict(inFile):
    dict = {}
    for line in inFile:
        if line[0] == '>':
            splitted = line.split(' ')
            seqID = splitted[0].strip()[1:]  # unique transcript ID
            dict[seqID] = ""
        else:
            dict[seqID] += line.strip().upper()
    return dict


def calculate_prob_MM(seq , transition_matrix):
    transition_dict ={  "AA":transition_matrix[0],"AC":transition_matrix[1],"AG":transition_matrix[2],"AT":transition_matrix[3],
                      "CA":transition_matrix[4],"CC":transition_matrix[5],"CG":transition_matrix[6],"CT":transition_matrix[7],
                      "GA":transition_matrix[8],"GC":transition_matrix[9],"GG":transition_matrix[10],"GT":transition_matrix[11],
                      "TA":transition_matrix[12],"TC":transition_matrix[13],"TG":transition_matrix[14],"TT":transition_matrix[15]}
    p = 0.25
    for i in range(1,len(seq)):
        p *=  transition_dict[seq[i-1] + seq[i]]
    return p

def firstOrderMMConstruction(sequence):
    seq_len=len(sequence)
    transition_counts ={  "AA":0,"AC":0,"AG":0,"AT":0,
                          "CA":0,"CC":0,"CG":0,"CT":0,
                          "GA":0,"GC":0,"GG":0,"GT":0,
                          "TA":0,"TC":0,"TG":0,"TT":0}
    for i in range(1,seq_len):
        if sequence[i] in [ 'A', 'C' ,'G', 'T'] and  sequence[i-1] in [ 'A' ,'C' ,'G' ,'T']:
            transition_counts[sequence[i-1]+sequence[i]] += 1.0
    matrix= OrderedDict(sorted(transition_counts.items(), key=lambda t: t[0]))   # to restore dictionary original order
    values =  np.fromiter(matrix.values(), dtype=int)
    values = values + 1.0    # adding pseudo_count
    new_values = values
    new_values[0:4] = values[0:4]/sum(values[0:4])
    new_values[4:8] = values[4:8]/sum(values[4:8])
    new_values[8:12] = values[8:12]/sum(values[8:12])
    new_values[12:16] = values[12:16]/sum(values[12:16])
    print(matrix)
    return new_values

def firstOrderMMGeneration(transitionMatrix,seq_len):
    alphabet = ['A','C','G','T']
    s= ''.join(np.random.choice(alphabet,1,[.25,.25,.25,.25]))  # generate random first letter
    for i in range(1,seq_len):
        if s[i-1] == 'A':
            s += ''.join(np.random.choice(alphabet,1,p = transitionMatrix[0:4]))
        elif s[i-1] == 'C':
            s += ''.join(np.random.choice(alphabet,1,p = transitionMatrix[4:8]))
        elif s[i-1] == 'G':
            s += ''.join(np.random.choice(alphabet,1,p = transitionMatrix[8:12]))
        else:
            s += ''.join(np.random.choice(alphabet,1,p = transitionMatrix[12:16]))
    return s


def encoding(seq):
    k=len(seq)
    str=''
    c=''
    for i in range(0,k):
        if(seq[i] == 'A' or seq[i] == 'a'):
            str='00'
        elif (seq[i] == 'C' or seq[i] == 'c'):
            str='01'
        elif (seq[i] == 'G' or seq[i] == 'g'):
            str='10'
        elif (seq[i] == 'T' or seq[i] == 't'):
            str='11'
        elif (seq[i] == 'N' or seq[i] == 'n'):
            return -1
        c=c+str
    return int(c,2)

def reverseComplement (s):
    complement = {'A':'T','C':'G','G':'C','T':'A'}
    return ''.join([ complement[s[i]] for i in range(len(s)-1,-1,-1)])

def indexing(subject,k): # building hash table for subject database (sequences to be searched)
    num_seqs=len(subject)
    hash_tbl= [[] for i in range(4**k)]
    keys = list(subject.keys())
    for i in range(0,num_seqs):
        seq=subject[keys[i]]
        for j in range(0,len(seq) -k + 1): # enumerating overlapping k-mers
            k_mer= seq[j:j+k]
            hash_tbl[encoding(k_mer)].append([i+1,j+1]) # 1_based indexing
    return hash_tbl,list(subject.keys())

def main(miR_file,k):
    utr_regions = fastaToDict(inFile)
    len_dict = {}
    firstOrderModels = {}
    utr_regions_keys = list(utr_regions.keys())
    for p in range(0,len(utr_regions)):
        key = utr_regions_keys[p]
        firstOrderModels[key] = firstOrderMMConstruction(utr_regions[key])
        len_dict [key] = float(len(utr_regions[key]))
    [hash_tbl, transcriptIDs]  = indexing(utr_regions,k)
    dict = {}
    miR_seed = {}
    for line in miR_file:
        splitted = line.split(" ")
        if line[0] == '>':
            if len(splitted[0]) > 1:
                id = splitted[0][1:].strip()
            else:
                id = splitted[1].strip()

        else :
            if k == 8:
                seq = splitted[0][0:8]
            elif k ==7:
                seq = splitted[0][1:8]               # seed region only
            elif k == 6:
                seq = splitted[0][1:7]
            else:
                seq = splitted[0][3:8]
            seq = seq.replace("U","T")
            seq = seq.replace("u","t")
            miR_seed [id] = seq.strip()
            dict[id] = {}
            for n in range(0,len(transcriptIDs)):
                dict[id] [transcriptIDs[n]] = 0.0
            for j in range(0,len(seq) - k + 1): # enumerating overlapping kmers
                k_mer= seq[j:j+k]
                kmer_rev = reverseComplement(k_mer)
                hits = hash_tbl[encoding(kmer_rev)]
                for i in range(0,len(hits)):
                    hits_file.write(k_mer +"\t"+ id +"\t"+ transcriptIDs[hits[i][0]-1] +"\t" + str(hits[i][1]) +"\t" + str(len_dict[transcriptIDs[hits[i][0]-1]]) +  "\n" )
                    dict[id] [transcriptIDs[hits[i][0]-1]] += 1
    return dict,len_dict,firstOrderModels ,miR_seed
#########################################################

hits_file.write("miRNA_seed" + "\t" +"miRNA_ID" + '\t' +"targetSeqID" + "\t" + "hit_start_position" + "\t" +
                "targetSeqLength"+ "\n")

[dict, len_dict,models,miR_seed] =  main(miR_file,k)

out.write('%-14s \t %20s \t %20s \t %25s \t %30s \t %5s\n'  %("miRNA_seed_" + str(k)+ "nt","miRNA_ID", "targetSeqID",
                                                                    "#_of_occurrences_of_"+str(k)+"mers","#_of_all_possible_" + str(k)+"mers_in_3UTR", "pValue"))

pValue = {}

dKeys = list(dict.keys())
for i in range(0,len(dict)):
    mir =  dKeys[i]
    keys = list(dict[mir].keys())  # transcripts IDs
    pValue[mir] ={}
    num_ocur = 0
    num_possible_kmers = 0
    for j in range(0,len(dict[mir])):
        pValue[mir][keys[j]] = 0.0
    for j in range(0,len(dict[mir])):
        num_ocur = dict[mir][keys[j]]
        num_possible_kmers = float(len_dict[keys[j]]-k+1)
        p = calculate_prob_MM(miR_seed[mir], models[keys[j]])
        transc_p = scipy.stats.binom.sf(num_ocur-1,num_possible_kmers,p )
        if transc_p >0:
            pValue[mir][keys[j]] += -2.0 * math.log(transc_p)   # gene-based p-values  !!!!!!!!!!!!
        else:
            transc_p = 1
            pValue[mir][keys[j]] += -2.0 * math.log(transc_p)
        out.write('%-14s \t %20s \t %20s \t %25s \t %30s \t %4f \n' %(miR_seed[mir],mir,keys[j], str(num_ocur), str(num_possible_kmers)
                                                   , transc_p ))


