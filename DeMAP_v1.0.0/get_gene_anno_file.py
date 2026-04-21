#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Zehong Wu wang_lab, Sun Yat-sen University"""
"""Feb, 2025"""
"""Email: wuzh93@mail2.sysu.edu.cn"""
"""Usage: This script extracts annotation information from the GTF file."""


#!/usr/bin/env python
import collections
import sys
import re

"""
argv:
    1, input file:  gtf file from gencode
"""
title = ["chrome","start","end","ensg","biotype","strand","genename","enst","trancript_type","exon_start","exon_end","cds_start","cds_end",
         "promoter","exon_len","exon_number","cds_len","UTR5_len","UTR3_len","trancript_level"]

gtffile = sys.argv[1]
ensg_enst = collections.OrderedDict()
with open(gtffile, "r") as input_file:
    for line in input_file:
        if line[0] == '#':
            continue
        a = line.split('\t')
        info = a[8].split('"')
        info = info[:-1]
        chrome = a[0]
        strand = a[6]
        start = int(a[3])-1
        end = int(a[4])
        
        if a[2] == "gene":
            ensg = info[1]
            biotype = info[3]
            genename = info[5]
            start = str(start)
            end = str(end)
            if strand == "+":
                promotor_start = int(start) - 1000
                if promotor_start <= 0:
                    promotor_start = 0
                promotor_start = str(promotor_start)
                promotor_end = start
            if strand == "-":
                promotor_start = end
                promotor_end = str(int(end) + 1000)
            ensg_enst[ensg] = [chrome,start,end,ensg,biotype,strand,genename,{},promotor_start+","+promotor_end]
            
        if a[2] == "transcript":
            ensg = info[1]
            enst = info[3]
            transcript_type = info[9]
            if "Ensembl_canonical" in line:
                ensg_enst[ensg][7][enst] = [transcript_type,[],[],"Ensembl_canonical"]
            elif "basic" in line:
                ensg_enst[ensg][7][enst] = [transcript_type,[],[],"Basic"]
            else:
                ensg_enst[ensg][7][enst] = [transcript_type,[],[],"Others"]
        
        if a[2] == "exon":
            ensg = info[1]
            enst = info[3]
            region = [start,end]
            ensg_enst[ensg][7][enst][1].append(region)
            
        if a[2] == "CDS":
            ensg = info[1]
            enst = info[3]
            region = [start,end]
            ensg_enst[ensg][7][enst][2].append(region)


priority = {'protein_coding':0,'lncRNA':1,'pseudogene':2}

def sort_key(item):
    gene_type = item[1][4] 
    return priority.get(gene_type,999)

ensg_enst = sorted(ensg_enst.items(), key=sort_key)
ensg_enst = dict(ensg_enst)

print("\t".join(title))
for gene in ensg_enst:
    strand = ensg_enst[gene][5]
    for transcript_id in ensg_enst[gene][7]:
        transcript_type = ensg_enst[gene][7][transcript_id][0]
        transcript_level = ensg_enst[gene][7][transcript_id][3]
        exonlist = ensg_enst[gene][7][transcript_id][1]
        exonlist.sort(key=lambda x:x[0])
        cdslist = ensg_enst[gene][7][transcript_id][2]
        cdslist.sort(key=lambda x:x[0])
        exon_length = str(sum([i[1]-i[0] for i in exonlist]))
        exonnum = str(len(exonlist))

        UTR3 = []
        UTR5 = []
        if not cdslist:
            cds_len = "0"
            utr3_len = "0"
            utr5_len = "0" 
        else:
            cds_start = cdslist[0][0]
            cds_end = cdslist[-1][1]
            if strand == "+":
                for i in range(0,len(exonlist)):
                    if exonlist[i][0] <= cds_end and cds_end <= exonlist[i][1]:
                        UTR3.append([cds_end,exonlist[i][1]])
                        for j in range(i+1,len(exonlist)):
                            UTR3.append(exonlist[j])
                        break
                    else:
                        continue
                    
                for i in range(0,len(exonlist)):
                    if not (exonlist[i][0] <= cds_start and cds_start <= exonlist[i][1]):
                        UTR5.append(exonlist[i])
                    else:
                        UTR5.append([exonlist[i][0],cds_start])
                        break
            else:
                for i in range(0,len(exonlist)):
                    if not (exonlist[i][0] <= cds_start and cds_start <= exonlist[i][1]):
                        UTR3.append(exonlist[i])
                    else:
                        UTR3.append([exonlist[i][0],cds_start])
                        break
                for i in range(0,len(exonlist)):
                    if exonlist[i][0] <= cds_end and cds_end <= exonlist[i][1]:
                        UTR5.append([cds_end,exonlist[i][1]])
                        for j in range(i+1,len(exonlist)):
                            UTR5.append(exonlist[j])
                        break
                    else:
                        continue
            cds_len = str(sum([i[1]-i[0] for i in cdslist]))
            utr3_len = str(sum([i[1]-i[0] for i in UTR3]))
            utr5_len = str(sum([i[1]-i[0] for i in UTR5]))        

        exonstart = ",".join([str(i[0]) for i in exonlist])
        exonend = ",".join([str(i[1]) for i in exonlist])
        cdstart = ",".join([str(i[0]) for i in cdslist])
        cdsend = ",".join([str(i[1]) for i in cdslist])
        promoter = ensg_enst[gene][8]
        if cdslist:
            print("\t".join(ensg_enst[gene][0:7])+"\t"+"\t".join([transcript_id,transcript_type,exonstart,exonend,cdstart,cdsend,promoter,exon_length,exonnum,cds_len,utr5_len,utr3_len,transcript_level]))
        else:
            print("\t".join(ensg_enst[gene][0:7])+"\t"+"\t".join([transcript_id,transcript_type,exonstart,exonend,"cmpl","cmpl",promoter,exon_length,exonnum,cds_len,utr5_len,utr3_len,transcript_level]))
            
            
