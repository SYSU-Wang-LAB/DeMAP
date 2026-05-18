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


def parse_attributes(attr_string):
    attr_dict = {}
    for attr in attr_string.split(';'):
        attr = attr.strip()
        if not attr:
            continue
        parts = re.split(r'[\s\t]+', attr, maxsplit=1)
        if len(parts) == 2:
            key = parts[0]
            val = parts[1].strip('"') 
            attr_dict[key] = val
    return attr_dict

with open(gtffile, "r") as input_file:
    for line in input_file:
        if line[0] == '#':
            continue
        a = line.strip().split('\t')
        if len(a) < 9:
            continue
            
        chrome = a[0]
        feature_type = a[2] # gene, transcript, exon, CDS
        strand = a[6]
        start = int(a[3]) - 1
        end = int(a[4])
        
        info = parse_attributes(a[8])
        
        ensg = info.get("gene_id", "")
        enst = info.get("transcript_id", "")
        
        if feature_type == "gene":
            biotype = info.get("gene_type", info.get("gene_biotype", ""))
            genename = info.get("gene_name", ensg)
            start_str = str(start)
            end_str = str(end)
            
            if strand == "+":
                promotor_start = start - 1000
                if promotor_start <= 0:
                    promotor_start = 0
                promotor_start = str(promotor_start)
                promotor_end = start_str
            else: # strand == "-"
                promotor_start = end_str
                promotor_end = str(int(end) + 1000)
                
            ensg_enst[ensg] = [chrome, start_str, end_str, ensg, biotype, strand, genename, collections.OrderedDict(), promotor_start+","+promotor_end]
            
        elif feature_type == "transcript":
            if ensg not in ensg_enst:
                continue
            transcript_type = info.get("transcript_type", info.get("transcript_biotype", ""))
            
            if "Ensembl_canonical" in line:
                ensg_enst[ensg][7][enst] = [transcript_type, [], [], "Ensembl_canonical"]
            elif "basic" in line:
                ensg_enst[ensg][7][enst] = [transcript_type, [], [], "Basic"]
            else:
                ensg_enst[ensg][7][enst] = [transcript_type, [], [], "Others"]
        
        elif feature_type == "exon":
            if ensg in ensg_enst and enst in ensg_enst[ensg][7]:
                region = [start, end]
                ensg_enst[ensg][7][enst][1].append(region)
            
        elif feature_type == "CDS":
            if ensg in ensg_enst and enst in ensg_enst[ensg][7]:
                region = [start, end]
                ensg_enst[ensg][7][enst][2].append(region)

priority = {'protein_coding':0, 'lncRNA':1, 'pseudogene':2}

def sort_key(item):
    gene_type = item[1][4] 
    return priority.get(gene_type, 999)

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
            cds_len = str(sum([i[1]-i[0] for i in cdslist]))
            utr3_len = str(sum([i[1]-i[0] for i in UTR3]))
            utr5_len = str(sum([i[1]-i[0] for i in UTR5]))        

        exonstart = ",".join([str(i[0]) for i in exonlist])
        exonend = ",".join([str(i[1]) for i in exonlist])
        cdstart = ",".join([str(i[0]) for i in cdslist])
        cdsend = ",".join([str(i[1]) for i in cdslist])
        promoter = ensg_enst[gene][8]
        
        out_base = "\t".join(ensg_enst[gene][0:7])
        if cdslist:
            out_tail = [transcript_id,transcript_type,exonstart,exonend,cdstart,cdsend,promoter,exon_length,exonnum,cds_len,utr5_len,utr3_len,transcript_level]
        else:
            out_tail = [transcript_id,transcript_type,exonstart,exonend,"cmpl","cmpl",promoter,exon_length,exonnum,cds_len,utr5_len,utr3_len,transcript_level]
            
        print(out_base + "\t" + "\t".join(out_tail))