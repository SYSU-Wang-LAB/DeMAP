#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Zehong Wu wang_lab, Sun Yat-sen University"""
"""Feb, 2025"""
"""Email: wuzh93@mail2.sysu.edu.cn"""
"""Usage: This script is used to annotate the identified methylation sites."""


import sys
import re
import math
""""""
#对bed格式位点注释
#extend 3bp for gene tss because 5' isoform
#usage：python site_annotation.py gene_anno_file peakfile > result
""""""

def load_anno_file(file):
    """
    chrome  start   end     ensg    biotype strand  genename        enst    trancript_type  exon_start      exon_end        cds_start       cds_end promoter        exon_len   exon_number     cds_len UTR5_len        UTR3_len        trancript_level
    chr1    65418   71585   ENSG00000186092.7       protein_coding  +       OR4F5   ENST00000641515.2       protein_coding  65418,65519,69036       65433,65573,71585 65564,69036      65573,70005     64418,65418     2618    3       978     60      1580    Ensembl_canonical
    chr1    450739  451678  ENSG00000284733.2       protein_coding  -       OR4F29  ENST00000426406.4       protein_coding  450739  451678  450742  451678  451678,452678      939     1       936     0       3       Ensembl_canonical
    """
    result_dic = {}
    infile = open(file,"r")
    for line in infile:
        info = line.rstrip().split("\t")
        chrome = info[0]
        try:
            start = int(info[1])
            end = int(info[2])
        except:
            continue
        ensg = info[3]
        ensg = re.sub(r"\.\d+","",ensg)
        biotype = info[4]
        strand = info[5]
        genename = info[6]
        enst = info[7]
        transcript_type = info[8]
        transcript_level = info[19]
        exonstart = list(map(int,info[9].split(",")))
        exonend = list(map(int,info[10].split(",")))
        if info[11] == "cmpl":
            cdsstart = "cmpl"
            cdsend = "cmpl"
        else:
            cdsstart = list(map(int,info[11].split(",")))
            cdsend = list(map(int,info[12].split(",")))
        promter_start = int(info[13].split(",")[0])
        promter_end = int(info[13].split(",")[1])

        if strand == "+":
            start = start - 10
        if strand == "-":
            end = end + 10

        if chrome+strand not in result_dic:
            result_dic[chrome+strand] = [[ensg,genename,biotype,start,end,strand,exonstart,exonend,cdsstart,cdsend,promter_start,promter_end,enst,transcript_type,transcript_level]]
        else:
            result_dic[chrome+strand].append([ensg,genename,biotype,start,end,strand,exonstart,exonend,cdsstart,cdsend,promter_start,promter_end,enst,transcript_type,transcript_level])
        
    infile.close()
    return(result_dic)

def calculate_distance(pos,exons,strand):
    if strand == "+":
        if pos < exons[0][0]: #upstream
            distance = pos - exons[0][0]
            return distance
        for exon_index in range(len(exons)):
            if exons[exon_index][0] <= pos and pos < exons[exon_index][1]:
                distance = sum([exon[1]- exon[0] for exon in exons[0:exon_index]]) + (pos - exons[exon_index][0])
                return distance
        #intron or downstream
        distance = pos - exons[0][0]
        return distance
    else:
        if pos > (exons[-1][1] - 1): #upstream
            distance = (exons[-1][1] - 1) - pos
            return distance
        for exon_index in range(len(exons)):
            if exons[exon_index][0] <= pos and pos < exons[exon_index][1]:
                distance = sum([exon[1]- exon[0] for exon in exons[(exon_index + 1):]]) + ((exons[exon_index][1] - 1) - pos)
                return distance
        distance = (exons[-1][1] - 1) - pos
        return distance         


def got_annotation(anno_dic,peak_loci):
    chrome1 = peak_loci[0]
    start1 = peak_loci[1]
    end1 = peak_loci[2]
    strand1 = peak_loci[3]
    location = math.ceil((start1+end1)/float(2))
    state_record = []  ### 考虑到一个位点可能overlap到不同的基因（少数基因在基因组上有重叠）or overlap 到同一基因多个转录本
    if chrome1+strand1 in anno_dic:
        for gene in anno_dic[chrome1+strand1]:
            state = "unfound"
            genetype = gene[2]
            start2 = gene[3]
            end2 = gene[4]
            transcript_id = gene[12]
            transcript_type = gene[13]
            transcript_level = gene[14]

            information = "\t".join(gene[0:3])
            if end2 >= start1 and end1 >= start2:  ### 在基因内
                bonafound = "NO"
                Exon_start = list(gene[6])
                Exon_end = list(gene[7])
                CDS_start = gene[8]
                CDS_end = gene[9]
                exonlist = list(zip(Exon_start,Exon_end))
                distance_to_tss = calculate_distance(location,exonlist,strand1)
                if strand1 == "+":
                    #distance_to_tss = abs(location - (Exon_start[0] + 1))
                    Exon_start[0] = Exon_start[0] - 10
                else:
                    #distance_to_tss = abs(location - Exon_end[-1])
                    Exon_end[-1] = Exon_end[-1] + 10

                if CDS_start == "cmpl":  ########## 在基因范围内，无CDS时不是Exon就是Intron
                    for i in range(len(Exon_start)):
                        if int(Exon_start[i]) <= end1 and start1 <= int(Exon_end[i]):
                            state = "Exon"
                            bonafound = "YES"
                    if state == "unfound":
                        #print("aaa")
                        state = "Intron"
                    
                else:   #### 蛋白编码基因
                    for i in range(len(CDS_start)):
                        if int(CDS_start[i]) <= end1 and start1 <= int(CDS_end[i]):
                            state = "CDS"
                            bonafound = "YES"
                    if state == "unfound": 
                        for i in range(len(Exon_start)):
                            if int(Exon_start[i]) <= end1 and start1 <= int(Exon_end[i]):
                                state = "Exon"
                                bonafound = "YES"
                        if state == "Exon":
                            if strand1 == "+":
                                if location > int(CDS_end[-1]):
                                    state = "3UTR"
                                if location < int(CDS_start[0]):
                                    state = "5UTR"
                            else:
                                if location > int(CDS_end[-1]):
                                    state = "5UTR"
                                if location < int(CDS_start[0]):
                                    state = "3UTR"
                    
                        if state == "unfound":   #不在exon上
                            if strand1 == "+":
                                if location <= int(Exon_end[-1]) and location >= int(Exon_start[0]):
                                    state = "Intron"
                                    bonafound = "YES"
                            else:
                                if location <= int(Exon_end[-1]) and location >= int(Exon_start[0]):
                                    state = "Intron"
                                    bonafound = "YES"        

                        if state == "unfound":  #不在intron上////即在基因内,但在转录本上游或者下游
                            if strand1 == "+":
                                if location < Exon_start[0]:
                                    state = "5UTR"
                                if location > Exon_end[-1]:
                                    state = "3UTR"
                            else:
                                if location < Exon_start[0]:
                                    state = "3UTR"
                                if location > Exon_end[-1]:
                                    state = "5UTR"            

                state_record.append([bonafound,genetype,state,distance_to_tss,transcript_type,transcript_level,information + "\t" + state+"\t"+transcript_id+"\t"+transcript_type])                     
            else:
                continue
    if not state_record:
        if chrome1+strand1 in anno_dic:
            for gene in anno_dic[chrome1+strand1]:
                genetype = gene[2]
                if genetype == "protein_coding":
                    bonafound = "YES"
                else:
                    bonafound = "NO"
                start2 = gene[10]
                end2 = gene[11]
                information = "\t".join(gene[0:3])
                transcript_id = gene[12]
                transcript_type = gene[13]
                transcript_level = gene[14]
                if end2 > start1 and end1 > start2:
                    if strand1 == "+":
                        gene_start = gene[3] + 10
                        #distance_to_tss = abs(start1 - gene_start)
                        distance_to_tss = start1 - gene_start
                    else:
                        gene_start = gene[4] - 10 
                        #distance_to_tss = abs(end1 - gene_start)  
                        distance_to_tss = gene_start - end1
                    state = "Promoter"
                    state_record.append([bonafound,genetype,state,distance_to_tss,transcript_type,transcript_level,information + "\t" + state+"\t"+transcript_id+"\t"+transcript_type])          
                else:
                    continue
    #print(state_record)
    return(state_record)

        
if __name__ == "__main__":
    geno_anno_file = sys.argv[1]
    peakfile = sys.argv[2]
    info_dic = load_anno_file(geno_anno_file)
    infile = open(peakfile,"r")

    priority_1 = {'YES': 0, 'NO': 1} #### 确切匹配
    priority_2 = {'protein_coding':0,'lncRNA':1,'pseudogene':2}
    def sortgenetype(item):
        gene_type = item
        for keyword, priority in priority_2.items():
            if keyword in gene_type:
                return priority
        return 999
    priority_3 = {'5UTR': 0, 'CDS': 0, '3UTR': 0, 'Exon': 3, 'Intron': 4, 'Promoter': 5,'unfound':6} #### 匹配到多个优先顺序5UTR CDS 3UTR Exon Intron Promoter
    priority_4 = {'Ensembl_canonical':0,'Basic':1,'Others':2} ### 转录本等级

    for line in infile:
        info = line.rstrip().split()
        chrome = info[0]
        try:
            start = int(info[1])
            end = int(info[2])
        except:
            print(line.rstrip()+"\t"+"\t".join(["disTSS","ensg","gene_name","gene_type","location","enst","transcript_type"]))
            continue
        strand = info[5]
        if chrome+strand not in info_dic:
            print(line.rstrip()+ "\t" + "999999" + "\t" +"\t".join(["Intergenic"]*6))
        else:
            annotations = got_annotation(info_dic,[chrome,start,end,strand])
            #print(annotations)
            if len(annotations) == 0:
                print(line.rstrip()+"\t" + "999999" + "\t" +"\t".join(["Intergenic"]*6))
            else:
                #[bonafound,genetype,state,distance_to_tss,transcript_type,transcript_level]
                final_anno_sorted = sorted(annotations, key=lambda x: (priority_1[x[0]],sortgenetype(x[1]), priority_3[x[2]],priority_4[x[5]],x[3]))
                print(line.rstrip()+"\t"+ str(final_anno_sorted[0][3]) + "\t" +final_anno_sorted[0][6])
        
    infile.close()
