#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Zehong Wu wang_lab, Sun Yat-sen University"""
"""Feb, 2025"""
"""Email: wuzh93@mail2.sysu.edu.cn"""
"""Usage: This program is used for calling methylation sites"""


from collections import OrderedDict
import scipy.stats
import argparse
from Bio import SeqIO
import math
import pandas as pd
from statsmodels.stats.multitest import multipletests
import os
import numpy as np
import gzip

parser = argparse.ArgumentParser(description = "Calculate methylation for all A sites",formatter_class=argparse.RawTextHelpFormatter)
group1 = parser.add_argument_group("Calling methylation")
group1.add_argument("-pileup", "--pileup", type=str, required=True, help = "Pileup file")
group1.add_argument("-fa","--fasta",type=str,help="Unconverted genome file in fasta format")
group1.add_argument("-AGcov","--AGCoverage",type=int, default=15,help="minimum A + G counts. (default:15)")
group1.add_argument("-Acov","--ACoverage",type=int, default=5,help="minimum A counts. (default:5)")
group1.add_argument("-ratio","--ARatio",type=float, default=0.05,help="minimum methylation level: A counts / AG counts. (default:0.05)")
group1.add_argument("-enrich","--foldenrich",type=float, default=3,help="minimum fold change of methylation level / background methylation level. (default:3)")
group1.add_argument("-AGRatio","--AGRatio",type=float, default=0.6,help="Filter sites with A+G counts/total counts < AGRatio. (default:0.6)\n"
                    "This parameter is used to filter out mutation sites originating from SNPs.")
group1.add_argument("-max_A_ratio","--max_A_ratio",type=float, default=0.6,help="Filter sites with max Acounts reads /total reads > max_A_ratio.(default:0.6)\n"
                    "This parameter is used to filter out modified sites resulting from incomplete conversion.")
group1.add_argument("-A_neighbor","--A_neighbor",type=int, default=50,help="Number of A bases in the vicinity on both sides as background. (default:50 indicating 25 on each side)")
group1.add_argument("-pre", "--outname_prefix", type=str,  help="Output name prefix")
group1.add_argument("-o", "--outputdir", type=str, help="Output directory")

group2 = parser.add_argument_group("thresholds for output m1A,m6Am,m6A")
group2.add_argument("-fdr","--fdr",type=float, default=0.005,help="FDR thresholds for a site to be methylated. (default:0.005)")
#m1A
group2.add_argument("-m1A_fdr","--m1A_fdr",type=float, default=0.005,help="For m1A calling, FDR thresholds for m1A. (default:0.005)")
group2.add_argument("-mut_ratio","--mut_ratio",type=float, default=0.03,help="For m1A calling, thresholds for other mutation ratio. T+C counts/total counts. (default:0.03)\n"
                    "During reverse transcription, m1A modifications cause sequencing error.")
group2.add_argument("-trunc_ratio","--trunc_ratio",type=float, default=0.05,help="For m1A calling, thresholds for truncation ratio:(Downstream - Upstream)/Downstream. (default:0.05)\n"
                    "During reverse transcription, m1A modifications cause sequencing truncations. 0.05 means 5%% truncations")
group2.add_argument("-read_drop","--read_drop",type=int, default=10,help="For m1A calling, Downstream coverage - Upstream coverage >= read_drop in a truncated site. (default:10)\n")

#m6Am
group2.add_argument("-softlen","--softlen",type=int, default=30,help="For m6Am calling, only using read with 5'end softclip length < softlen. (default:30)")
group2.add_argument("-m6Am_fdr","--m6Am_fdr",type=float, default=0.005,help="For m6Am calling, FDR thresholds for m6Am. (default:0.005)")
group2.add_argument("-start_percentile","--start_percentile",type=float, default=0.05,help="For m6Am calling, If a base is in the first 5%% of the read, it is considered to be at the start.(default:0.05)")
group2.add_argument("-A_start","--A_start",type=int, default=3,help="For m6Am calling, minimum unconverted A bases located in read start. (default:3)")
group2.add_argument("-A_start_ratio","--A_start_ratio",type=float, default= 0.6,help="For m6Am calling, minimum ratio of unconverted A in read start/ all unconverted A. (default:0.6)")

args = parser.parse_args()
pileupfile = args.pileup
genomefile = args.fasta
AG_Coverage = args.AGCoverage
A_Coverage = args.ACoverage
A_Ratio = args.ARatio
enrich_fold = args.foldenrich
AG_Ratio = args.AGRatio
max_A_ratio = args.max_A_ratio
A_neighbor = args.A_neighbor
flank_n = A_neighbor // 2

methy_fdr = args.fdr
m1A_fdr = args.m1A_fdr
mut_ratio = args.mut_ratio
trunc_ratio = args.trunc_ratio
read_drop = args.read_drop

softlen = args.softlen
m6Am_fdr = args.m6Am_fdr
start_percentile = args.start_percentile
A_start = args.A_start
A_start_ratio = args.A_start_ratio

outputdir = args.outputdir
outname_prx = args.outname_prefix

out_tmp = outputdir+"/"+outname_prx+"_methylation_tmp"
outfile = outputdir+"/"+outname_prx+"_methylation_result.bed"


def get_genomedic(filein):
    genomedic = {}
    for seq in SeqIO.parse(filein,"fasta"):
        if seq.id not in genomedic:
            genomedic[str(seq.id)] = str(seq.seq).upper()
    return genomedic

def get_motif(chrome,position,strand):
    convertdict = {"A":"T","C":"G","G":"C","T":"A","N":"N"}
    global genomedic
    if int(position)-3 < 0 or int(position) + 2 > len(genomedic[chrome]):
        motif = "NNANN"
    else:
        motif = genomedic[chrome][ (int(position)-3) : (int(position) + 2)]
    if strand == "-":
        convertmotif = [convertdict[base] for base in motif]
        motif = "".join(convertmotif[::-1])
    return motif

def get_nearest_Asite_index(Apos,maxindex,flank_n,all_Asite_position_list):
    current_index = all_Asite_position_list.index(Apos)
    if current_index < flank_n:
        min_index = 0
        max_index = current_index + flank_n + (flank_n - current_index)
    elif maxindex - current_index < flank_n:
        min_index = current_index - flank_n - (flank_n - (maxindex - current_index))
        max_index = maxindex
    else:
        min_index = current_index - flank_n
        max_index = current_index + flank_n

    outindex = list(range(min_index,max_index + 1))
    outindex.remove(current_index)
    return([all_Asite_position_list[i] for i in outindex])

def calculate_m6Am_significance(instring,cleavage_ratio,total_count):
    """
    2;68,0_G_62_80,0_G_62_62,0_G_46_55,0_G_56_66,0_G_67_86,0_G_76_87
    """
    alignments = instring.split(",")[1:]
    total_count = total_count - len([alignment for alignment in alignments if int(alignment.split("_")[0]) >= softlen ])

    alignments = [alignment for alignment in alignments if int(alignment.split("_")[0]) < softlen ] # only using read with 5'header softclip less than 3bp

    clean_instring = ",".join(alignments)
    Uncovrted_A = clean_instring.count('A')
    Uncovrted_A_at_readstart = 0
    base_at_readstart = 0
    
    # m6Am level: the number of Unconverted A happen in the first alignment base of read/the number of read alignment start at this base
    #A_in_first_base = clean_instring.count('A_1_')
    #read_start_at_first_base = clean_instring.count('_1_')
    #methy_ratio_at_firstbase = round(A_in_first_base/read_start_at_first_base,3) if read_start_at_first_base > 0 else 0
    #m6Amlevel = str(A_in_first_base) + "_" + str(read_start_at_first_base)  + "_" + str(methy_ratio_at_firstbase)

    # evaluate m6Am level: the number of Unconverted A happen in the first base of read/ (the number of read start at this base - read cleavage at this base)
    cleav_count = (total_count - Uncovrted_A)*cleavage_ratio
    all_read_start = len([alignment for alignment in alignments if 'A' in alignment or '_1_' in alignment])
    m6Amlevel = round(Uncovrted_A/all_read_start,3) if all_read_start > 0 else 0
    m6Amlevel = str(Uncovrted_A) + "_" + str(all_read_start) + "_" + str(m6Amlevel)
    evaluate_read_start_at_first_base = all_read_start - cleav_count if all_read_start - cleav_count > Uncovrted_A else Uncovrted_A
    evaluate_m6Amlevel = round(Uncovrted_A/evaluate_read_start_at_first_base,3) if evaluate_read_start_at_first_base > 0 else 0
    
    for alignment in alignments:
        aligninfo = alignment.split("_")
        base, position, readlength = aligninfo[1], int(aligninfo[2]), int(aligninfo[3])
        if position <= math.ceil(readlength*start_percentile):
            base_at_readstart += 1
            if base == "A":
                Uncovrted_A_at_readstart += 1
                
    if 'A' not in clean_instring:
        p_m6Am1 = 1
    else:
        Random_probability = start_percentile
        #Random_probability = Total_position_Uncovrted_A_at_readstart/Total_read_length_with_uncoverted_A
        p_m6Am1 = scipy.stats.binomtest(Uncovrted_A_at_readstart, n=Uncovrted_A, alternative='greater', p = Random_probability).pvalue
    return(p_m6Am1,Uncovrted_A_at_readstart,m6Amlevel,evaluate_m6Amlevel)

def calculate_m1A_significance(instring,strand,flank_n,p_mut):
    coveragelist = instring.split(",")
    coverage = [int(i) for i in coveragelist]
    adjacent_10 = coverage[(flank_n-5):(flank_n+6)]
    if strand == "-":
        adjacent_10.reverse()

    if min(adjacent_10) < 20: ## without sufficient coverage for truncation detection or mutation significance
        p_m1A = 1
        return("0",p_m1A)
    

    truncation_record = []
    max_stop_ratio = 0
    for i in range(0,len(adjacent_10) - 1):
        j = i + 1
        readdrop = adjacent_10[j] - adjacent_10[i]
        stop_ratio = round((adjacent_10[j] - adjacent_10[i])/adjacent_10[j],3)
        if stop_ratio > max_stop_ratio:
            max_stop_ratio = stop_ratio

        if readdrop >= read_drop and stop_ratio >= trunc_ratio and stop_ratio < 1 - trunc_ratio: # stop_ratio >= 1 - trunc_ratio may be a TSS
            truncation_record.append([readdrop,stop_ratio])

    truncation_record.sort(key=lambda x:x[1],reverse=True)
    if p_mut > m1A_fdr:
        p_m1A = 1
        if len(truncation_record) == 0:
            return(max_stop_ratio,p_m1A)
        else:
            return(truncation_record[0][1],p_m1A)
            
    if len(truncation_record) == 0:
        p_m1A = 1
        return(max_stop_ratio,p_m1A)
    else:
        stop_ratio = truncation_record[0][1]
        p_m1A = p_mut/(stop_ratio * 100)
        return(stop_ratio,p_m1A)


def A_site_pileup_state(filein,fileout):
    OUTPUT = open(fileout,"w")
    OUTPUT.write("\t".join(["chrome","start","end","sitename","motif","strand","Acount","AGcount","ratio","background_ratio","methy_level","TCcount","coverage","mut_ratio","A_in_readstart","m6Am_level","cleavage_ratio","estimated_m6Amlevel",
    "stop_ratio","p_methylation","p_m6Am","p_m1A"])+"\n")
    pileup_dic = OrderedDict()
    Asite_dic = OrderedDict()

    if filein.endswith(".gz"):
        infile = gzip.open(filein, "rt")
    else:
        infile = open(filein,"r")

    for line in infile:
        """
        chr1    360076  +       A       0       0       0       3       3       0.0     0;3,0_G_7_59,0_G_7_58,0_G_6_45
        chr1    360077  +       A       0       0       0       3       3       0.0     0;3,0_G_8_59,0_G_8_58,0_G_7_45
        chr1    360078  +       A       0       0       0       3       3       0.0     0;3,0_G_9_59,0_G_9_58,0_G_8_45
        """
        pileupinfo = line.rstrip().split()
        chrome = pileupinfo[0]
        position = pileupinfo[1]
        start = str(int(position) - 1)
        end = position
        strand = pileupinfo[2]
        ref = pileupinfo[3]
        chrome_strand = chrome+strand
        Acount,Gcount,Totalcount = int(pileupinfo[4]),int(pileupinfo[7]),int(pileupinfo[8])
        AGcount = Acount + Gcount
        Otherbase = int(pileupinfo[5])+int(pileupinfo[6])
        read_max_A_ratio = pileupinfo[9]
        base_position_in_read = pileupinfo[10]

        if (chrome_strand not in pileup_dic) and len(pileup_dic) == 0:  #### first line
            pileup_dic[chrome_strand] = OrderedDict()
            pileup_dic[chrome_strand][int(position)] = [chrome,position,Totalcount,base_position_in_read]
            Asite_dic[chrome_strand] = OrderedDict()
            if ref == "A" and AGcount/Totalcount >= AG_Ratio and AGcount >= AG_Coverage:
                motif = get_motif(chrome,position,strand)
                sitename = "_".join([chrome,start,end,strand,motif])
                Ratio = round(Acount/AGcount,3)
                Asite_dic[chrome_strand][int(position)] = [chrome,start,end,sitename,motif,strand,Acount,AGcount,Ratio,Otherbase,Totalcount,read_max_A_ratio,base_position_in_read]

        elif (chrome_strand in pileup_dic):
            pileup_dic[chrome_strand][int(position)] = [chrome,position,Totalcount,base_position_in_read]
            if ref == "A" and AGcount/Totalcount >= AG_Ratio and AGcount >= AG_Coverage:
                motif = get_motif(chrome,position,strand)
                sitename = "_".join([chrome,start,end,strand,motif])
                Ratio = round(Acount/AGcount,3)
                Asite_dic[chrome_strand][int(position)] = [chrome,start,end,sitename,motif,strand,Acount,AGcount,Ratio,Otherbase,Totalcount,read_max_A_ratio,base_position_in_read]

        elif (chrome_strand not in pileup_dic) and len(pileup_dic) != 0: # line meet a new chrome+stard,need to process sites in  previous chrome+stard
            previous_chrome = list(pileup_dic.keys())[0]
            all_Asite_position_list = list(Asite_dic[previous_chrome].keys())  #### position list of all A bases
            all_Asite_number = len(all_Asite_position_list)
            if all_Asite_number < A_neighbor + 1: # without sufficient background sites
                pileup_dic = OrderedDict()
                Asite_dic = OrderedDict()
                pileup_dic[chrome_strand] = OrderedDict()
                pileup_dic[chrome_strand][int(position)] = [chrome,position,Totalcount,base_position_in_read]
                Asite_dic[chrome_strand] = OrderedDict()
                if ref == "A" and AGcount/Totalcount >= AG_Ratio and AGcount >= AG_Coverage:
                    motif = get_motif(chrome,position,strand)
                    sitename = "_".join([chrome,start,end,strand,motif])
                    Ratio = round(Acount/AGcount,3)
                    Asite_dic[chrome_strand][int(position)] = [chrome,start,end,sitename,motif,strand,Acount,AGcount,Ratio,Otherbase,Totalcount,read_max_A_ratio,base_position_in_read]

            else:
                chrome_A_ratio = sum([Asite_dic[previous_chrome][pos][6] for pos in Asite_dic[previous_chrome]])/sum([Asite_dic[previous_chrome][pos][7] for pos in Asite_dic[previous_chrome]])
                chrome_mut_ratio = sum([Asite_dic[previous_chrome][pos][9] for pos in Asite_dic[previous_chrome]])/sum([Asite_dic[previous_chrome][pos][10] for pos in Asite_dic[previous_chrome]])
                for A_position,value in Asite_dic[previous_chrome].items():
                    strand = value[5]
                    a_count = value[6]
                    ag_count = value[7]
                    a_ratio = value[8]
                    total_count = value[10]
                    Unconverted_read_ratio = float(value[11])
                    base_align_position_to_read = value[12]
                    if a_count >= A_Coverage and ag_count >= AG_Coverage and a_ratio >= A_Ratio and Unconverted_read_ratio < max_A_ratio:  #### skip A without 100 adjacent A due to low sequnce depth

                        ##### call methylation site
                        nearest_keys = get_nearest_Asite_index(A_position,all_Asite_number - 1,flank_n,all_Asite_position_list)
                        nearest_keys_keep = [keepkey for keepkey in nearest_keys if not (Asite_dic[previous_chrome][keepkey][8] >= chrome_A_ratio*enrich_fold and Asite_dic[previous_chrome][keepkey][6] >= 3)]
                        A_in_neighbor = sum([Asite_dic[previous_chrome][Apos][6] for Apos in nearest_keys_keep])
                        AG_in_neighbor = sum([Asite_dic[previous_chrome][Apos][7] for Apos in nearest_keys_keep])
                        A_ratio_in_neighbor = round(A_in_neighbor/AG_in_neighbor,3) if A_in_neighbor > 0 else round(chrome_A_ratio,3)
                        p_methylation = scipy.stats.binomtest(a_count, n=ag_count, alternative='greater', p = A_ratio_in_neighbor).pvalue
                        methy_level = (a_ratio - A_ratio_in_neighbor) / (1 - A_ratio_in_neighbor) if 1 - A_ratio_in_neighbor > 0 else 0

                        ##### call m6Am
                        cleavages = [pileup_dic[previous_chrome][adjacentsite][3] if adjacentsite in pileup_dic[previous_chrome] else "0;0" for adjacentsite in range(A_position-10,A_position+11) if adjacentsite != A_position]
                        cleavages = [cleavage.split(",")[0] for cleavage in cleavages]
                        cleavages = [cleavage for cleavage in cleavages if int(cleavage.split(";")[1]) >= AG_Coverage]
                        if len(cleavages) < 5: # Insufficient number of sites for calculating the cleavage rate
                            cleavage_ratio = 0
                        else:
                            rates = np.array([int(n)/int(d) for n, d in [x.split(';') for x in cleavages]])
                            #calculate MAD and z-score
                            median = np.median(rates)
                            mad = scipy.stats.median_abs_deviation(rates)
                            modified_z = 0.6745 * (rates - median) / (mad + 1) if mad == 0 else 0.6745 * (rates - median) / mad
                            valid_mask = modified_z < 3
                            #calculate mean cleavage ratio
                            nums = np.array([int(x.split(';')[0]) for x in cleavages])
                            dens = np.array([int(x.split(';')[1]) for x in cleavages])
                            cleavage_ratio = nums[valid_mask].sum() / dens[valid_mask].sum() if dens[valid_mask].sum() > 0 else 0

                        p_m6Am,A_in_readstart,m6Am_level,evaluate_m6Amlevel = calculate_m6Am_significance(base_align_position_to_read,cleavage_ratio,total_count)

                        ##### call m1A
                        Othermut = value[9]
                        Coverage = value[10]
                        mutation_ratio = round(Othermut/Coverage,4)
                        Other_in_neighbor = sum([Asite_dic[previous_chrome][Apos][9] for Apos in nearest_keys_keep])
                        Coverage_in_neighbor = sum([Asite_dic[previous_chrome][Apos][10] for Apos in nearest_keys_keep])
                        mutation_ratio_in_neighbor = Other_in_neighbor/Coverage_in_neighbor if Other_in_neighbor > 0 else chrome_mut_ratio
                        p_mutation = scipy.stats.binomtest(Othermut, n=Coverage, alternative='greater', p = mutation_ratio_in_neighbor).pvalue

                        adjacentsite_coverage = ",".join([str(pileup_dic[previous_chrome][adjacentsite][2]) if adjacentsite in pileup_dic[previous_chrome] else "0" for adjacentsite in range(A_position-flank_n,A_position+flank_n+1)])
                        stop_ratio,p_m1A = calculate_m1A_significance(adjacentsite_coverage,strand,flank_n,p_mutation)

                        ###### output result
                        if A_ratio_in_neighbor > 0 and a_ratio/A_ratio_in_neighbor >= enrich_fold:
                            outlist = value[0:9] + [A_ratio_in_neighbor,methy_level] + value[9:11] + [mutation_ratio,A_in_readstart,m6Am_level,cleavage_ratio,evaluate_m6Amlevel,stop_ratio,p_methylation,p_m6Am,p_m1A]
                            outlist =  list(map(str,outlist))
                            OUTPUT.write("\t".join(outlist)+"\n")

                pileup_dic = OrderedDict()
                Asite_dic = OrderedDict()

                pileup_dic[chrome_strand] = OrderedDict()
                pileup_dic[chrome_strand][int(position)] = [chrome,position,Totalcount,base_position_in_read]
                Asite_dic[chrome_strand] = OrderedDict()
                if ref == "A" and AGcount/Totalcount >= AG_Ratio and AGcount >= AG_Coverage:
                    motif = get_motif(chrome,position,strand)
                    sitename = "_".join([chrome,start,end,strand,motif])
                    Ratio = round(Acount/AGcount,3)
                    Asite_dic[chrome_strand][int(position)] = [chrome,start,end,sitename,motif,strand,Acount,AGcount,Ratio,Otherbase,Totalcount,read_max_A_ratio,base_position_in_read]

    #### last chrome+strand
    lastchrome = list(pileup_dic.keys())[0]
    all_Asite_position_list = list(Asite_dic[lastchrome].keys())
    all_Asite_number = len(all_Asite_position_list)    
    if all_Asite_number >= A_neighbor + 1:
        chrome_A_ratio = sum([Asite_dic[lastchrome][pos][6] for pos in Asite_dic[lastchrome]])/sum([Asite_dic[lastchrome][pos][7] for pos in Asite_dic[lastchrome]])
        chrome_mut_ratio = sum([Asite_dic[lastchrome][pos][9] for pos in Asite_dic[lastchrome]])/sum([Asite_dic[lastchrome][pos][10] for pos in Asite_dic[lastchrome]])
        for A_position,value in Asite_dic[lastchrome].items():
            strand = value[5]
            a_count = value[6]
            ag_count = value[7]
            a_ratio = value[8]
            total_count = value[10]
            Unconverted_read_ratio = float(value[11])
            base_align_position_to_read = value[12]
            if a_count >= A_Coverage and ag_count >= AG_Coverage and a_ratio >= A_Ratio and Unconverted_read_ratio < max_A_ratio:  #### skip A without 100 adjacent A due to low sequnce depth

                ##### call methylation site
                nearest_keys = get_nearest_Asite_index(A_position,all_Asite_number - 1,flank_n,all_Asite_position_list)
                nearest_keys_keep = [keepkey for keepkey in nearest_keys if not (Asite_dic[lastchrome][keepkey][8] >= chrome_A_ratio*enrich_fold and Asite_dic[lastchrome][keepkey][6] >= 3)]
                A_in_neighbor = sum([Asite_dic[lastchrome][Apos][6] for Apos in nearest_keys_keep])
                AG_in_neighbor = sum([Asite_dic[lastchrome][Apos][7] for Apos in nearest_keys_keep])
                A_ratio_in_neighbor = round(A_in_neighbor/AG_in_neighbor,3) if A_in_neighbor > 0 else round(chrome_A_ratio,3)
                p_methylation = scipy.stats.binomtest(a_count, n=ag_count, alternative='greater', p = A_ratio_in_neighbor).pvalue
                methy_level = (a_ratio - A_ratio_in_neighbor) / (1 - A_ratio_in_neighbor) if 1 - A_ratio_in_neighbor > 0 else 0

                ##### call m6Am
                cleavages = [pileup_dic[lastchrome][adjacentsite][3] if adjacentsite in pileup_dic[lastchrome] else "0;0" for adjacentsite in range(A_position-10,A_position+11) if adjacentsite != A_position]
                cleavages = [cleavage.split(",")[0] for cleavage in cleavages]
                cleavages = [cleavage for cleavage in cleavages if int(cleavage.split(";")[1]) >= AG_Coverage]
                if len(cleavages) < 5:
                    cleavage_ratio = 0
                else:
                    rates = np.array([int(n)/int(d) for n, d in [x.split(';') for x in cleavages]])
                    #calculate MAD and z-score
                    median = np.median(rates)
                    mad = scipy.stats.median_abs_deviation(rates)
                    modified_z = 0.6745 * (rates - median) / (mad + 1) if mad == 0 else 0.6745 * (rates - median) / mad
                    valid_mask = modified_z < 3
                    #calculate mean cleavage ratio
                    nums = np.array([int(x.split(';')[0]) for x in cleavages])
                    dens = np.array([int(x.split(';')[1]) for x in cleavages])
                    cleavage_ratio = nums[valid_mask].sum() / dens[valid_mask].sum() if dens[valid_mask].sum() > 0 else 0

                p_m6Am,A_in_readstart,m6Am_level,evaluate_m6Amlevel = calculate_m6Am_significance(base_align_position_to_read,cleavage_ratio,total_count)

                ##### call m1A
                Othermut = value[9]
                Coverage = value[10]
                mutation_ratio = round(Othermut/Coverage,4)
                Other_in_neighbor = sum([Asite_dic[lastchrome][Apos][9] for Apos in nearest_keys_keep])
                Coverage_in_neighbor = sum([Asite_dic[lastchrome][Apos][10] for Apos in nearest_keys_keep])
                mutation_ratio_in_neighbor = Other_in_neighbor/Coverage_in_neighbor if Other_in_neighbor > 0 else chrome_mut_ratio
                p_mutation = scipy.stats.binomtest(Othermut, n=Coverage, alternative='greater', p = mutation_ratio_in_neighbor).pvalue

                adjacentsite_coverage = ",".join([str(pileup_dic[lastchrome][adjacentsite][2]) if adjacentsite in pileup_dic[lastchrome] else "0" for adjacentsite in range(A_position-flank_n,A_position+flank_n+1)])
                stop_ratio,p_m1A = calculate_m1A_significance(adjacentsite_coverage,strand,flank_n,p_mutation)

                ###### output result
                if A_ratio_in_neighbor > 0 and a_ratio/A_ratio_in_neighbor >= enrich_fold:
                    outlist = value[0:9] + [A_ratio_in_neighbor,methy_level] + value[9:11] + [mutation_ratio,A_in_readstart,m6Am_level,cleavage_ratio,evaluate_m6Amlevel,stop_ratio,p_methylation,p_m6Am,p_m1A]
                    outlist =  list(map(str,outlist))
                    OUTPUT.write("\t".join(outlist)+"\n")
    infile.close() 
    OUTPUT.close()



def output_result(filein,outfile):

    df = pd.read_csv(filein, sep="\t")

    df["p_methylation_fdr"] = multipletests(df["p_methylation"], method="fdr_bh")[1]
    df["p_m6Am_fdr"] = multipletests(df["p_m6Am"], method="fdr_bh")[1]
    df["p_m1A_fdr"] = multipletests(df["p_m1A"], method="fdr_bh")[1]

    # --------------------------------
    mask_precondition = (df["p_methylation_fdr"] <= methy_fdr)

    mask_m6Am = (
        mask_precondition &
        (df["p_m6Am_fdr"] <= m6Am_fdr) &
        (df["A_in_readstart"] >= A_start) &
        (df["A_in_readstart"] / df["Acount"] >= A_start_ratio)  
    )

    mask_m1A = (
        mask_precondition &
        ~mask_m6Am &
        (df["p_m1A_fdr"] <= m1A_fdr) &
        (df["mut_ratio"] >= mut_ratio)
    )

    mask_m6A = (
        mask_precondition &
        ~mask_m1A &
        ~mask_m6Am
    )

    df["mod_type"] = "none"  
    df.loc[mask_m1A, "mod_type"] = "m1A"
    df.loc[mask_m6Am, "mod_type"] = "m6Am"
    df.loc[mask_m6A, "mod_type"] = "m6A"

    df.loc[df['mod_type'] == "m6A", 'm6Am_level'] = "NA"
    df.loc[df['mod_type'] == "m1A", 'm6Am_level'] = "NA"
    df.loc[df['mod_type'] == "m6A", 'estimated_m6Amlevel'] = "NA"
    df.loc[df['mod_type'] == "m1A", 'estimated_m6Amlevel'] = "NA"

    df_filtered = df[df["mod_type"] != "none"]
    cols_to_drop = ['m6Am_level', 'cleavage_ratio', 'p_methylation', 'p_m6Am', 'p_m1A']
    df_filtered = df_filtered.drop(columns=cols_to_drop)
    df_filtered.to_csv(outfile, sep="\t", index=False)


if __name__ == "__main__":
    print(genomefile)
    genomedic = get_genomedic(genomefile)
    if os.path.exists(out_tmp):
        pass
    else:
        A_site_pileup_state(pileupfile,out_tmp)
    output_result(out_tmp,outfile)

        
