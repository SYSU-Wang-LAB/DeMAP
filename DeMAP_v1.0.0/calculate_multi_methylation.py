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
group1.add_argument("-max_A_ratio","--max_A_ratio",type=float, default=0.4,help="Filter sites with max Acounts reads /total reads > max_A_ratio.(default:0.4)\n"
                    "This parameter is used to filter out modified sites resulting from incomplete conversion.")
group1.add_argument("-A_neighbor","--A_neighbor",type=int, default=50,help="Number of A bases in the vicinity on both sides as background. (default:50 indicating 25 on each side)")
group1.add_argument("-pre", "--outname_prefix", type=str,  help="Output name prefix")
group1.add_argument("-o", "--outputdir", type=str, help="Output directory")

group2 = parser.add_argument_group("thresholds for output m1A,m6Am,m6A")
group2.add_argument("-fdr","--fdr",type=float, default=0.005,help="FDR thresholds for a site to be methylated. (default:0.005)")

#m1A
group2.add_argument("-m1A_fdr","--m1A_fdr",type=float, default=0.005,help="For m1A calling, FDR thresholds for mutation significance. (default:0.005)")
group2.add_argument("-mut_ratio","--mut_ratio",type=float, default=0.03,help="For m1A calling, thresholds for other mutation ratio. T+C counts/total counts. (default:0.03)\n"
                    "During reverse transcription, m1A modifications cause sequencing error.")
group2.add_argument("-trunc_ratio","--trunc_ratio",type=float, default=-1.0,help="For m1A calling, thresholds for truncation ratio:(Downstream - Upstream)/Downstream.\n"
                    "Default is -1.0 (truncation ignored). Set to >0 (e.g., 0.05) to enable truncation filtering.")
group2.add_argument("-read_drop","--read_drop",type=int, default=0,help="For m1A calling, Downstream coverage - Upstream coverage >= read_drop in a truncated site. (default:0)\n")
group2.add_argument("-min_C","--min_C",type=int, default=0,help="Minimum C mutation counts to filter out single-base SNPs for m1A. (default:0)")
group2.add_argument("-min_T","--min_T",type=int, default=0,help="Minimum T mutation counts to filter out single-base SNPs for m1A. (default:0)")
group2.add_argument("-min_t_pct","--min_t_pct",type=float, default=0.2,help="Minimum ratio T/(C+T). (default:0.2)")
# filter heterozygous SNP
#group2.add_argument("-min_mut_pct","--min_mut_pct",type=float, default=0.01,help="Minimum ratio of C/(C+T) and T/(C+T) to avoid high coverage sequencing sequencing error. (default:0.01)")


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

min_C_val = args.min_C
min_T_val = args.min_T
min_t_pct = args.min_t_pct
#min_mut_pct = args.min_mut_pct

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

def get_neighbor_indices(current_index, maxindex, flank_n):
    if current_index < flank_n:
        min_index = 0
        max_index = current_index + flank_n + (flank_n - current_index)
    elif maxindex - current_index < flank_n:
        min_index = current_index - flank_n - (flank_n - (maxindex - current_index))
        max_index = maxindex
    else:
        min_index = current_index - flank_n
        max_index = current_index + flank_n

    return [i for i in range(min_index, max_index + 1) if i != current_index]

def calculate_m6Am_significance(instring, cleavage_ratio, total_count):
    parts = instring.split(",")
    alignments = parts[1:] if len(parts) > 1 else []

    Uncovrted_A = 0
    Uncovrted_A_at_readstart = 0
    base_at_readstart = 0
    all_read_start = 0
    
    for alignment in alignments:
        info = alignment.split("_")
        if len(info) < 4:
            continue
            
        if int(info[0]) >= softlen:
            total_count -= 1
            continue
            
        base = info[1]
        position = int(info[2])
        readlength = int(info[3])
        
        is_A = (base == "A")
        if is_A:
            Uncovrted_A += 1
            
        if is_A or position == 1:
            all_read_start += 1
            
        if position <= math.ceil(readlength * start_percentile):
            base_at_readstart += 1
            if is_A:
                Uncovrted_A_at_readstart += 1

    cleav_count = (total_count - Uncovrted_A) * cleavage_ratio
    m6Amlevel = round(Uncovrted_A / all_read_start, 3) if all_read_start > 0 else 0
    m6Amlevel_str = f"{Uncovrted_A}_{all_read_start}_{m6Amlevel}"
    
    evaluate_read_start_at_first_base = all_read_start - cleav_count
    if evaluate_read_start_at_first_base < Uncovrted_A:
        evaluate_read_start_at_first_base = Uncovrted_A
    evaluate_m6Amlevel = round(Uncovrted_A / evaluate_read_start_at_first_base, 3) if evaluate_read_start_at_first_base > 0 else 0
    
    if Uncovrted_A == 0:
        p_m6Am1 = 1.0
    else:
        p_m6Am1 = scipy.stats.binomtest(Uncovrted_A_at_readstart, n=Uncovrted_A, alternative='greater', p=start_percentile).pvalue
        
    return p_m6Am1, Uncovrted_A_at_readstart, m6Amlevel_str, evaluate_m6Amlevel


def calculate_m1A_truncation(instring, strand, flank_n):
    """
    Returns the maximum termination ratio (stop_ratio) within the window and 
    the corresponding maximum read drop count (max_readdrop).
    """
    coveragelist = instring.split(",")
    coverage = [int(i) for i in coveragelist]
    adjacent_10 = coverage[(flank_n-5):(flank_n+6)]
    
    if strand == "-":
        adjacent_10.reverse()

    max_stop_ratio = 0.0
    max_readdrop = 0
    
    for i in range(0, len(adjacent_10) - 1):
        j = i + 1
        readdrop = adjacent_10[j] - adjacent_10[i]
        
        if adjacent_10[j] <= adjacent_10[i] or adjacent_10[j] == 0:
            stop_ratio = 0.0
        else:
            stop_ratio = round(readdrop / adjacent_10[j], 3)
            
        if stop_ratio > max_stop_ratio:
            max_stop_ratio = stop_ratio
        if readdrop > max_readdrop:
            max_readdrop = readdrop

    return max_stop_ratio, max_readdrop


def A_site_pileup_state(filein,fileout):
    OUTPUT = open(fileout,"w")
    OUTPUT.write("\t".join(["chrome","start","end","sitename","motif","strand","Acount","AGcount","ratio","background_ratio","methy_level","TCcount","coverage","mut_ratio","mut_ratio_background","pmut","A_in_readstart","m6Am_level","cleavage_ratio","estimated_m6Amlevel",
    "stop_ratio","max_readdrop","p_methylation","p_m6Am","Tcount","Ccount"])+"\n")
    pileup_dic = OrderedDict()
    Asite_dic = OrderedDict()

    if filein.endswith(".gz"):
        infile = gzip.open(filein, "rt")
    else:
        infile = open(filein,"r")

    for line in infile:
        pileupinfo = line.rstrip().split()
        chrome = pileupinfo[0]
        position = pileupinfo[1]
        start = str(int(position) - 1)
        end = position
        strand = pileupinfo[2]
        ref = pileupinfo[3]
        chrome_strand = chrome+strand
        
        Tcount, Ccount = int(pileupinfo[5]), int(pileupinfo[6])
        Acount,Gcount,Totalcount = int(pileupinfo[4]),int(pileupinfo[7]),int(pileupinfo[8])
        AGcount = Acount + Gcount
        Otherbase = Tcount + Ccount
        read_max_A_ratio = pileupinfo[9]
        base_position_in_read = pileupinfo[10]

        if (chrome_strand not in pileup_dic) and len(pileup_dic) == 0:  
            pileup_dic[chrome_strand] = OrderedDict()
            pileup_dic[chrome_strand][int(position)] = [chrome,position,Totalcount,base_position_in_read]
            Asite_dic[chrome_strand] = OrderedDict()
            if ref == "A" and AGcount/Totalcount >= AG_Ratio and AGcount >= AG_Coverage:
                motif = get_motif(chrome,position,strand)
                sitename = "_".join([chrome,start,end,strand,motif])
                Ratio = round(Acount/AGcount,3)
                Asite_dic[chrome_strand][int(position)] = [chrome,start,end,sitename,motif,strand,Acount,AGcount,Ratio,Otherbase,Totalcount,read_max_A_ratio,base_position_in_read,Tcount,Ccount]

        elif (chrome_strand in pileup_dic):
            pileup_dic[chrome_strand][int(position)] = [chrome,position,Totalcount,base_position_in_read]
            if ref == "A" and AGcount/Totalcount >= AG_Ratio and AGcount >= AG_Coverage:
                motif = get_motif(chrome,position,strand)
                sitename = "_".join([chrome,start,end,strand,motif])
                Ratio = round(Acount/AGcount,3)
                Asite_dic[chrome_strand][int(position)] = [chrome,start,end,sitename,motif,strand,Acount,AGcount,Ratio,Otherbase,Totalcount,read_max_A_ratio,base_position_in_read,Tcount,Ccount]

        elif (chrome_strand not in pileup_dic) and len(pileup_dic) != 0: 
            previous_chrome = list(pileup_dic.keys())[0]
            
            asite_items = list(Asite_dic[previous_chrome].items())
            all_Asite_number = len(asite_items)
            
            if all_Asite_number < A_neighbor + 1: 
                pileup_dic = OrderedDict()
                Asite_dic = OrderedDict()
                pileup_dic[chrome_strand] = OrderedDict()
                pileup_dic[chrome_strand][int(position)] = [chrome,position,Totalcount,base_position_in_read]
                Asite_dic[chrome_strand] = OrderedDict()
                if ref == "A" and AGcount/Totalcount >= AG_Ratio and AGcount >= AG_Coverage:
                    motif = get_motif(chrome,position,strand)
                    sitename = "_".join([chrome,start,end,strand,motif])
                    Ratio = round(Acount/AGcount,3)
                    Asite_dic[chrome_strand][int(position)] = [chrome,start,end,sitename,motif,strand,Acount,AGcount,Ratio,Otherbase,Totalcount,read_max_A_ratio,base_position_in_read,Tcount,Ccount]

            else:
                chrome_A_ratio = sum(v[6] for _, v in asite_items) / sum(v[7] for _, v in asite_items)
                chrome_mut_ratio = sum(v[9] for _, v in asite_items) / sum(v[10] for _, v in asite_items)
                
                pileup_sub = pileup_dic[previous_chrome]
                
                for current_index, (A_position, value) in enumerate(asite_items):
                    strand = value[5]
                    a_count = value[6]
                    ag_count = value[7]
                    a_ratio = value[8]
                    total_count = value[10]
                    Unconverted_read_ratio = float(value[11])
                    base_align_position_to_read = value[12]
                    t_count_val = value[13]
                    c_count_val = value[14]
                    
                    if a_count >= A_Coverage and ag_count >= AG_Coverage and a_ratio >= A_Ratio and Unconverted_read_ratio < max_A_ratio:  

                        ##### call methylation site
                        neighbor_indices = get_neighbor_indices(current_index, all_Asite_number - 1, flank_n)
                        
                        A_in_neighbor = 0
                        AG_in_neighbor = 0
                        Other_in_neighbor = 0
                        Coverage_in_neighbor = 0
                        
                        for idx in neighbor_indices:
                            v = asite_items[idx][1]
                            if not (v[8] >= chrome_A_ratio * enrich_fold and v[6] >= A_Coverage):
                                A_in_neighbor += v[6]
                                AG_in_neighbor += v[7]
                                Other_in_neighbor += v[9]
                                Coverage_in_neighbor += v[10]
                            #A_in_neighbor += v[6]
                            #AG_in_neighbor += v[7]
                            #Other_in_neighbor += v[9]
                            #Coverage_in_neighbor += v[10]

                        A_ratio_in_neighbor = round(A_in_neighbor/AG_in_neighbor,3) if A_in_neighbor > 0 else round(chrome_A_ratio,3)
                        p_methylation = scipy.stats.binomtest(a_count, n=ag_count, alternative='greater', p = A_ratio_in_neighbor).pvalue
                        methy_level = (a_ratio - A_ratio_in_neighbor) / (1 - A_ratio_in_neighbor) if 1 - A_ratio_in_neighbor > 0 else 0

                        ##### call m6Am
                        cleavages = []
                        for adjacentsite in range(A_position-10, A_position+11):
                            if adjacentsite == A_position:
                                continue
                            p_info = pileup_sub.get(adjacentsite)
                            if p_info:
                                cleavage = p_info[3].split(",")[0]
                                if ";" in cleavage and int(cleavage.split(";")[1]) >= AG_Coverage:
                                    cleavages.append(cleavage)
                                    
                        if len(cleavages) < 5: 
                            cleavage_ratio = 0
                        else:
                            rates = np.array([int(n)/int(d) for n, d in [x.split(';') for x in cleavages]])
                            median = np.median(rates)
                            mad = scipy.stats.median_abs_deviation(rates)
                            modified_z = 0.6745 * (rates - median) / (mad + 1) if mad == 0 else 0.6745 * (rates - median) / mad
                            valid_mask = modified_z < 3
                            nums = np.array([int(x.split(';')[0]) for x in cleavages])
                            dens = np.array([int(x.split(';')[1]) for x in cleavages])
                            cleavage_ratio = nums[valid_mask].sum() / dens[valid_mask].sum() if dens[valid_mask].sum() > 0 else 0

                        p_m6Am,A_in_readstart,m6Am_level,evaluate_m6Amlevel = calculate_m6Am_significance(base_align_position_to_read,cleavage_ratio,total_count)

                        ##### call m1A
                        Othermut = value[9]
                        Coverage = value[10]
                        mutation_ratio = round(Othermut/Coverage,4)
                        mutation_ratio_in_neighbor = Other_in_neighbor/Coverage_in_neighbor if Other_in_neighbor > 0 else chrome_mut_ratio
                        p_mutation = scipy.stats.binomtest(Othermut, n=Coverage, alternative='greater', p = mutation_ratio_in_neighbor).pvalue

                        adjacentsite_coverage = ",".join(str(pileup_sub[adj][2]) if adj in pileup_sub else "0" for adj in range(A_position-flank_n, A_position+flank_n+1))
                        # truncation
                        stop_ratio, max_readdrop = calculate_m1A_truncation(adjacentsite_coverage, strand, flank_n)

                        ###### output result
                        if A_ratio_in_neighbor > 0 and a_ratio/A_ratio_in_neighbor >= enrich_fold:
                            outlist = value[0:9] + [A_ratio_in_neighbor,methy_level] + value[9:11] + [mutation_ratio,mutation_ratio_in_neighbor,p_mutation,
                            A_in_readstart,m6Am_level,cleavage_ratio,evaluate_m6Amlevel,stop_ratio,max_readdrop,p_methylation,p_m6Am,t_count_val,c_count_val]
                            outlist = list(map(str,outlist))
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
                    Asite_dic[chrome_strand][int(position)] = [chrome,start,end,sitename,motif,strand,Acount,AGcount,Ratio,Otherbase,Totalcount,read_max_A_ratio,base_position_in_read,Tcount,Ccount]

    #### last chrome+strand
    if pileup_dic:
        lastchrome = list(pileup_dic.keys())[0]
        asite_items = list(Asite_dic[lastchrome].items())
        all_Asite_number = len(asite_items)    
        if all_Asite_number >= A_neighbor + 1:
            chrome_A_ratio = sum(v[6] for _, v in asite_items) / sum(v[7] for _, v in asite_items)
            chrome_mut_ratio = sum(v[9] for _, v in asite_items) / sum(v[10] for _, v in asite_items)
            
            pileup_sub = pileup_dic[lastchrome]
            
            for current_index, (A_position, value) in enumerate(asite_items):
                strand = value[5]
                a_count = value[6]
                ag_count = value[7]
                a_ratio = value[8]
                total_count = value[10]
                Unconverted_read_ratio = float(value[11])
                base_align_position_to_read = value[12]
                t_count_val = value[13]
                c_count_val = value[14]
                if a_count >= A_Coverage and ag_count >= AG_Coverage and a_ratio >= A_Ratio and Unconverted_read_ratio < max_A_ratio: 

                    ##### call methylation site
                    neighbor_indices = get_neighbor_indices(current_index, all_Asite_number - 1, flank_n)
                    
                    A_in_neighbor = 0
                    AG_in_neighbor = 0
                    Other_in_neighbor = 0
                    Coverage_in_neighbor = 0
                    
                    for idx in neighbor_indices:
                        v = asite_items[idx][1]
                        if not (v[8] >= chrome_A_ratio * enrich_fold and v[6] >= A_Coverage): 
                            A_in_neighbor += v[6]
                            AG_in_neighbor += v[7]
                            Other_in_neighbor += v[9]
                            Coverage_in_neighbor += v[10]
                        #A_in_neighbor += v[6]
                        #AG_in_neighbor += v[7]
                        #Other_in_neighbor += v[9]
                        #Coverage_in_neighbor += v[10]

                    A_ratio_in_neighbor = round(A_in_neighbor/AG_in_neighbor,3) if A_in_neighbor > 0 else round(chrome_A_ratio,3)
                    p_methylation = scipy.stats.binomtest(a_count, n=ag_count, alternative='greater', p = A_ratio_in_neighbor).pvalue
                    methy_level = (a_ratio - A_ratio_in_neighbor) / (1 - A_ratio_in_neighbor) if 1 - A_ratio_in_neighbor > 0 else 0

                    ##### call m6Am
                    cleavages = []
                    for adjacentsite in range(A_position-10, A_position+11):
                        if adjacentsite == A_position:
                            continue
                        p_info = pileup_sub.get(adjacentsite)
                        if p_info:
                            cleavage = p_info[3].split(",")[0]
                            if ";" in cleavage and int(cleavage.split(";")[1]) >= AG_Coverage:
                                cleavages.append(cleavage)
                                
                    if len(cleavages) < 5:
                        cleavage_ratio = 0
                    else:
                        rates = np.array([int(n)/int(d) for n, d in [x.split(';') for x in cleavages]])
                        median = np.median(rates)
                        mad = scipy.stats.median_abs_deviation(rates)
                        modified_z = 0.6745 * (rates - median) / (mad + 1) if mad == 0 else 0.6745 * (rates - median) / mad
                        valid_mask = modified_z < 3
                        nums = np.array([int(x.split(';')[0]) for x in cleavages])
                        dens = np.array([int(x.split(';')[1]) for x in cleavages])
                        cleavage_ratio = nums[valid_mask].sum() / dens[valid_mask].sum() if dens[valid_mask].sum() > 0 else 0

                    p_m6Am,A_in_readstart,m6Am_level,evaluate_m6Amlevel = calculate_m6Am_significance(base_align_position_to_read,cleavage_ratio,total_count)

                    ##### call m1A
                    Othermut = value[9]
                    Coverage = value[10]
                    mutation_ratio = round(Othermut/Coverage,4)
                    mutation_ratio_in_neighbor = Other_in_neighbor/Coverage_in_neighbor if Other_in_neighbor > 0 else chrome_mut_ratio
                    p_mutation = scipy.stats.binomtest(Othermut, n=Coverage, alternative='greater', p = mutation_ratio_in_neighbor).pvalue

                    adjacentsite_coverage = ",".join(str(pileup_sub[adj][2]) if adj in pileup_sub else "0" for adj in range(A_position-flank_n, A_position+flank_n+1))
                    stop_ratio, max_readdrop = calculate_m1A_truncation(adjacentsite_coverage, strand, flank_n)

                    ###### output result
                    if A_ratio_in_neighbor > 0 and a_ratio/A_ratio_in_neighbor >= enrich_fold:
                        outlist = value[0:9] + [A_ratio_in_neighbor,methy_level] + value[9:11] + [mutation_ratio,mutation_ratio_in_neighbor,p_mutation,A_in_readstart,m6Am_level,cleavage_ratio,evaluate_m6Amlevel,stop_ratio,max_readdrop,p_methylation,p_m6Am,t_count_val,c_count_val]
                        outlist = list(map(str,outlist))
                        OUTPUT.write("\t".join(outlist)+"\n")
                        
    infile.close() 
    OUTPUT.close()


def output_result(filein,outfile):

    df = pd.read_csv(filein, sep="\t")

    df["p_methylation_fdr"] = multipletests(df["p_methylation"], method="fdr_bh")[1]
    df["p_m6Am_fdr"] = multipletests(df["p_m6Am"], method="fdr_bh")[1]
    df["pmut_fdr"] = multipletests(df["pmut"], method="fdr_bh")[1]

    # --------------------------------
    mask_precondition = (df["p_methylation_fdr"] <= methy_fdr)

    mask_m6Am = (
        mask_precondition &
        (df["p_m6Am_fdr"] <= m6Am_fdr) &
        (df["A_in_readstart"] >= A_start) &
        (df["A_in_readstart"] / df["Acount"] >= A_start_ratio)  
    )

    total_ct = df["Tcount"] + df["Ccount"]
    valid_ct_mask = (total_ct > 0)
    
    mask_hetero_snp = (
        valid_ct_mask &
        (df["Ccount"] >= min_C_val) & 
        (df["Tcount"] >= min_T_val) & 
        ((df["Tcount"] / total_ct) >= min_t_pct) #& 
        #((df["Ccount"] / total_ct) > min_mut_pct) & 
        #((df["Tcount"] / total_ct) > min_mut_pct) 
    )

    mask_m1A = (
        mask_precondition &
        ~mask_m6Am &
        (df["pmut_fdr"] <= m1A_fdr) &
        (df["mut_ratio"] >= mut_ratio) &
        mask_hetero_snp &
        (df["stop_ratio"] >= trunc_ratio) &
        (df["max_readdrop"] >= read_drop)
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


    df["m6Am_level"] = df["m6Am_level"].astype(object)
    df["estimated_m6Amlevel"] = df["estimated_m6Amlevel"].astype(object)

    df.loc[df['mod_type'] == "m6A", 'm6Am_level'] = "NA"
    df.loc[df['mod_type'] == "m1A", 'm6Am_level'] = "NA"
    df.loc[df['mod_type'] == "m6A", 'estimated_m6Amlevel'] = "NA"
    df.loc[df['mod_type'] == "m1A", 'estimated_m6Amlevel'] = "NA"

    df_filtered = df[df["mod_type"] != "none"]
    cols_to_drop = ['m6Am_level', 'cleavage_ratio', 'p_methylation', 'p_m6Am', 'Tcount', 'Ccount']
    df_filtered = df_filtered.drop(columns=cols_to_drop)
    df_filtered.to_csv(outfile, sep="\t", index=False)


if __name__ == "__main__":
    print(genomefile)
    genomedic = get_genomedic(genomefile)
    #if os.path.exists(out_tmp):
    #    pass
    #else:
    A_site_pileup_state(pileupfile,out_tmp)
    output_result(out_tmp,outfile)
