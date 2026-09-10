#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Zehong Wu wang_lab, Sun Yat-sen University"""
"""Feb, 2025"""
"""Email: wuzh93@mail2.sysu.edu.cn"""
"""Usage: This script merges the alignments mapped to the A-to-G and T-to-C converted genomes."""


import os
import argparse
import sys
import time
from time import strftime
from Bio import SeqIO
import pysam


parser = argparse.ArgumentParser(description="Merge AG and TC samfile and split them to uniq-map and multi-map reads", formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-AG", "--samAG", nargs="?", type=str, default="none", help="samfile of reads mapping to AG converted genome")
parser.add_argument("-TC", "--samTC", nargs="?", type=str, default="none", help="samfile of reads mapping to TC converted genome")
parser.add_argument("-genome", "--genome", nargs="?", type=str, default="none", help="unconverted genome")
parser.add_argument("-pre", "--outname_prefix", nargs="?", type=str, default='default', help="outname_prefix")
parser.add_argument("-o", "--outputdir", nargs="?", type=str, default="none", help="outputdir")


parser.add_argument("-d", "--diff", type=int, default=4, help="score difference to filter secondary alignments (default: 4)")
parser.add_argument("-m", "--min_score", type=int, default=20, help="minimum max-AS score required for a read to be kept (default: 20)")

args = parser.parse_args()

AG_samfile = args.samAG
TC_samfile = args.samTC
outdir = args.outputdir
out_pre = args.outname_prefix
genomefile = args.genome
diff_score = args.diff          
min_as_score = args.min_score   


def get_genomedic(filein):
    """load genome sequence"""
    genomedic = {}
    for seq in SeqIO.parse(filein, "fasta"):
        if seq.id not in genomedic:
            genomedic[str(seq.id)] = str(seq.seq)
    return genomedic


def calculate_new_AS(read, genome_dict, sample_tag):
    """
    calculate AS
    """
    seq = read.query_sequence  
    is_reverse = read.is_reverse  
    
    try:
        AS = read.get_tag('AS')
    except KeyError:
        AS = 0

    if (not is_reverse and 'A' in seq) or (is_reverse and 'T' in seq):
        ref_seq = genome_dict[read.reference_name]
        
        for read_pos, ref_pos in read.get_aligned_pairs(matches_only=True):
            if read_pos is None or ref_pos is None:
                continue  
                
            if not is_reverse: 
                if seq[read_pos] == 'A' and ref_seq[ref_pos] == 'G':
                    AS -= 2
            else:  
                if seq[read_pos] == 'T' and ref_seq[ref_pos] == 'C':
                    AS -= 2
                    
        if AS < 0:
            AS = 0
    
    read.set_tag('AS', AS)
    read.set_tag('To', sample_tag)
    
    return read


def get_read_id_tuple(read):
    if read is None:
        return None
    try:
        return (int(read.query_name.split('.')[-1]),)
    except (ValueError, IndexError):
        import re
        return tuple(map(int, re.findall(r'\d+', read.query_name)))


def merge_samfile_by_readid_order(AGsam, TCsam, mergesam, genome_dict, diff_thresh, min_score_thresh):
    """
    merge A-to-G and T-to-C sam file
    """
    with pysam.AlignmentFile(AGsam, 'r') as f1, pysam.AlignmentFile(TCsam, 'r') as f2:
        with pysam.AlignmentFile(mergesam, 'wh', header=f1.header) as out:
            read1 = next(f1, None)
            read2 = next(f2, None)

            id1_tuple = get_read_id_tuple(read1)
            id2_tuple = get_read_id_tuple(read2)

            while read1 or read2:
                current_group = []
                
                if id1_tuple and id2_tuple:
                    target_id = id1_tuple if id1_tuple <= id2_tuple else id2_tuple
                elif id1_tuple:
                    target_id = id1_tuple
                else:
                    target_id = id2_tuple


                while read1 and id1_tuple == target_id:
                    processed_read1 = calculate_new_AS(read1, genome_dict, "AG")
                    current_group.append(processed_read1)
                    read1 = next(f1, None)
                    id1_tuple = get_read_id_tuple(read1)


                while read2 and id2_tuple == target_id:
                    processed_read2 = calculate_new_AS(read2, genome_dict, "TC")
                    current_group.append(processed_read2)
                    read2 = next(f2, None)
                    id2_tuple = get_read_id_tuple(read2)


                if current_group:
                    max_as = max(r.get_tag('AS') for r in current_group)                   
                    if max_as < min_score_thresh:
                        continue 
                    threshold = max_as - diff_thresh
                
                    for r in current_group:
                        if r.get_tag('AS') >= threshold:
                            out.write(r)

    return mergesam


if __name__ == "__main__":
    merge_sam = os.path.join(outdir, f"{out_pre}_merge_AG_TC.sam")

    #if os.path.exists(merge_sam):
    #    sys.stderr.write("[%s] Warning: merge file already exists, please make sure the input file is correct\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    #else:
    sys.stderr.write("[%s] loading genome...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    genomic_dic = get_genomedic(genomefile)
        
    sys.stderr.write("[%s] start merging AG and TC samfile with diff=%d and min_score=%d...\n" % 
                     (strftime("%Y-%m-%d %H:%M:%S", time.localtime()), diff_score, min_as_score))
        

    merge_samfile_by_readid_order(AG_samfile, TC_samfile, merge_sam, genomic_dic, diff_score, min_as_score)
    sys.stderr.write("[%s] merge finished.\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))