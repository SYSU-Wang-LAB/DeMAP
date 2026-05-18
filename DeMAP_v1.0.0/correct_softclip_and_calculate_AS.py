#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Zehong Wu Wang_lab, Sun Yat-sen University"""
"""Feb, 2025"""
"""Email: wuzh93@mail2.sysu.edu.cn"""
"""Usage: This program is used to correct soft-clipping."""

import argparse
import pysam
import numpy as np
from collections import defaultdict
import os
import multiprocessing as mp
import time  
import copy
import math
import subprocess
import re

parser = argparse.ArgumentParser(description = "Correct the softclip due to high C to U mutataion",formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-insam", "--insam", nargs="?", type=str, default="none", help = "samfile of read mapping to AG/TC converted genome")
parser.add_argument("-outname", "--outname", nargs="?", type=str, default="none", help = "output file name")
parser.add_argument("-genome", "--genome", nargs="?", type=str, default="none", help = "unconverted genome")
parser.add_argument("-gtf", "--gtf", nargs="?", type=str, default="none", help = "GTF file from gencode")
parser.add_argument("-p", "--Threads", nargs="?", type=str, default='1', help = "The number of threads")
parser.add_argument("-bz", "--batchsize", nargs="?", type=int, default=1000000, help = "The number of read processed in a batch. (default:1000000)")
parser.add_argument("-C2U", "--C2U", nargs="?", type=float, default=0.05, help = "If C2U mutation level higher than <C2U> in mapped region, treat two 2 C-to-U mismatch as a mismatch")
parser.add_argument("-SC_pval", "--SC_pval", nargs="?", type=float, default=0.001, help = "P-value for the cutoff of LLR test. (default:0.001)")
parser.add_argument("-match_ratio", "--match_ratio", nargs="?", type=float, default=0.8, help = "The minimum ratio of bases match to reference in softclip region. (default:0.8)")
parser.add_argument("--bam", action="store_true", help="If specified, output file as bam.")

args = parser.parse_args()

insam = args.insam
outname = args.outname
genomefile = args.genome
gtffile = args.gtf
Threads = int(args.Threads)
pvalue = args.SC_pval
match_ratio = args.match_ratio
chunk_size = args.batchsize
C2Ulevel = args.C2U
outbam = args.bam

BASE_TO_IDX = {'A': 0, 'T': 1, 'C': 2, 'G': 3}
IDX_TO_BASE = {0: 'A', 1: 'T', 2: 'C', 3: 'G'}
EPSILON = 1e-5  

def load_exons_and_splice_sites(gtf_file):
    """
    extract exon annotations and splice junctions of multi-exon transcripts from the GTF file
    return: (exons_by_transcript, splice_sites)
    """
    exons_by_transcript = {}
    splice_sites = {}
    

    tid_pattern = re.compile(r'transcript_id[\s\t]+(?:\"([^\"]+)\"|([^;\s\t]+))')
    
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t', 8)
            
            if len(fields) < 9 or fields[2] != 'exon':
                continue
                
            chrom = fields[0]
            start = int(fields[3]) - 1  # 1-based to 0-based
            end = int(fields[4]) - 1    # 1-based to 0-based
            strand = fields[6]
            info = fields[8]
            
            match = tid_pattern.search(info)
            if not match:
                continue

            transcript_id = match.group(1) if match.group(1) else match.group(2)
            
            if transcript_id not in exons_by_transcript:
                exons_by_transcript[transcript_id] = []
            exons_by_transcript[transcript_id].append((chrom, start, end, strand))
    
    # only keep multi-exon transcripts
    multi_exon_transcripts = {}
    for tid, exons in exons_by_transcript.items():
        if len(exons) > 1:
            exons_sorted = sorted(exons, key=lambda x: x[1])
            multi_exon_transcripts[tid] = exons_sorted
            
            chrom = exons_sorted[0][0]
            if chrom not in splice_sites:
                splice_sites[chrom] = {}
            strand = exons_sorted[0][3]
            
            for i in range(len(exons_sorted) - 1):
                # donor site: exon end + 1
                donor_pos = exons_sorted[i][2] + 1
                if donor_pos not in splice_sites[chrom]:
                    splice_sites[chrom][donor_pos] = []
                splice_sites[chrom][donor_pos].append((strand, tid, 'donor'))
                
                # acceptor site: next exon start
                acceptor_pos = exons_sorted[i+1][1]  # 0-based
                if acceptor_pos not in splice_sites[chrom]:
                    splice_sites[chrom][acceptor_pos] = []
                splice_sites[chrom][acceptor_pos].append((strand, tid, 'acceptor'))
                
    return multi_exon_transcripts, splice_sites



def compute_softclip_base_freq(sam_path, n_reads=5000000):
    """
    Calculate strand-specific nucleotide frequencies in soft-clipped regions
    """
    sc_pos_count = defaultdict(int)
    sc_neg_count = defaultdict(int)
    sc_pos_total = 0
    sc_neg_total = 0
    
    with pysam.AlignmentFile(sam_path, "r") as samfile:
        for i, read in enumerate(samfile):
            if i >= n_reads:
                break
                
            cigar = read.cigartuples
            read_seq = read.query_sequence
            read_pos = 0
            
            for op, length in cigar:
                if read_pos >= len(read_seq):
                    break
                    
                if op == 4:  # Softclip
                    end_pos = min(read_pos + length, len(read_seq))
                    sc_frag = read_seq[read_pos:end_pos]
                    
                    for base in sc_frag:
                        if base in 'ATCG':
                            if read.is_reverse:
                                sc_neg_count[base] += 1
                                sc_neg_total += 1
                            else:
                                sc_pos_count[base] += 1
                                sc_pos_total += 1
                    
                    read_pos = end_pos
                
                elif op in [0, 1, 7, 8]:  # M/I/=/X
                    read_pos = min(read_pos + length, len(read_seq))
    
    sc_pos_theta = {b: sc_pos_count[b] / sc_pos_total if sc_pos_count[b] > 0 else EPSILON for b in 'ATCG'}
    sc_neg_theta = {b: sc_neg_count[b] / sc_neg_total if sc_neg_count[b] > 0 else EPSILON for b in 'ATCG'}
    
    sc_pos_theta = {k: v / sum(sc_pos_theta.values()) for k, v in sc_pos_theta.items()}
    sc_neg_theta = {k: v / sum(sc_neg_theta.values()) for k, v in sc_neg_theta.items()}
    # 预计算对数版本
    sc_pos_log_theta = {b: np.log(sc_pos_theta[b]) for b in 'ATCG'}
    sc_neg_log_theta = {b: np.log(sc_neg_theta[b]) for b in 'ATCG'}
    
    return sc_pos_theta, sc_neg_theta, sc_pos_log_theta, sc_neg_log_theta


def build_strand_specific_transfer_matrices(sam_path, ref_fasta_path, n_reads=5000000):
    """
    Construct strand-specific nucleotide transition matrices for aligned matched regions
    """
    ref_fasta = pysam.FastaFile(ref_fasta_path)
    ref_bases = 'ATCG'
    
    pos_counts = np.zeros((4, 4), dtype=np.int64)
    neg_counts = np.zeros((4, 4), dtype=np.int64)
    
    with pysam.AlignmentFile(sam_path, "r") as samfile:
        for i, read in enumerate(samfile):
            if i >= n_reads:
                break
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
                
            chrom = read.reference_name
            aligned_pairs = read.get_aligned_pairs(matches_only=True)
            AS = read.get_tag('AS')

            if aligned_pairs is None or AS < 30:
                continue

            for read_pos, ref_pos in aligned_pairs:
                if ref_pos is None:
                    continue
                    
                try:
                    ref_base = ref_fasta.fetch(chrom, ref_pos, ref_pos + 1).upper()
                    read_base = read.query_sequence[read_pos]
                    
                    if ref_base in ref_bases and read_base in ref_bases:
                        ref_idx = BASE_TO_IDX[ref_base]
                        read_idx = BASE_TO_IDX[read_base]
                        
                        if read.is_reverse:
                            neg_counts[ref_idx, read_idx] += 1
                        else:
                            pos_counts[ref_idx, read_idx] += 1
                except:
                    continue
    
    ref_fasta.close()
    
    def build_prob_matrix(counts):
        row_sums = counts.sum(axis=1, keepdims=True)
        safe_sums = np.where(row_sums == 0, 1, row_sums)
        prob_matrix = counts / safe_sums

        log_prob_matrix = np.log(np.maximum(prob_matrix, EPSILON))
        return prob_matrix, log_prob_matrix
    
    pos_prob_matrix, pos_log_prob_matrix = build_prob_matrix(pos_counts)
    neg_prob_matrix, neg_log_prob_matrix = build_prob_matrix(neg_counts)
    
    return (pos_prob_matrix, neg_prob_matrix, pos_log_prob_matrix, neg_log_prob_matrix)

def generate_random_sequence(length, bg_freq):
    bases = list(bg_freq.keys())
    probs = list(bg_freq.values())
    return ''.join(np.random.choice(bases, size=length, p=probs))


def calculate_log_likelihood_ratio(seq, ref, bg_log_freq, log_trans_mat, 
                                   sc_type=None, return_best_subinterval=False):
    """
    calculate LLR
    sc_type: '5prime' or '3prime'
    return:best_subinterval
    """
    if not return_best_subinterval:
        log_p_h0 = 0.0
        log_p_h1 = 0.0
        
        for i, base in enumerate(seq):
            ref_base = ref[i]
            if ref_base not in 'ATCG' or base not in 'ATCG':
                continue        
            # H0 likelihood:log(P(bg))
            log_p_h0 += bg_log_freq[base]
            
            # H0 likelihood: log(P(seq_base|ref_base))
            ref_idx = BASE_TO_IDX[ref_base]
            seq_idx = BASE_TO_IDX[base]
            log_p_h1 += log_trans_mat[ref_idx, seq_idx]

        log_lr = log_p_h1 - log_p_h0 if log_p_h0 != 0 else 0
        return log_lr
    
    else:
        seq_len = len(seq)
        
        if seq_len < 5:
            return None
        
        best_log_lr = -float('inf')
        best_start = 0
        best_end = 0
        best_length = 0
        
        if sc_type == '5prime':

            for sub_len in range(5, seq_len + 1):
                start = seq_len - sub_len  
                end = seq_len
    
                log_p_h0 = 0.0
                log_p_h1 = 0.0
                for i in range(start, end):
                    base = seq[i]
                    ref_base = ref[i]
                    if ref_base not in 'ATCG' or base not in 'ATCG':
                        continue
                    
                    # H0 likelihood
                    log_p_h0 += bg_log_freq[base]
                    
                    # H1 likelihood
                    ref_idx = BASE_TO_IDX[ref_base]
                    seq_idx = BASE_TO_IDX[base]
                    log_p_h1 += log_trans_mat[ref_idx, seq_idx]
                
                log_lr = log_p_h1 - log_p_h0 if log_p_h0 != 0 else 0
                
                if log_lr > best_log_lr:
                    best_log_lr = log_lr
                    best_start = start
                    best_end = end
                    best_length = sub_len
        
        elif sc_type == '3prime':
            for sub_len in range(5, seq_len + 1):
                start = 0
                end = sub_len
                
                log_p_h0 = 0.0
                log_p_h1 = 0.0
                for i in range(start, end):
                    base = seq[i]
                    ref_base = ref[i]
                    if ref_base not in 'ATCG' or base not in 'ATCG':
                        continue
                    
                    # H0 likelihood
                    log_p_h0 += bg_log_freq[base]
                    
                    # H1 likelihood
                    ref_idx = BASE_TO_IDX[ref_base]
                    seq_idx = BASE_TO_IDX[base]
                    log_p_h1 += log_trans_mat[ref_idx, seq_idx]
                
                log_lr = log_p_h1 - log_p_h0 if log_p_h0 != 0 else 0
                
                if log_lr > best_log_lr:
                    best_log_lr = log_lr
                    best_start = start
                    best_end = end
                    best_length = sub_len
        
        if best_log_lr == -float('inf'):
            return None
        
        return {
            'log_lr': best_log_lr,
            'start': best_start,
            'end': best_end,
            'length': best_length,
            'sc_type': sc_type
        }


def simulate_h0_distribution(readlens, n_sim, bg_freq, bg_log_freq, log_trans_mat):
    """LLR distribution under H0"""
    results = defaultdict(list)
    for readlen in readlens:
        for _ in range(n_sim):
            seq = generate_random_sequence(readlen, bg_freq)
            ref = generate_random_sequence(readlen, bg_freq)
            log_lr = calculate_log_likelihood_ratio(seq, ref, bg_log_freq, log_trans_mat)
            results[readlen].append(log_lr)

    return results

def calculate_thresholds(sim_results, alpha=0.005):
    """LLR threthod (alpha分位数)"""
    thresholds = {}
    for length, lr_values in sim_results.items():
        threshold = np.percentile(lr_values, 100 * (1 - alpha))
        thresholds[length] = threshold
    
    return thresholds


def calculate_mismatch(seq1, seq2, is_reverse,tolerate_c2u=False):
    """
    calculate mismatch number
    '+' strand allow A→G mismatch; '-' allow T→C mismatch
    """
    if len(seq1) != len(seq2) or len(seq1) == 0:
        return 0
    
    mismatches = 0
    C2U_count = 0
    for a, b in zip(seq1, seq2):
        if a != b:
            if not is_reverse and a == 'A' and b == 'G':
                continue  
            elif is_reverse and a == 'T' and b == 'C':
                continue 
            else:
                mismatches += 1
                if not is_reverse and a == 'C' and b == 'T':
                    C2U_count += 1
                if is_reverse and a == 'G' and b == 'A':
                    C2U_count += 1
    if tolerate_c2u:
        return mismatches - (C2U_count/2)
    else:
        return mismatches

def create_shifted_read(read, is_5_end):
    """
    create a read with 1 base shift
    """
    shifted_read = copy.copy(read)
    cigar = list(read.cigartuples)
    
    if is_5_end and cigar[0][0] == 4:  # 5' softclip
        cigar[0] = (4, cigar[0][1] + 1)  # increase softclip length
        for i, op in enumerate(cigar):
            if op[0] == 0:  
                cigar[i] = (0, op[1] - 1)
                break
        shifted_read.reference_start += 1
    
    elif not is_5_end and cigar[-1][0] == 4:  # 3' softclip
        cigar[-1] = (4, cigar[-1][1] + 1) 
        for i in range(len(cigar)-1, -1, -1):
            if cigar[i][0] == 0:  
                cigar[i] = (0, cigar[i][1] - 1)
                break
    
    shifted_read.cigartuples = cigar
    return shifted_read

def process_softclip_from_splice(read, ref_fasta, splice_sites, exons_by_transcript, softclip_len, is_5_end, basetheta, transmatrix,tolerate_c2u):

    chrom = read.reference_name
    ref_start = read.reference_start
    ref_end = read.reference_end
    is_reverse = read.is_reverse
    seq = read.query_sequence
    
    junction_pos = ref_start if is_5_end else ref_end
    site_type = 'acceptor' if is_5_end else 'donor'
    

    if chrom not in splice_sites:
        return 0.0, None, 0, 0, 0, 0
    if junction_pos not in splice_sites[chrom]:
        return 0.0, None, 0, 0, 0, 0
    
    best_match = 0.0
    best_tid = None
    best_intron_len = 0
    best_exon_start = 0
    best_delta_AS = 0
    best_LLR = 0
    
    for site_info in splice_sites[chrom][junction_pos]:
        site_strand, tid, stype = site_info
        if stype != site_type:
            continue
        
        if (site_strand == '+' and not is_reverse) or (site_strand == '-' and is_reverse):
            exons = exons_by_transcript[tid]
            
            if is_5_end:
                prev_exon = None
                for i in range(len(exons)):
                    if exons[i][1] == junction_pos and i > 0:  # 找到acceptor位置
                        prev_exon = exons[i-1]
                        break
                if prev_exon is None:
                    continue
                
                exon_length = prev_exon[2] - prev_exon[1] + 1
                if softclip_len > exon_length:
                    continue
                
                exon_seq = ref_fasta.fetch(chrom, prev_exon[2] - softclip_len + 1, prev_exon[2] + 1)
                softclip_seq = seq[:softclip_len]
                tmp_LLR = calculate_log_likelihood_ratio(softclip_seq,exon_seq,basetheta,transmatrix)             

            else:
                next_exon = None
                for i in range(len(exons)):
                    if exons[i][2] + 1 == junction_pos and i < len(exons) - 1:  # 找到donor位置
                        next_exon = exons[i+1]
                        break
                if next_exon is None:
                    continue
                
                exon_length = next_exon[2] - next_exon[1] + 1
                if softclip_len > exon_length:
                    continue
                
                exon_seq = ref_fasta.fetch(chrom, next_exon[1], next_exon[1] + softclip_len)
                softclip_seq = seq[-softclip_len:]
                tmp_LLR = calculate_log_likelihood_ratio(softclip_seq,exon_seq,basetheta,transmatrix)    
            

            mismatch_count = calculate_mismatch(exon_seq, softclip_seq, is_reverse,tolerate_c2u)
            match_rate = (len(exon_seq) - mismatch_count) / len(exon_seq)

            if match_rate > best_match:
                best_match = match_rate
                best_tid = tid
                best_delta_AS = len(exon_seq) - int(mismatch_count + mismatch_count)
                best_LLR = tmp_LLR
                

                if is_5_end:
                    for i in range(len(exons)):
                        if exons[i][1] == junction_pos and i > 0:
                            best_intron_len = exons[i][1] - exons[i-1][2] - 1
                            best_exon_start = prev_exon[2] - softclip_len + 1
                            break
                else:
                    for i in range(len(exons)):
                        if exons[i][2] + 1 == junction_pos and i < len(exons) - 1:
                            best_intron_len = exons[i+1][1] - exons[i][2] - 1
                            best_exon_start = next_exon[1]
                            break
        else:
            continue

    return best_match, best_tid, best_intron_len, best_exon_start ,best_delta_AS,best_LLR


def process_softclip_from_genome(read, ref_fasta,tolerate_c2u):
    regions = []
    cigar_tuples = read.cigartuples
    read_seq = read.query_sequence
    chrom = read.reference_name
    is_reverse = read.is_reverse

    # 1. softclip (5'S)
    if cigar_tuples[0][0] == 4 and cigar_tuples[0][1] >= 5:  # 4 = BAM_CSOFT_CLIP
        length = cigar_tuples[0][1]
        ref_start = read.reference_start
        soft_ref_start = max(0, ref_start - length)
        soft_ref_end = ref_start
        mismatch_count = 0
        try:
            ref_seq = ref_fasta.fetch(chrom, soft_ref_start, soft_ref_end)
            read_seq_clip = read_seq[:length]
            
            if len(ref_seq) == length and len(read_seq_clip) == length:
                mismatch_count = calculate_mismatch(ref_seq, read_seq_clip, is_reverse,tolerate_c2u)
                matchratio = round(1 - (mismatch_count - 1)/(length - 1), 3) if mismatch_count - 1 > 0 else 1
                if matchratio >= 0.5: 
                    regions.append((read_seq_clip,ref_seq,mismatch_count,matchratio,'5prime'))
            
        except:
            pass
    
    # 2. softclip (3'S)
    if cigar_tuples[-1][0] == 4 and cigar_tuples[-1][1] >= 5:  # 4 = BAM_CSOFT_CLIP
        length = cigar_tuples[-1][1]
        ref_end = read.reference_end
        soft_ref_start = ref_end
        soft_ref_end = ref_end + length
        mismatch_count = 0
        try:
            ref_seq = ref_fasta.fetch(chrom, soft_ref_start, soft_ref_end)
            read_seq_clip = read_seq[-length:]
            
            if len(ref_seq) == length and len(read_seq_clip) == length:
                mismatch_count = calculate_mismatch(ref_seq, read_seq_clip, is_reverse,tolerate_c2u)
                matchratio = round(1 - (mismatch_count - 1)/(length - 1), 3) if mismatch_count - 1 > 0 else 1
                if matchratio >= 0.5: 
                    regions.append((read_seq_clip,ref_seq,mismatch_count,matchratio,'3prime'))
        except:
            pass
    
    return regions



def merge_cigar_tuples(cigar_tuples):
    """merge same cigar"""
    if not cigar_tuples:
        return []
    
    merged = []
    current_op, current_len = cigar_tuples[0]
    
    for op, length in cigar_tuples[1:]:
        if op == current_op:
            current_len += length
        else:
            merged.append((current_op, current_len))
            current_op = op
            current_len = length
    
    merged.append((current_op, current_len))
    return merged

def correct_softclip_chunk(reads, ref_fasta_path, 
                          pos_log_theta, neg_log_theta,
                          pos_log_transfer_matrix, neg_log_transfer_matrix,
                          thresholds_plus, thresholds_minus,exons_by_transcript, splice_sites,tolerate_c2u, alpha=0.005):
    ref_fasta = pysam.FastaFile(ref_fasta_path)
    total_processed_reads = []
    process_in_splice = []
    
    for read in reads:
        chrome = read.reference_name
        seq = read.query_sequence   
        cigar_tuples = read.cigartuples
        AS = read.get_tag('AS')
        is_reverse = read.is_reverse
        log_theta = neg_log_theta if is_reverse else pos_log_theta
        log_transfer_matrix = neg_log_transfer_matrix if is_reverse else pos_log_transfer_matrix
        thresholds = thresholds_minus if is_reverse else thresholds_plus


        ########## correct alignment score error during mapping strategy of converted genome ###########
        ########## correct alignment score error during mapping strategy of converted genome ###########

        if ('A' in seq and not is_reverse) or ('T' in seq and is_reverse):
            for read_pos, ref_pos in read.get_aligned_pairs(matches_only=True):
                if read_pos is None or ref_pos is None:
                    continue  
                read_base = seq[read_pos]
                ref_base = ref_fasta.fetch(chrome, ref_pos, ref_pos+1).upper()
                if not is_reverse: #"+" 
                    if tolerate_c2u:
                        if ref_base == 'C' and read_base == 'T': 
                            AS += 1
                else:   #"-" 
                    if tolerate_c2u:
                        if ref_base == 'G' and read_base == 'A': 
                            AS += 1
            if AS < 0:
                AS = 0

        read.set_tag('AS', AS, value_type='i')

        ################## correct fake softclip in splice junctions by match previous or next exon  ####################################### 
        ################## correct fake softclip in splice junctions by match previous or next exon   #######################################      
        # extract softclip
        softclip_5 = 0
        softclip_3 = 0
        if cigar_tuples[0][0] == 4:  # 5' softclip
            softclip_5 = cigar_tuples[0][1]
        if cigar_tuples[-1][0] == 4:  # 3' softclip
            softclip_3 = cigar_tuples[-1][1]
        
        # require softclip length >= 5
        if softclip_5 < 5 and softclip_3 < 5:
            read.set_tag('SC', "NO_short")
            total_processed_reads.append(read)
            continue
        
        # process 5' softclip
        best_5_match, best_5_tid, best_5_intron_len, best_5_exon_start, best_5_delta_AS, best_5_LLR = 0.0, None, 0, 0, 0, 0
        need_shift_5 = False
        
        if softclip_5 >= 5:
            best_5_match, best_5_tid, best_5_intron_len, best_5_exon_start, best_5_delta_AS, best_5_LLR = process_softclip_from_splice(
                read, ref_fasta, splice_sites, exons_by_transcript, softclip_5, True, log_theta, log_transfer_matrix,tolerate_c2u)
            
            if best_5_tid is None: #（junction_pos not in splice_sites）
                need_shift_5 = True
                shifted_read = create_shifted_read(read, True)
                shifted_softclip_5 = softclip_5 + 1
                    
                match_shift, tid_shift, intron_shift, exon_start_shift, delta_AS_shift, LLR_shift = process_softclip_from_splice(
                    shifted_read, ref_fasta, splice_sites, exons_by_transcript, shifted_softclip_5, True, log_theta, log_transfer_matrix,tolerate_c2u)
                    
                if match_shift > 0.5:
                    best_5_match, best_5_tid, best_5_intron_len, best_5_exon_start, best_5_delta_AS, best_5_LLR = \
                        match_shift, tid_shift, intron_shift, exon_start_shift, delta_AS_shift, LLR_shift
        
        #  3' softclip
        best_3_match, best_3_tid, best_3_intron_len, best_3_exon_start, best_3_delta_AS, best_3_LLR = 0.0, None, 0, 0, 0, 0
        need_shift_3 = False
        
        if softclip_3 >= 5:
            best_3_match, best_3_tid, best_3_intron_len, best_3_exon_start, best_3_delta_AS, best_3_LLR = process_softclip_from_splice(
                read, ref_fasta, splice_sites, exons_by_transcript, softclip_3, False, log_theta, log_transfer_matrix,tolerate_c2u)
            
            # try shift
            if best_3_tid is None:
                need_shift_3 = True
                shifted_read = create_shifted_read(read, False)
                shifted_softclip_3 = softclip_3 + 1
                    
                match_shift, tid_shift, intron_shift, exon_start_shift, delta_AS_shift, LLR_shift = process_softclip_from_splice(
                    shifted_read, ref_fasta, splice_sites, exons_by_transcript, shifted_softclip_3, False, log_theta, log_transfer_matrix,tolerate_c2u)
                    
                if match_shift > 0.5:
                    best_3_match, best_3_tid, best_3_intron_len, best_3_exon_start, best_3_delta_AS, best_3_LLR = \
                        match_shift, tid_shift, intron_shift, exon_start_shift, delta_AS_shift, LLR_shift
        
        # only correct read with soft-clipped bases >= 5bp
        if softclip_5 >= 5:
            threshold_len = min(softclip_5, 50)
            modify_5 = (best_5_match >= match_ratio and best_5_LLR >= thresholds[threshold_len]) or best_5_match == 1
        else:
            modify_5 = False

        if softclip_3 >= 5:
            threshold_len = min(softclip_3, 50)
            modify_3 = (best_3_match >= match_ratio and best_3_LLR >= thresholds[threshold_len]) or best_3_match == 1
        else:
            modify_3 = False

        max_str = str(round(max(best_5_match,best_3_match),3))

        if not modify_5 and not modify_3:
            read.set_tag('SC', "NO_splice_"+max_str)
            process_in_splice.append(read)
            continue
        
        if modify_5 and modify_3:
            if best_5_tid != best_3_tid:
                if best_5_match > best_3_match:
                    modify_3 = False
                elif best_3_match > best_5_match:
                    modify_5 = False
                else:
                    if softclip_5 > softclip_3:
                        modify_3 = False
                    elif softclip_3 > softclip_5:
                        modify_5 = False
                    else:
                        modify_5 = False
                        modify_3 = False

        if not modify_5 and not modify_3:
            read.set_tag('SC', "diff_transcript")
            total_processed_reads.append(read)
            continue

        # correct CIGAR
        new_cigar = list(cigar_tuples)
        new_pos = read.reference_start
        
        if modify_5:
            if new_cigar[0][0] == 4:
                new_cigar.pop(0)
            
            if need_shift_5:
                # S+1, M-1
                new_cigar.insert(0, (0, softclip_5 + 1))  # M 
                for i, op in enumerate(new_cigar):
                    if op[0] == 0 and i > 0:  
                        new_cigar[i] = (0, op[1] - 1)
                        break
            else:
                new_cigar.insert(0, (0, softclip_5))  # M
            
            new_cigar.insert(1, (3, best_5_intron_len))  # N
            
            new_pos = best_5_exon_start
            current_AS = read.get_tag('AS')
            read.set_tag('AS', int(current_AS + best_5_delta_AS), value_type='i')
        
        if modify_3:
            if new_cigar[-1][0] == 4:
                new_cigar.pop()
            
            if need_shift_3:
                new_cigar.append((3, best_3_intron_len))  # N
                new_cigar.append((0, softclip_3 + 1))  # M (增加1)
                
                for i in range(len(new_cigar)-2, -1, -1):
                    if new_cigar[i][0] == 0:  # M操作
                        new_cigar[i] = (0, new_cigar[i][1] - 1)
                        break
            else:
                new_cigar.append((3, best_3_intron_len))  # N
                new_cigar.append((0, softclip_3))  # M

            current_AS = read.get_tag('AS')
            read.set_tag('AS', int(current_AS + best_3_delta_AS), value_type='i')

        # update read information
        read.cigartuples = new_cigar
        read.reference_start = new_pos

        if modify_5 and modify_3:
            read.set_tag('SC', "YES_53_" + max_str)
            total_processed_reads.append(read)
            continue

        elif modify_5:
            if need_shift_5:
                read.set_tag('SC', "YES_5_shift_" + max_str)
            else:
                read.set_tag('SC', "YES_5_" + max_str)

        elif modify_3:
            if need_shift_3:
                read.set_tag('SC', "YES_3_shift" + max_str)     
            else:
                read.set_tag('SC', "YES_3_" + max_str)                    
        process_in_splice.append(read)


    ################## correct fake softclip by match previous or next genome region #######################################
    ################## correct fake softclip by match previous or next genome region #######################################
    for read in process_in_splice:
        chrome = read.reference_name
        seq = read.query_sequence   
        cigar_tuples = read.cigartuples
        is_reverse = read.is_reverse
        log_theta = neg_log_theta if is_reverse else pos_log_theta
        log_transfer_matrix = neg_log_transfer_matrix if is_reverse else pos_log_transfer_matrix
        thresholds = thresholds_minus if is_reverse else thresholds_plus

        soft_regions = process_softclip_from_genome(read, ref_fasta,tolerate_c2u)
        
        if not soft_regions:
            read.set_tag('SC', "NO")
            total_processed_reads.append(read)
            continue
        

        AS_changes = 0
        corrected = False
        new_cigar = list(cigar_tuples)  
        result_record = []

        for (soft_seq, ref_seq,mismatch_count,matchratio, sc_type) in soft_regions:
            soft_len = len(soft_seq)

            if soft_len < 5: 
                continue

            elif matchratio >= 0.95: #softclip region highly match to genome region
                corrected = True
                AS_delta = (soft_len - mismatch_count) - mismatch_count
                AS_changes += AS_delta
                    
                # change cigar from S to M
                if sc_type == '5prime':
                    new_cigar[0] = (0, new_cigar[0][1])  # (BAM_CMATCH, length)
                    read.reference_start = read.reference_start - soft_len
                else:
                    new_cigar[-1] = (0, new_cigar[-1][1])  # (BAM_CMATCH, length)                
                result_record.append([True, matchratio, sc_type, "none"])

            else:
                best_interval = calculate_log_likelihood_ratio(
                    soft_seq, ref_seq, log_theta, log_transfer_matrix, 
                    sc_type=sc_type, return_best_subinterval=True
                )
                if best_interval is not None:
                    sub_start = best_interval['start']
                    sub_end = best_interval['end']
                    sub_seq = soft_seq[sub_start:sub_end]
                    sub_ref = ref_seq[sub_start:sub_end]
                    
                    sub_mismatches = calculate_mismatch(sub_ref, sub_seq, is_reverse, tolerate_c2u)
                    sub_matchratio = round(1 - (sub_mismatches - 1) / (best_interval['length'] - 1),3) if sub_mismatches - 1 > 0 else 1
                    
                    threshold_len = min(best_interval['length'], 50)
                    threshold = thresholds[threshold_len]
                    if (best_interval['log_lr'] >= threshold and sub_matchratio >= match_ratio) or sub_matchratio == 1:
                        corrected = True
                        AS_delta = (best_interval['length'] - sub_mismatches) - sub_mismatches
                        AS_changes += AS_delta
                        if sc_type == '5prime':
                            left_soft_len = sub_start  
                            match_len = best_interval['length']  
                            
                            del new_cigar[0] 
                            if left_soft_len > 0:
                                new_cigar.insert(0, (4, left_soft_len))
                                new_cigar.insert(1, (0, match_len))
                            else:
                                new_cigar.insert(0, (0, match_len))
                            
                            read.reference_start = read.reference_start - match_len
                            
                        else:  # '3prime'
                            match_len = best_interval['length'] 
                            right_soft_len = soft_len - sub_end 
                            

                            del new_cigar[-1]
                            new_cigar.append((0, match_len))
                            if right_soft_len > 0:

                                new_cigar.append((4, right_soft_len))
                        if match_len < soft_len:
                            result_record.append([True, sub_matchratio, sc_type, f"sub_{best_interval['length']}bp_{ref_seq}_{soft_seq}"])
                        else:
                            result_record.append([True, matchratio, sc_type, "none"])
                    else:
                        result_record.append([False, matchratio, sc_type,"none"])

        if corrected:
            # rebuild cigar
            merged_cigar = merge_cigar_tuples(new_cigar)
            cigar_str_parts = []
            for op, length in merged_cigar:
                op_char = "MIDNSHP=X"[op]
                cigar_str_parts.append(f"{length}{op_char}")
            read.cigarstring = ''.join(cigar_str_parts)
            
            # re-score
            current_AS = read.get_tag('AS')
            read.set_tag('AS', int(current_AS + AS_changes), value_type='i')

        if result_record:
            splicetag = read.get_tag('SC')
            if any(r[0] for r in result_record): 
                all_match_ratio = [r[1] for r in result_record]
                best_match = max([r[1] for r in result_record if r[0]])
                best_index = all_match_ratio.index(best_match)
                match_tag = result_record[best_index][2] + "_" + result_record[best_index][3]
                if 'sub' in match_tag:
                    read.set_tag('SC', splicetag + "_YES_" + str(best_match) + "_" + match_tag)
                else:
                    read.set_tag('SC', splicetag + "_YES_" + str(best_match))
            else:
                best_match = max([r[1] for r in result_record])
                read.set_tag('SC', splicetag + "_NO_" +  str(best_match))
        total_processed_reads.append(read)

    ref_fasta.close()
    #sorted_read_list = sorted(total_processed_reads, key=lambda read: int(read.query_name.split('.')[1]))
    sorted_read_list = sorted(total_processed_reads, key=lambda read: (int(read.query_name[2:].split('_')[0]),int(read.query_name.split('_')[1]) if '_' in read.query_name else 0))
    return sorted_read_list

def correct_softclip_chunk_wrapper(args):
    chunk_id, read_strs, ref_fasta_path,  pos_log_theta, neg_log_theta, pos_log_transfer_matrix,\
    neg_log_transfer_matrix, thresholds_plus, thresholds_minus, alpha, samheader,exons_by_transcript, splice_sites,tolerate_c2u = args
    
    reads = []
    header = pysam.AlignmentHeader().from_dict(samheader)
    
    for read_str in read_strs:
        try:
            read =  pysam.AlignedSegment().fromstring(read_str, header)
            reads.append(read)
        except:
            dummy_read = pysam.AlignedSegment()
            reads.append(dummy_read)
    
    processed_reads = correct_softclip_chunk(
        reads, ref_fasta_path, 
        pos_log_theta, neg_log_theta,
        pos_log_transfer_matrix, neg_log_transfer_matrix,
        thresholds_plus, thresholds_minus,exons_by_transcript, splice_sites,tolerate_c2u, alpha
    )
    
    processed_strs = []
    for read in processed_reads:
        try:
            processed_strs.append(read.tostring())
        except:
            processed_strs.append("")
    
    return chunk_id, processed_strs

def write_completed_results(outfile, submitted_tasks, results_cache, next_write_id, header):
    """
    output results
    """
    i = 0
    while i < len(submitted_tasks):
        chunk_id, future = submitted_tasks[i]
        
        if future.ready():
            try:
                result_id, processed_strs = future.get()
                results_cache[result_id] = processed_strs
                submitted_tasks.pop(i)
                #print(f"Chunk {result_id} processed successfully")
            except Exception as e:
                print(f"Error processing chunk {chunk_id}: {e}")

                results_cache[chunk_id] = [""] * 100
                submitted_tasks.pop(i)
                continue
        else:
            i += 1
    
    while next_write_id in results_cache:
        processed_strs = results_cache.pop(next_write_id)
        for read_str in processed_strs:
            if read_str: 
                try:
                    read = pysam.AlignedSegment().fromstring(read_str, header)
                    outfile.write(read)
                except:
                    pass
        #print(f"Chunk {next_write_id} written to output")
        next_write_id += 1
    
    return next_write_id

def main(sam_input, reference_fasta, output_sam, chunk_size=chunk_size, alpha=pvalue, n_processes=Threads):

    exons_by_transcript, splice_sites = load_exons_and_splice_sites(gtffile)

    print("Precomputing strand-specific parameters from first 2,000,000 aligned reads...")
    
    pos_theta, neg_theta, pos_log_theta, neg_log_theta = compute_softclip_base_freq(sam_input, n_reads=2000000)
    print(f"\nPositive strand softclip region base frequencies (H0): {pos_theta}")
    print(f"Negative strand softclip region base frequencies (H0): {neg_theta}")

    pos_transfer_matrix, neg_transfer_matrix, pos_log_transfer_matrix, neg_log_transfer_matrix = build_strand_specific_transfer_matrices(
        sam_input, reference_fasta, n_reads=2000000
    )
    
    print("\nPositive strand transfer matrix (H1:ref → read):")
    print("       A      T      C      G")
    for i, ref_base in enumerate('ATCG'):
        row = pos_transfer_matrix[i]
        print(f"{ref_base}  {row[0]:.4f}  {row[1]:.4f}  {row[2]:.4f}  {row[3]:.4f}")
    
    print("\nNegative strand transfer matrix (H1:ref → read):")
    print("       A      T      C      G")
    for i, ref_base in enumerate('ATCG'):
        row = neg_transfer_matrix[i]
        print(f"{ref_base}  {row[0]:.4f}  {row[1]:.4f}  {row[2]:.4f}  {row[3]:.4f}")
    

    tolerate_c2u = pos_transfer_matrix[BASE_TO_IDX['C']][BASE_TO_IDX['T']] > C2Ulevel or neg_transfer_matrix[BASE_TO_IDX['G']][BASE_TO_IDX['A']] > C2Ulevel


    print("\nSimulating H0 distribution for LR thresholds:")    
    len_list = [i for i in range(5, 51)]  # if sotflen > 50 treated as 50
    N_SIM = 100000  # 1×10^-5
    np.random.seed(22)

    # simulate H0 distribution
    sim_results_plus = simulate_h0_distribution(len_list, N_SIM, pos_theta, pos_log_theta, pos_log_transfer_matrix)  # using strand '+' 
    sim_results_minus = simulate_h0_distribution(len_list, N_SIM, neg_theta, neg_log_theta, neg_log_transfer_matrix)  # using strand '-' 

    # Calculate thresholds for different lengths
    thresholds_plus = calculate_thresholds(sim_results_plus, alpha=alpha)
    thresholds_minus = calculate_thresholds(sim_results_minus, alpha=alpha)


    print(f"\nStreaming processing in chunks (size={chunk_size}) using {n_processes} processes...")
    
    chunk_id = 0
    total_reads = 0
    submitted_tasks = []
    results_cache = {}
    next_write_id = 0
    
    pool = mp.Pool(processes=n_processes)
    
    try:
        with pysam.AlignmentFile(sam_input, "rb") as samfile:
            header_dict = samfile.header.to_dict()
            
            with pysam.AlignmentFile(output_sam, "w", header=samfile.header) as outfile:
                buffer = []
                max_pending_tasks = n_processes * 4
                
                for read in samfile:
                    buffer.append(read)
                    total_reads += 1
                    
                    if len(buffer) >= chunk_size:
                        if len(submitted_tasks) >= max_pending_tasks:
                            time.sleep(0.01)
                            next_write_id = write_completed_results(
                                outfile, submitted_tasks, results_cache, next_write_id, samfile.header
                            )
                        
                        read_strs = [r.tostring() for r in buffer]
                        
                        task_args = (
                            chunk_id, read_strs, reference_fasta, 
                            pos_log_theta, neg_log_theta,
                            pos_log_transfer_matrix, neg_log_transfer_matrix,
                            thresholds_plus, thresholds_minus, alpha, header_dict, exons_by_transcript, splice_sites,tolerate_c2u
                        )
                        
                        future = pool.apply_async(correct_softclip_chunk_wrapper, (task_args,))
                        submitted_tasks.append((chunk_id, future))
                        
                        buffer = []
                        chunk_id += 1
                        
                        next_write_id = write_completed_results(
                            outfile, submitted_tasks, results_cache, next_write_id, samfile.header
                        )
                
                # Process the last batch
                if buffer:
                    read_strs = [r.tostring() for r in buffer]
                    task_args = (
                        chunk_id, read_strs, reference_fasta, 
                        pos_log_theta, neg_log_theta,
                        pos_log_transfer_matrix, neg_log_transfer_matrix,
                        thresholds_plus, thresholds_minus, alpha, header_dict,exons_by_transcript, splice_sites,tolerate_c2u
                    )
                    
                    future = pool.apply_async(correct_softclip_chunk_wrapper, (task_args,))
                    submitted_tasks.append((chunk_id, future))
                    chunk_id += 1
                
                # Process all remaining tasks
                while submitted_tasks:
                    next_write_id = write_completed_results(
                        outfile, submitted_tasks, results_cache, next_write_id, samfile.header)
                    if submitted_tasks:
                        time.sleep(0.01)
    
    except Exception as e:
        print(f"Main processing error: {e}")
    finally:
        pool.terminate()
        pool.join()
        print(f"Processed {total_reads} reads in total")
    
    print(f"\nDone! Corrected softclips written to {output_sam}")

    if outbam:
        pysam.sort("-o", output_sam[:-3] + "bam", "-@", "8", output_sam)
        pysam.index(output_sam[:-3] + "bam")
        subprocess.call("rm " + output_sam, shell=True)

if __name__ == "__main__":
    if outname == "none":
        output_sam = insam[:-4] + "_correct.sam"
    else:
        output_sam = outname
    
    start_time = time.time()
    main(insam, genomefile, output_sam)
    print(f"Total execution time: {time.time() - start_time:.2f} seconds")



