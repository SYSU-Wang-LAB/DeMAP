#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Zehong Wu wang_lab, Sun Yat-sen University"""
"""Feb, 2025"""
"""Email: wuzh93@mail2.sysu.edu.cn"""
"""Usage: This program is designed for the allocation of multi-mapped reads."""

import argparse
import pysam
import os
import multiprocessing
import gc
import numpy as np
from scipy.sparse import csr_matrix
import re
from functools import partial
import subprocess
import sys
import time
from time import strftime

parser = argparse.ArgumentParser(description = "Run EM model",formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-multi", "--multiAlign", nargs="?", type=str, default="none", help = "samfile for multi-alignment reads")
parser.add_argument("-uniq", "--uniqAlign", nargs="?", type=str, default="none", help = "samfile of uniq mapping reads or matched GLORI- reads")
parser.add_argument("-pre", "--outname_prefix", nargs="?", type=str, default='default',help = "--outname_prefix")
parser.add_argument("-o", "--outputdir", nargs="?", type=str, default="none", help="outputdir")
parser.add_argument("-max", "--maxmap", nargs="?", type=int, default=5, help = "Discard multi-mapping reads with more than maxmap mappings")
parser.add_argument("-w", "--winsize", nargs="?", type=int, default=150, help="window size (bp) for searching read location. (default:150)")
parser.add_argument("-maxcount", "--maxcount", nargs="?", type=int, default=100000, help="The max number of reads mapping to a single window. (default:100000)")
parser.add_argument("-p", "--Threads", nargs="?", type=int, default=1, help = "The number of alignment threads. (default:1)")
parser.add_argument("-max_iter", "--max_iter", nargs="?", type=int, default=50, help = "The max iteration times in EM model. (default:50)")
parser.add_argument("-tol", "--tol", nargs="?", type=float, default=1e-2, help = "A small number for testing convergence between iterations. (default:0.01)")
parser.add_argument("-SF", "--AS_factor", nargs="?", type=int, default=3, help = "Exponential parameter for amplifying the influence of AS, 0 for ignore alignemnt scores")
parser.add_argument("-cutoff", "--cutoff", nargs="?", type=float, default=0.8, help = "Probability cutoff for reassigning multi-mapping reads. (default:0.8 cutoff could be 0.60-0.95)\n" 
                    "choose region when region score/sum(score in all regions) >= cutoff")
parser.add_argument("-size", "--size", nargs="?", type=int, default=800, help = "The max number of multi-mapping reads used to estimate"
                    "the probability of the connected-region in a multi-mapping read. (default:600: to reduce memory consumption)")
parser.add_argument("-multisize", "--multi_size", nargs="?", type=int, default=20000000, help = "The max multi-mapping reads using in evaluate the region abundance. (default:20000000)\n"
                    "This parameter was used for large data")
parser.add_argument("-ms", "--min_score", nargs="?", type=int, default=28, help = "Only using high quality alignments (AS > min_score) for evaluating the region abundance. (default:28)")
parser.add_argument("-bz", "--batchsize", nargs="?", type=int, default=1000000, help = "The number of multi-mapping reads process in a batch. (default:1000000)")

args = parser.parse_args()

multi_sam = args.multiAlign
uniq_sam = args.uniqAlign
outdir = args.outputdir
out_pre = args.outname_prefix
winsize = args.winsize
maxcount = args.maxcount
threads = args.Threads
max_iter = args.max_iter
tol = args.tol
ScoreFactor = args.AS_factor
cutoff = args.cutoff
sizenum = args.size
batchsize = args.batchsize if args.batchsize >= 200000 else 200000
min_score = args.min_score
multi_size = args.multi_size
maxmap = args.maxmap

def convert_readid_to_float(input_string):
    """only for saving memory"""
    number_part = input_string[2:]
    if "_" in number_part:
        integer_part, decimal_part = number_part.split("_", 1)
        result = float(integer_part) + float("0." + decimal_part)
    else:
        result = float(number_part)
    return result


def create_genome_bins(chrom_lengths, bin_size=100):
    genome_bins = {}
    for chrom, length in chrom_lengths.items():
        for strand in ["+", "-"]:
            chrom_strand = chrom + ";" + strand
            bin_number = (length + bin_size) // bin_size
            for bin in range(0,bin_number):
                genome_bins[chrom_strand + "_" + str(bin)] = [0,0,[]]
    return genome_bins


def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join(complement.get(base, 'N') for base in reversed(seq))


def initial_bin_counts(sam_uniq,sam_multi,genome_bins,bin_size = 100,max_multi_reads=10000000):
    with pysam.AlignmentFile(sam_uniq, "r") as uniq_in:
        for read in uniq_in:
            Alignment_score = read.get_tag('AS')
            if Alignment_score < min_score:
                continue
            chrom = read.reference_name
            strand = "+" if not read.is_reverse else "-"
            chrom_strand = chrom + ";" + strand
            read_at_reference_position = read.get_reference_positions()
            center_index = round(len(read_at_reference_position)/2)
            read_center_pos = read_at_reference_position[center_index]
            bin_idx = read_center_pos // bin_size
            target_bin = chrom_strand + "_" + str(bin_idx)
            genome_bins[target_bin][0] += 1

    total_multi_reads = 0
    present_id = None
    read_to_location = {}
    with pysam.AlignmentFile(sam_multi, "r") as multi_in:
        for read in multi_in:
            Alignment_score = read.get_tag('AS')
            if Alignment_score < min_score:
                continue
            readid = convert_readid_to_float(read.query_name)
            if readid != present_id:
                total_multi_reads += 1
                present_id = readid

            if total_multi_reads <= max_multi_reads:
                chrom = read.reference_name
                strand = "+" if not read.is_reverse else "-"
                chrom_strand = chrom + ";" + strand
                read_at_reference_position = read.get_reference_positions()
                center_index = round(len(read_at_reference_position)/2)
                read_center_pos = read_at_reference_position[center_index]
                bin_idx = read_center_pos // bin_size
                
                target_bin = sys.intern(chrom_strand + "_" + str(bin_idx))

                if readid not in read_to_location:
                     read_to_location[readid] = [target_bin] 
                else:
                    read_to_location[readid].append(target_bin) 

                genome_bins[target_bin][1] += 1 
                genome_bins[target_bin][2].append(readid) 

    if total_multi_reads <= max_multi_reads:
        used_ratio = 1
    else:
        used_ratio = max_multi_reads/total_multi_reads


    for rid in read_to_location:
        read_to_location[rid] = tuple(read_to_location[rid])

    used_region_count = {}
    for target_bin,region_counts in genome_bins.items():
        uniqcounts = region_counts[0]
        multicounts = region_counts[1]
        bin_read_id = region_counts[2]
        
        if len(bin_read_id) >= maxcount:
            bin_read_id = np.random.choice(bin_read_id, maxcount, replace=False)
        else:
            bin_read_id = np.array(bin_read_id, dtype=np.float64)
            
        if uniqcounts != 0 or multicounts != 0: 
            used_region_count[target_bin] = [uniqcounts,bin_read_id] 

    genome_bins.clear()
    return used_region_count,read_to_location,used_ratio


def run_EM_model(connect_regions,maxsizenum=600,maxiter=50,mintol=1e-2):
    expression = {}

    reads_in_connect_region = set()
    for region in connect_regions:
        if region not in genome_bins_for_unsed:
            expression[region] = 1
        else:
            expression[region] = genome_bins_for_unsed[region][0] + 1
            reads_in_connect_region.update(genome_bins_for_unsed[region][1])

    if not reads_in_connect_region: 
        return None

    valid_reads = []
    connect_regions_set = set(connect_regions) 

    for read in reads_in_connect_region:
        read_mapped_regions = read_to_location.get(read, [])
        if all(r in connect_regions_set for r in read_mapped_regions):
            valid_reads.append(read)
            
    all_reads = valid_reads
    multi_readnum = len(all_reads)

    if multi_readnum == 0:
        return None
    
    if multi_readnum > maxsizenum:
        all_reads = np.random.choice(all_reads, maxsizenum, replace=False)
        scale_size = (multi_readnum/used_ratio) / maxsizenum
        for region in expression:
            expression[region] = expression[region] / scale_size
    else:
        for region in expression:
            expression[region] = expression[region] * used_ratio

    connect_regions = list(expression.keys())
    region_idx = {region: i for i, region in enumerate(connect_regions)}

    row_indices = []
    col_indices = []
    data = []
    for r_idx, read in enumerate(all_reads):
        for region in read_to_location.get(read, []):
            if region in connect_regions:
                row_indices.append(r_idx)
                col_indices.append(region_idx[region])
                data.append(1)

    read_region_matrix = csr_matrix((data, (row_indices, col_indices)), shape=(len(all_reads), len(connect_regions)))
    expression_array = np.array([expression[region] for region in connect_regions], dtype=np.float64)
    total_counts = expression_array.sum()
    prob_array = (expression_array / total_counts) 

    # start iteration
    for _ in range(maxiter):
        weighted_matrix = read_region_matrix.multiply(prob_array) 
        total_p = weighted_matrix.sum(axis=1).A1 
        weights_matrix = weighted_matrix.multiply(1 / total_p[:, None])
        new_expression_array = expression_array + weights_matrix.sum(axis=0).A1
        total_counts = new_expression_array.sum()
        new_prob_array = new_expression_array / total_counts

        max_diff = np.max(np.abs(new_prob_array - prob_array))
        if max_diff < mintol:
            break
        prob_array = new_prob_array

    final_expression_ratio = list(new_prob_array)
    return final_expression_ratio


def process_chunk(chunk, sizenum, max_iter, tol):
    result = {}
    for connect_regions, _ in chunk:
        connect_bins = connect_regions.split(";;")
        ratio = run_EM_model(connect_bins, sizenum, max_iter, tol)
        if ratio is not None:
            result[connect_regions] = {connect_bins[i]: ratio[i] for i in range(len(connect_bins))}
    return result


def get_relative_abundance_in_connected_regions(read_to_location, num_workers=4, chunk_size=100):
    relative_region_abundance = {}
    connect_region_counts = {}
    
    for read in read_to_location:
        if len(set(read_to_location[read])) != len(read_to_location[read]):
            continue
        connect_regions = ";;".join(sorted(read_to_location[read]))
        connect_region_counts[connect_regions] = connect_region_counts.get(connect_regions, 0) + 1

    sorted_regions = sorted(connect_region_counts.items(), key=lambda x: x[1], reverse=True)

    worker_func = partial(process_chunk, sizenum=sizenum, max_iter=max_iter, tol=tol)
    chunks = [sorted_regions[i:i+chunk_size] for i in range(0, len(sorted_regions), chunk_size)]
    n_process = num_workers if num_workers <= 4 else 4

    with multiprocessing.Pool(processes=n_process, maxtasksperchild=50) as pool:
        chunk_results = pool.map(worker_func, chunks)
    
    for res in chunk_results:
        relative_region_abundance.update(res)
    
    return relative_region_abundance


def rank_list_unique(lst):
    sorted_lst = sorted(set(lst)) 
    rank_dict = {val: rank + 1 for rank, val in enumerate(sorted_lst)} 
    return [rank_dict[x] for x in lst]


def reassign_reads(reads,maxsizenum,maxiter=50, mintol=1e-2):
    outlist_success = []
    outlist_fail = []
    for single_read in reads:
        connect_regions = single_read[0]
        Alignments = single_read[1]
        if len(set(connect_regions)) != len(Alignments):
            continue

        connect_regions_key = ";;".join(sorted(connect_regions))
        if connect_regions_key in high_repeat_region_abundance:
            final_expression_ratio = [high_repeat_region_abundance[connect_regions_key][region] for region in connect_regions]
        else:
            final_expression_ratio = run_EM_model(connect_regions,maxsizenum,maxiter,mintol)
            if final_expression_ratio is None:
                continue

        Mapping_Score = [int(re.search(r'AS:i:(\d+)', alignment).group(1)) for alignment in Alignments]
        Mapping_Score_ranking = rank_list_unique(Mapping_Score) 
        ranking_Score = [score ** ScoreFactor for score in Mapping_Score_ranking] 

        final_weights = [ranking_Score[i]*final_expression_ratio[i] for i in range(0,len(ranking_Score))]
        final_weights = [final_weight/sum(final_weights) for final_weight in final_weights]
        final_weights_sorted = sorted(final_weights,reverse=True)

        max_weights_index = final_weights.index(final_weights_sorted[0])

        if final_weights[max_weights_index] >= cutoff:
            tmp_str = re.sub(r'(NH|HI):i:\d+', r'\1:i:1', Alignments[max_weights_index])
            outlist_success.append(tmp_str)
        else:
            outlist_fail.append(Alignments[max_weights_index])
    return outlist_success,outlist_fail


def process_batch(onebatch, batch_id, pool):
    info_record = {}
    for read_id, target_bin, sam_text in onebatch:
        if read_id not in info_record:
            info_record[read_id] = [[target_bin],[sam_text]]
        else:
            info_record[read_id][0].append(target_bin)
            info_record[read_id][1].append(sam_text)

    onebatch_info = list(info_record.values())
    del info_record


    chunk_size = 1000 
    chunks = [onebatch_info[i:i+chunk_size] for i in range(0, len(onebatch_info), chunk_size)]

    results_success = []
    results_fail = []
    args = [(chunk, sizenum,max_iter,tol) for chunk in chunks]
    

    for result in pool.starmap(reassign_reads, args):
        results_success.extend(result[0])
        results_fail.extend(result[1])
    
    tmp_file1 = os.path.join(outdir, f"{out_pre}_tmp_batch_{batch_id}_success.sam")
    tmp_file2 = os.path.join(outdir, f"{out_pre}_tmp_batch_{batch_id}_fail.sam")
    with open(tmp_file1, "w") as f1, open(tmp_file2, "w") as f2:
        f1.write("\n".join(results_success) + "\n")
        f2.write("\n".join(results_fail) + "\n")


def read_and_process_file(multi_file, max_readnum=500000):
    batch = []  
    unique_readid = set() 
    batch_id = 1
    

    with multiprocessing.Pool(processes=threads) as pool:
        with pysam.AlignmentFile(multi_file, "rb") as multi_align:
            for alignment in multi_align:
                read_id = alignment.query_name

                if read_id not in unique_readid:
                    if len(unique_readid) >= max_readnum:
                        process_batch(batch, batch_id, pool)  
                        gc.collect() 
                        batch_id += 1
                        batch = []  
                        unique_readid.clear()

                unique_readid.add(read_id)
                

                chrom = alignment.reference_name
                strand = "+" if not alignment.is_reverse else "-"
                chrom_strand = chrom + ";" + strand
                read_at_reference_position = alignment.get_reference_positions()
                center_index = round(len(read_at_reference_position)/2)
                read_center_pos = read_at_reference_position[center_index]
                bin_index = read_center_pos // winsize
                target_bin = chrom_strand + "_" + str(bin_index)
                sam_text = alignment.to_string()
                
                batch.append((read_id, target_bin, sam_text))

            if batch:
                process_batch(batch, batch_id, pool)


if __name__ == "__main__":
    parent_pid = os.getpid()
    with pysam.AlignmentFile(uniq_sam, "r") as sam_in:
        chrom_lengths = {ref["SN"]: ref["LN"] for ref in sam_in.header["SQ"]}

    genome_bins = create_genome_bins(chrom_lengths, winsize)

    mergesam = outdir + "/" + out_pre + "_final.sam"
    mergesam_sorted = outdir + "/" + out_pre + "_final_sorted.bam"

    if maxmap == 1:
        sys.stderr.write("[%s] Only using uniq-mapping reads, skip EM model.\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
        subprocess.call("samtools view -h -bS " + uniq_sam + " | samtools sort -@ " + str(threads) + " -o " + mergesam_sorted, shell=True)
        subprocess.call("samtools index " + mergesam_sorted, shell=True)
    else:
        sys.stderr.write("[%s] Start running EM to assign multi-mapping reads.\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))

        genome_bins_for_unsed,read_to_location,used_ratio = initial_bin_counts(uniq_sam,multi_sam,genome_bins,winsize,multi_size)
        

        gc.collect() 
        
        high_repeat_region_abundance = get_relative_abundance_in_connected_regions(read_to_location,threads,100)
        read_and_process_file(multi_sam,batchsize)

        tmp_names1 = outdir + "/" + out_pre + "_tmp_batch_*_success.sam"
        tmp_names2 = outdir + "/" + out_pre + "_tmp_batch_*_fail.sam"
        tmp_merge_name1 = outdir + "/" + out_pre + "_reassign_success"
        tmp_merge_name2 = outdir + "/" + out_pre + "_reassign_fail"

        os.system("cat %s > %s" % (tmp_names1,tmp_merge_name1))
        os.system("cat %s > %s" % (tmp_names2,tmp_merge_name2))
        all_tmp_names = outdir + "/" + out_pre + "_tmp_batch_*.sam "
        print("Removing tmp files...")
        os.system('rm %s' % all_tmp_names)

        mergesam = outdir + "/" + out_pre + "_final.sam"
        mergesam_sorted = outdir + "/" + out_pre + "_final_sorted.bam"

        subprocess.call("cat " + uniq_sam + " " + tmp_merge_name1 + " > " + mergesam, shell=True)
        subprocess.call("samtools view -h -bS " + mergesam + " | samtools sort -@ " + str(threads) + " -o " + mergesam_sorted, shell=True)
        subprocess.call("samtools index " + mergesam_sorted, shell=True)
        subprocess.call("rm " + mergesam, shell=True)