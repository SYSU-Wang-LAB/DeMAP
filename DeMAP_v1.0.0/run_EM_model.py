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
    """
    创建基因组分区，每个分区的大小为 bin_size。
    """
    genome_bins = {}
    for chrom, length in chrom_lengths.items():
        for strand in ["+", "-"]:
            chrom_strand = chrom + ";" + strand
            bin_number = (length + bin_size) // bin_size
            for bin in range(0,bin_number):
                genome_bins[chrom_strand + "_" + str(bin)] = [0,0,[]]

    return genome_bins

def reverse_complement(seq):
    """生成反向互补序列（仅处理大写ATCGN）"""
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
            bin = read_center_pos // bin_size
            target_bin = chrom_strand + "_" + str(bin)
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
                bin = read_center_pos // bin_size
                target_bin = chrom_strand + "_" + str(bin)

                if readid not in read_to_location:
                     read_to_location[readid] = [target_bin]  # 初始化 read id 的记录
                else:
                    read_to_location[readid].append(target_bin)  # 追加比对信息

                genome_bins[target_bin][1] += 1  # 增加比对计数
                genome_bins[target_bin][2].append(readid)  # 记录 read id

    if total_multi_reads <= max_multi_reads:
        used_ratio = 1
    else:
        used_ratio = max_multi_reads/total_multi_reads

    used_region_count = {}
    for target_bin,region_counts in genome_bins.items():
        uniqcounts = region_counts[0]
        multicounts = region_counts[1]
        bin_read_id = region_counts[2]
        if len(bin_read_id) >= maxcount:
            bin_read_id = np.random.choice(bin_read_id, maxcount, replace=False)
        if uniqcounts != 0 or multicounts != 0:  # count be changes
            used_region_count[target_bin] = [uniqcounts,bin_read_id] #[] for record readname of multi-mapping reads

    genome_bins.clear()
    return used_region_count,read_to_location,used_ratio

def run_EM_model(connect_regions,maxsizenum=600,maxiter=50,mintol=1e-2):
    expression = {}
    reads_in_connect_region = []
    for region in connect_regions:
        if region not in genome_bins_for_unsed:
            expression[region] = 1
            reads_in_connect_region.extend([])
        else:
            expression[region] = genome_bins_for_unsed[region][0] + 1
            reads_in_connect_region.extend(genome_bins_for_unsed[region][1])

    if not reads_in_connect_region: # unreliable multi-mapping reads
        return None


    valid_reads = []
    connect_regions_set = set(connect_regions) # 转换为 set 提高查找速度

    # 严格过滤：只保留其所有映射区域都在当前 connect_regions 内的 reads
    for read in reads_in_connect_region:
        read_mapped_regions = set(read_to_location.get(read, []))
        # 如果该 read 映射的区域是当前 connect_regions 的子集，才保留
        if read_mapped_regions.issubset(connect_regions_set):
            valid_reads.append(read)
            
    all_reads = list(set(valid_reads))
    multi_readnum = len(all_reads)

    # 如果过滤后没有 reads 了，直接返回
    if multi_readnum == 0:
        return None
    
    #all_reads = list(set(reads_in_connect_region))
    #multi_readnum = len(all_reads)
    if multi_readnum > maxsizenum:
        all_reads = np.random.choice(all_reads, maxsizenum, replace=False)
        scale_size = (multi_readnum/used_ratio) / maxsizenum
        for region in expression:
            expression[region] = expression[region] / scale_size
    else:
        for region in expression:
            expression[region] = expression[region] * used_ratio

    # region index :{'region1': 0, 'region2': 1, 'region3': 2}
    connect_regions = list(expression.keys())
    region_idx = {region: i for i, region in enumerate(connect_regions)}

    # 构建 sparse matrix (R × C)
    row_indices = []
    col_indices = []
    data = []
    for r_idx, read in enumerate(all_reads):
        for region in read_to_location.get(read, []):
            if region in connect_regions:
                row_indices.append(r_idx)
                col_indices.append(region_idx[region])
                data.append(1)

    # 用稀疏矩阵存储 reads 与 regions 的关系
    read_region_matrix = csr_matrix((data, (row_indices, col_indices)), shape=(len(all_reads), len(connect_regions)))

    # 初始化表达量数组 & 初始化概率数组
    expression_array = np.array([expression[region] for region in connect_regions], dtype=np.float64)
    total_counts = expression_array.sum()
    prob_array = (expression_array / total_counts) 

    ######### ******** EM-MODEL start ********** ################ 

    # start iteration
    for _ in range(maxiter):
        # **E步: 计算所有 reads 分配在各个 region 上的概率**
        # 先构造加权矩阵：将 read-region 的非零值乘以对应的 p_j
        weighted_matrix = read_region_matrix.multiply(prob_array)  # 每个非零条目变为 p_j
        # 计算每个 read 的分母：映射区域 p 值之和
        total_p = weighted_matrix.sum(axis=1).A1 
        # 计算权重矩阵: 每个 read 在每个映射区域的权重为 p_j / sum(p_k)
        weights_matrix = weighted_matrix.multiply(1 / total_p[:, None])

        # **M步: 更新表达量**
        # 唯一比对表达量不变，加上所有 read 分配的权重之和
        new_expression_array = expression_array + weights_matrix.sum(axis=0).A1
        
        # **计算新概率**
        total_counts = new_expression_array.sum()
        new_prob_array = new_expression_array / total_counts

        # check convergence
        max_diff = np.max(np.abs(new_prob_array - prob_array))
        if max_diff < mintol:
            break

        # updata values
        prob_array = new_prob_array

    final_expression_ratio = list(new_prob_array)
    return final_expression_ratio
    ######### ******** EM-MODEL end ********** ################     


def process_chunk(chunk, sizenum, max_iter, tol):
    result = {}
    for connect_regions, _ in chunk:
        connect_bins = connect_regions.split(";;")
        ratio = run_EM_model(connect_bins, sizenum, max_iter, tol)
        if ratio is not None:
            result[connect_regions] = {connect_bins[i]: ratio[i] for i in range(len(connect_bins))}
    return result


def get_relative_abundance_in_connected_regions(read_to_location, num_workers=4, chunk_size=100):
    """
    get relative abundance of connected-regions using EM model
    """
    relative_region_abundance = {}
    connect_region_counts = {}
    
    for read in read_to_location:
        if len(set(read_to_location[read])) != len(read_to_location[read]):
            continue
        connect_regions = ";;".join(sorted(read_to_location[read]))
        connect_region_counts[connect_regions] = connect_region_counts.get(connect_regions, 0) + 1

    sorted_regions = sorted(connect_region_counts.items(), key=lambda x: x[1], reverse=True)#[:50000] # top 50000 connected regions

    worker_func = partial(process_chunk, sizenum=sizenum, max_iter=max_iter, tol=tol)
    chunks = [sorted_regions[i:i+chunk_size] for i in range(0, len(sorted_regions), chunk_size)]
    n_process = num_workers if num_workers <= 4 else 4

    with multiprocessing.Pool(processes=n_process, maxtasksperchild=50) as pool:
        chunk_results = pool.map(worker_func, chunks)
    
    for res in chunk_results:
        relative_region_abundance.update(res)
    
    return relative_region_abundance

def rank_list_unique(lst):
    sorted_lst = sorted(set(lst))  # 去重并排序
    rank_dict = {val: rank + 1 for rank, val in enumerate(sorted_lst)}  # 构造排名字典
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

        ######### ******** reassign multi-mapping reads accroding to region expression level and mappings Score ****** ---------- 
        Mapping_Score = [int(re.search(r'AS:i:(\d+)', alignment).group(1)) for alignment in Alignments]
        Mapping_Score_ranking = rank_list_unique(Mapping_Score) #[AS:33,29,29,28] -> [3,2,2,1]
        ranking_Score = [score ** ScoreFactor for score in Mapping_Score_ranking] #ScoreFactor = 3: [3,2,2,1] -> [27,8,8,1]

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


def process_batch(onebatch,batch_id):
    """
    处理一个批次的数据
    """ 
    info_record = {}
    for alignment in onebatch:
        read_id = alignment.query_name
        sam_text = alignment.to_string()
        chrom = alignment.reference_name
        strand = "+" if not alignment.is_reverse else "-"
        chrom_strand = chrom + ";" + strand
        read_at_reference_position = alignment.get_reference_positions()
        center_index = round(len(read_at_reference_position)/2)
        read_center_pos = read_at_reference_position[center_index]
        bin_index = read_center_pos // winsize
        target_bin = chrom_strand + "_" + str(bin_index)
        if read_id not in info_record:
            info_record[read_id] = [[target_bin],[sam_text]]
        else:
            info_record[read_id][0].append(target_bin)
            info_record[read_id][1].append(sam_text)

    onebatch_info = list(info_record.values())
    del info_record

    chunk_size = 1000 
    chunks = [onebatch_info[i:i+chunk_size] for i in range(0, len(onebatch_info), chunk_size)]

    with multiprocessing.Pool(processes=threads) as pool:
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
    """
    逐行读取文件，并按第一列的值累积批次
    """
    batch = []  
    unique_readid = set() 
    batch_id = 1
    with pysam.AlignmentFile(multi_file, "rb") as multi_align:
        for alignment in multi_align:
            read_id = alignment.query_name

            if read_id not in unique_readid:
                if len(unique_readid) >= max_readnum:

                    process_batch(batch,batch_id)  
                    gc.collect() 
                    batch_id += 1
                    batch = []  
                    unique_readid = set()  

                unique_readid.add(read_id)
            batch.append(alignment)

        if batch:
            process_batch(batch,batch_id)


if __name__ == "__main__":
    parent_pid = os.getpid()
    # extract chromesome length
    with pysam.AlignmentFile(uniq_sam, "r") as sam_in:
        chrom_lengths = {ref["SN"]: ref["LN"] for ref in sam_in.header["SQ"]}

    # generate genomic bins
    genome_bins = create_genome_bins(chrom_lengths, winsize)

    mergesam = outdir + "/" + out_pre + "_final.sam"
    mergesam_sorted = outdir + "/" + out_pre + "_final_sorted.bam"

    if maxmap == 1:
        sys.stderr.write("[%s] Only using uniq-mapping reads, skip EM model.\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
        subprocess.call("samtools view -h -bS " + uniq_sam + " | samtools sort -@ " + str(threads) + " -o " + mergesam_sorted, shell=True)
        subprocess.call("samtools index " + mergesam_sorted, shell=True)
    else:
        # generate region counts for uniq_bam and filter out unsed region
        sys.stderr.write("[%s] Start running EM to assign multi-mapping reads.\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))

        genome_bins_for_unsed,read_to_location,used_ratio = initial_bin_counts(uniq_sam,multi_sam,genome_bins,winsize,multi_size)
        high_repeat_region_abundance = get_relative_abundance_in_connected_regions(read_to_location,threads,100)
        read_and_process_file(multi_sam,batchsize)

        tmp_names1 = outdir + "/" + out_pre + "_tmp_batch_*_success.sam"
        tmp_names2 = outdir + "/" + out_pre + "_tmp_batch_*_fail.sam"
        tmp_merge_name1 = outdir + "/" + out_pre + "_reassign_success"
        tmp_merge_name2 = outdir + "/" + out_pre + "_reassign_fail"

        #sys.stderr.write("[%s] Merging TEMPs\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
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
    