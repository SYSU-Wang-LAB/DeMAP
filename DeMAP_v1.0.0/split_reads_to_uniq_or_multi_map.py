#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Zehong Wu wang_lab, Sun Yat-sen University"""
"""Feb, 2025"""
"""Email: wuzh93@mail2.sysu.edu.cn"""
"""Usage: This script separates the reads in the merged SAM file into multi-mapped and uniquely mapped reads."""


import argparse
import re


parser = argparse.ArgumentParser(description = "Merge AG and TC samfile and split them to uniq-map and multi-map reads",formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-pre", "--outname_prefix", nargs="?", type=str, default='default',help = "--outname_prefix")
parser.add_argument("-o", "--outputdir", nargs="?", type=str, default="none", help="outputdir")
parser.add_argument("-p", "--Threads", nargs="?", type=str, default='1', help = "The number of alignment threads")
parser.add_argument("-max", "--maxmap", nargs="?", type=int, default=5, help = "Discard multi-mapping reads with more than maxmap mappings")
parser.add_argument("-s1", "--score1", nargs="?", type=int, default=2, help = "Discard low quality alignment with AS (alignment score) < max(AS) - score1 in a multi-mapping read" + "\n" + \
                    "(2 in default: Best alignment has more than 2 bases match to reference)")
parser.add_argument("-s2", "--minScore", nargs="?", type=int, default=25, help = "Discard read if AS < minScore in all alignments of that read. (defalut:25)")

args = parser.parse_args()


outdir = args.outputdir
out_pre = args.outname_prefix
Threads = args.Threads
maxmap = args.maxmap # if args.maxmap > 1 else 2
score1 = args.score1
minScore = args.minScore


def split_merge_sam(mergesam,uniqsam,multisam):
    with open(mergesam, 'r') as f1, open(uniqsam, 'w') as uniq_out, open(multisam, 'w') as multi_out:
        record_dic = {}
        for line in f1:
            if line.startswith('@'):
                uniq_out.write(line)
                multi_out.write(line)
            else:
                alignment_info = line.split("\t")
                read_id = alignment_info[0]
                if read_id not in record_dic and len(record_dic) == 0: # first line
                    record_dic[read_id] = [line]
                elif read_id in record_dic:
                    record_dic[read_id].append(line)
                elif read_id not in record_dic and len(record_dic) != 0: # meet a new read id
                    previous_read_id = list(record_dic.keys())[0]
                    if len(record_dic[previous_read_id]) >= 2: # multi-alignment
                        #print(previous_read_id)
                        mapping_score = [int(re.search(r"AS:i:(\d+)", alignment).group(1)) for alignment in record_dic[previous_read_id]]
                        if max(mapping_score) < minScore:
                            record_dic.clear()
                            record_dic[read_id] = [line]
                            continue
                        # filter alignments with extremely low mapping score compare to best alignment
                        filtered_alignments = [alignment for AS, alignment in zip(mapping_score, record_dic[previous_read_id]) if max(mapping_score) - AS <=  score1]
                        cigars = ";".join([alignment.split("\t")[5] for alignment in filtered_alignments])
                        scaffolds = [alignment.split("\t")[2] for alignment in filtered_alignments if 'chr' not in alignment]
                        if "N" in cigars or len(filtered_alignments) > maxmap or len(scaffolds) > 0:  # filter read cigar with N and reads mapped to too many locations
                            record_dic.clear()
                            record_dic[read_id] = [line]
                            continue

                        if len(filtered_alignments) < 2:
                            tmpstr = re.sub(r'(NH|HI):i:\d+', r'\1:i:1', filtered_alignments[0])
                            mapping_score = int(re.search(r"AS:i:(\d+)", tmpstr).group(1))
                            if mapping_score >= minScore:
                                uniq_out.write(tmpstr)
                            record_dic.clear()
                            record_dic[read_id] = [line]
                            continue
                        else:
                            mapping_num = len(filtered_alignments)
                            for multi_num_index in range(0,mapping_num):
                                tmpstr = re.sub(r'NH:i:(\d+)', f'NH:i:{mapping_num}', filtered_alignments[multi_num_index])
                                tmpstr1 = re.sub(r'HI:i:(\d+)', f'NI:i:{multi_num_index + 1}', tmpstr)
                                multi_out.write(tmpstr1)
                            record_dic.clear()
                            record_dic[read_id] = [line]
                    else:
                        mapping_score = int(re.search(r"AS:i:(\d+)", record_dic[previous_read_id][0]).group(1))
                        if mapping_score >= minScore and 'chr' in record_dic[previous_read_id][0]:
                            uniq_out.write(record_dic[previous_read_id][0])
                        record_dic.clear()
                        record_dic[read_id] = [line]
        # output last read            
        last_read_id = list(record_dic.keys())[0]
        if len(record_dic[last_read_id]) >= 2: # multi-alignment
            mapping_score = [int(re.search(r"AS:i:(\d+)", alignment).group(1)) for alignment in record_dic[last_read_id]]
            cigars = ";".join([alignment.split("\t")[5] for alignment in record_dic[last_read_id]])
            scaffolds = [alignment.split("\t")[2] for alignment in record_dic[last_read_id] if 'chr' not in alignment]
            if "N" not in cigars and max(mapping_score) >= minScore and len(scaffolds) == 0:  # read cigar with N 
                # filter alignments with extremely low mapping score compare to best alignment
                filtered_alignments = [alignment for AS, alignment in zip(mapping_score, record_dic[last_read_id]) if max(mapping_score) - AS <=  score1]
                if len(filtered_alignments) < 2:
                    tmpstr = re.sub(r'(NH|HI):i:\d+', r'\1:i:1', filtered_alignments[0])
                    mapping_score = int(re.search(r"AS:i:(\d+)", tmpstr).group(1))
                    if mapping_score >= minScore:
                        uniq_out.write(tmpstr)
                else:
                    if len(filtered_alignments) <= maxmap:
                        mapping_num = len(filtered_alignments)
                        for multi_num_index in range(0,mapping_num):
                            tmpstr = re.sub(r'NH:i:(\d+)', f'NH:i:{mapping_num}', filtered_alignments[multi_num_index])
                            tmpstr1 = re.sub(r'HI:i:(\d+)', f'NI:i:{multi_num_index + 1}', tmpstr)
                            multi_out.write(tmpstr1)
        else:
            mapping_score = int(re.search(r"AS:i:(\d+)", record_dic[last_read_id][0]).group(1))
            if mapping_score >= minScore and 'chr' in record_dic[last_read_id][0]:
                uniq_out.write(record_dic[last_read_id][0])


if __name__ == "__main__":
    ###### merge AG and TC alinement reads
    merge_sam = outdir + "/" + out_pre + "_merge_AG_TC_correct.sam"
    ##### split merge sam to uniq or multi alignment reads
    uniq_sam = outdir + "/" + out_pre + "_" + "uniq_map.sam"
    multi_sam = outdir + "/" + out_pre + "_" + "multi_map.sam"
    split_merge_sam(merge_sam,uniq_sam,multi_sam)



    