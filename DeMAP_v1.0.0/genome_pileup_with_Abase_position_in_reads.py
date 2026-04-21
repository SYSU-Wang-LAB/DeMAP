#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Zehong Wu wang_lab, Sun Yat-sen University"""
"""Feb, 2025"""
"""Email: wuzh93@mail2.sysu.edu.cn"""
"""Usage: This program is used for genome pileup"""

import glob
from heapq import merge
import subprocess
import sys
import os
import argparse
import pysam
from Bio.Seq import reverse_complement
from Bio import SeqIO
import multiprocessing
import time
from time import strftime

parser = argparse.ArgumentParser(prog="pileup_genome_multiprocessing",fromfile_prefix_chars='@',description="Gonome pileup",formatter_class=argparse.RawTextHelpFormatter)
#Required
group_required = parser.add_argument_group("Required")
group_required.add_argument("-bam",dest="input",help="Sorted bam and indexed bam input")
group_required.add_argument("-o",dest="output",default="genome_merge.pileup.txt",help="Output text")
group_required.add_argument("-fa",dest="fasta",help="Unconverted genome file in fasta format")
#Filter
group_optional = parser.add_argument_group("Optional")
group_optional.add_argument("-p",dest="process",type=int,default=1,help="Process number, default=1")
group_optional.add_argument("-q",dest="qual",type=int,default=10,help="Quality filter, default=10")
group_optional.add_argument("-s",dest="step",type=int,default=50000,help="Steps for pileup, default=50000")
group_optional.add_argument("-m","--max-depth",dest="max_depth",type=int,default=1000000,help="Max depth for pileup, default=1000000")
group_optional.add_argument("--trim-head",dest="trim_head",type=int,default=0,help="Trim the head of reads, default=0")
group_optional.add_argument("--trim-tail",dest="trim_tail",type=int,default=0,help="Trim the tail of reads, default=0")
group_optional.add_argument("--rm_5end_softclip", action="store_true", help="If specified, Remove reads with soft-clipped starts, as alignment accuracy drops after soft-clipping.\n"
                "This is recommended only when soft-clipped reads are rare (e.g., <20%%) and for m6Am-specific analysis (e.g., CROWN-seq).\n")
group_optional.add_argument("-max_A",dest="max_A",type=int,default=3,help="Filter reads with more than max_A unconverted-A in a read.\n"
                            "If a read contains more than max_A unconverted-A, this read may not be fully converted. default=3")
group_optional.add_argument("-min_convert", dest="min_convert", type=float, default=0.8, help="If A->G conversion rate >= this value, keep the read even if unconverted A > max_A. default=0.8")
group_other = parser.add_argument_group("Other")
group_other.add_argument("--version",action="version",version="%(prog)s 1.3")
options = parser.parse_args()

rm_5end_softclip = options.rm_5end_softclip

def worker(conting,start,stop):
    global parent_pid,dictRefSeq,options
    pid = os.getpid()
    tmp_file_name = 'tmp_' + str(parent_pid) + '_' + str(pid)
    max_depth = options.max_depth

    max_depth = options.max_depth
    with pysam.AlignmentFile(options.input, 'rb') as samfile, open(tmp_file_name,'a') as tmp_file:
        try:
            for pileupcolumn in samfile.pileup(conting,start=start,stop=stop,max_depth=max_depth,ignore_orphans=False,ignore_overlaps=False,min_base_quality=0,truncate=True):
                chr = pileupcolumn.reference_name
                ref_seq = dictRefSeq.get(pileupcolumn.reference_name)
                ref_base = ref_seq[pileupcolumn.pos].upper()
                pos = pileupcolumn.pos + 1
                PF_positive_base = {}
                PF_negative_base = {}
                PF_positive_A = []
                PF_negative_A = []
                positive_base_position_to_read = {}
                negative_base_position_to_read = {}

                for pileupread in pileupcolumn.pileups:
                    if pileupread.query_position is not None and pileupread.alignment.query_qualities[pileupread.query_position] >= options.qual:
                        not_at_end = True
                        if pileupread.alignment.is_reverse:
                            if options.trim_tail and pileupread.query_position < options.trim_tail:
                                not_at_end = False
                            if not_at_end == True and options.trim_head and pileupread.alignment.query_length - options.trim_head <= pileupread.query_position:
                                not_at_end = False
                        else:
                            if options.trim_tail and pileupread.alignment.query_length - options.trim_tail < pileupread.query_position:
                                not_at_end = False
                            if not_at_end == True and pileupread.query_position < options.trim_head:
                                not_at_end = False

                        if not_at_end:
                            query_name,query_base,query_position,align_len,align_start,align_end = pileupread.alignment.query_name,\
                                pileupread.alignment.query_sequence[pileupread.query_position],\
                                pileupread.query_position,pileupread.alignment.query_alignment_length,\
                                pileupread.alignment.query_alignment_start,pileupread.alignment.query_alignment_end
                            cigartuples = pileupread.alignment.cigartuples
                            has_5prime_softclip = (cigartuples[0][0] == 4)
                            has_3prime_softclip = (cigartuples[-1][0] == 4)


                            ref_slice = str(ref_seq[align_start:align_end]).upper()
                            total_ref_A = ref_slice.count('A')
                            total_ref_T = ref_slice.count('T')


                            if rm_5end_softclip:
                                if pileupread.alignment.is_reverse and (not has_3prime_softclip):
                                    #header_softlen_negative.append("0")
                                    A_counts = pileupread.alignment.query_sequence.count('T')
                                    PF_negative_A.append(A_counts)

                                    conversion_rate = (total_ref_T - A_counts) / total_ref_T if total_ref_T > 0 else 0
                                    if A_counts <= options.max_A or conversion_rate >= options.min_convert:
                                        base_position_and_alignment_length = str(align_end - query_position)+ "_" + str(align_len)
                                        query_base = reverse_complement(query_base)
                                        if PF_negative_base.get(query_name) is None:
                                            PF_negative_base[query_name] = query_base
                                            negative_base_position_to_read[query_name] = "0_"+query_base+"_"+base_position_and_alignment_length

                                if (not pileupread.alignment.is_reverse) and (not has_5prime_softclip):
                                    #header_softlen_positive.append("0")
                                    A_counts = pileupread.alignment.query_sequence.count('A')
                                    PF_positive_A.append(A_counts)
                                    conversion_rate = (total_ref_A - A_counts) / total_ref_A if total_ref_A > 0 else 0
                                    if A_counts <= options.max_A or conversion_rate >= options.min_convert:
                                        base_position_and_alignment_length = str(query_position - align_start + 1) + "_" + str(align_len)
                                        if PF_positive_base.get(query_name) is None:
                                            PF_positive_base[query_name] = query_base
                                            positive_base_position_to_read[query_name] = "0_"+query_base+"_"+base_position_and_alignment_length

                            else:
                                if pileupread.alignment.is_reverse:
                                    A_counts = pileupread.alignment.query_sequence.count('T')
                                    PF_negative_A.append(A_counts)
                                    conversion_rate = (total_ref_T - A_counts) / total_ref_T if total_ref_T > 0 else 0
                                    if A_counts <= options.max_A or conversion_rate >= options.min_convert:
                                        base_position_and_alignment_length = str(align_end - query_position)+ "_" + str(align_len)
                                        query_base = reverse_complement(query_base)
                                        if PF_negative_base.get(query_name) is None:
                                            PF_negative_base[query_name] = query_base
                                            if has_3prime_softclip:
                                                negative_base_position_to_read[query_name] = str(cigartuples[-1][1]) + "_" + query_base+"_"+base_position_and_alignment_length
                                            else:
                                                negative_base_position_to_read[query_name] = "0_" + query_base+"_"+base_position_and_alignment_length

                                else:
                                    A_counts = pileupread.alignment.query_sequence.count('A')
                                    PF_positive_A.append(A_counts)
                                    conversion_rate = (total_ref_A - A_counts) / total_ref_A if total_ref_A > 0 else 0
                                    if A_counts <= options.max_A or conversion_rate >= options.min_convert:
                                        base_position_and_alignment_length = str(query_position - align_start + 1) + "_" + str(align_len)
                                        if PF_positive_base.get(query_name) is None:
                                            PF_positive_base[query_name] = query_base
                                            if has_5prime_softclip:
                                                positive_base_position_to_read[query_name] = str(cigartuples[0][1]) + "_" + query_base+"_"+base_position_and_alignment_length
                                            else:
                                                positive_base_position_to_read[query_name] = "0_" + query_base+"_"+base_position_and_alignment_length

                list_tmp_positive = PF_positive_base.values()
                list_tmp_negative = PF_negative_base.values()
                if len(list_tmp_positive) >= 1:
                    positive_seq = PF_positive_base.values()
                    seq_base = ','.join(positive_seq)
                    read_with_maxA_ratio = str(round(len([A_in_read for A_in_read in PF_positive_A if A_in_read >= options.max_A])/len(PF_positive_A),3)) 
                    A = str(seq_base.count('A'))
                    T = str(seq_base.count('T'))
                    C = str(seq_base.count('C'))
                    G = str(seq_base.count('G'))
                    relativeposition = ",".join(positive_base_position_to_read.values())   
                    cleavage_ratio = str(relativeposition.count('_1_')) + ";" + str(len(positive_base_position_to_read.values()))                 
                    if ref_base == "A":                        
                        tmp_file.write("\t".join([chr,str(pos),'+',ref_base,A,T,C,G,str(len(positive_seq)),read_with_maxA_ratio,cleavage_ratio+","+relativeposition]))
                    else:
                        tmp_file.write("\t".join([chr,str(pos),'+',ref_base,A,T,C,G,str(len(positive_seq)),"NA",cleavage_ratio]))
                    tmp_file.write("\n")
                    
                if len(list_tmp_negative) >= 1:
                    negative_seq = PF_negative_base.values()
                    seq_base = ','.join(negative_seq)
                    read_with_maxA_ratio = str(round(len([A_in_read for A_in_read in PF_negative_A if A_in_read >= options.max_A])/len(PF_negative_A),3)) 
                    A = str(seq_base.count('A'))
                    T = str(seq_base.count('T'))
                    C = str(seq_base.count('C'))
                    G = str(seq_base.count('G'))                    
                    ref_base = reverse_complement(ref_base)
                    relativeposition = ",".join(negative_base_position_to_read.values())
                    cleavage_ratio = str(relativeposition.count('_1_')) + ";" + str(len(negative_base_position_to_read.values()))
                    if ref_base == "A":                        
                        tmp_file.write("\t".join([chr,str(pos),'-',ref_base,A,T,C,G,str(len(negative_seq)),read_with_maxA_ratio,cleavage_ratio+","+relativeposition]))
                    else:
                        tmp_file.write("\t".join([chr,str(pos),'-',ref_base,A,T,C,G,str(len(negative_seq)),"NA",cleavage_ratio]))
                    tmp_file.write("\n")
        except ValueError:
            print("Conting [%s] does not exist in @SQ header, pass" % conting)

def sort_columns(inline):
    info_for_sort = inline.strip().split("\t")
    info_for_sort = tuple([int(info_for_sort[x]) if info_for_sort[x].isdigit() else info_for_sort[x] for x in [0,2,1]]) #sort by chrome,strand,position
    return info_for_sort

def sort_files(file_for_sort,file_sorted,chunksize):
    chunkfiles = file_for_sort+"_chunk_*"
    chunkid = 1
    lines = []
    with open(file_for_sort, 'r') as f_in:
        f_out = open(file_for_sort+'_chunk_{}'.format(chunkid), 'w')
        for line_num, line in enumerate(f_in, 1):
            lines.append(line)
            if not line_num % chunksize:
                lines = sorted(lines, key=sort_columns)  
                f_out.writelines(lines)
                f_out.close()
                lines = []
                chunkid += 1
                f_out = open(file_for_sort+'_chunk_{}'.format(chunkid), 'w')
        # last chunk
        if lines:
            lines = sorted(lines, key=sort_columns)
            f_out.writelines(lines)
            f_out.close()
            lines = []

    chunks = []
    for filename in glob.glob(chunkfiles):
        chunks += [open(filename, 'r')]

    with open(file_sorted, 'w') as f_out:
        f_out.writelines(merge(*chunks, key=sort_columns))
    subprocess.call("rm -f " + file_for_sort + "_chunk_*", shell=True)
    subprocess.call("rm -f " + file_for_sort, shell=True)


if __name__ == "__main__":
    parent_pid = os.getpid()
    if os.path.exists(options.output):
        sys.stderr.write("[%s] Pileup file already exists, skip this precess\n" % (strftime("%Y-%m-%d %H:%M:%S", time.localtime())))
        pass
    else:
        sys.stderr.write("[%s] Pileup genome, processers: %d, pid: %d\n" % (strftime("%Y-%m-%d %H:%M:%S", time.localtime()),options.process ,parent_pid))
        dictRefSeq = {}
        RefBins = []
        for seq in SeqIO.parse(options.fasta,"fasta"):
            if seq.id not in dictRefSeq:
                dictRefSeq[str(seq.id)] = seq.seq
                if "chrM" in seq.id:
                    for start in range(0,len(seq.seq),200):
                        RefBins.append((seq.id,start,start+200))
                else:
                    for start in range(0,len(seq.seq),options.step):   # split genome into bins
                        RefBins.append((seq.id,start,start+options.step))
                        
        pool = multiprocessing.Pool(options.process)
        try:
            for item in RefBins:
                chr,start,stop = item
                pool.apply_async(func=worker,args=(chr,start,stop,))
            pool.close()
            print("Merging tmp files222...")
            pool.join()
            print("Merging tmp files333...")
            tmp_names = "tmp_" + str(parent_pid) + "*" 

            print("Merging tmp files...")
            sys.stderr.write("[%s] Merging TEMPs\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
            os.system("cat %s > %s" % (tmp_names,options.output))
            print("Removing pileup tmp files...")
            os.system('rm %s' % tmp_names)
        finally:
            pool.terminate()
        sys.stderr.write("[%s] Genome pileup finished.\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))

    sys.stderr.write("[%s] start sort pielup file \n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))

    sort_files(options.output,options.output+"_sorted",10000000)
    sys.stderr.write("[%s] Done! \n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))