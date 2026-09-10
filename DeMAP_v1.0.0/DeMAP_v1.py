#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Zehong Wu wang_lab, Sun Yat-sen University"""
"""Feb, 2025"""
"""Email: wuzh93@mail2.sysu.edu.cn"""
"""Usage: This is the main script of the DeMAP pipeline."""

import argparse
import subprocess
import sys
import time
from time import strftime
import os


parser = argparse.ArgumentParser(description="Mapping reads of GLORI data",formatter_class=argparse.RawTextHelpFormatter)


Required = parser.add_argument_group("Required Parameters")
Required.add_argument("-software", "--software", type=str, help="The path to the tool directory.")
Required.add_argument("-genomedir", "--genomedir", type=str, help = "The path to the genome index directory.")
Required.add_argument("-fa","--fasta",type=str,help="Unconverted genome file in FASTA format.")
Required.add_argument("-gtf", "--gtf", nargs="?", type=str, default="none", help = "GTF file from gencode")
Required.add_argument("-q", "--fastq", type=str, help="Fastq file for single-end data; read should be same direction to sense strand.\n"
                    "For most widely used stranded library, dUTP-based method and SMART mothod: R2. Input files in gzip format (.gz) are supported.")
Required.add_argument("-pre", "--outname_prefix", type=str,  help="Output name prefix")
Required.add_argument("-o", "--outputdir", type=str,  help="Output directory")

# pre-alignment
mapping = parser.add_argument_group("Mapping Parameters (Optional)")
mapping.add_argument("-p", "--Threads", type=str, default='1', 
                     help="The number of threads. (default:1)\n"
                           "In the EM model, more than 20 threads may consume over 100 GB of memory.")
mapping.add_argument("-m", "--mismatch", type=str, default='2', help="Maximum mapping mismatch. (default:2)")
mapping.add_argument("-max", "--maxmap", type=str, default='10', help="Only output the best top <maxmap> alignments for a multi-mapping reads. (default:10)\n"
                     "EM model will select a best alignment for the multi-mapping reads")
mapping.add_argument("--re_align", action="store_true", help="If specified, split and re-align reads that fail to map to both AG- and TC-converted genomes.")


# correct softclip due to C-to-U mutation
correct_softclip = parser.add_argument_group("Correct softclip due to C-to-U mutation")
correct_softclip.add_argument("-SC_pval", "--SC_pval", nargs="?", type=str, default="0.001", help="P-value for the cutoff of LLR test. (default:0.001)")
correct_softclip.add_argument("-match_ratio", "--match_ratio", nargs="?", type=str, default="0.8", help="The minimum ratio of bases match to reference in softclip region. (default:0.8)")
correct_softclip.add_argument("-minSC", "--min_softclip", nargs="?", type=str, default='5', help = "The minimum softclip length included in correction. (default:5)")


# merge and split reads
merge_alignemts = parser.add_argument_group("Merge AG,TC alignemts and reclassify reads into uniq/multi-mapping reads")
merge_alignemts.add_argument("-s1", "--score1", type=str, default='2', 
                    help="Discard low-quality alignments with AS (alignment score) < max(AS) - score1 in a multi-mapping read. (default:2)\n"
                         "This parameter removes low-quality secondary alignments from multiply aligned reads.")
merge_alignemts.add_argument("-s2", "--minScore", type=str, default='22', help="Discard a read if its AS < minScore in all alignments of that read. (default:22)\n"
                    "For STAR aligner, The AS equals matched bases minus penalty terms (mismatches, insertions, deletions, soft clipping, and splice junctions).")

# run EM model
runEM = parser.add_argument_group("EM model for reassigning multi-mapping reads")
runEM.add_argument("-w", "--winsize", type=str, default='150', help="window size (bp) for searching read alignment region. (default:150)")
runEM.add_argument("-size", "--size", type=str, default='800', help = "The max number of multi-mapping reads used to estimate the probability of\n"
                    "the connected-region for a multi-mapping read. (defalut:800: to reduce memory consumption)")
runEM.add_argument("-bz", "--batchsize", type=str, default='1000000', help = "The number of multi-mapping reads process in a batch. (defalut:1000000)")
runEM.add_argument("-cutoff", "--cutoff", type=str, default='0.8', help = "Probability cutoff for reassigning multi-mapping reads. (default:0.8; cutoff range: 0.60-0.95)\n" 
                    "choose region when region score/sum(score in all regions) >= cutoff")

# calling methylation 
calling_methy = parser.add_argument_group("Calling methylation.\n"
                                          "When methylated sites meet the thresholds simultaneously, the priority is: m6Am > m1A >  m6A")
calling_methy.add_argument("-max_A",dest="max_A",type=str,default='3',help="Filter reads containing more than max_A unconverted A bases. (default=3)\n"
                            "Reads with more than max_A unconverted A bases may not be fully converted.")
calling_methy.add_argument("-max_A_ratio","--max_A_ratio",type=str, default='0.4',help="Filter sites with max_A reads /total reads > max_A_ratio.(default:0.4)\n"
                            "This parameter filters out modified sites resulting from incomplete conversion.")
calling_methy.add_argument("-min_convert", dest="min_convert", type=str, default='0.8', help="Keep reads with conversion rate >= this value even if unconverted A > max_A. default=0.8")

calling_methy.add_argument("--rm_5end_softclip", action="store_true", help="If specified, Remove reads with soft-clipped starts, as alignment accuracy drops after soft-clipping.\n"
                "This is recommended only when soft-clipped reads are rare (e.g., <20%%) and for m6Am-specific analyses (e.g., CROWN-seq).\n")
calling_methy.add_argument("-AGcov","--AGCoverage",type=str, default='15',help="Minimum A + G counts. (default:15)")
calling_methy.add_argument("-Acov","--ACoverage",type=str, default='5',help="Minimum A counts. (default:5)")
calling_methy.add_argument("-ratio","--ARatio",type=str, default='0.05',help="Minimum methylation level: A counts / AG counts. (default:0.05)")
calling_methy.add_argument("-enrich","--foldenrich",type=str, default='5',help="Minimum fold change between methylation level and background methylation level. (default:5)")
calling_methy.add_argument("-AGRatio","--AGRatio",type=str, default='0.6',help="Filter sites with A+G counts/total counts < AGRatio. (default:0.6)\n"
                            "This parameter filters out mutation sites caused by SNPs.")
calling_methy.add_argument("-A_neighbor","--A_neighbor",type=str, default='50',help="Number of A bases in the vicinity (on both sides) used as background.\n"
                           "(default:50 indicating 25 A bases on each side)")
calling_methy.add_argument("-fdr","--fdr",type=str, default='0.005',help="FDR thresholds for a site to be considered methylated. (default:0.005)")

calling_methy.add_argument("-m1A_fdr","--m1A_fdr",type=str, default='0.005',help="For m1A calling, FDR thresholds for m1A. (default:0.005)")
calling_methy.add_argument("-mut_ratio","--mut_ratio",type=str, default='0.03',help="For m1A calling, thresholds for other mutation ratio: T+C counts/total counts. (default:0.03)\n"
                            "During reverse transcription, m1A modifications can introduce sequencing errors.")
calling_methy.add_argument("-trunc_ratio","--trunc_ratio",type=str, default='0.05',help="For m1A calling, thresholds for truncation ratio:(Downstream - Upstream)/Downstream. (default:0.05)\n"
                            "During reverse transcription, m1A modifications cause sequencing truncations. 0.05 means 5%% truncations")
calling_methy.add_argument("-read_drop","--read_drop",type=str, default='10',help="For m1A calling, Downstream coverage - Upstream coverage >= read_drop in a truncated site. (default:10)\n")

calling_methy.add_argument("-softlen","--softlen",type=str, default='3',help="For m6Am calling, only using read with 5'end softclip length < softlen. (default:3)")
calling_methy.add_argument("-m6Am_fdr","--m6Am_fdr",type=str, default='0.005',help="For m6Am calling, FDR thresholds for m6Am. (default:0.005)")
calling_methy.add_argument("-start_percentile","--start_percentile",type=str, default='0.05',help="For m6Am calling, If a base is in the first 5%% of the read, it is considered to be at the start. (default:0.05)")
calling_methy.add_argument("-A_start","--A_start",type=str, default='3',help="For m6Am calling, minimum unconverted A bases located at the read start. (default:3)")
calling_methy.add_argument("-A_start_ratio","--A_start_ratio",type=str, default= '0.6',help="For m6Am calling, minimum ratio of unconverted A in read start/ all unconverted A. (default:0.6)")


args = parser.parse_args()
########### pre-alignment
tooldir = args.software
genomedir = args.genomedir
genome = args.fasta #genomedir+"/"+[f for f in os.listdir(genomedir) if f.endswith('genome.fa')][0]
gtf_file = args.gtf
fastq = args.fastq
Threads = args.Threads
mismatch = args.mismatch
maxmap = args.maxmap
re_align = args.re_align
outputdir = args.outputdir
outname_prx = args.outname_prefix

####### correct softclip due to C-to-U mutation 
SC_pval = args.SC_pval
match_ratio = args.match_ratio
min_softclip = args.min_softclip

########### # merge and split reads
score1 = args.score1
score2 = args.minScore


########### # run EM model
#dedup = args.dedup
winsize = args.winsize
maxsize = args.size
batchsize = args.batchsize
cutoff = args.cutoff

#childnum = args.child

################ calling methylation 
max_A = args.max_A
max_A_ratio = args.max_A_ratio
min_convert = args.min_convert
rm_5end_softclip = args.rm_5end_softclip

AG_Coverage = args.AGCoverage
A_Coverage = args.ACoverage
A_Ratio = args.ARatio
enrich_fold = args.foldenrich
AGRatio = args.AGRatio
A_neighbor = args.A_neighbor

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

################ methylation reslut filter

if __name__ == "__main__":

    ################# *** Mapping reads *** #####################
    finalbam = outputdir + "/" + outname_prx + "_final_sorted.bam"
    if os.path.exists(finalbam) and os.path.exists(finalbam + ".bai"):
        sys.stderr.write("[%s] Warning：bam files already exists, please make sure the bam file is correct. skip mapping step...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
        pass
    else:
        sys.stderr.write("[%s] Start mapping reads\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
        if re_align:
            command_mappingreads = "python " + tooldir+"/" + "mapping_reads_and_correct_softclip.py -software " + tooldir + " -q " + fastq + " -p " + Threads + " -genomedir " + genomedir + \
            " -fa " + genome + " -gtf " + gtf_file + " -m " + mismatch + " -max " + maxmap + " -SC_pval " + SC_pval + " -match_ratio " + match_ratio + " -minSC " + min_softclip + " -pre " + outname_prx + " -o " + outputdir + " -s1 " + score1 + " -s2 " + score2 + " --re_align "
        else:
            command_mappingreads = "python " + tooldir+"/" + "mapping_reads_and_correct_softclip.py -software " + tooldir + " -q " + fastq + " -p " + Threads + " -genomedir " + genomedir + \
            " -fa " + genome + " -gtf " + gtf_file + " -m " + mismatch + " -max " + maxmap + " -SC_pval " + SC_pval + " -match_ratio " + match_ratio + " -minSC " + min_softclip + " -pre " + outname_prx + " -o " + outputdir + " -s1 " + score1 + " -s2 " + score2
        subprocess.call(command_mappingreads, shell=True)
        sys.stderr.write("[%s] mapping reads done.\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))



        ################# *** Run EM model *** #####################
        uniqmap = outputdir + "/" + outname_prx + "_" + "uniq_map.sam"
        multimap = outputdir + "/" + outname_prx + "_" + "multi_map.sam"
        command_EM = "python " + tooldir+"/" + "run_EM_model.py -multi " + multimap + " -uniq " + uniqmap + " -pre " + outname_prx + " -o " + outputdir + " -w " + winsize + \
        " -size " + maxsize + " -bz " + batchsize +  " -cutoff " + cutoff  + " -p " + Threads + " -max " + maxmap
        subprocess.call(command_EM, shell=True)



    ################# *** Genome pileup *** #####################
    pileupfile = outputdir + "/" + outname_prx + ".pileup"
    if rm_5end_softclip:
        command_genome_pileup = "python " + tooldir + "/" + "genome_pileup_with_Abase_position_in_reads.py -bam " + finalbam + " -o " + pileupfile + " -fa " + genome + \
        " -p " + Threads + " -max_A " + max_A + " -min_convert " + min_convert +  " --rm_5end_softclip "
    else:
        command_genome_pileup = "python " + tooldir + "/" + "genome_pileup_with_Abase_position_in_reads.py -bam " + finalbam + " -o " + pileupfile + " -fa " + genome + \
        " -p " + Threads + " -max_A " + max_A + " -min_convert " + min_convert
    subprocess.call(command_genome_pileup, shell=True)


    ################# *** Calling methylation *** #####################
    command_calling_methylation = "python " + tooldir + "/" + "calculate_multi_methylation.py -pileup " + pileupfile + "_sorted" + " -fa " + genome + \
    " -pre " + outname_prx + " -o " + outputdir + " -AGcov " + AG_Coverage + " -Acov " + A_Coverage + " -ratio " + A_Ratio + " -max_A_ratio " + max_A_ratio + " -enrich " + enrich_fold + \
    " -AGRatio " + AGRatio + " -A_neighbor " + A_neighbor + " -fdr " + methy_fdr + " -m1A_fdr " + m1A_fdr + " -mut_ratio " + mut_ratio + " -read_drop " + read_drop + " -trunc_ratio " + trunc_ratio + \
    " -softlen " + softlen + " -m6Am_fdr " + m6Am_fdr + " -start_percentile " + start_percentile + " -A_start " + A_start + " -A_start_ratio " + A_start_ratio
    subprocess.call(command_calling_methylation, shell=True)


    ################# *** Methylation result annotation *** #####################
    site_result = outputdir+"/"+outname_prx+"_methylation_result.bed"
    site_result_anno = outputdir+"/"+outname_prx+"_methylation_result_anno.bed"
    annotation_file = tooldir + "/" + os.path.basename(gtf_file)[:-3]+"txt"

    if os.path.exists(annotation_file):
        command_site_annotation = "python " + tooldir + "/" + "site_annotation.py " + annotation_file + " " + site_result + " > " + site_result_anno
        subprocess.call(command_site_annotation, shell=True)
    else:
        command_get_annotation_file = "python " + tooldir + "/" + "get_gene_anno_file.py " +  gtf_file  + " > " + annotation_file
        subprocess.call(command_get_annotation_file, shell=True)
        command_site_annotation = "python " + tooldir + "/" + "site_annotation.py " + annotation_file + " " + site_result + " > " + site_result_anno
        subprocess.call(command_site_annotation, shell=True)
