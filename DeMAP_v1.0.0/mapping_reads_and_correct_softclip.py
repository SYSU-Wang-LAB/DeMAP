#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Zehong Wu Wang_lab, Sun Yat-sen University"""
"""Feb, 2025"""
"""Email: wuzh93@mail2.sysu.edu.cn"""
"""Usage: This program is used for mapping reads and correct soft-clipping"""

import sys
import os
import re
import argparse
import multiprocessing
import subprocess
import itertools
import time
from heapq import merge
import glob
from time import strftime
from Bio.Seq import reverse_complement
from Bio import SeqIO
import gzip

parser = argparse.ArgumentParser(description="Mapping reads of GLORI data")

# pre-alignment
group1 = parser.add_argument_group("Alignment Parameters")
group1.add_argument("-software", "--software", type=str, default="none", help="The path to the tool directory")
group1.add_argument("-genomedir", "--genomedir", nargs="?", type=str, default="none", help = "The path to the genome index directory")
group1.add_argument("-fa","--fasta",type=str,help="Unconverted genome file in fasta format")
group1.add_argument("-q", "--fastq", type=str, default="none", help="Fastq file for single end data,read should be sense strand.\n"
                    "For most widely used standard library: dUTP-based method: R2")
group1.add_argument("-p", "--Threads", type=str, default='1', help="The number of alignment threads. (defalut:1)")
group1.add_argument("-t", "--tools", type=str, default="STAR", help="Alignment tools: STAR,bowtie2. (defalut:STAR)")
group1.add_argument("-m", "--mismatch", type=int, default=2, help="Mapping mismatch. (defalut:2)")
group1.add_argument("-F", "--FilterN", type=str, default='0.4', help="MinOverLread. (defalut:0.4)")
group1.add_argument("-max", "--maxmap", type=str, default='5', help="Discard multi-mapping reads with more than maxmap mappings. (defalut:5)")
group1.add_argument("-pre", "--outname_prefix", type=str, default='default', help="Output name prefix")
group1.add_argument("-o", "--outputdir", type=str, default="none", help="Output directory")

# correct softclip due to C-to-U mutation
group2 = parser.add_argument_group("Correct softclip due to C-to-U mutation")
group2.add_argument("-gtf", "--gtf", nargs="?", type=str, default="none", help = "GTF file from gencode")
group2.add_argument("-SC_pval", "--SC_pval", nargs="?", type=str, default="0.001", help="P-value for the cutoff of LLR test. (default:0.001)")
group2.add_argument("-match_ratio", "--match_ratio", nargs="?", type=str, default="0.80", help="The minimum ratio of bases match to reference in softclip region. (default:0.80)")
group2.add_argument("-minSC", "--min_softclip", nargs="?", type=str, default='5', help = "The minimum softclip length included in correction. (default:5)")

# merge and split reads
group3 = parser.add_argument_group("Merge AG,TC alignemts and reclassify reads into uniq/multi-mapping reads")
group3.add_argument("-s1", "--compareScore", type=str, default='2', 
                    help="Discard alignment with AS (alignment score) < max(AS) - compareScore in a multi-mapping read\n"
                         "(2 in default: Remove alignments with AS at least 2 points lower than the best alignment of that read)")
group3.add_argument("-s2", "--minScore", type=str, default='22', 
                    help="Discard read if AS < minScore in all alignments of that read (defalut:22).\n"
                    "For STAR aligner,The AS equals matched bases minus penalty terms (mismatches, insertions, deletions, soft clipping, and splice junctions)")

# re-alignment of unmap reads
group4 = parser.add_argument_group("Re-alignment of unmap reads")
group4.add_argument("--re_align", action="store_true", help="If specified, split and re-align the unmapped reads")
group4.add_argument("-len", "--unmaplen", type=str, default='60', 
                    help="Discard reads if read length < unmaplen; the remaining reads will be split into 2 equal parts for re-alignment")


args = parser.parse_args()
tooldir = args.software
genomedir = args.genomedir
AGgenome = genomedir+"/"+[f for f in os.listdir(genomedir) if f.endswith('AG_conversion.fa')][0]
TCgenome = genomedir+"/"+[f for f in os.listdir(genomedir) if f.endswith('TC_conversion.fa')][0]
rRNAref = genomedir+"/rRNA/"+[f for f in os.listdir(genomedir+"/rRNA/") if f.endswith('rRNA.AG_conversion.fa')][0]
#rRNAref = genomedir+"/rRNA.AG_conversion.fa"
genome = args.fasta #genomedir+"/"+[f for f in os.listdir(genomedir) if f.endswith('genome.fa')][0]
fastq = args.fastq
Threads = args.Threads

tools = args.tools
FilterN = args.FilterN

mismatch = args.mismatch
maxmap = args.maxmap
outputdir = args.outputdir
outname_prx = args.outname_prefix
re_align = args.re_align

SC_pval = args.SC_pval
match_ratio = args.match_ratio
gtf_file = args.gtf
min_softclip = args.min_softclip

score1 = args.compareScore
score2 = args.minScore
unmaplen = args.unmaplen


re_digits = re.compile(r'(\d+)')

def embedded_numbers(s):
    s2=s.strip().split("\t")
    pieces = re_digits.split(s2[0])
    pieces[1::2] = map(int, pieces[1::2])
    return pieces


def sort_bedfiles(bedfiles,outputfiles):
    prx2 = bedfiles[:-4]
    path = prx2+"_chunk_*.bed"
    chunksize = 5000000
    fid = 1
    lines = []
    with open(bedfiles, 'r') as f_in:
        f_out = open(prx2+'_chunk_{}.bed'.format(fid), 'w')
        for line_num, line in enumerate(f_in, 1):
            lines.append(line)
            if not line_num % chunksize:
                lines = sorted(lines, key=embedded_numbers)  ###将readname 切割成字符和数值的列表，再排序
                f_out.writelines(lines)
                f_out.close()
                lines = []
                fid += 1
                f_out = open(prx2+'_chunk_{}.bed'.format(fid), 'w')
        # last chunk
        if lines:
            lines = sorted(lines, key=embedded_numbers)
            f_out.writelines(lines)
            f_out.close()
            lines = []

    chunks = []
    for filename in glob.glob(path):
        chunks += [open(filename, 'r')]

    with open(outputfiles, 'w') as f_out:
        print('merging bedfiles')
        f_out.writelines(merge(*chunks, key=embedded_numbers))
    subprocess.call("rm -f " + prx2 + "_chunk_*.bed", shell=True)

def change_reads(fastq, changename, output_bed, outputdir, rename, step=40000):

    def smart_open(filename, mode="rt"):
        if filename.endswith(".gz"):
            return gzip.open(filename, mode)
        else:
            return open(filename, mode)

    os.makedirs(outputdir, exist_ok=True)
    with smart_open(fastq) as file, \
         open(changename, 'w') as file_change, \
         open(output_bed, 'w') as file_bed:
        
        read_id = 1
        while True:
            fr = list(itertools.islice(file, step))
            if not fr:
                break
            batch_change = []
            batch_bed = []
            
            for i in range(0, len(fr), 4):
                header_line = fr[i].strip()
                seq_line = fr[i+1].strip().upper()
                plus_line = fr[i+2].strip()
                qual_line = fr[i+3].strip()
                
                if rename:
                    reads_name = f"@R.{read_id}"
                else:
                    reads_name = header_line.split(" ")[0]
                
                A_sites = [str(pos) for pos, char in enumerate(seq_line) if char == 'A']
                new_seq = seq_line.replace('A', 'G')
                
                batch_change.extend([reads_name, new_seq, plus_line, qual_line])
                
                if A_sites:
                    bed_line = f"{reads_name[1:]}\t{'_'.join(A_sites)}\n"
                else:
                    bed_line = f"{reads_name[1:]}\tNA\n"
                batch_bed.append(bed_line)
                read_id += 1
            
            file_change.write("\n".join(batch_change) + "\n")
            file_bed.writelines(batch_bed)

    sort_bedfiles(output_bed, output_bed + "_sorted")
    subprocess.call("rm -f " + output_bed, shell=True)


def mapping_files(tool,fastq,reference,Threads,muta_N,fqname,outputdir,maxmap,flag):
    tmpfile = outputdir +fqname+"tmp.sam"
    outputfile = outputdir +fqname+".sam"  ## test.sam
    unmapfastq = outputdir +fqname+"_unmap.fq"   ### test_un_2.fq
    strand_filter = "-F 16 " if flag == "AG" else "-f 16 "

    # Keep all possible alignments in the pre-alignment stage to enable downstream SAM merging and subsequent score-based filtering
    out_maxmap = "10" if int(maxmap) <= 5 else "20"


    if tool == "bowtie2":  
        para_0 = 'bowtie2 -p ' + str(Threads)
        para_A = ' -x ' + reference + ' -U ' + fastq
        para_B = ' --no-unal -S ' + tmpfile
        para_unmap = ' --un ' + unmapfastq
        para_end = ' 2>' + outputfile + '.output'

        command = para_0 + para_A + para_B + para_unmap + para_end
        print(command)
        subprocess.call(command, shell=True)
        
        print("samtools view -F " + flag + " -@ " + Threads + " -h " + tmpfile +" | samtools sort -@ " + Threads + " -n -O SAM > " + outputfile)
        subprocess.call("samtools view -F " + flag + " -@ " + Threads + " -h " + tmpfile + " | samtools sort -@ " + Threads + " -n -O SAM > " + outputfile, shell=True)
        subprocess.call("rm -f " + tmpfile, shell=True)


    elif tool == "STAR":
        para_0 = "STAR --runThreadN "+ Threads
        para_g = " --genomeDir "+ reference[:-3]
        para_A = " --limitOutSJcollapsed 5000000 --outFilterMultimapScoreRange 4 " #--outFilterScoreMin "  + score2
        para_B = " --outFilterMismatchNmax " + str(muta_N)  ### mismatch number
        para_B_1=' --outFilterScoreMinOverLread '+FilterN+' --outFilterMatchNminOverLread '+FilterN+' --seedSearchStartLmax 30 '# mean mapping score/mean mapping base ratio/increase overall mapping sensitivity
        para_C = " --outSAMattributes All --outSAMprimaryFlag AllBestScore --outMultimapperOrder Random --outSAMmultNmax " + out_maxmap + " --outSAMtype BAM Unsorted"
        para_D = " --outFilterMultimapNmax 100 "
        para_E = " --outFileNamePrefix " + outputfile[:-3] + " --readFilesIn " + fastq
        para_unmap = " --outReadsUnmapped Fastx"
        line_command = para_0+para_g+para_A+para_B+para_B_1+para_C+para_D+para_E + para_unmap
        print(line_command)
        subprocess.call(line_command, shell=True)
        print("samtools view " + strand_filter + " -@ " + Threads+" -h " + outputfile[:-3] + "Aligned.out.bam | samtools sort -@ " + Threads + " -n -O SAM > " + outputfile)
        subprocess.call("samtools view " + strand_filter + " -@ " + Threads+" -h " + outputfile[:-3] + "Aligned.out.bam | samtools sort -@ " + Threads + " -n -O SAM > " + outputfile, shell=True) #test.sam
        subprocess.call("mv " + outputfile[:-3] + 'Unmapped.out.mate1 ' + unmapfastq, shell=True) #test_un_2.fq
        subprocess.call("rm -f " + outputfile[:-3] + 'Aligned.out.bam', shell=True)
    return outputfile,unmapfastq


def return_A_to_samfile(samfile,output_bed,base_to_change): # could be used to multi-mapping reads
    reverse_sam = samfile[:-4] + "_r.sam"
    with open(reverse_sam, 'w') as OUTPUT, open(samfile, 'r') as file1, open(output_bed, 'r') as file2:
        sam_line = file1.readline().rstrip()
        read_line = file2.readline().rstrip() 
        read_id = read_line.split()[0] 
        read_A_pos = read_line.split()[1] 

        while sam_line and sam_line[0] == "@":
            OUTPUT.write(sam_line+"\n")
            sam_line = file1.readline().rstrip()

        while sam_line:
            sam_info = sam_line.split()
            sam_readid = sam_info[0]
            sam_readseq = sam_info[9]
            if sam_readid == read_id:
                if read_A_pos == "NA":  
                    OUTPUT.write(sam_line+"\n")
                    sam_line = file1.readline().rstrip()
                else:
                    Alist = map(int, read_A_pos.split("_"))
                    if base_to_change == "TC":
                        old_seq = reverse_complement(sam_readseq)
                        old_seq_list = list(old_seq)
                        for index in Alist:
                            old_seq_list[index] = "A"
                        newseq = ''.join(old_seq_list)
                        newseq = reverse_complement(newseq)
                    else:
                        old_seq = sam_readseq
                        old_seq_list = list(old_seq)
                        for index in Alist:
                            old_seq_list[index] = "A"
                        newseq = ''.join(old_seq_list)
                    sam_info[9] = newseq

                    OUTPUT.write("\t".join(sam_info)+"\n")
                    sam_line = file1.readline().rstrip()
            else:
                read_line = file2.readline().rstrip() # 读取bed 下一行
                read_id = read_line.split()[0]
                read_A_pos = read_line.split()[1]
    subprocess.call("rm -f " + samfile, shell=True)
    return reverse_sam

def getbamfiles(outputfile,fac,Threads,flag):
    output_bam = outputfile[:-4] + fac
    print("samtools view -F " + flag + " -bS -@ " + Threads+" -h "+outputfile + \
                    " | samtools sort > " + output_bam)
    subprocess.call("samtools view -F " + flag + " -bS -@ " + Threads+" -h "+outputfile + \
                    " | samtools sort > " + output_bam,shell=True)
    subprocess.call("samtools index " + output_bam, shell=True)
    subprocess.call("rm -f " + outputfile, shell=True)
    return output_bam


def get_both_unmap(unmapAG,unmapTC, output_file):
    sequences_id = set()
    with open(output_file, 'w') as out_file:
        for record in SeqIO.parse(unmapAG, 'fastq'):
            sequences_id.add(record.id)
        for record in SeqIO.parse(unmapTC, 'fastq'):
            if record.id in sequences_id:
                SeqIO.write(record, out_file, 'fastq')

                

if __name__ == "__main__":
    global step
    step = 10000
    if outname_prx == 'default':
        outname_prx = "_".join(os.path.basename(fastq).split(".")[:-1])

    if os.path.exists(outputdir):
        pass
    else:
        os.makedirs(outputdir)

    outputdir2 = outputdir+"/"
    final_out_prx = outputdir2 + "/" + outname_prx
    changefastq = final_out_prx + "_" + "AGchanged_2.fq"  
    output_bed = final_out_prx + "_A.bed"   
    if os.path.isfile(changefastq) or os.path.isfile(changefastq+'.gz'):
        sys.stderr.write("[%s] Warning：the changed files already exists, skip this precess\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
        pass
    else:
        sys.stderr.write("[%s] change to A>G...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
        start_time = time.time()
        change_reads(fastq, changefastq,output_bed,outputdir2,True)
        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"Elapsed time: {elapsed_time:.2f} seconds")




    ############ **********  remove reads mapping to rRNA  ********** ############
    sys.stderr.write("[%s] map to rRNA...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    rRNA_unmap_fq = final_out_prx + "_" + "rRNA_unmap.fq"
    if os.path.isfile(rRNA_unmap_fq) and os.path.getsize(rRNA_unmap_fq) > 100000:
        sys.stderr.write("[%s] Warning：the clean fastq file already exists, skip this precess\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
        pass
    else:        
        sam_rRNA, rRNA_unmap_fq, = mapping_files('bowtie2', changefastq, rRNAref, Threads, mismatch,outname_prx + "_rRNA", outputdir2, maxmap,'20')
        reverse_sam_rRNA = return_A_to_samfile(sam_rRNA, output_bed+"_sorted", "AG")
        rRNA_bam = getbamfiles(reverse_sam_rRNA, 's.bam', Threads,'20')   # test_rRNA_rs.bam
    



    ############ **********  mapping to AG/TC converted genome  ********** ############
    merge_ag_tc = final_out_prx +  "_" + "merge_AG_TC.sam"

    unmapfastqAG = final_out_prx + "_AG_unmap.fq"
    reverse_samAG = final_out_prx + "_AG_r.sam"
    unmapfastqTC = final_out_prx + "_TC_unmap.fq"
    reverse_samTC = final_out_prx + "_TC_r.sam"

    if os.path.isfile(merge_ag_tc) and os.path.getsize(merge_ag_tc) > 100000:
        sys.stderr.write("[%s] Warning：the merge sam file already exists, skip the process of mapping to A-G/T-C converted genome\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))

    else:
        def run_AG():
            sys.stderr.write("[%s] map to A-to-G converted genome...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
            samAG, _ = mapping_files(tools, rRNA_unmap_fq, AGgenome, Threads, mismatch, outname_prx + "_AG", outputdir2, maxmap, "AG")
            return_A_to_samfile(samAG, output_bed+"_sorted", "AG")

        def run_TC():
            sys.stderr.write("[%s] map to T-to-C converted genome...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
            samTC, _ = mapping_files(tools, rRNA_unmap_fq, TCgenome, Threads, mismatch, outname_prx + "_TC", outputdir2, maxmap, "TC")
            return_A_to_samfile(samTC, output_bed+"_sorted", "TC")

        p1 = multiprocessing.Process(target=run_AG)
        p2 = multiprocessing.Process(target=run_TC)

        p1.start()
        p2.start()

        p1.join() 
        p2.join()

        ############ **********  merge reads mapping to A-to-G converted genome and reads mapping to T-to-C converted genome  ********** ############ 
        sys.stderr.write("[%s] merge AG and TC alignments...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
        command_pre_merge = "python " + tooldir + "/pre_merge_sam.py"+ " -AG " + reverse_samAG  + " -TC " + reverse_samTC  + " -genome " + genome + \
        " -pre " + outname_prx  + " -o " + outputdir2 + " -m " + score2
        subprocess.call(command_pre_merge, shell=True)
        subprocess.call("rm -rf " + reverse_samAG + " " + reverse_samTC, shell=True)

    ############ **********  correct sotfclip due to C-to-U mutation  ********** ############    
    sys.stderr.write("[%s] Start correct sotfclip due to C-to-U mutation \n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    merge_sam = final_out_prx + "_merge_AG_TC.sam"
    correct_merge = final_out_prx + "_merge_AG_TC_correct.sam"
    command_correct1 = "python " + tooldir + "/correct_softclip_and_calculate_AS.py" + " -insam " + merge_sam  + \
      " -genome " + genome + " --Threads " + Threads +  " -outname " + correct_merge + " -gtf " + gtf_file + " -SC_pval " + SC_pval + " -match_ratio " + match_ratio + " -minSC " + min_softclip
    print(command_correct1)
    subprocess.call(command_correct1, shell=True)
    subprocess.call("rm  " + merge_sam, shell=True)

    ############ **********  split reads to multi or uniq mapping ********** ############   
    sys.stderr.write("[%s] split merge reads it uniq-mapping and multi-mapping reads...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    command_merge_and_split = "python " + tooldir + "/split_reads_to_uniq_or_multi_map.py" +  \
      " -pre " + outname_prx +  " -o " + outputdir2 + " -s1 " + score1 + " -s2 " + score2 + " -max " + maxmap
    print(command_merge_and_split)
    subprocess.call(command_merge_and_split, shell=True)
    uniqmap = final_out_prx + "_" + "uniq_map.sam"
    multimap = final_out_prx + "_" + "multi_map.sam"

    subprocess.call("rm  " + correct_merge, shell=True)

    if re_align:
        mismatch = mismatch if mismatch <= 2 else 2
        total_both_unmap = final_out_prx + "_" + "genome_both_ummap.fq"
        get_both_unmap(unmapfastqAG,unmapfastqTC,total_both_unmap)
        command_split_unmap = "python " + tooldir + "/sort_and_reverse_A.py" + " -bed " + output_bed+"_sorted" + " -pre " + outname_prx  + " -o " + outputdir2 + \
        " --unmaplen " + unmaplen
        print(command_split_unmap)
        subprocess.call(command_split_unmap, shell=True)

        unmap_fq = final_out_prx + "_realign.fq"
        unmap_fq_AGchanged = final_out_prx + "_realign_" + "AGchanged_2.fq"
        output_bed_unmap = final_out_prx + "_realign_A.bed"
        change_reads(unmap_fq, unmap_fq_AGchanged,output_bed_unmap,outputdir2,False)

        sys.stderr.write("[%s] unmap reads map to rRNA...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
        sam_rRNA, rRNA_unmap_fq, = mapping_files('bowtie2', unmap_fq_AGchanged, rRNAref, Threads, mismatch,outname_prx + "_realign_rRNA", outputdir2, maxmap,'20')


        reverse_samAG_realign = outputdir2 + outname_prx + "_realign_AG_r.sam"
        reverse_samTC_realign = outputdir2 + outname_prx + "_realign_TC_r.sam"


        def run_realign_AG():
            sys.stderr.write("[%s] unmap reads map to A-to-G converted genome...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
            samAG, _ = mapping_files(tools, rRNA_unmap_fq, AGgenome, Threads, mismatch,outname_prx + "_realign_AG", outputdir2,maxmap,'AG')
            return_A_to_samfile(samAG, output_bed_unmap+"_sorted", "AG")


        def run_realign_TC():
            sys.stderr.write("[%s] unmap reads map to T-to-C converted genome...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
            samTC, _ = mapping_files(tools, rRNA_unmap_fq, TCgenome, Threads, mismatch,outname_prx + "_realign_TC", outputdir2,maxmap,'TC')
            return_A_to_samfile(samTC, output_bed_unmap+"_sorted", "TC")


        p3 = multiprocessing.Process(target=run_realign_AG)
        p4 = multiprocessing.Process(target=run_realign_TC)

        p3.start()
        p4.start()

        p3.join()
        p4.join()


        # merge 
        sys.stderr.write("[%s] merge AG and TC alignments...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
        command_pre_merge = "python " + tooldir + "/pre_merge_sam.py"+ " -AG " + reverse_samAG_realign  + " -TC " + reverse_samTC_realign  + " -genome " + genome + \
        " -pre " + outname_prx + "_realign"  + " -o " + outputdir2 + " -m " + score2
        subprocess.call(command_pre_merge, shell=True)


        # correct
        merge_sam = final_out_prx + "_realign_merge_AG_TC.sam"
        correct_merge = final_out_prx + "_realign_merge_AG_TC_correct.sam"
        command_correct1 = "python " + tooldir + "/correct_softclip_and_calculate_AS.py" + " -insam " + merge_sam  + \
        " -genome " + genome + " --Threads " + Threads +  " -outname " + correct_merge + " -gtf " + gtf_file + " -SC_pval " + SC_pval + " -match_ratio " + match_ratio + " -minSC " + min_softclip
        print(command_correct1)
        subprocess.call(command_correct1, shell=True)



        # split
        score2 = str(int(score2) + 3) if int(score2) >= 24 else '27'
        command_merge_and_split = "python " + tooldir + "/split_reads_to_uniq_or_multi_map.py" +  \
        " -pre " + outname_prx  + "_realign" + " -o " + outputdir2 + " -s1 " + score1 + " -s2 " + score2 + " -max " + maxmap
        print(command_merge_and_split)
        subprocess.call(command_merge_and_split, shell=True)


        uniqmap_realign = final_out_prx + "_realign_" + "uniq_map.sam"
        multimap_realign = final_out_prx + "_realign_" + "multi_map.sam"

        subprocess.call("grep -v '^@' " + uniqmap_realign + " >> " + uniqmap, shell=True)
        subprocess.call("grep -v '^@' " + multimap_realign + " >> " + multimap, shell=True)
        subprocess.call("rm -rf " + " ".join([total_both_unmap,unmap_fq,unmap_fq_AGchanged,sam_rRNA,unmapfastqAG,unmapfastqTC,reverse_samAG_realign,reverse_samTC_realign]), shell=True)
        subprocess.call("rm -rf " + final_out_prx + "_realign* " , shell=True)

    ### mapping information
    subprocess.call("mkdir -p "+outputdir2+"/mapping-info_"+outname_prx, shell=True)
    subprocess.call("mv "+final_out_prx+"_AG*out "+outputdir2+"/mapping-info_"+outname_prx, shell=True)
    #subprocess.call("mv "+final_out_prx+"*put "+outputdir2+"/mapping-info_"+outname_prx, shell=True)
    subprocess.call("mv "+final_out_prx+"_AG*tab "+outputdir2+"/mapping-info_"+outname_prx, shell=True)
    subprocess.call("mv "+final_out_prx+"_TC*out "+outputdir2+"/mapping-info_"+outname_prx, shell=True)
    subprocess.call("mv "+final_out_prx+"_TC*tab "+outputdir2+"/mapping-info_"+outname_prx, shell=True)

    ### remove tmp file
    subprocess.call("rm -rf " + final_out_prx + "_AG*_STARtmp", shell=True)
    subprocess.call("rm -rf " + final_out_prx + "_TC*_STARtmp", shell=True)
    subprocess.call("rm -rf " + final_out_prx+"_AG_unmap.fq " + final_out_prx+ "_TC_unmap.fq " , shell=True)  #+ final_out_prx + "_merge_AG_TC_correct.sam "
    subprocess.call("rm -rf " + final_out_prx+"*_r.sam " + final_out_prx+"*_AG_unmap.fq " + final_out_prx+ "*_TC_unmap.fq ", shell=True)

    sys.stderr.write("[%s]Mapping reads done !!!...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
