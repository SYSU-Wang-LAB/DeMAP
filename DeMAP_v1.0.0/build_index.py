#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Zehong Wu wang_lab, Sun Yat-sen University"""
"""Feb, 2025"""
"""Email: wuzh93@mail2.sysu.edu.cn"""
"""Usage: This program is used to build index for DeMAP"""

import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
import re
import subprocess

parser = argparse.ArgumentParser(description = "building index for genome and transcriptome")
parser.add_argument("-g", "--genome", nargs="?", type=str, default='none', help="genome fasta")
parser.add_argument("-r", "--rRNA", nargs="?", type=str, default='none', help="rRNA fasta")
parser.add_argument("-p", "--Threads", nargs="?", type=str, default='1', help = "the number of threads")
parser.add_argument("-prefix", "--outfile_prefix", nargs="?", type=str, default='none',help = "outname prefix")
parser.add_argument("-o", "--outputdir", nargs="?", type=str, default='./', help="index directory")

args = parser.parse_args()

genome = args.genome
rRNA = args.rRNA
Threads = args.Threads
prefix = args.outfile_prefix
outputdir = args.outputdir

def converted_fasta_sequence(genome,rRNA,prefix,outputdir):

    if not os.path.exists(outputdir):
        os.makedirs(outputdir)

    rRNA_dir = os.path.join(outputdir, "rRNA")
    if not os.path.exists(rRNA_dir):
        os.makedirs(rRNA_dir)

    AG_genome_fasta = outputdir+"/" + prefix + ".AG_conversion.fa"
    TC_genome_fasta = outputdir+"/" + prefix + ".TC_conversion.fa"

    AG_rRNA_fasta = os.path.join(rRNA_dir, "rRNA.AG_conversion.fa")
    #AG_rRNA_fasta = outputdir+"/rRNA.AG_conversion.fa"

    AGfastafile = open(AG_genome_fasta,'w')
    for record in SeqIO.parse(genome, "fasta"):
        record.seq = Seq(re.sub('A','G',str(record.seq).upper()))
        SeqIO.write(record, AGfastafile, "fasta")
    AGfastafile.close()
    print("********** Genome A-to-G converted done ************")

    TCfastafile = open(TC_genome_fasta,'w')
    for record in SeqIO.parse(genome, "fasta"):
        record.seq = Seq(re.sub('T','C',str(record.seq).upper()))
        SeqIO.write(record, TCfastafile, "fasta")
    TCfastafile.close()
    print("********** Genome T-to-C converted done ************")

    AGrRNA_fastafile = open(AG_rRNA_fasta,'w')
    for record in SeqIO.parse(rRNA, "fasta"):
        record.seq = Seq(re.sub('A','G',str(record.seq).upper()))
        SeqIO.write(record, AGrRNA_fastafile, "fasta")
    AGrRNA_fastafile.close()
    print("********** rRNA and tRNA A-to-G converted done ************")

    return AG_genome_fasta,TC_genome_fasta,AG_rRNA_fasta

def build_STAR_index(AG_convert_fasta,TC_convert_fasta,Threads):
    AG_index_dir = AG_convert_fasta[:-3]
    TC_index_dir = TC_convert_fasta[:-3]

    command_AG = 'STAR --runMode genomeGenerate --runThreadN ' + str(Threads) + ' --genomeDir ' + AG_index_dir + ' --genomeFastaFiles ' + AG_convert_fasta + \
         ' --genomeSAindexNbases 14 --limitGenomeGenerateRAM 128000000000 '
    print(command_AG)
    subprocess.call(command_AG,shell=True)

    command_TC = 'STAR --runMode genomeGenerate --runThreadN ' + str(Threads) + ' --genomeDir ' + TC_index_dir + ' --genomeFastaFiles ' + TC_convert_fasta + \
         ' --genomeSAindexNbases 14 --limitGenomeGenerateRAM 128000000000 '
    print(command_TC)
    subprocess.call(command_TC,shell=True)


def build_bowtie2_index(fasta_in, Threads):

    command_rRNA = 'bowtie2-build --threads ' + str(Threads) + ' ' + fasta_in + ' ' + fasta_in
    print(command_rRNA)
    subprocess.call(command_rRNA, shell=True)

if __name__ == "__main__":
    genomename = os.path.basename(genome)
    if prefix == "none":
        prefix = genomename.split(".")[0]
        
    print("********** generate A-to-G and T-to-C converted fasta... ************")
    AG_convert_fasta,TC_convert_fasta,AG_convert_rRNA_fasta = converted_fasta_sequence(genome,rRNA,prefix,outputdir)

    print("********** Building genome index using STAR.... *************")
    build_STAR_index(AG_convert_fasta,TC_convert_fasta,Threads)

    print("********** Building rRNA index using bowtie2.... ************")
    build_bowtie2_index(AG_convert_rRNA_fasta, Threads)

    print("********** All indexes have been built ************")