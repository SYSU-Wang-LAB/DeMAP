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
parser.add_argument("-r", "--rRNA", nargs="?", type=str, default='none', help="rRNA tasta")
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
    if os.path.exists(outputdir):
        pass
    else:
        os.mkdir(outputdir)
    AG_genome_fasta = outputdir+"/" + prefix + ".AG_conversion.fa"
    TC_genome_fasta = outputdir+"/" + prefix + ".TC_conversion.fa"
    AG_rRNA_fasta = outputdir+"/rRNA.AG_conversion.fa"

    AGfastafile = open(AG_genome_fasta,'a')
    for record in SeqIO.parse(genome, "fasta"):
        record.seq = Seq(re.sub('A','G',str(record.seq).upper()))
        SeqIO.write(record, AGfastafile, "fasta")
    print("**********Genome AG converted done ************")

    TCfastafile = open(TC_genome_fasta,'a')
    for record in SeqIO.parse(genome, "fasta"):
        record.seq = Seq(re.sub('T','C',str(record.seq).upper()))
        SeqIO.write(record, TCfastafile, "fasta")
    print("**********Genome TC converted done ************")

    AGrRNA_fastafile = open(AG_rRNA_fasta,'a')
    for record in SeqIO.parse(rRNA, "fasta"):
        record.seq = Seq(re.sub('A','G',str(record.seq).upper()))
        SeqIO.write(record, AGrRNA_fastafile, "fasta")
    print("**********rRNA and tRNA AG converted done ************")

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


def build_bowtie_index(fasta_in):

    command_trancript_TC = 'bowtie-build ' + fasta_in + ' ' +fasta_in
    print(command_trancript_TC)
    subprocess.call('bowtie-build ' + fasta_in + ' ' + fasta_in,shell=True)


if __name__ == "__main__":
    genomename = os.path.basename(genome)
    if prefix == "none":
        prefix = genomename.split(".")[0]
    #subprocess.call("cp " + genome + " " + outputdir+"/" + prefix + ".genome.fa" ,shell=True)
    print("********** generate A-to-G and T-to-C converted fasta... ************")
    AG_convert_fasta,TC_convert_fasta,AG_convert_rRNA_fasta = converted_fasta_sequence(genome,rRNA,prefix,outputdir)

    print("********** Building genome index using STAR.... *************")
    build_STAR_index(AG_convert_fasta,TC_convert_fasta,Threads)

    print("**********Building rRNA index using bowtie.... ************")
    build_bowtie_index(AG_convert_rRNA_fasta)

    print("********** All indexes have been built ************")
