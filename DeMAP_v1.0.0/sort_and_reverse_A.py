#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Zehong Wu wang_lab, Sun Yat-sen University"""
"""Feb, 2025"""
"""Email: wuzh93@mail2.sysu.edu.cn"""
"""Usage: This script reverts the temporarily substituted G bases back to A bases."""



import argparse

parser = argparse.ArgumentParser(description = "Run EM model")
#parser.add_argument("-unmap", "--unmap", nargs="?", type=str, default="none", help = "unmap reads")
parser.add_argument("-bed", "--bed", nargs="?", type=str, default="none", help = "bed file with information of A bases")
parser.add_argument("-pre", "--outname_prefix", nargs="?", type=str, default='default',help = "--outname_prefix")
parser.add_argument("-o", "--outputdir", nargs="?", type=str, default="none", help="outputdir")
parser.add_argument("-len", "--unmaplen", nargs="?", type=int, default=50, help="Discard reads if read length < unmaplen, the remains will be split to 2 reads equally for re-alignemt")

args = parser.parse_args()


#unmap_fq = args.unmap
bedfile = args.bed
outdic = args.outputdir
prefix = args.outname_prefix
minlen = args.unmaplen
unmap_fq = outdic + "/" + prefix + "_genome_both_ummap.fq"
outfile = outdic + "/" + prefix + "_realign.fq"

def load_unmap(filein):
    unmap_dic = {}
    infile = open(filein,"r")
    for line in infile:
        line = line.rstrip()
        if line[0] == "@":
            readid = line.split()[0]
            unmap_dic[readid]=[]
        else:
            unmap_dic[readid].append(line)
    infile.close()
    return unmap_dic

def return_A_bases(unmap_dic,bedfile):
    infile = open(bedfile)
    for line in infile:
        info = line.rstrip().split()
        readnum = info[0]
        Apos = info[1]
        if "@"+readnum not in unmap_dic:
            continue
        else:
            if info[1] == "NA":
                continue
            else:
                Alist = map(int, Apos.split("_"))
                oldstr = unmap_dic["@"+readnum][0]
                oldstr = list(oldstr)
                for index in Alist:
                    oldstr[index] = "A"
                newstr = ''.join(oldstr)
                unmap_dic["@"+readnum][0] = newstr

    infile.close()

def split_fq(unmap_dic,minlen):
    fileout = open(outfile,'w')
    for read_id in unmap_dic:
        sequence = unmap_dic[read_id][0]
        Quality = unmap_dic[read_id][2]
        if len(sequence) < minlen:
            continue
        else:
            mid_pos = len(sequence) // 2
            seq1, seq2 = sequence[:mid_pos], sequence[mid_pos:]
            qual1, qual2 = Quality[:mid_pos], Quality[mid_pos:]
            read_id1 = f"{read_id}_1"
            read_id2 = f"{read_id}_2"
            outstr = "\n".join([read_id1,seq1,"+",qual1,read_id2,seq2,"+",qual2])
            fileout.write(outstr + "\n")
    fileout.close()

if __name__ == "__main__":
    unmap_dic = load_unmap(unmap_fq)
    return_A_bases(unmap_dic,bedfile)
    split_fq(unmap_dic,minlen)
