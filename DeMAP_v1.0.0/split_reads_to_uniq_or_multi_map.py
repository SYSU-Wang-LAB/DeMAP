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
parser.add_argument("-s2", "--minScore", nargs="?", type=int, default=24, help = "Discard read if AS < minScore in all alignments of that read. (defalut:24)")

args = parser.parse_args()


outdir = args.outputdir
out_pre = args.outname_prefix
Threads = args.Threads
maxmap = args.maxmap # if args.maxmap > 1 else 2
score1 = args.score1
minScore = args.minScore


min_overlap = 15

def get_alignment_score(alignment):
    match = re.search(r'(?:^|\t)AS:i:(-?\d+)(?:\t|$)', alignment)
    if match is None:
        raise ValueError("AS tag was not found:\n" + alignment.rstrip("\n"))
    return int(match.group(1))


def get_extended_reference_interval(alignment):
    """
    Return:
        chromosome, strand, extended_start, extended_end

    Reference-consuming CIGAR operations:
        M, D, N, =, X

    Query-only insertion I does not consume the reference.

    Terminal S operations are used to extend both sides of the
    genomic interval. Terminal H operations are ignored.
    """
    fields = alignment.rstrip("\n").split("\t")
    flag = int(fields[1])
    chrom = fields[2]
    pos = int(fields[3])
    cigar = fields[5]
    strand = flag & 16

    cigar_operations = re.findall(r'(\d+)([MIDNSHP=X])', cigar)

    if not cigar_operations or ''.join(length + operation for length, operation in cigar_operations) != cigar:
        raise ValueError("Invalid or unsupported CIGAR:\n" + alignment.rstrip("\n"))

    left_softclip = 0
    right_softclip = 0

    # Skip terminal hard clipping and read the left terminal soft clip.
    for length, operation in cigar_operations:
        if operation == 'H':
            continue
        if operation == 'S':
            left_softclip = int(length)
        break

    # Skip terminal hard clipping and read the right terminal soft clip.
    for length, operation in reversed(cigar_operations):
        if operation == 'H':
            continue
        if operation == 'S':
            right_softclip = int(length)
        break

    reference_length = sum(int(length) for length, operation in cigar_operations if operation in 'MDN=X')

    extended_start = pos - left_softclip
    extended_end = pos + reference_length - 1 + right_softclip

    return chrom, strand, extended_start, extended_end



def merge_overlapping_alignments(alignments, minimum_overlap=15):
    """
    For alignments from the same read:

    1. Calculate the reference interval after extending both terminal
       soft-clipped regions.
    2. If two intervals are on the same chromosome and strand and
       overlap by at least minimum_overlap bases, treat them as the
       same local mapping.
    3. Keep the alignment with the highest AS.
    4. If AS is tied, preferentially keep the primary alignment.

    Alignments are processed from highest to lowest priority. This
    avoids retaining a lower-scoring alignment before a better
    overlapping alignment is examined.
    """
    ranked_alignments = []

    for original_index, alignment in enumerate(alignments):
        fields = alignment.rstrip("\n").split("\t")
        flag = int(fields[1])
        score = get_alignment_score(alignment)
        is_secondary = bool(flag & 256)
        chrom, strand, start, end = get_extended_reference_interval(alignment)

        # Higher AS first; if tied, primary first; then original order.
        ranked_alignments.append((-score, is_secondary, original_index, chrom, strand, start, end, alignment))

    ranked_alignments.sort(key=lambda item: (item[0], item[1], item[2]))

    retained_alignments = []
    retained_intervals = []

    for negative_score, is_secondary, original_index, chrom, strand, start, end, alignment in ranked_alignments:
        merge_with_existing = False

        for retained_chrom, retained_strand, retained_start, retained_end in retained_intervals:
            if chrom != retained_chrom or strand != retained_strand:
                continue

            overlap_length = min(end, retained_end) - max(start, retained_start) + 1

            if overlap_length >= minimum_overlap:
                merge_with_existing = True
                break

        if not merge_with_existing:
            retained_alignments.append((original_index, alignment))
            retained_intervals.append((chrom, strand, start, end))

    # Restore the relative order from the original SAM file.
    retained_alignments.sort(key=lambda item: item[0])

    return [alignment for original_index, alignment in retained_alignments]



def mark_alignment_as_unique(alignment):
    fields = alignment.rstrip("\n").split("\t")

    # Clear the secondary-alignment bit 0x100.
    fields[1] = str(int(fields[1]) & ~256)

    alignment = "\t".join(fields) + "\n"
    alignment = re.sub(r'NH:i:\d+', 'NH:i:1', alignment)
    alignment = re.sub(r'HI:i:\d+', 'HI:i:1', alignment)

    return alignment


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
                        mapping_score = [get_alignment_score(alignment) for alignment in record_dic[previous_read_id]]

                        if max(mapping_score) < minScore:
                            record_dic.clear()
                            record_dic[read_id] = [line]
                            continue

                        # filter alignments with extremely low mapping score compare to best alignment
                        filtered_alignments = [alignment for AS, alignment in zip(mapping_score, record_dic[previous_read_id]) if max(mapping_score) - AS <= score1]

                        filtered_alignments = merge_overlapping_alignments(filtered_alignments, min_overlap)

                        if len(filtered_alignments) > maxmap:
                            record_dic.clear()
                            record_dic[read_id] = [line]
                            continue

                        if len(filtered_alignments) < 2:
                            tmpstr = mark_alignment_as_unique(filtered_alignments[0])

                            mapping_score = get_alignment_score(tmpstr)

                            if mapping_score >= minScore:
                                uniq_out.write(tmpstr)

                            record_dic.clear()
                            record_dic[read_id] = [line]
                            continue

                        else:
                            mapping_num = len(filtered_alignments)

                            for multi_num_index in range(0,mapping_num):
                                tmpstr = re.sub(r'NH:i:(\d+)', f'NH:i:{mapping_num}', filtered_alignments[multi_num_index])
                                tmpstr1 = re.sub(r'HI:i:(\d+)', f'HI:i:{multi_num_index + 1}', tmpstr)
                                multi_out.write(tmpstr1)

                            record_dic.clear()
                            record_dic[read_id] = [line]

                    else:
                        mapping_score = get_alignment_score(record_dic[previous_read_id][0])

                        if mapping_score >= minScore:
                            uniq_out.write(record_dic[previous_read_id][0])

                        record_dic.clear()
                        record_dic[read_id] = [line]

        # output last read
        if len(record_dic) == 0:
            return

        last_read_id = list(record_dic.keys())[0]

        if len(record_dic[last_read_id]) >= 2: # multi-alignment
            mapping_score = [get_alignment_score(alignment) for alignment in record_dic[last_read_id]]

            if max(mapping_score) >= minScore:
                # filter alignments with extremely low mapping score compare to best alignment
                filtered_alignments = [alignment for AS, alignment in zip(mapping_score, record_dic[last_read_id]) if max(mapping_score) - AS <= score1]

                filtered_alignments = merge_overlapping_alignments(filtered_alignments, min_overlap)

                if len(filtered_alignments) < 2:
                    tmpstr = mark_alignment_as_unique(filtered_alignments[0])

                    mapping_score = get_alignment_score(tmpstr)

                    if mapping_score >= minScore:
                        uniq_out.write(tmpstr)

                else:
                    if len(filtered_alignments) <= maxmap:
                        mapping_num = len(filtered_alignments)

                        for multi_num_index in range(0,mapping_num):
                            tmpstr = re.sub(r'NH:i:(\d+)', f'NH:i:{mapping_num}', filtered_alignments[multi_num_index])
                            tmpstr1 = re.sub(r'HI:i:(\d+)', f'HI:i:{multi_num_index + 1}', tmpstr)
                            multi_out.write(tmpstr1)

        else:
            mapping_score = get_alignment_score(record_dic[last_read_id][0])

            if mapping_score >= minScore:
                uniq_out.write(record_dic[last_read_id][0])


if __name__ == "__main__":
    ###### merge AG and TC alinement reads
    merge_sam = outdir + "/" + out_pre + "_merge_AG_TC_correct.sam"
    ##### split merge sam to uniq or multi alignment reads
    uniq_sam = outdir + "/" + out_pre + "_uniq_map.sam"
    multi_sam = outdir + "/" + out_pre + "_multi_map.sam"
    split_merge_sam(merge_sam,uniq_sam,multi_sam)