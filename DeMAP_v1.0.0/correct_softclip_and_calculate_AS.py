#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Zehong Wu Wang_lab, Sun Yat-sen University"""
"""Feb, 2025"""
"""Email: wuzh93@mail2.sysu.edu.cn"""
"""Usage: This program is used to correct soft-clipping."""

import argparse
import pysam
import numpy as np
from collections import defaultdict
import os
import multiprocessing as mp
import copy
import time
import re

parser = argparse.ArgumentParser(description = "Correct the softclip due to high C to U mutataion",formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-insam", "--insam", nargs="?", type=str, default="none", help = "samfile of read mapping to AG/TC converted genome")
parser.add_argument("-outname", "--outname", nargs="?", type=str, default="none", help = "output file name")
parser.add_argument("-genome", "--genome", nargs="?", type=str, default="none", help = "unconverted genome")
parser.add_argument("-gtf", "--gtf", nargs="?", type=str, default="none", help = "GTF file from gencode")
parser.add_argument("-p", "--Threads", nargs="?", type=str, default='1', help = "The number of threads")
parser.add_argument("-bz", "--batchsize", nargs="?", type=int, default=100000, help = "The number of reads processed in a batch. (default:100000)")
parser.add_argument("--c2u_rate_threshold", type=float, default=0.05,help="Background C-to-U rate threshold for half-weighted mismatches (default: 0.05).")
parser.add_argument("-SC_pval", "--SC_pval", nargs="?", type=float, default=0.001, help = "P-value for the cutoff of LLR test. (default:0.001)")
parser.add_argument("-match_ratio", "--match_ratio", nargs="?", type=float, default=0.7, help = "The minimum ratio of bases match to reference in softclip region. (default:0.75)")
parser.add_argument("-minSC", "--min_softclip", nargs="?", type=int, default=5, help = "The minimum softclip length included in correction. (default:5)")
parser.add_argument("--bootstrap_n", type=int, default=100000, help="H0 simulations per reference condition. (default:100000)")
parser.add_argument("--bam", action="store_true", help="If specified, output file as bam.")

args = parser.parse_args()


insam = args.insam
outname = args.outname
genomefile = args.genome
gtffile = args.gtf
Threads = int(args.Threads)
pvalue = args.SC_pval
match_ratio = args.match_ratio
min_softclip = args.min_softclip
chunk_size = args.batchsize
C2Ulevel = args.c2u_rate_threshold
outbam = args.bam

BASE_TO_IDX = {'A': 0, 'T': 1, 'C': 2, 'G': 3}
EPSILON = 1e-5
BOOTSTRAP_N_SIM = args.bootstrap_n  # Set --bootstrap_n 10000 to reduce simulation count.
BOOTSTRAP_BATCH_SIZE = 2000  # Simulate in small batches to limit temporary memory.
BOOTSTRAP_CACHE_SIZE = 1000  # Maximum cached reference conditions per chunk and strand.

if not 1 / (BOOTSTRAP_N_SIM + 1) <= pvalue < 1:
    parser.error(f"--SC_pval must be >= 1/({BOOTSTRAP_N_SIM}+1) and < 1")

def load_exons_and_splice_sites(gtf_file):
    """
    extract exon annotations and splice junctions of multi-exon transcripts from the GTF file
    return: (exons_by_transcript, splice_sites)
    """
    exons_by_transcript = {}
    splice_sites = {}
    # only keep multi-exon transcripts
    multi_exon_transcripts = {}

    if gtf_file == "none":
        return multi_exon_transcripts, splice_sites


    tid_pattern = re.compile(r'transcript_id[\s\t]+(?:\"([^\"]+)\"|([^;\s\t]+))')

    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t', 8)

            if len(fields) < 9 or fields[2] != 'exon':
                continue

            chrom = fields[0]
            start = int(fields[3]) - 1  # 1-based to 0-based
            end = int(fields[4]) - 1    # 1-based to 0-based
            strand = fields[6]
            info = fields[8]

            match = tid_pattern.search(info)
            if not match:
                continue

            transcript_id = match.group(1) if match.group(1) else match.group(2)

            if transcript_id not in exons_by_transcript:
                exons_by_transcript[transcript_id] = []
            exons_by_transcript[transcript_id].append((chrom, start, end, strand))


    for tid, exons in exons_by_transcript.items():
        if len(exons) > 1:
            exons_sorted = sorted(exons, key=lambda x: x[1])
            multi_exon_transcripts[tid] = exons_sorted

            chrom = exons_sorted[0][0]
            if chrom not in splice_sites:
                splice_sites[chrom] = {}
            strand = exons_sorted[0][3]

            for i in range(len(exons_sorted) - 1):
                # donor site: exon end + 1
                donor_pos = exons_sorted[i][2] + 1
                if donor_pos not in splice_sites[chrom]:
                    splice_sites[chrom][donor_pos] = []
                splice_sites[chrom][donor_pos].append((strand, tid, 'donor'))

                # acceptor site: next exon start
                acceptor_pos = exons_sorted[i+1][1]  # 0-based
                if acceptor_pos not in splice_sites[chrom]:
                    splice_sites[chrom][acceptor_pos] = []
                splice_sites[chrom][acceptor_pos].append((strand, tid, 'acceptor'))

    return multi_exon_transcripts, splice_sites


def compute_softclip_base_freq(sam_path, n_reads=5000000):
    """
    Calculate strand-specific nucleotide frequencies in soft-clipped regions
    """
    sc_pos_count = defaultdict(int)
    sc_neg_count = defaultdict(int)
    sc_pos_total = 0
    sc_neg_total = 0

    with pysam.AlignmentFile(sam_path, "r") as samfile:
        for i, read in enumerate(samfile):
            if i >= n_reads:
                break

            cigar = read.cigartuples
            read_seq = read.query_sequence
            read_pos = 0

            for op, length in cigar:
                if read_pos >= len(read_seq):
                    break

                if op == 4:  # Softclip
                    end_pos = min(read_pos + length, len(read_seq))
                    sc_frag = read_seq[read_pos:end_pos]

                    for base in sc_frag:
                        if base in 'ATCG':
                            if read.is_reverse:
                                sc_neg_count[base] += 1
                                sc_neg_total += 1
                            else:
                                sc_pos_count[base] += 1
                                sc_pos_total += 1

                    read_pos = end_pos

                elif op in [0, 1, 7, 8]:  # M/I/=/X
                    read_pos = min(read_pos + length, len(read_seq))

    sc_pos_theta = {b: sc_pos_count[b] / sc_pos_total if sc_pos_count[b] > 0 else EPSILON for b in 'ATCG'}
    sc_neg_theta = {b: sc_neg_count[b] / sc_neg_total if sc_neg_count[b] > 0 else EPSILON for b in 'ATCG'}

    sc_pos_theta = {k: v / sum(sc_pos_theta.values()) for k, v in sc_pos_theta.items()}
    sc_neg_theta = {k: v / sum(sc_neg_theta.values()) for k, v in sc_neg_theta.items()}
    # Precompute log frequencies.
    sc_pos_log_theta = {b: np.log(sc_pos_theta[b]) for b in 'ATCG'}
    sc_neg_log_theta = {b: np.log(sc_neg_theta[b]) for b in 'ATCG'}

    return sc_pos_theta, sc_neg_theta, sc_pos_log_theta, sc_neg_log_theta


def fetch_aligned_reference(read, ref_fasta):
    """Fetch reference bases by alignment block, skipping intervening introns."""
    ref_by_pos = {}
    for start, end in read.get_blocks():
        try:
            ref_seq = ref_fasta.fetch(read.reference_name, start, end).upper()
        except Exception:
            continue
        for offset, base in enumerate(ref_seq):
            ref_by_pos[start + offset] = base
    return ref_by_pos


def build_strand_specific_transfer_matrices(sam_path, ref_fasta_path, n_reads=5000000):
    """
    Construct strand-specific nucleotide transition matrices for aligned matched regions
    """
    ref_fasta = pysam.FastaFile(ref_fasta_path)
    ref_bases = 'ATCG'

    pos_counts = np.zeros((4, 4), dtype=np.int64)
    neg_counts = np.zeros((4, 4), dtype=np.int64)

    with pysam.AlignmentFile(sam_path, "r") as samfile:
        for i, read in enumerate(samfile):
            if i >= n_reads:
                break
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue

            chrom = read.reference_name
            aligned_pairs = read.get_aligned_pairs(matches_only=True)
            AS = read.get_tag('AS')

            if aligned_pairs is None or AS < 40:
                continue

            ref_by_pos = fetch_aligned_reference(read, ref_fasta)
            for read_pos, ref_pos in aligned_pairs:
                if ref_pos is None:
                    continue

                try:
                    ref_base = ref_by_pos[ref_pos]
                    read_base = read.query_sequence[read_pos]

                    if ref_base in ref_bases and read_base in ref_bases:
                        ref_idx = BASE_TO_IDX[ref_base]
                        read_idx = BASE_TO_IDX[read_base]

                        if read.is_reverse:
                            neg_counts[ref_idx, read_idx] += 1
                        else:
                            pos_counts[ref_idx, read_idx] += 1
                except:
                    continue

    ref_fasta.close()

    def build_prob_matrix(counts):
        row_sums = counts.sum(axis=1, keepdims=True)
        safe_sums = np.where(row_sums == 0, 1, row_sums)
        prob_matrix = counts / safe_sums

        log_prob_matrix = np.log(np.maximum(prob_matrix, EPSILON))
        return prob_matrix, log_prob_matrix

    pos_prob_matrix, pos_log_prob_matrix = build_prob_matrix(pos_counts)
    neg_prob_matrix, neg_log_prob_matrix = build_prob_matrix(neg_counts)

    return (pos_prob_matrix, neg_prob_matrix, pos_log_prob_matrix, neg_log_prob_matrix)

def calculate_log_likelihood_ratio(seq, ref, bg_log_freq, log_trans_mat):
    """Calculate the LLR for the complete soft-clipped fragment."""
    if len(seq) != len(ref):
        raise ValueError("Read and reference fragments must have the same length")

    seq = seq.upper()
    ref = ref.upper()
    log_lr = 0.0
    for base, ref_base in zip(seq, ref):
        if base in BASE_TO_IDX and ref_base in BASE_TO_IDX:
            log_lr += (log_trans_mat[BASE_TO_IDX[ref_base], BASE_TO_IDX[base]]
                       - bg_log_freq[base])
    return float(log_lr)


def simulate_h0_distribution(ref, n_sim, bg_freq, bg_log_freq, log_trans_mat):
    """Simulate whole-fragment LLRs under H0, keeping the reference fixed."""
    ref = ref.upper()
    if not ref:
        raise ValueError("Reference fragment must not be empty")

    probs = np.array([bg_freq[b] for b in 'ATCG'], dtype=float)
    probs /= probs.sum()
    bg_logs = np.array([bg_log_freq[b] for b in 'ATCG'])
    score_table = np.vstack((log_trans_mat - bg_logs, np.zeros((1, 4))))
    ref_idx = np.array([BASE_TO_IDX.get(b, 4) for b in ref])

    # Each row is a simulated read; ref is the fixed candidate reference.
    seq_idx = np.random.choice(4, size=(n_sim, len(ref)), p=probs).astype(np.uint8)
    scores = score_table[ref_idx[None, :], seq_idx]
    del seq_idx
    return scores.sum(axis=1)


def bootstrap_supports_alignment(seq, ref, bg_log_freq, log_trans_mat, cache,
                                 alpha, observed_llr=None):
    """Test the whole fragment; cache one cutoff per reference composition."""
    if len(seq) != len(ref) or not seq:
        return False
    seq = seq.upper()
    ref = ref.upper()
    # Ignore the same non-ATCG positions in both observation and simulation.
    fixed_ref = ''.join(r if r in BASE_TO_IDX and b in BASE_TO_IDX else 'N'
                        for b, r in zip(seq, ref))
    if all(b == 'N' for b in fixed_ref):
        return False
    if observed_llr is None:
        observed_llr = calculate_log_likelihood_ratio(
            seq, ref, bg_log_freq, log_trans_mat)

    # Whole-fragment sums have the same null distribution for the same base counts.
    ref_key = tuple(fixed_ref.count(b) for b in 'ATCG')
    key = (ref_key, alpha, BOOTSTRAP_N_SIM)

    if key not in cache:
        bg_freq = {b: float(np.exp(bg_log_freq[b])) for b in 'ATCG'}
        sim_scores = np.empty(BOOTSTRAP_N_SIM)
        for start in range(0, BOOTSTRAP_N_SIM, BOOTSTRAP_BATCH_SIZE):
            end = min(start + BOOTSTRAP_BATCH_SIZE, BOOTSTRAP_N_SIM)
            sim_scores[start:end] = simulate_h0_distribution(
                fixed_ref, end - start, bg_freq, bg_log_freq, log_trans_mat)

        # For p = (1 + exceedances) / (N + 1), compare with the rank-th largest score.
        rank = int(alpha * (BOOTSTRAP_N_SIM + 1))
        sim_scores.sort()
        cutoff = float(sim_scores[-rank])
        if len(cache) >= BOOTSTRAP_CACHE_SIZE:
            cache.clear()
        cache[key] = cutoff

    return observed_llr > cache[key]


def calculate_mismatch(seq1, seq2, is_reverse,tolerate_c2u=False):
    """
    calculate mismatch number
    '+' strand allow A→G mismatch; '-' allow T→C mismatch
    """
    if len(seq1) != len(seq2) or len(seq1) == 0:
        return 0

    mismatches = 0
    C2U_count = 0
    for a, b in zip(seq1, seq2):
        if a != b:
            if not is_reverse and a == 'A' and b == 'G':
                continue
            elif is_reverse and a == 'T' and b == 'C':
                continue
            else:
                mismatches += 1
                if not is_reverse and a == 'C' and b == 'T':
                    C2U_count += 1
                if is_reverse and a == 'G' and b == 'A':
                    C2U_count += 1
    if tolerate_c2u:
        return mismatches - (C2U_count/2)
    else:
        return mismatches


def terminal_clip_length(seq, ref, is_reverse, is_5_end, tolerate_c2u=False):
    """Count outermost consecutive mismatches, retaining permitted conversions."""
    seq, ref = seq.upper(), ref.upper()
    positions = range(len(seq)) if is_5_end else range(len(seq) - 1, -1, -1)
    clipped = 0
    for i in positions:
        if seq[i] in BASE_TO_IDX and ref[i] in BASE_TO_IDX:
            if calculate_mismatch(ref[i], seq[i], is_reverse, False) < 1:
                break
        clipped += 1
    return clipped


def create_shifted_read(read, is_5_end, ref_fasta, tolerate_c2u):
    """Move one boundary match into S; its previous AS contribution is +1."""
    cigar = list(read.cigartuples)
    s_idx, m_idx = (0, 1) if is_5_end else (len(cigar) - 1, len(cigar) - 2)

    cigar[s_idx] = (4, cigar[s_idx][1] + 1)
    cigar[m_idx] = (0, cigar[m_idx][1] - 1)

    shifted = copy.copy(read)
    shifted.cigartuples = cigar
    if is_5_end:
        shifted.reference_start += 1

    return shifted, 1

def process_softclip_from_splice(read, ref_fasta, splice_sites, exons_by_transcript, softclip_len, is_5_end, basetheta, transmatrix,tolerate_c2u):

    chrom = read.reference_name
    ref_start = read.reference_start
    ref_end = read.reference_end
    is_reverse = read.is_reverse
    seq = read.query_sequence

    junction_pos = ref_start if is_5_end else ref_end
    site_type = 'acceptor' if is_5_end else 'donor'


    if chrom not in splice_sites:
        return 0.0, None, 0, 0, 0, 0, ''
    if junction_pos not in splice_sites[chrom]:
        return 0.0, None, 0, 0, 0, 0, ''

    best_match = 0.0
    best_tid = None
    best_intron_len = 0
    best_exon_start = 0
    best_delta_AS = 0
    best_LLR = 0
    best_ref_seq = ''

    for site_info in splice_sites[chrom][junction_pos]:
        site_strand, tid, stype = site_info
        if stype != site_type:
            continue

        if (site_strand == '+' and not is_reverse) or (site_strand == '-' and is_reverse):
            exons = exons_by_transcript[tid]

            if is_5_end:
                prev_exon = None
                for i in range(len(exons)):
                    if exons[i][1] == junction_pos and i > 0:  # Find the acceptor position.
                        prev_exon = exons[i-1]
                        break
                if prev_exon is None:
                    continue

                exon_length = prev_exon[2] - prev_exon[1] + 1
                if softclip_len > exon_length:
                    continue

                exon_seq = ref_fasta.fetch(chrom, prev_exon[2] - softclip_len + 1, prev_exon[2] + 1).upper()
                softclip_seq = seq[:softclip_len]
                tmp_LLR = calculate_log_likelihood_ratio(softclip_seq,exon_seq,basetheta,transmatrix)

            else:
                next_exon = None
                for i in range(len(exons)):
                    if exons[i][2] + 1 == junction_pos and i < len(exons) - 1:  # Find the donor position.
                        next_exon = exons[i+1]
                        break
                if next_exon is None:
                    continue

                exon_length = next_exon[2] - next_exon[1] + 1
                if softclip_len > exon_length:
                    continue

                exon_seq = ref_fasta.fetch(chrom, next_exon[1], next_exon[1] + softclip_len).upper()
                softclip_seq = seq[-softclip_len:]
                tmp_LLR = calculate_log_likelihood_ratio(softclip_seq,exon_seq,basetheta,transmatrix)


            mismatch_count = calculate_mismatch(exon_seq, softclip_seq, is_reverse,tolerate_c2u)
            match_rate = (len(exon_seq) - mismatch_count) / len(exon_seq)

            if match_rate > best_match:
                best_match = match_rate
                best_tid = tid
                best_delta_AS = len(exon_seq) - int(mismatch_count + mismatch_count)
                best_LLR = tmp_LLR
                best_ref_seq = exon_seq


                if is_5_end:
                    for i in range(len(exons)):
                        if exons[i][1] == junction_pos and i > 0:
                            best_intron_len = exons[i][1] - exons[i-1][2] - 1
                            best_exon_start = prev_exon[2] - softclip_len + 1
                            break
                else:
                    for i in range(len(exons)):
                        if exons[i][2] + 1 == junction_pos and i < len(exons) - 1:
                            best_intron_len = exons[i+1][1] - exons[i][2] - 1
                            best_exon_start = next_exon[1]
                            break
        else:
            continue

    return best_match, best_tid, best_intron_len, best_exon_start, best_delta_AS, best_LLR, best_ref_seq


def process_softclip_from_genome(read, ref_fasta,tolerate_c2u):
    regions = []
    cigar_tuples = read.cigartuples
    read_seq = read.query_sequence
    chrom = read.reference_name
    is_reverse = read.is_reverse

    # 1. softclip (5'S)
    if cigar_tuples[0][0] == 4 and cigar_tuples[0][1] >= min_softclip:  # 4 = BAM_CSOFT_CLIP
        length = cigar_tuples[0][1]
        ref_start = read.reference_start
        soft_ref_start = max(0, ref_start - length)
        soft_ref_end = ref_start
        mismatch_count = 0
        try:
            ref_seq = ref_fasta.fetch(chrom, soft_ref_start, soft_ref_end).upper()
            read_seq_clip = read_seq[:length]

            if len(ref_seq) == length and len(read_seq_clip) == length:
                mismatch_count = calculate_mismatch(ref_seq, read_seq_clip, is_reverse,tolerate_c2u)
                matchratio = round(1 - (mismatch_count - 1)/(length - 1), 3) if mismatch_count - 1 > 0 else 1
                if matchratio >= 0.5: 
                    regions.append((read_seq_clip,ref_seq,mismatch_count,matchratio,'5prime'))

        except:
            pass

    # 2. softclip (3'S)
    if cigar_tuples[-1][0] == 4 and cigar_tuples[-1][1] >= min_softclip:  # 4 = BAM_CSOFT_CLIP
        length = cigar_tuples[-1][1]
        ref_end = read.reference_end
        soft_ref_start = ref_end
        soft_ref_end = ref_end + length
        mismatch_count = 0
        try:
            ref_seq = ref_fasta.fetch(chrom, soft_ref_start, soft_ref_end).upper()
            read_seq_clip = read_seq[-length:]

            if len(ref_seq) == length and len(read_seq_clip) == length:
                mismatch_count = calculate_mismatch(ref_seq, read_seq_clip, is_reverse,tolerate_c2u)
                matchratio = round(1 - (mismatch_count - 1)/(length - 1), 3) if mismatch_count - 1 > 0 else 1
                if matchratio >= 0.5:
                    regions.append((read_seq_clip,ref_seq,mismatch_count,matchratio,'3prime'))
        except:
            pass

    return regions


def merge_cigar_tuples(cigar_tuples):
    """Merge adjacent operations of the same type."""
    merged = []
    for op, length in cigar_tuples:
        if merged and merged[-1][0] == op:
            merged[-1] = (op, merged[-1][1] + length)
        else:
            merged.append((op, length))
    return merged

def correct_softclip_chunk(reads, ref_fasta_path,
                          pos_log_theta, neg_log_theta,
                          pos_log_transfer_matrix, neg_log_transfer_matrix,
                          bootstrap_cache_plus, bootstrap_cache_minus,exons_by_transcript, splice_sites,tolerate_c2u, alpha=0.005):
    ref_fasta = pysam.FastaFile(ref_fasta_path)
    total_processed_reads = []
    process_in_splice = []

    for read in reads:
        chrome = read.reference_name
        seq = read.query_sequence
        cigar_tuples = read.cigartuples
        AS = read.get_tag('AS')
        is_reverse = read.is_reverse
        log_theta = neg_log_theta if is_reverse else pos_log_theta
        log_transfer_matrix = neg_log_transfer_matrix if is_reverse else pos_log_transfer_matrix
        bootstrap_cache = bootstrap_cache_minus if is_reverse else bootstrap_cache_plus


        # Adjust AS for the converted-genome mapping strategy.

        if ('A' in seq and not is_reverse) or ('T' in seq and is_reverse):
            if tolerate_c2u:
                ref_by_pos = fetch_aligned_reference(read, ref_fasta)
                for read_pos, ref_pos in read.get_aligned_pairs(matches_only=True):
                    if read_pos is None or ref_pos is None:
                        continue
                    read_base = seq[read_pos]
                    ref_base = ref_by_pos[ref_pos]
                    if not is_reverse and ref_base == 'C' and read_base == 'T':
                        AS += 1
                    elif is_reverse and ref_base == 'G' and read_base == 'A':
                        AS += 1
            if AS < 0:
                AS = 0

        read.set_tag('AS', AS, value_type='i')

        ################## correct fake softclip in splice junctions by match previous or next exon  #######################################
        ################## correct fake softclip in splice junctions by match previous or next exon   #######################################
        # extract softclip
        read.set_tag('SC', "NO")
        softclip_5 = 0
        softclip_3 = 0
        if cigar_tuples[0][0] == 4:  # 5' softclip
            softclip_5 = cigar_tuples[0][1]
        if cigar_tuples[-1][0] == 4:  # 3' softclip
            softclip_3 = cigar_tuples[-1][1]

        # require softclip length >= min_softclip
        if softclip_5 < min_softclip and softclip_3 < min_softclip:
            read.set_tag('SC', "NO_short")
            total_processed_reads.append(read)
            continue

        shift_5 = shift_3 = 0
        old_score_5 = old_score_3 = 0
        trim_5 = trim_3 = 0
        sites = splice_sites.get(chrome, {})

        # process 5' softclip
        best_5_match, best_5_tid, best_5_intron_len, best_5_exon_start, best_5_delta_AS, best_5_LLR, best_5_ref_seq = 0.0, None, 0, 0, 0, 0, ''


        if softclip_5 >= min_softclip:
            best_5_match, best_5_tid, best_5_intron_len, best_5_exon_start, best_5_delta_AS, best_5_LLR, best_5_ref_seq = process_softclip_from_splice(
                read, ref_fasta, splice_sites, exons_by_transcript, softclip_5, True, log_theta, log_transfer_matrix,tolerate_c2u)

            # Try the annotated +1 boundary only when the original boundary is unannotated.
            if best_5_tid is None and read.reference_start not in sites and read.reference_start + 1 in sites:
                shifted, moved_score = create_shifted_read(read, True, ref_fasta, tolerate_c2u)
                if shifted is not None:
                    candidate = process_softclip_from_splice(
                        shifted, ref_fasta, splice_sites, exons_by_transcript,
                        softclip_5 + 1, True, log_theta, log_transfer_matrix, tolerate_c2u)
                    if candidate[1] is not None and candidate[0] > 0.5:
                        best_5_match, best_5_tid, best_5_intron_len, best_5_exon_start, best_5_delta_AS, best_5_LLR, best_5_ref_seq = candidate
                        shift_5, old_score_5 = 1, moved_score


        #  3' softclip
        best_3_match, best_3_tid, best_3_intron_len, best_3_exon_start, best_3_delta_AS, best_3_LLR, best_3_ref_seq = 0.0, None, 0, 0, 0, 0, ''

        if softclip_3 >= min_softclip:
            best_3_match, best_3_tid, best_3_intron_len, best_3_exon_start, best_3_delta_AS, best_3_LLR, best_3_ref_seq = process_softclip_from_splice(
                read, ref_fasta, splice_sites, exons_by_transcript, softclip_3, False, log_theta, log_transfer_matrix,tolerate_c2u)

            if best_3_tid is None and read.reference_end not in sites and read.reference_end - 1 in sites:
                shifted, moved_score = create_shifted_read(read, False, ref_fasta, tolerate_c2u)
                if shifted is not None:
                    candidate = process_softclip_from_splice(
                        shifted, ref_fasta, splice_sites, exons_by_transcript,
                        softclip_3 + 1, False, log_theta, log_transfer_matrix, tolerate_c2u)
                    if candidate[1] is not None and candidate[0] > 0.5:
                        best_3_match, best_3_tid, best_3_intron_len, best_3_exon_start, best_3_delta_AS, best_3_LLR, best_3_ref_seq = candidate
                        shift_3, old_score_3 = 1, moved_score


        # only correct read with soft-clipped bases >= min_softclip
        candidate_len_5 = softclip_5 + shift_5
        candidate_len_3 = softclip_3 + shift_3
        if softclip_5 >= min_softclip:
            modify_5 = best_5_match == 1 or (
                best_5_match >= match_ratio and bootstrap_supports_alignment(
                    seq[:candidate_len_5], best_5_ref_seq, log_theta, log_transfer_matrix,
                    bootstrap_cache, alpha, observed_llr=best_5_LLR))
        else:
            modify_5 = False

        if softclip_3 >= min_softclip:
            modify_3 = best_3_match == 1 or (
                best_3_match >= match_ratio and bootstrap_supports_alignment(
                    seq[-candidate_len_3:], best_3_ref_seq, log_theta, log_transfer_matrix,
                    bootstrap_cache, alpha, observed_llr=best_3_LLR))
        else:
            modify_3 = False

        # Clip terminal errors only AFTER the complete candidate has passed its test.
        if modify_5:
            trim_5 = terminal_clip_length(seq[:candidate_len_5], best_5_ref_seq,
                                          is_reverse, True, tolerate_c2u=False)
            kept_5 = candidate_len_5 - trim_5
            modify_5 = kept_5 >= min_softclip and kept_5 > shift_5
            if modify_5:
                mm = calculate_mismatch(best_5_ref_seq[trim_5:], seq[trim_5:candidate_len_5],
                                        is_reverse, tolerate_c2u)
                best_5_delta_AS = kept_5 - 2 * mm - old_score_5
        if modify_3:
            trim_3 = terminal_clip_length(seq[-candidate_len_3:], best_3_ref_seq,
                                          is_reverse, False, tolerate_c2u=False)
            kept_3 = candidate_len_3 - trim_3
            modify_3 = kept_3 >= min_softclip and kept_3 > shift_3
            if modify_3:
                mm = calculate_mismatch(best_3_ref_seq[:kept_3], seq[-candidate_len_3:][:kept_3],
                                        is_reverse, tolerate_c2u)
                best_3_delta_AS = kept_3 - 2 * mm - old_score_3

        max_str = str(round(max(best_5_match,best_3_match),3))

        if not modify_5 and not modify_3:
            read.set_tag('SC', "NO_splice_"+max_str)
            process_in_splice.append(read)
            continue

        if modify_5 and modify_3:
            if best_5_tid != best_3_tid:
                if best_5_match > best_3_match:
                    modify_3 = False
                elif best_3_match > best_5_match:
                    modify_5 = False
                else:
                    if softclip_5 > softclip_3:
                        modify_3 = False
                    elif softclip_3 > softclip_5:
                        modify_5 = False
                    else:
                        modify_5 = False
                        modify_3 = False

        if not modify_5 and not modify_3:
            read.set_tag('SC', "diff_transcript")
            total_processed_reads.append(read)
            continue

        # Both shifts may consume the same short M block; keep at least one base in it.
        if modify_5 and modify_3 and shift_5 and shift_3 and len(cigar_tuples) == 3:
            if cigar_tuples[1][1] <= 2:
                if best_3_match > best_5_match:
                    modify_5 = False
                else:
                    modify_3 = False

        corrected_5_bp = softclip_5 - trim_5 if modify_5 else 0
        corrected_3_bp = softclip_3 - trim_3 if modify_3 else 0

        # correct CIGAR
        new_cigar = list(cigar_tuples)
        new_pos = read.reference_start

        if modify_5:
            new_cigar.pop(0)
            if shift_5:
                new_cigar[0] = (0, new_cigar[0][1] - 1)
            new_cigar[0:0] = [(0, kept_5), (3, best_5_intron_len)]
            if trim_5:
                new_cigar.insert(0, (4, trim_5))
            new_pos = best_5_exon_start + trim_5
            current_AS = read.get_tag('AS')
            read.set_tag('AS', int(current_AS + best_5_delta_AS) - 8, value_type='i')

        if modify_3:
            new_cigar.pop()
            if shift_3:
                new_cigar[-1] = (0, new_cigar[-1][1] - 1)
            new_cigar.extend([(3, best_3_intron_len), (0, kept_3)])
            if trim_3:
                new_cigar.append((4, trim_3))
            current_AS = read.get_tag('AS')
            read.set_tag('AS', int(current_AS + best_3_delta_AS) - 8, value_type='i')

        # update read information
        read.cigartuples = new_cigar
        read.reference_start = new_pos

        shift_tag = "shift_" if (modify_5 and shift_5) or (modify_3 and shift_3) else ""
        if modify_5 and modify_3:
            read.set_tag('SC', "YES_53_" + shift_tag + max_str + "_CORR_" + str(corrected_5_bp + corrected_3_bp) + "bp")
            total_processed_reads.append(read)
            continue

        elif modify_5:
            read.set_tag('SC', "YES_5_" + shift_tag + max_str + "_CORR_" + str(corrected_5_bp) + "bp")

        elif modify_3:
            read.set_tag('SC', "YES_3_" + shift_tag + max_str + "_CORR_" + str(corrected_3_bp) + "bp")
        process_in_splice.append(read)


    # Rescue soft clips in the adjacent genomic sequence.
    for read in process_in_splice:
        chrome = read.reference_name
        seq = read.query_sequence
        cigar_tuples = read.cigartuples
        is_reverse = read.is_reverse
        log_theta = neg_log_theta if is_reverse else pos_log_theta
        log_transfer_matrix = neg_log_transfer_matrix if is_reverse else pos_log_transfer_matrix
        bootstrap_cache = bootstrap_cache_minus if is_reverse else bootstrap_cache_plus

        soft_regions = process_softclip_from_genome(read, ref_fasta,tolerate_c2u)

        if not soft_regions:
            sctag = read.get_tag('SC')
            if 'YES' in sctag:
                total_processed_reads.append(read)
                continue
            else:
                read.set_tag('SC', "NO")
                total_processed_reads.append(read)
                continue


        AS_changes = 0
        corrected = False
        corrected_bp = 0
        new_cigar = list(cigar_tuples)
        result_record = []

        for (soft_seq, ref_seq,mismatch_count,matchratio, sc_type) in soft_regions:
            soft_len = len(soft_seq)

            if soft_len < min_softclip:
                continue

            # Keep the existing full-fragment acceptance rules.
            accepted = matchratio >= 0.95
            if not accepted and matchratio >= match_ratio:
                accepted = bootstrap_supports_alignment(
                    soft_seq, ref_seq, log_theta, log_transfer_matrix, bootstrap_cache, alpha)
            if not accepted:
                result_record.append((False, matchratio))
                continue

            trim = terminal_clip_length(soft_seq, ref_seq, is_reverse,
                                        sc_type == '5prime', tolerate_c2u=False)
            kept = soft_len - trim
            if kept < min_softclip:
                result_record.append((False, matchratio))
                continue
            if sc_type == '5prime':
                kept_seq, kept_ref = soft_seq[trim:], ref_seq[trim:]
                new_cigar[:1] = ([(4, trim)] if trim else []) + [(0, kept)]
                read.reference_start -= kept
            else:
                kept_seq, kept_ref = soft_seq[:kept], ref_seq[:kept]
                new_cigar[-1:] = [(0, kept)] + ([(4, trim)] if trim else [])
            kept_mm = calculate_mismatch(kept_ref, kept_seq, is_reverse, tolerate_c2u)
            AS_changes += kept - 2 * kept_mm
            corrected = True
            corrected_bp += kept
            result_record.append((True, matchratio))

        if corrected:
            # rebuild cigar
            merged_cigar = merge_cigar_tuples(new_cigar)
            cigar_str_parts = []
            for op, length in merged_cigar:
                op_char = "MIDNSHP=X"[op]
                cigar_str_parts.append(f"{length}{op_char}")
            read.cigarstring = ''.join(cigar_str_parts)

            # re-score
            current_AS = read.get_tag('AS')
            read.set_tag('AS', int(current_AS + AS_changes), value_type='i')

        if result_record:
            splicetag = read.get_tag('SC')
            if any(r[0] for r in result_record):
                best_match = max([r[1] for r in result_record if r[0]])
                read.set_tag('SC', splicetag + "_YES_" + str(best_match) + "_CORR_" + str(corrected_bp) + "bp")
            else:
                best_match = max([r[1] for r in result_record])
                read.set_tag('SC', splicetag + "_NO_" +  str(best_match))
        total_processed_reads.append(read)

    ref_fasta.close()
    sorted_read_list = sorted(total_processed_reads, key=lambda read: (int(read.query_name[2:].split('_')[0]),int(read.query_name.split('_')[1]) if '_' in read.query_name else 0))
    return sorted_read_list

def correct_softclip_chunk_wrapper(args):
    chunk_id, read_strs, ref_fasta_path,  pos_log_theta, neg_log_theta, pos_log_transfer_matrix,\
    neg_log_transfer_matrix, bootstrap_cache_plus, bootstrap_cache_minus, alpha, samheader,exons_by_transcript, splice_sites,tolerate_c2u = args

    # Seed by chunk, so process scheduling does not change the random stream.
    np.random.seed((22 + chunk_id) % (2 ** 32))

    reads = []
    header = pysam.AlignmentHeader().from_dict(samheader)

    for read_str in read_strs:
        read = pysam.AlignedSegment().fromstring(read_str, header)
        reads.append(read)

    read_strs.clear()  # Read objects now hold the data; release the SAM input strings.

    processed_reads = correct_softclip_chunk(
        reads, ref_fasta_path,
        pos_log_theta, neg_log_theta,
        pos_log_transfer_matrix, neg_log_transfer_matrix,
        bootstrap_cache_plus, bootstrap_cache_minus,exons_by_transcript, splice_sites,tolerate_c2u, alpha
    )

    processed_strs = [read.tostring() for read in processed_reads]

    return chunk_id, processed_strs

def write_completed_results(outfile, submitted_tasks, results_cache, next_write_id, header):
    """
    output results
    """
    i = 0
    while i < len(submitted_tasks):
        chunk_id, future = submitted_tasks[i]

        if future.ready():
            result_id, processed_strs = future.get()
            results_cache[result_id] = processed_strs
            submitted_tasks.pop(i)
        else:
            i += 1

    while next_write_id in results_cache:
        processed_strs = results_cache.pop(next_write_id)
        for read_str in processed_strs:
            read = pysam.AlignedSegment().fromstring(read_str, header)
            outfile.write(read)
        next_write_id += 1

    return next_write_id

def main(sam_input, reference_fasta, output_sam, chunk_size=chunk_size, alpha=pvalue, n_processes=Threads):

    exons_by_transcript, splice_sites = load_exons_and_splice_sites(gtffile)

    print("Precomputing strand-specific parameters from first 10,000,000 aligned reads...")

    pos_theta, neg_theta, pos_log_theta, neg_log_theta = compute_softclip_base_freq(sam_input, n_reads=10000000)
    print(f"\nPositive strand softclip region base frequencies (H0): {pos_theta}")
    print(f"Negative strand softclip region base frequencies (H0): {neg_theta}")

    pos_transfer_matrix, neg_transfer_matrix, pos_log_transfer_matrix, neg_log_transfer_matrix = build_strand_specific_transfer_matrices(
        sam_input, reference_fasta, n_reads=10000000
    )

    print("\nPositive strand transfer matrix (H1:ref → read):")
    print("       A      T      C      G")
    for i, ref_base in enumerate('ATCG'):
        row = pos_transfer_matrix[i]
        print(f"{ref_base}  {row[0]:.4f}  {row[1]:.4f}  {row[2]:.4f}  {row[3]:.4f}")

    print("\nNegative strand transfer matrix (H1:ref → read):")
    print("       A      T      C      G")
    for i, ref_base in enumerate('ATCG'):
        row = neg_transfer_matrix[i]
        print(f"{ref_base}  {row[0]:.4f}  {row[1]:.4f}  {row[2]:.4f}  {row[3]:.4f}")


    tolerate_c2u = pos_transfer_matrix[BASE_TO_IDX['C']][BASE_TO_IDX['T']] > C2Ulevel or neg_transfer_matrix[BASE_TO_IDX['G']][BASE_TO_IDX['A']] > C2Ulevel


    print(f"\nWhole-fragment conditional H0 bootstrap: {BOOTSTRAP_N_SIM} simulations per reference condition")
    # Reference-dependent cutoffs are computed only when a candidate needs them.
    # Each chunk receives separate caches for the two strands.
    bootstrap_cache_plus = {}
    bootstrap_cache_minus = {}


    print(f"\nStreaming processing in chunks (size={chunk_size}) using {n_processes} processes...")

    chunk_id = 0
    total_reads = 0
    submitted_tasks = []
    results_cache = {}
    next_write_id = 0

    pool = mp.Pool(processes=n_processes)

    try:
        with pysam.AlignmentFile(sam_input, "rb") as samfile:
            header_dict = samfile.header.to_dict()

            with pysam.AlignmentFile(output_sam, "w", header=samfile.header) as outfile:
                buffer = []
                max_pending_tasks = n_processes * 2

                for read in samfile:
                    buffer.append(read)
                    total_reads += 1

                    if len(buffer) >= chunk_size:
                        # Count unfinished tasks and results waiting for ordered output.
                        while len(submitted_tasks) + len(results_cache) >= max_pending_tasks:
                            next_write_id = write_completed_results(
                                outfile, submitted_tasks, results_cache, next_write_id, samfile.header
                            )
                            if len(submitted_tasks) + len(results_cache) >= max_pending_tasks:
                                time.sleep(0.01)

                        read_strs = [r.tostring() for r in buffer]

                        task_args = (
                            chunk_id, read_strs, reference_fasta,
                            pos_log_theta, neg_log_theta,
                            pos_log_transfer_matrix, neg_log_transfer_matrix,
                            bootstrap_cache_plus, bootstrap_cache_minus, alpha, header_dict, exons_by_transcript, splice_sites,tolerate_c2u
                        )

                        future = pool.apply_async(correct_softclip_chunk_wrapper, (task_args,))
                        submitted_tasks.append((chunk_id, future))

                        buffer = []
                        chunk_id += 1

                        next_write_id = write_completed_results(
                            outfile, submitted_tasks, results_cache, next_write_id, samfile.header
                        )

                # Process the last batch
                if buffer:
                    while len(submitted_tasks) + len(results_cache) >= max_pending_tasks:
                        next_write_id = write_completed_results(
                            outfile, submitted_tasks, results_cache, next_write_id, samfile.header
                        )
                        if len(submitted_tasks) + len(results_cache) >= max_pending_tasks:
                            time.sleep(0.01)
                    read_strs = [r.tostring() for r in buffer]
                    task_args = (
                        chunk_id, read_strs, reference_fasta,
                        pos_log_theta, neg_log_theta,
                        pos_log_transfer_matrix, neg_log_transfer_matrix,
                        bootstrap_cache_plus, bootstrap_cache_minus, alpha, header_dict,exons_by_transcript, splice_sites,tolerate_c2u
                    )

                    future = pool.apply_async(correct_softclip_chunk_wrapper, (task_args,))
                    submitted_tasks.append((chunk_id, future))
                    chunk_id += 1

                # Process all remaining tasks
                while submitted_tasks:
                    next_write_id = write_completed_results(
                        outfile, submitted_tasks, results_cache, next_write_id, samfile.header)
                    if submitted_tasks:
                        time.sleep(0.01)

    finally:
        pool.terminate()
        pool.join()
        print(f"Processed {total_reads} reads in total")

    print(f"\nDone! Corrected softclips written to {output_sam}")

    if outbam:
        pysam.sort("-o", output_sam[:-3] + "bam", "-@", "8", output_sam)
        pysam.index(output_sam[:-3] + "bam")
        os.remove(output_sam)

if __name__ == "__main__":
    if outname == "none":
        output_sam = insam[:-4] + "_correct.sam"
    else:
        output_sam = outname

    start_time = time.time()
    main(insam, genomefile, output_sam)
    print(f"Total execution time: {time.time() - start_time:.2f} seconds")
