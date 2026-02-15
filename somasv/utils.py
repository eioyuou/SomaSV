import os
import pysam
import numpy as np
from numba import jit, njit
from math import floor, ceil


@njit
def calculate_statistics_numba(starts, mapqs, event_sizes):
    starts_array = np.array(starts, dtype=np.float64)
    mapqs_array = np.array(mapqs, dtype=np.float64)
    event_sizes_array = np.array(event_sizes, dtype=np.float64)

    starts_std = np.std(starts_array)
    event_size_std = np.std(event_sizes_array)

    mapq_mean = np.mean(mapqs_array)
    event_size_mean = np.mean(event_sizes_array)
    event_size_median = np.median(event_sizes_array)

    return starts_std, mapq_mean, event_size_std, event_size_median, event_size_mean


@njit
def calculate_median_numba(values):
    if len(values) == 0:
        return 0.0
    sorted_values = np.sort(np.array(values, dtype=np.float64))
    n = len(sorted_values)
    if n % 2 == 0:
        return (sorted_values[n // 2 - 1] + sorted_values[n // 2]) / 2.0
    else:
        return sorted_values[n // 2]


@njit
def calculate_mean_numba(values):
    if len(values) == 0:
        return 0.0
    return np.mean(np.array(values, dtype=np.float64))


def split_intervals(chrom_length, interval_size):
    """ Split chromosome length into multiple intervals, return (start, end) pairs """
    intervals = []
    for start in range(0, chrom_length, interval_size):
        end = min(start + interval_size, chrom_length)
        intervals.append((start, end))
    return intervals


def get_chrom_lengths(fasta_file):
    fasta = pysam.FastaFile(fasta_file)
    chrom_lengths = {ref: fasta.get_reference_length(ref) for ref in fasta.references}
    return chrom_lengths


def get_chrom_lengths_from_bam(bam_file):
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        chrom_lengths = {ref: bam.get_reference_length(ref) for ref in bam.references}
    return chrom_lengths


def get_first_24_chrom_lengths(chrom_lengths):
    chrom_list = ['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY']
    return {chrom: chrom_lengths[chrom] for chrom in chrom_list if chrom in chrom_lengths}


def validate_coordinates(start, end):
    if start > end:
        return end, start
    return start, end


def ensure_float(val):
    if isinstance(val, str):
        try:
            return float(val)
        except:
            return 0
    return val if val is not None else 0


def ensure_int(val):
    if isinstance(val, str):
        try:
            return int(float(val))
        except:
            return 0
    return int(val) if val is not None else 0


def get_breakpoint_af(bp):
    if 'allele_fractions' not in bp:
        return 0.0, 0.0

    tumor_af = max(ensure_float(bp['allele_fractions']['tumor'][0]),
                   ensure_float(bp['allele_fractions']['tumor'][1]))
    normal_af = max(ensure_float(bp['allele_fractions']['normal'][0]),
                    ensure_float(bp['allele_fractions']['normal'][1]))

    return tumor_af, normal_af


def get_breakpoint_support(bp):
    tumor_support = ensure_float(bp['read_support_counts'].get('tumor', 0))
    normal_support = ensure_float(bp['read_support_counts'].get('normal', 0))
    support_ratio = tumor_support / (normal_support + 1)

    return tumor_support, normal_support, support_ratio


def get_breakpoint_mapq(bp):
    origin_mapq = ensure_float(bp.get('start_cluster', {}).get('stats', {}).get('mapq_mean', 0))
    end_mapq = ensure_float(bp.get('end_cluster', {}).get('stats', {}).get('mapq_mean', 0))

    return origin_mapq, end_mapq


def get_breakpoint_clustered_reads(bp):
    clustered_tumor = ensure_float(bp['total_read_counts'].get('tumor', 0))
    clustered_normal = ensure_float(bp['total_read_counts'].get('normal', 0))

    return clustered_tumor, clustered_normal


def get_breakpoint_std_dev(bp):
    origin_starts_std_dev = ensure_float(bp.get('start_cluster', {}).get('stats', {}).get('starts_std_dev', 0))
    end_starts_std_dev = ensure_float(bp.get('end_cluster', {}).get('stats', {}).get('starts_std_dev', 0))

    return origin_starts_std_dev, end_starts_std_dev