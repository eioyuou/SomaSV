import time
import multiprocessing
import numpy as np

from somasv.utils import split_intervals
from somasv.read_processing import get_potential_breakpoints_task
from somasv.coverage import merge_coverage_arrays_efficient


def merge_breakpoint_dicts(dict1, dict2):
    merged_dict = {}
    for chr_key, breakpoints in dict1.items():
        if chr_key not in merged_dict:
            merged_dict[chr_key] = breakpoints
        else:
            merged_dict[chr_key].extend(breakpoints)
    for chr_key, breakpoints in dict2.items():
        if chr_key not in merged_dict:
            merged_dict[chr_key] = breakpoints
        else:
            merged_dict[chr_key].extend(breakpoints)
    return merged_dict


def run_multiprocess_for_all_chromosomes(aln_filename_tumor, aln_filename_normal, length, mapq,
                                         contig_order, chrom_lengths, internal_num, coverage_binsize, platform,
                                         num_processes=None):

    all_tasks_tumor = []
    all_tasks_normal = []

    print(f"Processing {len(chrom_lengths)} chromosomes")

    for contig, chrom_length in chrom_lengths.items():
        intervals = split_intervals(chrom_length, internal_num)

        for (start, end) in intervals:
            task_args_tumor = (
                aln_filename_tumor, length, mapq, 'tumor', contig_order,
                contig, start, end, coverage_binsize, platform
            )
            task_args_normal = (
                aln_filename_normal, length, mapq, 'normal', contig_order,
                contig, start, end, coverage_binsize, platform
            )
            all_tasks_tumor.append(task_args_tumor)
            all_tasks_normal.append(task_args_normal)

    print(f"Created {len(all_tasks_tumor)} tasks")

    with multiprocessing.Pool(processes=num_processes) as pool:
        results_tumor = pool.map(get_potential_breakpoints_task, all_tasks_tumor)
        results_normal = pool.map(get_potential_breakpoints_task, all_tasks_normal)

    print(f"Finished processing {len(chrom_lengths)} chromosomes")
    time_count = time.time()

    breakpoints_tumor = [result[0] for result in results_tumor]
    coverage_tumor = [result[1] for result in results_tumor]
    read_lengths_tumor = [result[2] for result in results_tumor]

    breakpoints_normal = [result[0] for result in results_normal]
    coverage_normal = [result[1] for result in results_normal]
    read_lengths_normal = [result[2] for result in results_normal]

    total_read_length_tumor = sum(read_lengths_tumor)
    total_read_length_normal = sum(read_lengths_normal)

    print(f"Time taken for processing: {time.time() - time_count:.2f} seconds")
    time_count = time.time()

    merged_result_tumor = {}
    merged_result_normal = {}
    print(f"Merging tumor breakpoint results")
    for result in breakpoints_tumor:
        for chr_key, breakpoints in result.items():
            merged_result_tumor.setdefault(chr_key, []).extend(breakpoints)
    print(f"Time taken for merging tumor breakpoints: {time.time() - time_count:.2f} seconds")
    print(f"Merging normal breakpoint results")
    time_count = time.time()
    for result in breakpoints_normal:
        for chr_key, breakpoints in result.items():
            merged_result_normal.setdefault(chr_key, []).extend(breakpoints)
    print(f"Time taken for merging normal breakpoints: {time.time() - time_count:.2f} seconds")
    print(f"Finished merging breakpoint results")
    time_count = time.time()
    merged_result = merge_breakpoint_dicts(merged_result_tumor, merged_result_normal)
    print(f"Time taken for merging all breakpoints: {time.time() - time_count:.2f} seconds")
    print(f"Total merged breakpoints: {len(merged_result)}")

    merged_coverage_arrays = {
        'tumor': {},
        'normal': {}
    }
    print(f"Merging tumor coverage arrays")
    time_count = time.time()
    for contig, chrom_length in chrom_lengths.items():
        merged_coverage_arrays['tumor'][contig] = merge_coverage_arrays_efficient(
            coverage_tumor, contig, chrom_length, coverage_binsize
        )
        merged_coverage_arrays['normal'][contig] = merge_coverage_arrays_efficient(
            coverage_normal, contig, chrom_length, coverage_binsize
        )

    print(f"Time taken for merging coverage arrays: {time.time() - time_count:.2f} seconds")
    print(f"Finished merging coverage arrays")

    genome_size = sum(chrom_lengths.values())
    tumor_global_coverage = total_read_length_tumor / genome_size if genome_size > 0 else 0
    normal_global_coverage = total_read_length_normal / genome_size if genome_size > 0 else 0

    print(f"Global coverage - Tumor: {tumor_global_coverage:.2f}X, Normal: {normal_global_coverage:.2f}X")

    return merged_result, merged_coverage_arrays, {
        'tumor_total_read_length': total_read_length_tumor,
        'normal_total_read_length': total_read_length_normal,
        'genome_size': genome_size,
        'tumor_global_coverage': tumor_global_coverage,
        'normal_global_coverage': normal_global_coverage
    }